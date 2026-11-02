"""Add 220-degree coplanar electrodes to the rings in the SiN mask.

The source GDS is imported with Layer 10 removed.  Every isolated ring and
only the bus-remote ring of each double-ring device receives two concentric
220-degree electrodes.  The top single-ring and lower double-ring groups use
independent 0.1 um closest-edge gap sweeps.  All selected rings use the same
one-side circular-pad geometry.
"""

from pathlib import Path
import csv

import gdsfactory as gf
import numpy as np


SOURCE_GDS = Path("/Users/yuyue/Downloads/123024_sin_250nm_final.gds")
OUTPUT_GDS = Path(
    "/Users/yuyue/Documents/Github/layout/gds/"
    "123024_sin_250nm_final_with_220deg_electrodes.gds"
)
RING_MAP_CSV = OUTPUT_GDS.with_name("123024_sin_250nm_ring_map.csv")

OPTICAL_LAYER = (1, 0)
IGNORED_LAYER = (10, 0)
ELECTRODE_LAYER = (20, 0)

ELECTRODE_ANGLE = 220.0
ELECTRODE_GAP = 7.2
SINGLE_RING_GAP_START = 7.2
DOUBLE_RING_GAP_START = 7.2
ELECTRODE_GAP_STEP = 0.1
DOUBLE_RING_DEVICES_PER_ROW = 8
ELECTRODE_WIDTH = 10.0
JUNCTION_SMOOTHING = 2.0
OPTICAL_CLEARANCE = 2.0
OPTICAL_RING_WIDTH = 1.6
PAD_SIZE = (80.0, 80.0)
PAD_CENTER_PITCH = 150.0
SINGLE_RING_PAD_DIAMETER = 80.0
SINGLE_RING_NEIGHBOR_DISTANCE = 300.0
PAIR_MAX_X_SEPARATION = 260.0
PAIR_MAX_Y_MISMATCH = 30.0

# Connected Layer-1 polygons in this size/area range are the optical rings.
RING_BBOX_MIN = 180.0
RING_BBOX_MAX = 240.0
RING_AREA_MIN = 500.0
RING_AREA_MAX = 2000.0


def find_rings(component: gf.Component) -> list[tuple[float, float, float]]:
    """Return ``(center_x, center_y, outer_radius)`` for Layer-1 rings."""
    region = component.get_region(OPTICAL_LAYER)
    region.merge()
    rings = []
    for polygon in region.each():
        bbox = polygon.bbox()
        width = bbox.width() / 1000
        height = bbox.height() / 1000
        area = polygon.area() / 1e6
        if not (
            RING_BBOX_MIN <= width <= RING_BBOX_MAX
            and RING_BBOX_MIN <= height <= RING_BBOX_MAX
            and abs(width - height) <= 1.0
            and RING_AREA_MIN <= area <= RING_AREA_MAX
        ):
            continue

        outer_radius = (width + height) / 4
        rings.append(
            (
                (bbox.left + bbox.right) / 2000,
                (bbox.bottom + bbox.top) / 2000,
                outer_radius,
            )
        )
    return sorted(rings, key=lambda ring: (ring[1], ring[0]))


def electrode_orientations(
    rings: list[tuple[float, float, float]],
) -> list[float]:
    """Place both concentric electrodes over the upper half of every ring."""
    return [90.0] * len(rings)


def electrode_radii(
    ring_outer_radius: float,
    electrode_gap: float,
    electrode_width: float,
) -> tuple[float, float]:
    """Return electrode radii centered on the optical ring waveguide.

    The optical waveguide centerline lies exactly halfway between the facing
    edges of the inner and outer electrodes.  The pad pair is positioned from
    these radii separately, which preserves this centering and the same swept
    gap through the straight connection region.
    """
    ring_center_radius = ring_outer_radius - OPTICAL_RING_WIDTH / 2
    center_offset = (electrode_gap + electrode_width) / 2
    return (
        ring_center_radius - center_offset,
        ring_center_radius + center_offset,
    )


def isolated_ring_flags(
    rings: list[tuple[float, float, float]],
) -> list[bool]:
    """Identify rings that are not part of a close two-ring device.

    The coupled rings in this mask are about 221 um apart, whereas the
    single-ring devices in the top row are about 400 um apart.
    """
    centers = np.asarray([(ring[0], ring[1]) for ring in rings], dtype=float)
    flags = []
    for index, center in enumerate(centers):
        distances = np.linalg.norm(centers - center, axis=1)
        distances[index] = np.inf
        flags.append(float(np.min(distances)) > SINGLE_RING_NEIGHBOR_DISTANCE)
    return flags


def double_ring_targets(
    rings: list[tuple[float, float, float]],
    single_flags: list[bool],
) -> dict[int, float]:
    """Return remote-ring indices and the pad side for every double ring.

    In this source mask the bus-coupled member of each pair has the larger
    radius (112.28 um), and the ring farther from the bus has the smaller
    radius (108.71 um).  The pad side is chosen away from the coupled ring:
    -1 for left and +1 for right.
    """
    remaining = {index for index, is_single in enumerate(single_flags) if not is_single}
    targets: dict[int, float] = {}
    while remaining:
        first = min(remaining)
        candidates = remaining - {first}
        if not candidates:
            raise RuntimeError("Unpaired ring remains in the double-ring map")

        first_xy = np.asarray(rings[first][:2], dtype=float)
        second = min(
            candidates,
            key=lambda index: float(
                np.linalg.norm(np.asarray(rings[index][:2], dtype=float) - first_xy)
            ),
        )
        dx = abs(rings[first][0] - rings[second][0])
        dy = abs(rings[first][1] - rings[second][1])
        if dx > PAIR_MAX_X_SEPARATION or dy > PAIR_MAX_Y_MISMATCH:
            raise RuntimeError(f"Could not form double-ring pair for ring {first}")

        remote = min((first, second), key=lambda index: rings[index][2])
        coupled = second if remote == first else first
        pad_side = 1.0 if rings[remote][0] > rings[coupled][0] else -1.0
        targets[remote] = pad_side
        remaining.remove(first)
        remaining.remove(second)

    return targets


def electrode_gap_assignments(
    rings: list[tuple[float, float, float]],
    single_flags: list[bool],
    double_targets: dict[int, float],
) -> dict[int, float]:
    """Assign the electrode-gap sweep across each device set.

    Double-ring targets are grouped into physical rows and swept left-to-right.
    The ten top single-ring devices use 7.2, 7.3, ..., 8.1 um.  Each
    eight-device double-ring row restarts at 7.2 um and ends at 7.9 um.
    """
    assignments: dict[int, float] = {}

    single_indices = sorted(
        (index for index, is_single in enumerate(single_flags) if is_single),
        key=lambda index: rings[index][0],
    )
    for position, index in enumerate(single_indices):
        assignments[index] = round(
            SINGLE_RING_GAP_START + position * ELECTRODE_GAP_STEP,
            10,
        )

    # The source arranges the lower double-ring array as four sets of eight
    # devices.  Within each set, the first four form the upper sub-row and the
    # next four form the lower sub-row.  Sorting by y then x preserves that
    # physical order and lets every eight-device set restart its sweep.
    double_indices = sorted(
        double_targets,
        key=lambda index: (rings[index][1], rings[index][0]),
    )
    rows = [
        double_indices[start : start + DOUBLE_RING_DEVICES_PER_ROW]
        for start in range(0, len(double_indices), DOUBLE_RING_DEVICES_PER_ROW)
    ]
    if any(len(row) != DOUBLE_RING_DEVICES_PER_ROW for row in rows):
        raise RuntimeError(
            "Double-ring target count is not divisible into eight-device rows"
        )

    for row in rows:
        for position, index in enumerate(row):
            assignments[index] = round(
                DOUBLE_RING_GAP_START + position * ELECTRODE_GAP_STEP,
                10,
            )

    return assignments


@gf.cell
def ring_coplanar_electrode(
    ring_radius: float,
    outward_angle: float = 0.0,
    electrode_gap: float = ELECTRODE_GAP,
    electrode_width: float = ELECTRODE_WIDTH,
    pad_center_pitch: float = PAD_CENTER_PITCH,
    layer=ELECTRODE_LAYER,
) -> gf.Component:
    """Return two concentric curved electrodes and two probe pads."""
    if electrode_gap <= 0 or electrode_width <= 0:
        raise ValueError("Electrode gap and width must be positive")
    if pad_center_pitch < PAD_SIZE[0]:
        raise ValueError("PAD_CENTER_PITCH must be at least the pad width")

    component = gf.Component()
    # Layer 1 lies between the two curved electrodes.  Their closest metal
    # edges are ELECTRODE_GAP apart, centered on the optical ring.
    inner_electrode_radius, outer_electrode_radius = electrode_radii(
        ring_outer_radius=ring_radius,
        electrode_gap=electrode_gap,
        electrode_width=electrode_width,
    )

    cross_section = gf.cross_section.cross_section(
        width=electrode_width,
        layer=layer,
    )
    for radius in (inner_electrode_radius, outer_electrode_radius):
        # Define the curved path around the ring center explicitly. The generic
        # gf.path.arc() origin is its path start, not its circle center, which
        # previously shifted every electrode by approximately one radius.
        theta = np.radians(
            np.linspace(
                outward_angle - ELECTRODE_ANGLE / 2,
                outward_angle + ELECTRODE_ANGLE / 2,
                257,
            )
        )
        arc = gf.Path(
            np.column_stack((radius * np.cos(theta), radius * np.sin(theta)))
        )
        component.add_ref(gf.path.extrude(arc, cross_section=cross_section))

    # Put one pad radially inside the optical ring and one radially outside.
    # Their centers lie on the ring's horizontal diameter and are separated
    # by exactly PAD_CENTER_PITCH.
    angle = np.radians(outward_angle)
    outward = np.array([np.cos(angle), np.sin(angle)])
    inner_pad_center = outward * (ring_radius - pad_center_pitch / 2)
    outer_pad_center = outward * (ring_radius + pad_center_pitch / 2)

    # Short radial connections meet the midpoint of each curved electrode. The
    # inner pad connects to the inner electrode and the outside pad connects
    # to the outer electrode, without crossing the 6 um coplanar gap.
    route_points = (
        [
            outward * inner_electrode_radius,
            inner_pad_center + outward * PAD_SIZE[0] / 2,
        ],
        [
            outward * outer_electrode_radius,
            outer_pad_center - outward * PAD_SIZE[0] / 2,
        ],
    )
    for points in route_points:
        path = gf.Path(np.asarray(points, dtype=float))
        component.add_ref(gf.path.extrude(path, cross_section=cross_section))

    for pad_center in (inner_pad_center, outer_pad_center):
        pad = component.add_ref(
            gf.components.rectangle(
                size=PAD_SIZE,
                centered=True,
                layer=layer,
            )
        )
        pad.move(tuple(pad_center))

    component.info["ring_radius"] = ring_radius
    component.info["electrode_gap"] = electrode_gap
    component.info["electrode_width"] = electrode_width
    component.info["pad_size"] = list(PAD_SIZE)
    component.info["pad_center_pitch"] = pad_center_pitch
    component.info["outward_angle"] = outward_angle
    return component


@gf.cell
def ring_same_side_pad_electrode(
    ring_radius: float,
    pad_side: float = -1.0,
    arc_style: str = "single_upper",
    electrode_gap: float = ELECTRODE_GAP,
    electrode_width: float = ELECTRODE_WIDTH,
    pad_center_pitch: float = PAD_CENTER_PITCH,
    pad_diameter: float = SINGLE_RING_PAD_DIAMETER,
    layer=ELECTRODE_LAYER,
) -> gf.Component:
    """Curved electrodes with one-side circular pads.

    ``single_upper`` draws the top arc from 200 to -20 degrees.
    ``double_side`` draws 110 to -110 degrees for right-side pads and mirrors
    that arc for
    left-side pads so the existing pad/feed position remains unchanged.
    """
    if electrode_gap <= 0 or electrode_width <= 0 or pad_diameter <= 0:
        raise ValueError("Electrode gap, width, and pad diameter must be positive")
    if pad_center_pitch < pad_diameter:
        raise ValueError("PAD_CENTER_PITCH must be at least the pad diameter")
    if pad_side not in (-1.0, 1.0):
        raise ValueError("pad_side must be -1.0 (left) or +1.0 (right)")
    if arc_style not in ("single_upper", "double_side"):
        raise ValueError("arc_style must be 'single_upper' or 'double_side'")

    component = gf.Component()
    inner_electrode_radius, outer_electrode_radius = electrode_radii(
        ring_outer_radius=ring_radius,
        electrode_gap=electrode_gap,
        electrode_width=electrode_width,
    )
    cross_section = gf.cross_section.cross_section(
        width=electrode_width,
        layer=layer,
    )

    # The top single-ring arc runs from 200 to -20 degrees.  For double rings,
    # the right-side geometry runs from 110 to -110 degrees and the left-side
    # geometry is mirrored about the y axis.  The feed always retains its
    # original horizontal connection coordinate.
    if arc_style == "single_upper":
        start_angle = 200.0
        end_angle = -20.0
        connection_angle = 180.0
    elif pad_side > 0:
        start_angle = 110.0
        end_angle = -110.0
        connection_angle = 0.0
    else:
        start_angle = 70.0
        end_angle = 290.0
        connection_angle = 180.0
    theta = np.radians(
        np.linspace(start_angle, end_angle, 257)
    )
    for radius in (inner_electrode_radius, outer_electrode_radius):
        arc = gf.Path(
            np.column_stack((radius * np.cos(theta), radius * np.sin(theta)))
        )
        component.add_ref(gf.path.extrude(arc, cross_section=cross_section))

        # Round both ends of the curved electrode.  The pad-side cap also
        # smooths the otherwise sharp 90-degree arc-to-feed junction.
        for angle in (start_angle, end_angle):
            angle_rad = np.radians(angle)
            cap_center = radius * np.array(
                [np.cos(angle_rad), np.sin(angle_rad)]
            )
            cap = component.add_ref(
                gf.components.circle(
                    radius=electrode_width / 2,
                    layer=layer,
                )
            )
            cap.move(tuple(cap_center))

    pad_radius = pad_diameter / 2
    # Put both pads on the left side of the ring and align their centers to the
    # ring's horizontal centerline.  The pair is shifted together toward the
    # ring so the outside pad clears the neighboring vertical bus waveguide.
    # Place the outer pad's ring-facing edge exactly on the inner edge of the
    # outer electrode.  This keeps the connection gap equal to the curved gap
    # while the optical waveguide remains centered between both electrodes.
    outside_pad_center = np.array(
        [
            pad_side
            * (outer_electrode_radius - electrode_width / 2 + pad_radius),
            0.0,
        ]
    )
    inside_pad_center = outside_pad_center + np.array(
        [-pad_side * pad_center_pitch, 0.0]
    )

    # Keep the pad-side feed at its original horizontal position.  For the
    # double-ring arc this is the midpoint of the 180-degree curved electrode.
    connection_angle_rad = np.radians(connection_angle)
    endpoint_direction = np.array(
        [np.cos(connection_angle_rad), np.sin(connection_angle_rad)]
    )
    outer_arc_end = outer_electrode_radius * endpoint_direction
    inner_arc_end = inner_electrode_radius * endpoint_direction

    def pad_edge_toward(pad_center: np.ndarray, target: np.ndarray) -> np.ndarray:
        direction = target - pad_center
        return pad_center + pad_radius * direction / np.linalg.norm(direction)

    connector_points = (
        [
            outer_arc_end,
            pad_edge_toward(outside_pad_center, outer_arc_end),
        ],
        [
            inner_arc_end,
            pad_edge_toward(inside_pad_center, inner_arc_end),
        ],
    )
    for points in connector_points:
        path = gf.Path(np.asarray(points, dtype=float))
        component.add_ref(gf.path.extrude(path, cross_section=cross_section))

    for pad_center in (outside_pad_center, inside_pad_center):
        pad = component.add_ref(
            gf.components.circle(radius=pad_radius, layer=layer)
        )
        pad.move(tuple(pad_center))

    component.info["ring_radius"] = ring_radius
    component.info["electrode_gap"] = electrode_gap
    component.info["electrode_width"] = electrode_width
    component.info["pad_diameter"] = pad_diameter
    component.info["pad_center_pitch"] = pad_center_pitch
    return component


def build_layout() -> gf.Component:
    """Import the source and add electrodes to selected rings."""
    gf.gpdk.PDK.activate()
    source = gf.import_gds(SOURCE_GDS)
    rings = find_rings(source)
    if not rings:
        raise RuntimeError("No Layer-1 rings were detected")

    source.remove_layers([IGNORED_LAYER], recursive=True)
    top = gf.Component("sin_250nm_with_220deg_coplanar_electrodes")
    top.add_ref(source)
    electrode_metal = gf.Component("selected_ring_same_side_pad_electrodes")
    single_flags = isolated_ring_flags(rings)
    double_targets = double_ring_targets(rings, single_flags)
    gap_by_index = electrode_gap_assignments(
        rings,
        single_flags,
        double_targets,
    )
    for index, ((center_x, center_y, ring_radius), is_single) in enumerate(
        zip(rings, single_flags, strict=True)
    ):
        if is_single:
            pad_side = -1.0
        elif index in double_targets:
            pad_side = double_targets[index]
        else:
            # This is the bus-coupled member of a double-ring device.
            continue

        electrode = electrode_metal.add_ref(
            ring_same_side_pad_electrode(
                ring_radius=ring_radius,
                pad_side=pad_side,
                arc_style="single_upper" if is_single else "double_side",
                electrode_gap=gap_by_index[index],
            )
        )
        electrode.move((center_x, center_y))

    # Round the concave feed-to-arc corners.  This closing distance stays below
    # half of the minimum 7.2 um coplanar gap, so neighboring electrodes remain
    # separate while the abrupt right-angle junctions are removed.
    electrode_metal.flatten()
    electrode_metal.over_under(
        layer=ELECTRODE_LAYER,
        distance=JUNCTION_SMOOTHING,
        corner_mode=2,
    )

    # Enforce the same optical keep-out for both single- and double-ring
    # devices and provides a final safeguard around every optical feature.
    optical_keepout = gf.Component("optical_keepout")
    optical_keepout.add_ref(source)
    optical_keepout.flatten()
    optical_keepout.remove_layers(
        [layer for layer in optical_keepout.layers if tuple(layer) != OPTICAL_LAYER],
        recursive=False,
    )
    optical_keepout.offset(
        layer=OPTICAL_LAYER,
        distance=OPTICAL_CLEARANCE,
        corner_mode=2,
    )
    safe_metal = gf.boolean(
        A=electrode_metal,
        B=optical_keepout,
        operation="not",
        layer=ELECTRODE_LAYER,
        layer1=ELECTRODE_LAYER,
        layer2=OPTICAL_LAYER,
    )
    top.add_ref(safe_metal)

    top.info["source_gds"] = str(SOURCE_GDS)
    top.info["ignored_layer"] = list(IGNORED_LAYER)
    top.info["electrode_layer"] = list(ELECTRODE_LAYER)
    top.info["ring_count"] = len(rings)
    top.info["single_ring_count"] = int(sum(single_flags))
    top.info["double_ring_device_count"] = len(double_targets)
    top.info["electrode_device_count"] = int(sum(single_flags)) + len(double_targets)
    top.info["ring_map"] = [list(ring) for ring in rings]
    top.info["electrode_gap"] = ELECTRODE_GAP
    top.info["single_ring_gap_sweep"] = [
        round(SINGLE_RING_GAP_START + index * ELECTRODE_GAP_STEP, 10)
        for index in range(int(sum(single_flags)))
    ]
    top.info["double_ring_gap_sweep_per_row"] = [
        round(DOUBLE_RING_GAP_START + index * ELECTRODE_GAP_STEP, 10)
        for index in range(DOUBLE_RING_DEVICES_PER_ROW)
    ]
    top.info["electrode_gap_by_ring_index"] = [
        [index, gap_by_index[index]] for index in sorted(gap_by_index)
    ]
    top.info["electrode_width"] = ELECTRODE_WIDTH
    top.info["junction_smoothing"] = JUNCTION_SMOOTHING
    top.info["optical_clearance"] = OPTICAL_CLEARANCE
    top.info["pad_size"] = list(PAD_SIZE)
    top.info["pad_center_pitch"] = PAD_CENTER_PITCH
    return top


if __name__ == "__main__":
    layout = build_layout()
    OUTPUT_GDS.parent.mkdir(parents=True, exist_ok=True)
    layout.write_gds(OUTPUT_GDS, with_metadata=False)
    mapped_rings = layout.info["ring_map"]
    with RING_MAP_CSV.open("w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(("ring_index", "center_x_um", "center_y_um", "outer_radius_um"))
        for index, (center_x, center_y, radius) in enumerate(mapped_rings):
            writer.writerow((index, center_x, center_y, radius))
    print(f"Detected rings: {layout.info['ring_count']}")
    print(f"Single-ring electrodes: {layout.info['single_ring_count']}")
    print(f"Double-ring electrodes: {layout.info['double_ring_device_count']}")
    print(f"Total electrode devices: {layout.info['electrode_device_count']}")
    print(f"Single-ring gap sweep: {layout.info['single_ring_gap_sweep']} um")
    print(
        "Double-ring gap sweep per row: "
        f"{layout.info['double_ring_gap_sweep_per_row']} um"
    )
    print(f"Saved: {OUTPUT_GDS}")
    print(f"Ring map: {RING_MAP_CSV}")
