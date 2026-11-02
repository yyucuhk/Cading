"""Directional tapered-beam stress-test patterns.

Each device contains a circular holding pad and one 100 um beam that tapers
from 0.8 um at the pad to 0.2 um at its outer end.  The first row has free
outer ends; the second row terminates each beam in a clamping slab.

Layer 10 contains the retained mechanical structure.  Layer 12 contains an
independent positive-resist release window: a halo around the complete Layer
10 structure.  No Boolean operation is applied between different devices.
"""

import os
import sys
from math import cos, radians, sin
from pathlib import Path

import gdsfactory as gf

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from layout.core import GDS
from layout.core.layers import Layers
from layout.core.markers import calipermarker_array


# ---------------------------- User parameters -----------------------------
STRUCTURE_LAYER = Layers.wg             # GDS 10/0
OPEN_WINDOW_LAYER = Layers.grating      # GDS 12/0, positive-resist opening
LABEL_LAYER = Layers.wg_doubleetch      # GDS 11/0
CENTER_PROTECTION_LAYER = Layers.ring   # GDS 13/0

# Full 360-degree sweep; 360 itself is omitted because it duplicates 0 degrees.
DIRECTION_STEP_DEG = 10
DIRECTIONS_DEG = tuple(range(0, 360, DIRECTION_STEP_DEG))

HOLDING_SPOT_RADIUS = 30.0
# Additional radius of the Layer-13 protection disk beyond the holding spot.
# Set to 0.0 for an exact duplicate of the support circle.
CENTER_PROTECTION_MARGIN = 10.0
BEAM_LENGTH = 80.0                     # measured from circle edge
BEAM_INNER_WIDTH = 0.8
FREE_BEAM_OUTER_WIDTH = 0.2
CLAMPED_BEAM_OUTER_WIDTHS = (0.2, 0.25, 0.3)
BEAM_CIRCLE_OVERLAP = 0.5

CLAMP_SLAB_SIZE = (30.0, 30.0)          # radial length, transverse width
BEAM_SLAB_OVERLAP = 0.5

OPEN_WINDOW_MARGIN = 0.8
PATTERN_X_PITCH = 300.0
PATTERN_Y_PITCH = 300.0
LABEL_RADIUS = 120.0
LABEL_SIZE = 3.0
# FREE_W... and CLAMPED_W... title settings. Offsets are measured from
# the center of each corresponding star pattern.
TITLE_OFFSET_X = -20.0
TITLE_OFFSET_Y = -150.0
TITLE_SIZE = 8.0
CHIP_MARGIN = 150.0


def _copy_layer(source: gf.Component, source_layer, destination_layer):
    """Return a flat copy of one layer remapped to another layer."""
    copied = gf.Component()
    copied.add_ref(source)
    copied.flatten()
    source_tuple = tuple(source_layer)
    copied.remove_layers(
        [layer for layer in copied.layers if tuple(layer) != source_tuple]
    )
    copied.remap_layers(
        {source_tuple: tuple(destination_layer)},
        recursive=False,
    )
    return copied


@gf.cell
def radial_stress_pattern(
    clamped: bool = False,
    beam_outer_width: float = 0.2,
) -> gf.Component:
    """Create one center spot with beams spanning all requested directions."""
    if BEAM_LENGTH <= 0:
        raise ValueError("BEAM_LENGTH must be positive")
    if not 0 < beam_outer_width <= BEAM_INNER_WIDTH:
        raise ValueError(
            "Widths must satisfy 0 < beam_outer_width <= BEAM_INNER_WIDTH"
        )
    if OPEN_WINDOW_MARGIN <= 0:
        raise ValueError("OPEN_WINDOW_MARGIN must be positive")

    structure = gf.Component()
    release_seed = gf.Component()
    holding_spot = gf.components.circle(
        radius=HOLDING_SPOT_RADIUS,
        layer=STRUCTURE_LAYER,
    )
    structure.add_ref(holding_spot)
    release_seed.add_ref(holding_spot)

    beam_arm = gf.Component()
    beam_start = HOLDING_SPOT_RADIUS - BEAM_CIRCLE_OVERLAP
    beam = beam_arm.add_ref(
        gf.components.taper(
            length=BEAM_LENGTH + BEAM_CIRCLE_OVERLAP,
            width1=BEAM_INNER_WIDTH,
            width2=beam_outer_width,
            layer=STRUCTURE_LAYER,
        )
    )
    beam.move((beam_start, 0))

    clamp_arm = None
    if clamped:
        clamp_arm = gf.Component()
        slab_length, slab_width = CLAMP_SLAB_SIZE
        slab_start = HOLDING_SPOT_RADIUS + BEAM_LENGTH - BEAM_SLAB_OVERLAP
        slab = clamp_arm.add_ref(
            gf.components.rectangle(
                size=(slab_length, slab_width),
                centered=False,
                layer=STRUCTURE_LAYER,
            )
        )
        slab.move((slab_start, -slab_width / 2))

    for angle in DIRECTIONS_DEG:
        structure.add_ref(beam_arm).rotate(angle)
        release_seed.add_ref(beam_arm).rotate(angle)
        if clamp_arm is not None:
            structure.add_ref(clamp_arm).rotate(angle)

    # Construct the positive-resist release halo from the center and beams.
    # For the clamped variant, the outer slabs are deliberately excluded from
    # the offset seed, then subtracted from the halo. Thus Layer 12 stops at
    # the inner clamp faces and never cuts around their outer edges; each beam
    # remains connected to the surrounding substrate at its outer end.
    release_seed_for_window = _copy_layer(
        release_seed,
        STRUCTURE_LAYER,
        OPEN_WINDOW_LAYER,
    )
    full_structure_for_window = _copy_layer(
        structure,
        STRUCTURE_LAYER,
        OPEN_WINDOW_LAYER,
    )
    expanded = gf.Component()
    expanded.add_ref(release_seed_for_window)
    expanded.flatten()
    expanded.offset(
        layer=OPEN_WINDOW_LAYER,
        distance=OPEN_WINDOW_MARGIN,
        corner_mode=2,
    )
    window = gf.boolean(
        A=expanded,
        B=full_structure_for_window,
        operation="A-B",
        layer=OPEN_WINDOW_LAYER,
    )

    pattern = gf.Component()
    # The right-hand doubly clamped test is fabricated directly from the
    # positive-resist openings on Layer 12. Its beam and clamps are the
    # unexposed material between those openings, so no Layer-10 copy is needed.
    # Keep Layer 10 for the left free-end reference pattern only.
    if not clamped:
        pattern.add_ref(structure)
    pattern.add_ref(window)
    pattern.add_ref(
        gf.components.circle(
            radius=HOLDING_SPOT_RADIUS + CENTER_PROTECTION_MARGIN,
            layer=CENTER_PROTECTION_LAYER,
        )
    )
    for angle in DIRECTIONS_DEG:
        theta = radians(angle)
        label = pattern.add_ref(
            gf.components.text(
                text=f"{angle}",
                size=LABEL_SIZE,
                layer=LABEL_LAYER,
            )
        )
        target_x = LABEL_RADIUS * cos(theta)
        target_y = LABEL_RADIUS * sin(theta)
        # Center the complete text bounding box on the radial target. This
        # keeps opposite labels symmetric even when one has more digits.
        label.move(
            (
                target_x - label.center[0],
                target_y - label.center[1],
            )
        )

    pattern.info["directions_deg"] = list(DIRECTIONS_DEG)
    pattern.info["direction_step_deg"] = DIRECTION_STEP_DEG
    pattern.info["clamped"] = bool(clamped)
    pattern.info["structure_layer_drawn"] = not clamped
    pattern.info["beam_length"] = BEAM_LENGTH
    pattern.info["beam_widths"] = [BEAM_INNER_WIDTH, beam_outer_width]
    pattern.info["open_window_margin"] = OPEN_WINDOW_MARGIN
    pattern.info["center_protection_layer"] = tuple(
        CENTER_PROTECTION_LAYER
    )
    pattern.info["center_protection_margin"] = CENTER_PROTECTION_MARGIN
    return pattern


def build_layout() -> gf.Component:
    """Place one free-end and three doubly clamped tests in a 2 x 2 array."""
    gf.gpdk.PDK.activate()
    top = gf.Component(Path(__file__).stem)

    half_x_pitch = PATTERN_X_PITCH / 2
    half_y_pitch = PATTERN_Y_PITCH / 2
    pattern_configs = (
        (False, FREE_BEAM_OUTER_WIDTH, -half_x_pitch, half_y_pitch, "FREE"),
        (True, CLAMPED_BEAM_OUTER_WIDTHS[0], half_x_pitch, half_y_pitch, "CLAMPED"),
        (True, CLAMPED_BEAM_OUTER_WIDTHS[1], -half_x_pitch, -half_y_pitch, "CLAMPED"),
        (True, CLAMPED_BEAM_OUTER_WIDTHS[2], half_x_pitch, -half_y_pitch, "CLAMPED"),
    )

    for clamped, beam_outer_width, x_position, y_position, name in pattern_configs:
        top.add_ref(
            radial_stress_pattern(
                clamped=clamped,
                beam_outer_width=beam_outer_width,
            )
        ).move((x_position, y_position))
        title = top.add_ref(
            gf.components.text(
                text=f"{name}_W{beam_outer_width:g}",
                size=TITLE_SIZE,
                layer=LABEL_LAYER,
            )
        )
        title.move(
            (
                x_position + TITLE_OFFSET_X,
                y_position + TITLE_OFFSET_Y,
            )
        )

    chip_width = 2 * PATTERN_X_PITCH + 2 * CHIP_MARGIN
    chip_height = PATTERN_Y_PITCH + 2 * (LABEL_RADIUS + 40) + 2 * CHIP_MARGIN
    top.add_ref(
        gf.components.rectangle(
            size=(chip_width, chip_height),
            centered=True,
            layer=Layers.outline,
        )
    )
    top.add_ref(
        calipermarker_array(
            chip_shape=(chip_width - 200, chip_height - 100),
            arm_length=30,
            layer1=LABEL_LAYER,
            layer2=None,
        )
    )

    top.info["directions_deg"] = list(DIRECTIONS_DEG)
    top.info["beam_length"] = BEAM_LENGTH
    top.info["beam_inner_width"] = BEAM_INNER_WIDTH
    top.info["free_beam_outer_width"] = FREE_BEAM_OUTER_WIDTH
    top.info["clamped_beam_outer_widths"] = list(
        CLAMPED_BEAM_OUTER_WIDTHS
    )
    top.info["holding_spot_radius"] = HOLDING_SPOT_RADIUS
    top.info["open_window_margin"] = OPEN_WINDOW_MARGIN
    top.info["variants"] = ["free", "clamped", "clamped", "clamped"]
    return top


if __name__ == "__main__":
    component = build_layout()
    save_options = gf.kdb.SaveLayoutOptions()
    save_options.gds2_max_vertex_count = 2000
    component.write_gds(
        os.path.join(GDS, f"{component.name}.gds"),
        save_options=save_options,
        with_metadata=False,
    )
