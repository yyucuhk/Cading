"""Generate a standalone geometry sweep of grating-coupled ring devices."""

from __future__ import annotations

import argparse
import sys
from functools import partial
from pathlib import Path

import numpy as np

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from layout.core.layers import Layers
from layout.core.layers import layer as resolve_layer
from layout.core.markers import marker_cross_array
from layout.core import GDS
from layout.core.qrcodes import add_coordinate_qrcode
from layout.resonators.racetracks import ring_point
from layout.coupling.grating import grating, grating_loopback


import gdsfactory as gf

STRAIGHT_WG_LENGTH = 5000.0
STRAIGHT_WG_WIDTHS = (1.0, 1.5, 2.0, 2.5, 3.0)
STRAIGHT_WG_Y_PITCH = 100.0
# Place the aligned straight-waveguide group near the lower chip edge.
STRAIGHT_WG_CENTER_Y = -1500.0
SCRIPT_NAME = Path(__file__).stem
OUTPUT_FILENAME = "2608_bto_Ian_yue_grating_ring_sweep.gds"
CHIP_SIZE = (5000.0, 5000.0)
CHIP_OUTLINE_WIDTH = 1.0
MARKER_EDGE_INSET = 200.0
MARKER_ARM_LENGTH = 150.0
MARKER_ARM_WIDTH = 2.0
GRATING_TOOTH_BIASES = tuple(np.linspace(-0.1, 0.1, 5, endpoint=True))
GRATING_GAP_BIASES = tuple(np.linspace(0, 0, 1, endpoint=True))
RING_GAPS = tuple(np.linspace(0.6, 1.0, 5, endpoint=True))
GRATING_PERIOD = 1.05
GRATING_N_PERIODS = 20
GRATING_TAPER_LENGTH = 20
GRATING_TAPER_ANGLE = 40
WAVEGUIDE_X_SHIFT = -160.0


@gf.cell
def chip_outline(
    size: tuple[float, float] = CHIP_SIZE,
    width: float = CHIP_OUTLINE_WIDTH,
    layer=Layers.outline,
) -> gf.Component:
    """Return a four-rail chip outline centered at the origin."""
    chip_width, chip_height = (float(value) for value in size)
    width = float(width)
    if chip_width <= 0 or chip_height <= 0:
        raise ValueError("chip dimensions must be positive")
    if width <= 0 or 2 * width >= min(chip_width, chip_height):
        raise ValueError(
            "outline width must be positive and less than half the chip size"
        )

    c = gf.Component()
    x0 = -chip_width / 2
    y0 = -chip_height / 2
    horizontal = gf.components.rectangle(
        size=(chip_width, width),
        layer=layer,
    )
    vertical = gf.components.rectangle(
        size=(width, chip_height - 2 * width),
        layer=layer,
    )
    c.add_ref(horizontal).move((x0, y0))
    c.add_ref(horizontal).move((x0, y0 + chip_height - width))
    c.add_ref(vertical).move((x0, y0 + width))
    c.add_ref(vertical).move((x0 + chip_width - width, y0 + width))
    return c

# QR-code locations expressed directly as (i, j, k), where:
#   i = grating tooth-bias index
#   j = ring-gap index
#   k = grating gap-bias index
QR_DEVICE_INDICES = {
    (0, 0, 0),
    (0, len(RING_GAPS) - 1, 0),
    (
        len(GRATING_TOOTH_BIASES) - 1,
        0,
        len(GRATING_GAP_BIASES) - 1,
    ),
    (
        len(GRATING_TOOTH_BIASES) - 1,
        len(RING_GAPS) - 1,
        len(GRATING_GAP_BIASES) - 1,
    ),
}

@gf.cell
def grating_ring_sweep():
    """Return the grating-loopback ring sweep with coordinate QR codes."""
    c = gf.Component()
    pattern = gf.Component()
    qr_code_count = 0

    dx = 600
    dy = 500
    dd = 500
    ddx = 10

    for i, tooth_bias in enumerate(GRATING_TOOTH_BIASES):
        for j, ring_gap in enumerate(RING_GAPS):
            for k, gap_bias in enumerate(GRATING_GAP_BIASES):
                dev = "dev " + str(i) + ", " + str(j) + ", " + str(k)

                x = i*dx + k*dd + j*ddx
                y = j*dy + 300

                place_qr = (i, j, k) in QR_DEVICE_INDICES
                if place_qr:
                    add_coordinate_qrcode(
                        pattern,
                        (x, y),
                        pattern_name="BTO_v1",
                        scale=4.5,
                    )
                    qr_code_count += 1

                dev = pattern.add_ref(
                    grating_loopback(
                        grating_spec=partial(
                            grating,
                            period=GRATING_PERIOD,
                            bias_tooth=tooth_bias,
                            bias_gap=gap_bias,
                            n_periods=GRATING_N_PERIODS,
                            taper_length=GRATING_TAPER_LENGTH,
                            taper_angle=GRATING_TAPER_ANGLE,
                            layer_slab=None,
                        ),
                        device_spec=partial(
                            ring_point,
                            metal=False,
                            radius=180,
                            width = 3,
                            bus_length=60,
                            bus_x_offset=WAVEGUIDE_X_SHIFT,
                            gap=ring_gap,
                        ),
                        reduce_etaloning=False,
                        label=(
                            f"{dev}\n\n"
                            f"grating\nbias_tooth={tooth_bias:.2f}\n"
                            f"bias_gap={gap_bias:.2f}\n\n"
                            f"ring\ngap={ring_gap:.2f}"
                        ),
                    )
                )
                dev.move((x, y))

    # Center all devices, labels, and QR codes at the chip origin. The chip
    # boundary is also centered at (0, 0), so the pattern has an unambiguous
    # origin for placement in a larger mask.
    pattern_reference = c.add_ref(pattern)
    pattern_center_before_centering = tuple(
        float(coordinate) for coordinate in pattern_reference.center
    )
    pattern_reference.move(pattern_reference.center, (0, 0))

    c.info["device_count"] = (
        len(GRATING_TOOTH_BIASES) * len(RING_GAPS) * len(GRATING_GAP_BIASES)
    )
    c.info["pattern_center"] = (0.0, 0.0)
    c.info["pattern_center_before_centering"] = pattern_center_before_centering

    # Draw a thin perimeter, not a filled 25 mm^2 polygon. A filled chip box
    # can be hidden in KLayout yet fractured into misleading shapes by BEAMER.
    c.add_ref(
        chip_outline(
            size=CHIP_SIZE,
            width=CHIP_OUTLINE_WIDTH,
            layer=Layers.outline,
        )
    )

    # Reuse the library's four-corner global alignment-marker array. Marker
    # centers sit 200 um inside the chip edge, leaving 125 um clearance from
    # the ends of the 150 um crosses to the physical boundary.
    c.add_ref(
        marker_cross_array(
            marker_x=CHIP_SIZE[0] / 2 - MARKER_EDGE_INSET,
            marker_y=CHIP_SIZE[1] / 2 - MARKER_EDGE_INSET,
            arm_length=MARKER_ARM_LENGTH,
            arm_width=MARKER_ARM_WIDTH,
            layer=Layers.globalmarker,
        )
    )

    c.info["qr_code_count"] = qr_code_count
    return c

@gf.cell
def straight_waveguide_sweep():
    c = gf.Component()

    for index, width in enumerate(STRAIGHT_WG_WIDTHS):
        cross_section = gf.cross_section.strip(
            width=width,
            layer=Layers.wg,
        )

        waveguide = c.add_ref(
            gf.components.straight(
                length=STRAIGHT_WG_LENGTH,
                cross_section=cross_section,
            )
        )

        # Every waveguide spans the same x coordinates. Only y changes, so the
        # different widths remain aligned for a common cleave while not
        # overlapping one another.
        y = (
            index - (len(STRAIGHT_WG_WIDTHS) - 1) / 2
        ) * STRAIGHT_WG_Y_PITCH

        waveguide.move((-STRAIGHT_WG_LENGTH / 2, y))
        width_tag = f"{width:g}".replace(".", "p")
        c.add_port(f"w{width_tag}_o1", port=waveguide.ports["o1"])
        c.add_port(f"w{width_tag}_o2", port=waveguide.ports["o2"])
    c.info["str_count"] = (
            len(STRAIGHT_WG_WIDTHS)
        )
    c.info["x_center"] = 0.0
    c.info["y_center"] = 0.0

    return c


def build_bto_layout(
    output_gds: str | Path = Path(GDS) / OUTPUT_FILENAME,
) -> Path:
    """Build the standalone sweep and return the written GDS path."""
    output_gds = Path(output_gds).expanduser().resolve()

    gf.gpdk.PDK.activate()
    top = gf.Component("BTO_GRATING_RING_SWEEP")
    sweep = top.add_ref(grating_ring_sweep())
    straight_sweep = top.add_ref(straight_waveguide_sweep())
    straight_sweep.move((0, STRAIGHT_WG_CENTER_Y))
    top.info["device_count"] = sweep.cell.info["device_count"]
    top.info["qr_code_count"] = sweep.cell.info["qr_code_count"]
    top.info["str_count"] = straight_sweep.cell.info["str_count"]
    top.info["pattern_center"] = (0.0, 0.0)
    top.info["straight_waveguide_center"] = (0.0, STRAIGHT_WG_CENTER_Y)

    # BEAMER-safe export: collapse the hierarchy and remove GDS TEXT metadata.
    # QR pixels remain as real waveguide-layer polygons; only their hidden text
    # payload records and the device-description text records are removed.
    # Keep the original polygons separate. ``merge=True`` can create GDS
    # boundary records larger than 32 KiB, which KLayout accepts but some
    # BEAMER versions parse incorrectly as apparently random geometry.
    top.flatten(merge=False)
    top.remove_layers([(resolve_layer(Layers.wg)[0], 1)], recursive=False)

    allowed_layers = {
        resolve_layer(Layers.outline),
        resolve_layer(Layers.globalmarker),
        resolve_layer(Layers.wg),
    }
    unexpected_layers = set(top.layers) - allowed_layers
    if unexpected_layers:
        raise ValueError(
            "unexpected fabrication layers in BTO export: "
            f"{sorted(unexpected_layers)}"
        )

    output_gds.parent.mkdir(parents=True, exist_ok=True)
    return top.write_gds(output_gds)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(GDS) / OUTPUT_FILENAME,
        help="Output GDS path.",
    )
    args = parser.parse_args()

    written = build_bto_layout(args.output)
    print(f"{SCRIPT_NAME}: wrote {written}")


if __name__ == "__main__":
    main()
