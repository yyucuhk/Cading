import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from functools import partial
import numpy as np

import gdsfactory as gf

from layout.core import GDS
from layout.core.layers import Layers
from layout.core.qrcodes import add_coordinate_qrcode
from layout.core.markers import calipermarker_array
from layout.core.diagnostics import lines_n_gaps
from layout.core.logos import get_cal, get_mj

from layout.resonators.racetracks import ring_point
from layout.devices.bragg import bragg
from layout.core.sweep import sweep
from layout.coupling.testing import connect_grating_couplers
from layout.coupling.grating import grating


RING_GAPS = (0.6, 0.7, 0.8, 0.9, 1.0, 1.1)
GRATING_PERIODS = (0.9, 1.0, 1.1, 1.2)
RING_WIDTH = 3.0
RING_RADIUS = 240.0
RING_BUS_LENGTH = 50.0
GRATING_N_PERIODS = 20
GRATING_TAPER_LENGTH = 20.0
GRATING_TAPER_ANGLE = 40.0
ROUTE_BEND_RADIUS = 50.0
ROUTE_TAPER_LENGTH = 40.0
DEVICE_PITCH = (600.0, 600.0)
CHIP_SHAPE = (5000.0, 5000.0)
INCLUDE_BRAGG_DIAGNOSTICS = False
STRAIGHT_WG_LENGTH = 6000.0
STRAIGHT_WG_WIDTHS = (1.0, 1.5, 2.0, 2.5, 3.0)
STRAIGHT_WG_Y_PITCH = 50.0
STRAIGHT_WG_CENTER_Y = -1550.0

name = __file__.split(os.sep)[-1].split(".")[0]
print(name)

gf.gpdk.PDK.activate()


@gf.cell
def bto_ring_with_gratings(
    gap: float,
    grating_period: float,
    label: str = "",
) -> gf.Component:
    """Ian-style grating-coupled ring for the BTO geometry sweep."""
    return connect_grating_couplers(
        device_spec=partial(
            ring_point,
            width=RING_WIDTH,
            radius=RING_RADIUS,
            bus_length=RING_BUS_LENGTH,
            gap=gap,
        ),
        grating_spec=partial(
            grating,
            period=grating_period,
            n_periods=GRATING_N_PERIODS,
            taper_length=GRATING_TAPER_LENGTH,
            taper_angle=GRATING_TAPER_ANGLE,
            layer_slab=None,
        ),
        ports=("o2", "o1"),
        bend_kwargs={
            "radius": ROUTE_BEND_RADIUS,
            "taper_length": ROUTE_TAPER_LENGTH,
        },
        label=label,
    )


@gf.cell
def straight_waveguide_sweep() -> gf.Component:
    """Five 5 mm waveguides centered on the ring array's x coordinate."""
    c = gf.Component()

    for index, width in enumerate(STRAIGHT_WG_WIDTHS):
        cross_section = gf.cross_section.strip(width=width, layer=Layers.wg)
        waveguide = c.add_ref(
            gf.components.straight(
                length=STRAIGHT_WG_LENGTH,
                cross_section=cross_section,
            )
        )
        y = (
            index - (len(STRAIGHT_WG_WIDTHS) - 1) / 2
        ) * STRAIGHT_WG_Y_PITCH
        waveguide.move((-STRAIGHT_WG_LENGTH / 2, y))

        width_tag = f"{width:g}".replace(".", "p")
        c.add_port(f"w{width_tag}_o1", port=waveguide.ports["o1"])
        c.add_port(f"w{width_tag}_o2", port=waveguide.ports["o2"])

    c.info["waveguide_count"] = len(STRAIGHT_WG_WIDTHS)
    c.info["length"] = STRAIGHT_WG_LENGTH
    c.info["widths"] = list(STRAIGHT_WG_WIDTHS)
    return c

for polarity in ["neg", "pos"]:
    name_ = name
    if polarity == "pos":
        name_ += "_pos"

    # The positive mask is produced by a boolean below, so its source cell
    # needs a distinct name from the final boolean result.
    source_name = f"{name_}_source" if polarity == "pos" else name_
    fin = gf.Component(source_name)

    # Chip outline
    chip_shape = CHIP_SHAPE
    fin.add_ref(
        gf.components.rectangle(size=chip_shape, centered=True, layer=(0, 0))
    )

    # Negative and positive masks use their respective QR fabrication layers.
    label_layer = Layers.label if polarity == "neg" else Layers.wg

    # Alignment marks in case there needs to be an MLA write or something
    fin.add_ref(
        calipermarker_array(
            chip_shape=(chip_shape[0], chip_shape[1]),
            arm_length=100,
            layer1=Layers.wg,
            layer2=None,    # If this layer is filled in caliper verniers will also be drawn
        )
    )

    # Ring gaps run along x and grating periods run along y.
    ring_sweep = fin.add_ref(
        sweep(
            bto_ring_with_gratings,
            pitch=DEVICE_PITCH,
            sweep=(
                {"gap": RING_GAPS},
                {"grating_period": GRATING_PERIODS},
            ),
            label=name_ + "\nBTO_ring",
        )
    )

    # Put one QR code inside each of the four corner rings. The sweep grid and
    # every ring are centered on their placement coordinates.
    corner_x = (len(RING_GAPS) - 1) * DEVICE_PITCH[0] / 2
    corner_y = (len(GRATING_PERIODS) - 1) * DEVICE_PITCH[1] / 2
    for qr_x, qr_y in (
        (-corner_x, -corner_y),
        (-corner_x, +corner_y),
        (+corner_x, -corner_y),
        (+corner_x, +corner_y),
    ):
        add_coordinate_qrcode(
            fin,
            (qr_x, qr_y),
            pattern_name=name_,
            scale=4.5,
            layer=label_layer,
        )

    # Five x-aligned calibration waveguides in the clear area below the rings.
    straight_sweep = fin.add_ref(straight_waveguide_sweep())
    straight_sweep.move((0, STRAIGHT_WG_CENTER_Y))

    # Logos and fiducials
    #fin.add_ref(get_cal(layer=label_layer)).move((1e3/6,0))
    #fin.add_ref(get_mj(width=150, layer=label_layer)).move((-1e3/6,0))

    #fin.add_ref(lines_n_gaps(support=+1, length=50, gap=.6)).rotate(90).move(( 5*1e3/6, 50))
    #fin.add_ref(lines_n_gaps(support=-1, length=50, gap=.6)).rotate(90).move((-5*1e3/6, 50))

   # fin.add_ref(lines_n_gaps(support=+1, length=50)).rotate(90).move(( 5*1e3/6, -55))
   # fin.add_ref(lines_n_gaps(support=-1, length=50)).rotate(90).move((-5*1e3/6, -55))

    if INCLUDE_BRAGG_DIAGNOSTICS:
        for i in [-1, 1]:
            l = -3.5
            for pert in (np.arange(3)*0.1 + 0.2):
                for d_c in [0.25, 0.5]:
                    fin.add_ref(
                        bragg(
                            perturbation=pert,
                            duty_cycle=d_c,
                            length=50,
                            markers=False,
                        )
                    ).move((i*5*1e3/6-25, -l*5))
                    l += 1

    # For this .gds, we also want to export a version compatible with positive resist.
    if polarity == "pos":
        subtract_from = gf.components.rectangle(
            size=chip_shape, centered=True, layer=Layers.wg
        )
        subtract_with = fin

        fin = gf.boolean(
            A=subtract_from,
            B=subtract_with,
            operation="not",
            layer=Layers.wg
        )
        fin.name = name_

    gds_save_options = gf.kdb.SaveLayoutOptions()
    gds_save_options.gds2_max_vertex_count = 2000
    fin.write_gds(
        os.path.join(GDS, f"{name_}.gds"),
        save_options=gds_save_options,
    )
