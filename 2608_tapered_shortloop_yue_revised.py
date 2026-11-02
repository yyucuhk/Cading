"""Dual-resist taper/tether sweep for a shared final etch.

Layer 10 is the negative-resist fine pattern: complete S-bend, tapers, and
tethers. Layer 12 is the positive-resist opening: a 5 um expansion around the
S-bend and optical tapers only; it stops at the taper/tether junction and does
not follow the tether geometry. Layer 12 also contains the original long
tapered slots, support surround slots, and glue holes. Layer 11 holds alignment
markers, taper-end ticks, and labels only. No Boolean merge is applied.
"""

import os
import sys
from math import cos, radians, sin
from functools import partial
from itertools import product
from pathlib import Path

import gdsfactory as gf

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from layout.core import GDS
from layout.core.layers import Layers
from layout.core.markers import calipermarker_array
from tapeout.yue.tapered_local import (
    edge_coupler,
    tapered_fiber_3d,
    tapered_fiber_coupler,
    tapered_fiber_inout,
    tapered_fiber_linear,
)
from tapeout.yue.tether_variants import TETHER_FACTORIES


NEGATIVE_RESIST_LAYER = Layers.wg       # GDS 10/0: negative-resist fine mask
MARKER_LAYER = Layers.wg_doubleetch     # GDS 11/0: markers, ticks, and labels
POSITIVE_RESIST_LAYER = Layers.grating  # GDS 12/0: positive-resist open mask

# Use the editable parametric generator. Set True only when an exact, read-only
# copy of REFERENCE_CAD_GDS is needed; that mode bypasses all parameters below.
RESTORE_REFERENCE_CAD = False
REFERENCE_CAD_GDS = (
    REPOSITORY_ROOT
    / "gds"
    / "2608_tapered_shortloop_yue_revised_nested_2um_3um_5um.gds"
)
POSITIVE_WINDOW_MARGIN = 5.0
# Exact Layer-12 window endpoint measured beyond the physical taper end.
# Keep this at least as large as POSITIVE_WINDOW_MARGIN.
POSITIVE_WINDOW_END_EXTENSION = 10.0
TETHER_VARIANTS = (1, 2, 3, 4)
TETHER_CONNECTION_WIDTHS = (0.25, 0.30, 0.35)
# The reference CAD used one 40 um taper length, giving 3 widths x 4 tether
# shapes = 12 devices. Add 50.0 or 60.0 here when a length sweep is wanted.
TAPER_LENGTHS = (40.0,)
TETHER_CONNECTION_BLEND_LENGTH = 2.0
TETHER_PAD_SCALE = 1.0
# Independent full taper angles for the two long slot layers.
LAYER12_CHANNEL_ANGLE = 3.0
LAYER14_CHANNEL_ANGLE = 1.5
LAYER12_CHANNEL_START_WIDTH = 2.0
LAYER14_CHANNEL_START_WIDTH = 1.0
# Each alignment tip is 0.5 um wide in x, centered at x=0. Layer 14 starts
# after its outward edge, therefore the offset is half the tip width.
# Move the Layer-14 channel start from x = +/-180 um to the centers of the
# two 0.5 x 10 um alignment tips at x = +/-200 um. The offset follows each
# coupler's outward direction, so the same positive value works on both sides.
LAYER14_CHANNEL_START_OFFSET = 25.0
# Layer 12 and the restored Layer 14 tapered slots are both 2000 um long.
POSITIVE_CHANNEL_LENGTH = 2000.0
LAYER14_CHANNEL_LENGTH = 2000.0
S_BEND_WIDTH = 0.8
S_BEND_TAPER_LENGTH = 10.0
# Historical placement values recovered from the aligned reference GDS.
# These three dimensions jointly determine the x position of the taper and
# glue circles. Keep them unchanged when aligning Layer 14 to old Layer 15.
S_BEND_HALF_LENGTH = 200.0
TAPER_REGION_LENGTH = 65.0
TAPER_INPUT_EXTENSION_LENGTH = 20.0
# Two Layer-12 side windows along each existing 0.5 x 10 um alignment tip.
ALIGNMENT_TIP_WINDOW_WIDTH = 2.0
ALIGNMENT_TIP_WINDOW_LENGTH = 8.0
DEVICE_Y_PITCH = 600.0
TETHER_LABEL_SIZE = 2.0

# The support island remains; only this border around it is exposed on L12.
SUPPORT_RECTANGLE_SIZE = (50.0, 40.0)
SUPPORT_SURROUND_SLOT_WIDTH = 5.0
EXTRA_RECTANGLE_GAP = 100.0
# Start the long Layer-12 tapered slot at the end of the short rectangular
# taper-window extension, not at the distant support rectangle.
LAYER12_CHANNEL_START_OFFSET = POSITIVE_WINDOW_END_EXTENSION


def rectangular_surround_slot() -> gf.Component:
    """Return the positive 5 um frame around the support island."""
    length, width = SUPPORT_RECTANGLE_SIZE
    margin = SUPPORT_SURROUND_SLOT_WIDTH
    component = gf.Component()
    for size, center in (
        ((length + 2 * margin, margin), (0, width / 2 + margin / 2)),
        ((length + 2 * margin, margin), (0, -width / 2 - margin / 2)),
        ((margin, width), (length / 2 + margin / 2, 0)),
        ((margin, width), (-length / 2 - margin / 2, 0)),
    ):
        component.add_ref(
            gf.components.rectangle(
                size=size,
                centered=True,
                layer=POSITIVE_RESIST_LAYER,
            )
        ).move(center)
    component.add_ref(
        gf.components.rectangle(
            size=SUPPORT_RECTANGLE_SIZE,
            centered=True,
            layer=Layers.ring,
        )
    )
    return component


def make_tapered_coupler(
    tether_connection_width: float,
    taper_length: float,
):
    """Return a coupler whose narrow taper end matches its tether."""
    if not 0 < tether_connection_width < S_BEND_WIDTH:
        raise ValueError(
            "tether_connection_width must be positive and narrower than "
            "S_BEND_WIDTH"
        )

    taper_spec = partial(
        tapered_fiber_linear,
        width1=S_BEND_WIDTH,
        width2=tether_connection_width,
        # Match the historical CAD: the swept optical taper occupies part of
        # a fixed 65 um region, leaving a straight section after the taper.
        length=TAPER_REGION_LENGTH,
        taper_length=taper_length,
        smooth=True,
    )
    return partial(
        tapered_fiber_coupler,
        tapered_fiber_spec=partial(
            tapered_fiber_3d,
            taper_spec=taper_spec,
            add_doubleetch_strip=False,
            input_extension_length=TAPER_INPUT_EXTENSION_LENGTH,
        ),
        output_bend=False,
        center_taper_length=0,
        center_straight_length=0,
        label=f"{tether_connection_width:g}",
        waveguide_label_size=TETHER_LABEL_SIZE,
        doubleetch_label_size=TETHER_LABEL_SIZE,
        waveguide_label_layer=MARKER_LAYER,
        doubleetch_label_layer=MARKER_LAYER,
        alignment_tick_layer=MARKER_LAYER,
        alignment_tick_window_layer=POSITIVE_RESIST_LAYER,
        alignment_tick_window_width=ALIGNMENT_TIP_WINDOW_WIDTH,
        alignment_tick_window_length=ALIGNMENT_TIP_WINDOW_LENGTH,
        glue_layers=(POSITIVE_RESIST_LAYER, Layers.via),
        oxide_open_layer=POSITIVE_RESIST_LAYER,
        oxide_open_angle=LAYER12_CHANNEL_ANGLE,
        oxide_open_length=POSITIVE_CHANNEL_LENGTH,
        oxide_open_start_width=LAYER12_CHANNEL_START_WIDTH,
        oxide_open_start_at_taper_end=True,
        oxide_open_start_offset=LAYER12_CHANNEL_START_OFFSET,
        oxide_open_duplicate_layer=Layers.via,
        oxide_open_duplicate_angle=LAYER14_CHANNEL_ANGLE,
        oxide_open_duplicate_length=LAYER14_CHANNEL_LENGTH,
        oxide_open_duplicate_start_width=LAYER14_CHANNEL_START_WIDTH,
        oxide_open_duplicate_start_at_taper_end=False,
        oxide_open_duplicate_start_offset=LAYER14_CHANNEL_START_OFFSET,
    )


def add_positive_window(
    device: gf.Component,
    window_source: gf.Component,
    margin: float = POSITIVE_WINDOW_MARGIN,
    end_extension: float = POSITIVE_WINDOW_END_EXTENSION,
) -> gf.Component:
    """Add Layer 12 around the S-bend/tapers, excluding the tethers."""
    if margin <= 0:
        raise ValueError("POSITIVE_WINDOW_MARGIN must be positive")
    if end_extension < margin:
        raise ValueError(
            "POSITIVE_WINDOW_END_EXTENSION must be at least "
            "POSITIVE_WINDOW_MARGIN"
        )

    negative_layer = tuple(NEGATIVE_RESIST_LAYER)
    positive_layer = tuple(POSITIVE_RESIST_LAYER)

    # Flatten only a temporary copy for the offset calculation.
    positive_window = gf.Component()
    positive_window.add_ref(window_source)
    positive_window.flatten()
    positive_window.remove_layers(
        [layer for layer in positive_window.layers if tuple(layer) != negative_layer]
    )
    positive_window.offset(
        layer=negative_layer,
        distance=margin,
        corner_mode=2,
    )
    positive_window.remap_layers(
        {negative_layer: positive_layer},
        recursive=False,
    )

    result = gf.Component()
    result.add_ref(device)
    result.add_ref(positive_window)
    for center, orientation, taper_width in zip(
        window_source.info["taper_end_centers"],
        window_source.info["taper_end_orientations"],
        window_source.info["taper_end_widths"],
        strict=True,
    ):
        endpoint_extension = result.add_ref(
            gf.components.rectangle(
                size=(end_extension, taper_width + 2 * margin),
                centered=True,
                layer=POSITIVE_RESIST_LAYER,
            )
        )
        endpoint_extension.rotate(orientation)
        angle = radians(orientation)
        endpoint_extension.move(
            (
                center[0] + end_extension / 2 * cos(angle),
                center[1] + end_extension / 2 * sin(angle),
            )
        )
    result.info["mask_mode"] = "dual_resist_shared_etch"
    result.info["negative_resist_layer"] = negative_layer
    result.info["positive_resist_layer"] = positive_layer
    result.info["positive_window_margin"] = float(margin)
    result.info["positive_window_end_extension"] = float(end_extension)
    return result


def build_layout() -> gf.Component:
    """Build all tether-width and tether-shape combinations."""
    if RESTORE_REFERENCE_CAD:
        if not REFERENCE_CAD_GDS.is_file():
            raise FileNotFoundError(
                f"Reference GDS not found: {REFERENCE_CAD_GDS}"
            )
        gf.gpdk.PDK.activate()
        return gf.import_gds(
            REFERENCE_CAD_GDS,
            rename_duplicated_cells=True,
        )

    name = Path(__file__).stem
    top_cell_name = f"{name}_L10_L12"

    gf.gpdk.PDK.activate()
    fin = gf.Component(top_cell_name)
    sweep_parameters = tuple(
        product(TETHER_CONNECTION_WIDTHS, TETHER_VARIANTS, TAPER_LENGTHS)
    )
    device_count = len(sweep_parameters)
    chip_height = max(10e3, (device_count - 1) * DEVICE_Y_PITCH + 1000)
    chip_shape = (5e3, chip_height)

    fin.add_ref(
        gf.components.rectangle(size=chip_shape, centered=True, layer=(0, 0))
    )
    fin.add_ref(
        rectangular_surround_slot()
    )

    for marker_chip_shape, arm_length in (
        ((chip_shape[0], chip_shape[1] - 2500), 50),
        ((chip_shape[0] - 1000, chip_shape[1] - 2500), 80),
    ):
        fin.add_ref(
            calipermarker_array(
                chip_shape=marker_chip_shape,
                arm_length=arm_length,
                layer1=MARKER_LAYER,
                layer2=None,
            )
        )

    y_positions = tuple(
        (index - (device_count - 1) / 2) * DEVICE_Y_PITCH
        for index in range(device_count)
    )

    for y_position, (connection_width, tether_variant, taper_length) in zip(
        y_positions,
        sweep_parameters,
        strict=True,
    ):
        tether_spec = partial(
            TETHER_FACTORIES[tether_variant],
            layer=NEGATIVE_RESIST_LAYER,
            pad_scale=TETHER_PAD_SCALE,
            connection_width=connection_width,
            connection_taper_length=TETHER_CONNECTION_BLEND_LENGTH,
        )
        common_device_kwargs = dict(
            chip_length=6000,
            tapered_fiber_angle=0,
            edge_coupler_spec=partial(
                edge_coupler,
                length=300,
                taper_length=100,
            ),
            tapered_fiber_spec=make_tapered_coupler(
                connection_width,
                taper_length,
            ),
            include_edge_couplers=False,
            flatten=False,
            s_bend_width=S_BEND_WIDTH,
            s_bend_taper_length=S_BEND_TAPER_LENGTH,
            s_bend_half_length=S_BEND_HALF_LENGTH,
        )
        device = tapered_fiber_inout(
            **common_device_kwargs,
            extra_rectangle_size=SUPPORT_RECTANGLE_SIZE,
            extra_rectangle_gap=EXTRA_RECTANGLE_GAP,
            extra_rectangle_layer=POSITIVE_RESIST_LAYER,
            extra_rectangle_fill_layer=Layers.ring,
            extra_rectangle_margin=0.0,
            extra_rectangle_surround_slot_width=(
                SUPPORT_SURROUND_SLOT_WIDTH
            ),
            tether_spec=tether_spec,
        )
        window_source = tapered_fiber_inout(**common_device_kwargs)
        fin.add_ref(
            add_positive_window(device, window_source)
        ).move((0, y_position))

    fin.info["device_count"] = device_count
    fin.info["tether_connection_widths"] = list(TETHER_CONNECTION_WIDTHS)
    fin.info["tether_variants"] = list(TETHER_VARIANTS)
    fin.info["taper_lengths"] = list(TAPER_LENGTHS)
    fin.info["s_bend_width"] = S_BEND_WIDTH
    fin.info["s_bend_half_length"] = S_BEND_HALF_LENGTH
    fin.info["alignment_tip_window_size"] = [
        ALIGNMENT_TIP_WINDOW_WIDTH,
        ALIGNMENT_TIP_WINDOW_LENGTH,
    ]
    fin.info["sweep_parameters"] = [
        [float(width), int(variant), float(taper_length)]
        for width, variant, taper_length in sweep_parameters
    ]
    fin.info["mask_mode"] = "dual_resist_shared_etch"
    fin.info["interlayer_boolean_applied"] = False
    fin.info["negative_resist_layer"] = tuple(NEGATIVE_RESIST_LAYER)
    fin.info["positive_resist_layer"] = tuple(POSITIVE_RESIST_LAYER)
    fin.info["positive_window_margin"] = POSITIVE_WINDOW_MARGIN
    fin.info["positive_window_end_extension"] = (
        POSITIVE_WINDOW_END_EXTENSION
    )
    fin.info["marker_layer"] = tuple(MARKER_LAYER)
    fin.info["tether_pad_scale"] = TETHER_PAD_SCALE
    fin.info["support_rectangle_size"] = list(SUPPORT_RECTANGLE_SIZE)
    fin.info["support_surround_slot_width"] = SUPPORT_SURROUND_SLOT_WIDTH
    fin.info["layer13_support_rectangles_only"] = True
    fin.info["layer12_channel_angle"] = LAYER12_CHANNEL_ANGLE
    fin.info["layer14_channel_angle"] = LAYER14_CHANNEL_ANGLE
    fin.info["layer12_channel_start_width"] = (
        LAYER12_CHANNEL_START_WIDTH
    )
    fin.info["layer14_channel_start_width"] = (
        LAYER14_CHANNEL_START_WIDTH
    )
    fin.info["layer14_channel_start_offset"] = (
        LAYER14_CHANNEL_START_OFFSET
    )
    fin.info["layer14_channel_length"] = LAYER14_CHANNEL_LENGTH
    fin.info["layer12_channel_start_offset"] = (
        LAYER12_CHANNEL_START_OFFSET
    )
    fin.info["positive_channel_length"] = POSITIVE_CHANNEL_LENGTH
    fin.info["tether_rectangle_slot"] = False
    return fin


if __name__ == "__main__":
    fin = build_layout()
    gds_save_options = gf.kdb.SaveLayoutOptions()
    gds_save_options.gds2_max_vertex_count = 2000
    fin.write_gds(
        os.path.join(
            GDS,
            f"{Path(__file__).stem}_L10_L12.gds",
        ),
        save_options=gds_save_options,
        with_metadata=False,
    )
