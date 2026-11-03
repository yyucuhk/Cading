"""Design 2 - offset two-stage ridge: editable, all-layer GDS generator.

QUICK START
    1. Edit Section 01 (all lengths/widths in um; angles in degrees).
    2. Run the WHOLE file in a fresh Python process after changing settings.
    3. The output path is printed at the end. Existing GDS is overwritten.

SECTION MAP
    00 imports/paths              06 coupler packaging and long channels
    01 USER PARAMETERS           07 Layer 2 / 14 / 16 exposure masks
    02 layers/derived values     08 move alignment tips to Layer 15
    03 taper polygon primitives  09 make_device(): one COMPLETE device
    04 PNC/meander support masks  10 chip placement and global markers
    05 two-layer taper assembly  11 settings checks and GDS export

DRAWING FLOW
    build_layout() -> make_device(L) -> taper/meander/support geometry
                   -> release/protection masks -> labels/tip-layer move
                   -> whole-device placement -> export_gds()

MEANDER_STRAIGHT_LENGTHS has one entry per device, ordered bottom to top.
Each entry is the COMMON length of all THREE straights; other settings are shared.
Meander centerline length = 3 * L + 2 * actual_180_degree_Euler_bend_length.
Complete end assemblies move with the ports; their local geometry is unchanged.
Both ends are constructed together, not copied or clipped from another device.
Use CAD_EDITING_GUIDE.md for the parameter and layer reference.
"""

# %% [00] IMPORTS AND REPOSITORY PATHS
# =============================================================================
from functools import partial
from math import cos, radians, sin
from pathlib import Path
import sys

import gdsfactory as gf
import numpy as np

REPOSITORY_ROOT = Path(__file__).resolve().parent
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))
GDS = REPOSITORY_ROOT / "gds"

from layers import Layers
from markers import calipermarker_array
from phononics import phononic_tether
from tapered_local import (
    edge_coupler,
    tapered_fiber_3d,
    tapered_fiber_coupler,
    tapered_fiber_inout,
)
from tether_variants import TETHER_FACTORIES


# %% [01] USER PARAMETERS - EDIT HERE
# =============================================================================
# Units: um for dimensions, degrees for angles, dimensionless for scales.
# Change values here, then rerun the entire file; decorated-cell defaults and
# the derived values in Section 02 are evaluated when the module is loaded.

# --- 1A. Device count, straight-length sweep, and whole-device placement ---
# One entry = one device. Example: (200.0, 800.0, 1400.0).
MEANDER_STRAIGHT_LENGTHS = (1400.0,)  # EACH of three equal straights; bottom to top
DEVICE_Y_PITCH = 600.0  # origin-to-origin vertical pitch
BOTTOM_DEVICE_Y = -700.0  # first device origin, NOT its lower taper tip
LEFT_BEND_ALIGNMENT_X = -900.0  # common leftmost Layer-11 PNC bend edge

# --- 1B. Optical taper stages (outward from the 3 um PNC endpoint) ---
#                   Stage 1          Stage 2           Stage 3
# Layer 10:       3 -> 1.65         1.65 -> 0.980      0.980 -> 0.200
# Layer 15:       1 -> 0.500        0.500 -> 0.200     NO third stage
# Lengths:          20                 100                 10
# L10 and L15 stage 1 share the first length; stage 2 shares the second.
TAPER_STAGE1_LENGTH = 20.0
TAPER_STAGE2_LENGTH = 100.0
TAPER_STAGE3_LENGTH = 10.0  # Layer 10 only
SLAB_STAGE1_END_WIDTH = 1.65
SLAB_STAGE2_END_WIDTH = 0.980
RIDGE_MEANDER_WIDTH = 1.0  # Layer-15 meander and stage-1 input width
RIDGE_EXTENSION_START_WIDTH = 0.500  # Layer-15 stage-1 END width
RIDGE_CONVERTER_END_WIDTH = 0.2  # Layer-15 terminal width
RIDGE_CONVERTER_Y_SHIFT = 0.840  # final CENTERLINE offset, not edge gap
REVIEW_TETHER_CONNECTION_WIDTH = 0.2  # Layer-10 terminal width AND tether neck
PNC_TAPER_SUPPORT_CUTOFF_WIDTH = 3.0  # Layer-10 input width after the PNC region
TAPER_ANGLE = 0.0  # paired taper placement angle; validated design uses 0
ADIABATIC_TAPER_POINTS = 201  # polygon samples per quintic transition

# --- 1C. Euler meander ---
S_BEND_COMPACT_180_RADIUS = 100.0  # effective radius of each 180-degree bend
S_BEND_EULER_P = 0.2  # Euler fraction: 0 circular, 1 full Euler
# Horizontal half-port separation is derived as L / 2 in make_device().
# Do not set a fixed S_BEND_HALF_LENGTH: it would make the outer leads unequal.
S_BEND_WIDTH = 15.0  # construction/support envelope; NOT optical ridge width
SLAB_EXTENSION_END_WIDTH = 15.0  # width at the start of the 15 -> 3 PNC transition
S_BEND_TAPER_LENGTH = 10.0  # legacy noncompact-route setting; compact path ignores it
# Optical width is RIDGE_MEANDER_WIDTH in 1B.

# --- 1D. Phononic-crystal (PNC) support pattern ---
MEANDER_STRAIGHT_PNC_PERIOD = 5.0  # exact center-to-center pitch on all three straights
REVIEW_SUPPORT_ANCHOR_PERIOD = 15.0  # PNC TAPER pitch only; preserves its 45 um length
REVIEW_SUPPORT_BEND_PERIOD = 3.0  # requested PNC pitch on bends
PHONONIC_SUPPORT_SEGMENTS = (
    0.626,
    0.600,
    0.560,
    0.540,
    0.516,
    0.504,
)  # six square-mass sizes in the PNC library cell
PHONONIC_SUPPORT_GAPS = 0.25  # gap between PNC masses
PNC_NARROW_END_EXTENSION = 0.5  # narrow neck extension at both radial ends
# The PNC taper length is three straight periods; see Section 02.

# --- 1E. Butterfly tethers and support islands ---
TETHER_VARIANT = 2  # 2 = split/forked butterfly in tether_variants.py
TETHER_PAD_SCALE = 1.0  # library pad scale
TETHER_CONNECTION_BLEND_LENGTH = 2.0  # neck blending length
SUPPORT_RECTANGLE_SIZE = (50.0, 40.0)  # (longitudinal length, transverse width)
EXTRA_RECTANGLE_GAP = 100.0  # edge-to-edge gap: taper endpoint to support island
SUPPORT_SURROUND_SLOT_WIDTH = 5.0  # Layer-2 slot around each support island
# Optional extra side-support tips are disabled in make_device(). They
# are NOT the fiber-alignment ticks (which are kept on Layer 15).
SUPPORT_TIP_POSITION = 30.0
SUPPORT_TIP_WIDTH = 0.5
SUPPORT_TIP_LENGTH = 10.0
SUPPORT_TIP_OVERLAP = 0.2

# --- 1F. Layer 2 / 14 release windows and Layer 16 protection ---
POSITIVE_WINDOW_MARGIN = 3.0  # Layer-2 taper-window expansion
POSITIVE_WINDOW_END_EXTENSION = 15.0  # must be >= POSITIVE_WINDOW_MARGIN
LAYER2_LAYER10_COPY_DELAY = 0.5  # inset of L10 copy from the PNC boundary
RIDGE_PROTECTION_MARGIN = 0.2  # Layer-16 offset around ridge; excludes alignment ticks
MLA_CHANNEL_ANGLE = 2.0  # full opening angle on Layer 2
LAYER14_CHANNEL_ANGLE = 1.5  # full opening angle on Layer 14
MLA_CHANNEL_START_WIDTH = 6.0
LAYER14_CHANNEL_START_WIDTH = 1.0
LAYER14_CHANNEL_START_OFFSET = 0.0  # from stage-1 end, NOT final optical endpoint
POSITIVE_CHANNEL_LENGTH = 2000.0  # Layer-2 long-channel length
LAYER14_CHANNEL_LENGTH = 2000.0  # Layer-14 long-channel length
# Dimensions below define REMOVED companion windows, not alignment ticks.
# Actual alignment ticks are fixed at 0.5 x 10 um in tapered_local.py and
# matched by alignment_tip_region(); edit both together to change tick shape.
ALIGNMENT_TIP_WINDOW_WIDTH = 2.0
ALIGNMENT_TIP_WINDOW_LENGTH = 8.0
ALIGNMENT_TIP_WINDOW_INWARD_SHIFT = 0.2

# --- 1G. Labels, chip outline, global overlay marks, and output ---
TETHER_LABEL_SIZE = 2.0
DETAIL_LABEL_SIZE = 3.0
DETAIL_LABEL_TIP_X_OFFSET = 5.0
DETAIL_LABEL_TIP_Y_OFFSET = 10.0
LAYER10_MARKER_ARM_WIDTH = 2.0
LAYER10_MARKER_CENTER_SQUARE_SIZE = 3.0
MLA_MARKER_SIZE_INCREASE = 10.0
MLA_MARKER_ARM_WIDTH = 3.0
MLA_MARKER_CENTER_SQUARE_SIZE = 5.0
CHIP_WIDTH = 5000.0
CHIP_MIN_HEIGHT = 10000.0
CHIP_Y_ALLOWANCE = 1000.0
INNER_MARKER_CHIP_WIDTH_REDUCTION = 1000.0
GLOBAL_MARKER_ARM_LENGTHS = (50.0, 80.0)  # outer row, inner row
GDS_MAX_VERTEX_COUNT = 2000
TOP_CELL_NAME = "Design2_L1400"
OUTPUT_GDS_FILENAME = "2611_tapered_yue_design2.gds"
# TOP_CELL_NAME is a naming label, not a device-count or geometry control.


# %% [02] MASK LAYERS AND DERIVED VALUES - NORMALLY DO NOT EDIT
# =============================================================================
NEGATIVE_RESIST_LAYER = Layers.wg  # 10/0: slab/tethers; alignment ticks move off this layer
FULL_ETCH_LAYER = Layers.globalmarker  # 2/0: continuous positive-resist MLA mask
POSITIVE_RESIST_LAYER = FULL_ETCH_LAYER  # alias for the shared Layer-2 mask
MARKER_LAYER = FULL_ETCH_LAYER  # 2/0: labels and coarse alignment markers
PNC_EBL_LAYER = Layers.wg_doubleetch  # 11/0: solid PNC structures
RIDGE_LAYER = Layers.mc  # 15/0: optical ridge AND isolated alignment ticks
RIDGE_PROTECTION_LAYER = Layers.mc2  # 16/0: optical-ridge protection
MEANDER_CONSTRUCTION_LAYER = (999, 0)  # temporary; never exported
# Other inherited layers: Layers.via = 14/0 (release), Layers.ring =
# 13/0 (island fill), Layers.fiber = 51/0 (fiber guide), 0/0 = chip outline.

DESIGN_ID = "Design 2"
RIDGE_SEPARATION_DEFINITION = "centerline_to_centerline"
RIDGE_TAPER_LENGTH = TAPER_STAGE1_LENGTH
SLAB_STAGE2_LENGTH = TAPER_STAGE2_LENGTH
SLAB_TAPER_LENGTH = TAPER_STAGE3_LENGTH
SLAB_EXTENSION_START_WIDTH = SLAB_STAGE1_END_WIDTH
RIDGE_CONVERTER_LENGTH = TAPER_STAGE2_LENGTH
S_BEND_VERTICAL_HALF_SPAN = 2 * S_BEND_COMPACT_180_RADIUS
PNC_TAPER_LENGTH = 3 * REVIEW_SUPPORT_ANCHOR_PERIOD
POST_PNC_TAPER_LENGTH = RIDGE_TAPER_LENGTH
SLAB_EXTENSION_LENGTH = (
    PNC_TAPER_LENGTH + POST_PNC_TAPER_LENGTH
)  # stage-1 endpoint / alignment landmark, not total taper length
PNC_TAPER_END_WIDTH = PNC_TAPER_SUPPORT_CUTOFF_WIDTH
RIDGE_EXTENSION_END_WIDTH = RIDGE_MEANDER_WIDTH
PHONONIC_SUPPORT_OVERSHOOT = PNC_NARROW_END_EXTENSION
MLA_CHANNEL_START_OFFSET = POSITIVE_WINDOW_END_EXTENSION


# %% [03] TAPER POLYGON PRIMITIVES
# =============================================================================
# Edit profiles here only when changing the taper law, not dimensions.
@gf.cell
def design2_tapered_fiber_adiabatic(
    length: float = 50.0,
    width1: float = 0.8,
    width2: float = 0.2,
    taper_length: float | None = None,
    npoints: int = ADIABATIC_TAPER_POINTS,
    layer=NEGATIVE_RESIST_LAYER,
) -> gf.Component:
    """Quintic adiabatic taper with zero slope/curvature at both ends."""
    if taper_length is None:
        taper_length = length
    if length <= 0 or taper_length <= 0 or taper_length > length:
        raise ValueError("Require 0 < taper_length <= length")
    if width1 <= 0 or width2 <= 0:
        raise ValueError("Taper widths must be positive")
    if npoints < 8:
        raise ValueError("npoints must be at least 8")

    straight_length = length - taper_length
    transition_x = np.linspace(0.0, taper_length, npoints)
    u = transition_x / taper_length
    # Quintic smootherstep: first and second derivatives vanish at u=0,1.
    weight = 6 * u**5 - 15 * u**4 + 10 * u**3
    transition_width = width1 + (width2 - width1) * weight
    x = transition_x + straight_length
    widths = transition_width
    if straight_length > 0:
        x = np.concatenate(([0.0, straight_length], x[1:]))
        widths = np.concatenate(([width1, width1], widths[1:]))

    upper = [(float(xi), float(wi / 2)) for xi, wi in zip(x, widths)]
    lower = [(float(xi), float(-wi / 2)) for xi, wi in zip(x[::-1], widths[::-1])]
    component = gf.Component()
    component.add_polygon(upper + lower, layer=layer)
    component.add_port(
        name="o1",
        center=(0.0, 0.0),
        width=width1,
        orientation=180,
        layer=layer,
    )
    component.add_port(
        name="o2",
        center=(length, 0.0),
        width=width2,
        orientation=0,
        layer=layer,
    )
    component.info["profile"] = "quintic_adiabatic"
    component.info["npoints"] = npoints
    return component


@gf.cell
def design2_offset_ridge_converter(
    length: float = RIDGE_CONVERTER_LENGTH,
    width1: float = RIDGE_EXTENSION_START_WIDTH,
    width2: float = RIDGE_CONVERTER_END_WIDTH,
    y_shift: float = RIDGE_CONVERTER_Y_SHIFT,
    layer=RIDGE_LAYER,
) -> gf.Component:
    """Linear width taper and linear lateral shift, matching the sketch.

    The input cross-section is at x=0, y=0. The output is at x=length,
    y=y_shift. Widths are vertical cross-section widths, not perpendicular
    to the oblique centerline. The far end is intentionally a taper tip.
    """
    if not all(np.isfinite(v) for v in (length, width1, width2, y_shift)):
        raise ValueError("Converter dimensions must be finite")
    if min(length, width1, width2) <= 0:
        raise ValueError("Converter length and widths must be positive")
    component = gf.Component()
    component.add_polygon(
        [
            (0.0, -width1 / 2),
            (length, y_shift - width2 / 2),
            (length, y_shift + width2 / 2),
            (0.0, width1 / 2),
        ],
        layer=layer,
    )
    component.add_port("o1", center=(0.0, 0.0), width=width1, orientation=180, layer=layer)
    component.add_port("o2", center=(length, y_shift), width=width2, orientation=0, layer=layer)
    component.info["length_um"] = float(length)
    component.info["widths_um"] = [float(width1), float(width2)]
    component.info["centerline_offset_um"] = float(y_shift)
    component.info["profile"] = "linear_width_and_centerline"
    return component


# %% [04] PNC SUPPORT AND MEANDER MASK GEOMETRY
# =============================================================================
# Writes Layer 11 PNCs, Layer 2 bands, and optional solid Layer 14 fill.
# The connected optical meander path itself comes from tapered_local.py.
@gf.cell
def design2_phononic_meander_shield(
    path,
    width: float = REVIEW_TETHER_CONNECTION_WIDTH,
    segments: tuple[float, ...] = PHONONIC_SUPPORT_SEGMENTS,
    gaps: float = PHONONIC_SUPPORT_GAPS,
    overshoot: float = PHONONIC_SUPPORT_OVERSHOOT,
    pitch: float = MEANDER_STRAIGHT_PNC_PERIOD,
    bend_pitch: float = REVIEW_SUPPORT_BEND_PERIOD,
    slab_width: float = S_BEND_WIDTH,
    layer=FULL_ETCH_LAYER,
    pnc_layer=PNC_EBL_LAYER,
    layer14_fill_layer=None,
    end_extension: float = PNC_NARROW_END_EXTENSION,
) -> gf.Component:
    """Draw the PNC masks and optional solid Layer-14 meander opening."""
    if pitch <= 0 or bend_pitch <= 0:
        raise ValueError("Straight and bend PNC pitches must be positive")
    if end_extension < 0:
        raise ValueError("PNC narrow-end extension must be non-negative")

    points = np.asarray(path.points, dtype=float)
    vectors = np.diff(points, axis=0)
    lengths = np.linalg.norm(vectors, axis=1)
    if len(lengths) == 0:
        raise ValueError("PNC shield path must contain at least one segment")

    # The constructed meander uses one long segment for every straight and
    # finely sampled short segments for each Euler bend. Grouping consecutive
    # segment types preserves the exact straight/Euler boundaries.
    bend_step_limit = max(1.0, 4.0 * float(np.median(lengths)))
    section_types = lengths <= bend_step_limit
    section_ranges = []
    section_start = 0
    for index in range(1, len(section_types)):
        if section_types[index] != section_types[section_start]:
            section_ranges.append((section_start, index, bool(section_types[section_start])))
            section_start = index
    section_ranges.append((section_start, len(section_types), bool(section_types[section_start])))

    tether = phononic_tether(
        width=width,
        segments=segments,
        gaps=gaps,
        overshoot=overshoot,
    )
    tether_length = float(tether.ports["a2"].center[0] - tether.ports["a1"].center[0])
    cell_extent = max(float(max(segments)), width)
    cell_half_length = cell_extent / 2 + 1e-3

    # Draw continuous Layer-2 side bands and separate solid Layer-11
    # PNC cells. Keep whole cells at straight/Euler boundaries; no PNC-shaped
    # holes are subtracted from the positive mask.
    band_offset = (slab_width + tether_length) / 2
    band_cross_section = gf.CrossSection(
        sections=tuple(
            gf.Section(
                width=tether_length,
                offset=side * band_offset,
                layer=Layers.wg_doubleetch,
                **({"port_names": ("o1", "o2")} if index == 0 else {}),
            )
            for index, side in enumerate((-1, 1))
        )
    )
    outer_bands = gf.path.extrude(path, cross_section=band_cross_section)

    cumulative_lengths = np.concatenate(([0.0], np.cumsum(lengths)))
    full_cells = gf.Component()
    support_positions = []

    def add_pnc_pair(path_position, tangent_override=None):
        """Place one PNC pair at one unique path coordinate."""
        segment_index = min(
            int(
                np.searchsorted(
                    cumulative_lengths,
                    path_position,
                    side="right",
                )
                - 1
            ),
            len(lengths) - 1,
        )
        segment_index = max(segment_index, 0)
        segment_fraction = (path_position - cumulative_lengths[segment_index]) / lengths[
            segment_index
        ]
        center = points[segment_index] + segment_fraction * vectors[segment_index]
        tangent = (
            vectors[segment_index] / lengths[segment_index]
            if tangent_override is None
            else tangent_override / np.linalg.norm(tangent_override)
        )
        normal = np.array([-tangent[1], tangent[0]])
        tangent_angle = float(np.degrees(np.arctan2(tangent[1], tangent[0])))
        for side in (-1, 1):
            anchor = center + side * normal * slab_width / 2
            full_cells.add_ref(tether).rotate(tangent_angle + side * 90).move(tuple(anchor))
        support_positions.append(float(path_position))

    # A single shared PNC is placed exactly at every straight/bend boundary.
    # Keep the first interior PNC one complete local pitch away from it, rather
    # than only one cell length away, to avoid a dense three-PNC cluster.
    for start, stop, is_bend in section_ranges:
        section_start_length = cumulative_lengths[start]
        section_end_length = cumulative_lengths[stop]
        section_length = section_end_length - section_start_length
        requested_pitch = bend_pitch if is_bend else pitch
        start_clearance = requested_pitch if start > 0 else cell_half_length
        end_clearance = requested_pitch if stop < len(lengths) else cell_half_length
        usable_start = section_start_length + start_clearance
        usable_end = section_end_length - end_clearance
        if usable_end < usable_start:
            continue
        # Straights use a fixed pitch, NOT a linspace fitted to the span.
        # Anchor the first straight backwards from its bend. Anchor the other
        # straights forwards from their starting bend. Shared boundary pairs
        # are added below exactly once. Any residual stays at an end clearance.
        if not is_bend:
            if start == 0:
                count = max(0, int(np.floor(
                    (section_length - cell_half_length) / requested_pitch + 1e-9
                )))
                path_positions = section_end_length - requested_pitch * np.arange(1, count + 1)
            else:
                clearance = requested_pitch if stop < len(lengths) else cell_half_length
                count = max(0, int(np.floor(
                    (section_length - clearance) / requested_pitch + 1e-9
                )))
                path_positions = section_start_length + requested_pitch * np.arange(1, count + 1)
            for path_position in np.sort(path_positions):
                add_pnc_pair(float(path_position))
            continue

        # Preserve the existing fitted-to-arc BEND sampler without changes.
        usable_length = usable_end - usable_start
        count = max(
            1,
            int(np.floor(usable_length / requested_pitch)) + 1,
        )
        path_positions = (
            np.array([(usable_start + usable_end) / 2])
            if count == 1
            else np.linspace(
                usable_start,
                usable_end,
                count,
            )
        )
        for path_position in path_positions:
            add_pnc_pair(float(path_position))

    # Exactly one PNC pair at each connection between a straight section and
    # an Euler-bend section. Use the average tangent for a smooth orientation.
    connection_positions = []
    for _, stop, _ in section_ranges[:-1]:
        boundary_segment = stop
        before = vectors[boundary_segment - 1] / lengths[boundary_segment - 1]
        after = vectors[boundary_segment] / lengths[boundary_segment]
        connection_tangent = before + after
        connection_position = float(cumulative_lengths[boundary_segment])
        add_pnc_pair(connection_position, connection_tangent)
        connection_positions.append(connection_position)

    component = gf.Component()
    # Layer 2 is one continuous expanded MLA region. It deliberately contains
    # no PNC-shaped exclusions; the solid PNC is defined separately on L11.
    mla_band = outer_bands.copy()
    mla_band.flatten()
    mla_band.remap_layers(
        {Layers.wg_doubleetch: layer},
        recursive=False,
    )
    component.add_ref(mla_band)

    # Optional Layer 14: fill the complete area enclosed by the OUTER
    # boundaries of the two Layer-2 meander bands. Each Layer-2 band has
    # width tether_length and begins at +/- slab_width/2, so one centered
    # band of width slab_width + 2*tether_length reaches exactly to their
    # outer edges. This opens the whole meander region in Layer 14 without
    # expanding beyond the Layer-2 envelope.
    if layer14_fill_layer is not None:
        layer14_full_width = slab_width + 2 * tether_length
        layer14_cross_section = gf.CrossSection(
            sections=(
                gf.Section(
                    width=layer14_full_width,
                    offset=0.0,
                    layer=layer14_fill_layer,
                    port_names=("o1", "o2"),
                ),
            )
        )
        layer14_full_open = gf.path.extrude(
            path,
            cross_section=layer14_cross_section,
        )
        component.add_ref(layer14_full_open)

    # Layer 11 is the negative-resist EBL mask: draw the complete PNC itself,
    # including the narrow beam extended at both radial ends.
    pnc_cells = full_cells.copy()
    pnc_cells.flatten()
    pnc_cells.remap_layers(
        {Layers.wg_doubleetch: pnc_layer},
        recursive=False,
    )
    component.add_ref(pnc_cells)
    component.info["support_type"] = "phononic_crystal_shield"
    component.info["pnc_ebl_layer"] = tuple(pnc_layer)
    component.info["pnc_mla_layer"] = tuple(layer)
    component.info["pnc_ebl_tone"] = "negative_solid_pnc"
    component.info["pnc_mla_tone"] = "positive_continuous_expanded_area"
    component.info["layer14_full_meander_open"] = layer14_fill_layer is not None
    if layer14_fill_layer is not None:
        component.info["layer14_fill_layer"] = tuple(layer14_fill_layer)
        component.info["layer14_full_open_width"] = slab_width + 2 * tether_length
    component.info["pnc_narrow_end_extension"] = end_extension
    component.info["full_pnc_cells"] = True
    component.info["straight_pnc_pitch_um"] = float(pitch)
    component.info["straight_pnc_sampling"] = "exact_pitch_bend_anchored"
    component.info["support_positions"] = support_positions
    component.info["connection_pnc_positions"] = connection_positions
    return component


@gf.cell
def design2_phononic_taper_shield(
    length: float = SLAB_EXTENSION_LENGTH,
    width1: float = SLAB_EXTENSION_END_WIDTH,
    width2: float = SLAB_EXTENSION_START_WIDTH,
    cutoff_width: float = PNC_TAPER_SUPPORT_CUTOFF_WIDTH,
    pitch: float = REVIEW_SUPPORT_ANCHOR_PERIOD,
    tether_width: float = REVIEW_TETHER_CONNECTION_WIDTH,
    segments: tuple[float, ...] = PHONONIC_SUPPORT_SEGMENTS,
    gaps: float = PHONONIC_SUPPORT_GAPS,
    overshoot: float = PHONONIC_SUPPORT_OVERSHOOT,
    layer=FULL_ETCH_LAYER,
    pnc_layer=PNC_EBL_LAYER,
    layer14_fill_layer=None,
    end_extension: float = PNC_NARROW_END_EXTENSION,
    npoints: int = ADIABATIC_TAPER_POINTS,
) -> gf.Component:
    """Draw PNC masks plus an optional solid Layer-14 taper opening."""
    if not width2 <= cutoff_width < width1:
        raise ValueError("Require width2 <= cutoff_width < width1")
    if length <= 0 or pitch <= 0:
        raise ValueError("Taper length and PNC pitch must be positive")
    if end_extension < 0:
        raise ValueError("PNC narrow-end extension must be non-negative")

    # Solve the same quintic width profile used by design2_tapered_fiber_adiabatic.
    # When the taper ends at the cutoff width, the PNC mask reaches the exact
    # taper endpoint and the Layer-15 taper can start at that same coordinate.
    if np.isclose(width2, cutoff_width):
        cutoff_x = length
    else:
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = (lo + hi) / 2
            weight = 6 * mid**5 - 15 * mid**4 + 10 * mid**3
            local_width = width1 + (width2 - width1) * weight
            if local_width > cutoff_width:
                lo = mid
            else:
                hi = mid
        cutoff_x = length * (lo + hi) / 2

    x = np.linspace(0.0, cutoff_x, npoints)
    u = x / length
    weight = 6 * u**5 - 15 * u**4 + 10 * u**3
    widths = width1 + (width2 - width1) * weight

    tether = phononic_tether(
        width=tether_width,
        segments=segments,
        gaps=gaps,
        overshoot=overshoot,
    )
    tether_length = float(tether.ports["a2"].center[0] - tether.ports["a1"].center[0])
    cell_extent = max(float(max(segments)), tether_width)
    cell_half_length = cell_extent / 2 + 1e-3

    outer_bands = gf.Component()
    inner_edge = widths / 2
    outer_edge = inner_edge + tether_length
    # Extend the PNC mask by half a cell at both longitudinal ends. The last
    # cell remains centered exactly at the 3 um taper endpoint but is no
    # longer cut in half by the end of the surrounding positive-resist band.
    band_end_extension = cell_half_length
    band_x = np.concatenate(([-band_end_extension], x, [cutoff_x + band_end_extension]))
    band_inner_edge = np.concatenate(([inner_edge[0]], inner_edge, [inner_edge[-1]]))
    band_outer_edge = np.concatenate(([outer_edge[0]], outer_edge, [outer_edge[-1]]))
    # Layer 14 should open the COMPLETE 15->3 um taper region, not just
    # duplicate the two Layer-2 side bands.  Fill from the upper OUTER
    # boundary of Layer 2 all the way to the lower OUTER boundary using the
    # exact same longitudinal extension as the Layer-2 mask.
    if layer14_fill_layer is not None:
        layer14_upper = list(zip(band_x, band_outer_edge, strict=True))
        layer14_lower = list(zip(band_x[::-1], -band_outer_edge[::-1], strict=True))
        layer14_full_open = gf.Component()
        layer14_full_open.add_polygon(
            layer14_upper + layer14_lower,
            layer=layer14_fill_layer,
        )
        outer_bands.add_ref(layer14_full_open)

    upper_polygon = list(zip(band_x, band_inner_edge, strict=True)) + list(
        zip(band_x[::-1], band_outer_edge[::-1], strict=True)
    )
    lower_polygon = [(xi, -yi) for xi, yi in upper_polygon]
    outer_bands.add_polygon(upper_polygon, layer=Layers.wg_doubleetch)
    outer_bands.add_polygon(lower_polygon, layer=Layers.wg_doubleetch)

    cuts = gf.Component()
    support_positions = np.arange(pitch, cutoff_x + 1e-9, pitch)
    for support_x in support_positions:
        support_u = support_x / length
        support_weight = 6 * support_u**5 - 15 * support_u**4 + 10 * support_u**3
        support_width = width1 + (width2 - width1) * support_weight
        cuts.add_ref(tether).rotate(90).move((support_x, support_width / 2))
        cuts.add_ref(tether).rotate(-90).move((support_x, -support_width / 2))

    component = gf.Component()
    # Layer 2 remains continuous across the taper and does not reveal PNC
    # detail. Layer 11 independently carries the solid negative-resist cells.
    mla_band = outer_bands.copy()
    mla_band.flatten()
    mla_band.remap_layers(
        {Layers.wg_doubleetch: layer},
        recursive=False,
    )
    component.add_ref(mla_band)
    pnc_cells = cuts.copy()
    pnc_cells.flatten()
    pnc_cells.remap_layers(
        {Layers.wg_doubleetch: pnc_layer},
        recursive=False,
    )
    component.add_ref(pnc_cells)
    component.info["cutoff_width"] = cutoff_width
    component.info["cutoff_x"] = float(cutoff_x)
    component.info["support_positions"] = [float(value) for value in support_positions]
    component.info["pnc_ebl_layer"] = tuple(pnc_layer)
    component.info["pnc_mla_layer"] = tuple(layer)
    component.info["pnc_ebl_tone"] = "negative_solid_pnc"
    component.info["pnc_mla_tone"] = "positive_continuous_expanded_area"
    component.info["pnc_narrow_end_extension"] = end_extension
    return component


# %% [05] TWO-LAYER TAPER ASSEMBLY
# =============================================================================
# Local +x points outward from the meander. Widths are along local y.
# The second device end is placed by rotation; never hand-patch one end.
@gf.cell
def design2_two_layer_tapered_fiber(
    tether_connection_width: float = 0.2,
    ridge_taper_length: float = RIDGE_TAPER_LENGTH,
    add_support_tips: bool = True,
    add_phononic_taper_support: bool = False,
    slab_layer=NEGATIVE_RESIST_LAYER,
) -> gf.Component:
    """Design-2 two-layer coupler after the unchanged PNC region.

    Layer 10: 3 -> 1.65 -> 0.980 -> 0.200 um over 20/100/10 um.
    Layer 15: 1 -> 0.500 -> 0.200 um over 20/100 um; no third stage.
    The preceding 45 um PNC transition and its 1 um ridge stay unchanged.
    """
    if not 0 < ridge_taper_length <= POST_PNC_TAPER_LENGTH:
        raise ValueError("ridge_taper_length must satisfy 0 < length <= " "POST_PNC_TAPER_LENGTH")
    component = gf.Component()

    # The first transition is intentionally not drawn on Layer 10 because its
    # width is greater than 3 um. Keep its longitudinal distance so all PNC,
    # ridge, and downstream taper coordinates remain unchanged.
    # --- 5A. Layer 10: first post-PNC slab stage ---
    post_pnc_slab_taper = component.add_ref(
        design2_tapered_fiber_adiabatic(
            length=POST_PNC_TAPER_LENGTH,
            width1=PNC_TAPER_END_WIDTH,
            width2=SLAB_EXTENSION_START_WIDTH,
            taper_length=POST_PNC_TAPER_LENGTH,
            layer=slab_layer,
        )
    )
    post_pnc_slab_taper.move((PNC_TAPER_LENGTH, 0.0))

    # Keep a standard 1 um Layer-15 optical guide across the complete PNC
    # section. This explicitly continues the 1 um ridge from the meander to
    # the exact point where the narrow ridge taper begins.
    # --- 5B. Layer 15: constant ridge through PNC, then first taper ---
    ridge_pnc_section = component.add_ref(
        design2_tapered_fiber_adiabatic(
            length=PNC_TAPER_LENGTH,
            width1=RIDGE_MEANDER_WIDTH,
            width2=RIDGE_MEANDER_WIDTH,
            taper_length=PNC_TAPER_LENGTH,
            layer=RIDGE_LAYER,
        )
    )

    ridge_extension = component.add_ref(
        design2_tapered_fiber_adiabatic(
            length=ridge_taper_length,
            width1=RIDGE_EXTENSION_END_WIDTH,
            width2=RIDGE_EXTENSION_START_WIDTH,
            taper_length=ridge_taper_length,
            layer=RIDGE_LAYER,
        )
    )
    ridge_extension.connect("o1", ridge_pnc_section.ports["o2"])

    # The 500 -> 200 nm, 100 um displaced converter begins at the first
    # ridge stage endpoint. It terminates without a third stage.
    # --- 5C. Layer 15: shifted second stage; no third ridge stage ---
    ridge_converter = component.add_ref(design2_offset_ridge_converter())
    ridge_converter.connect("o1", ridge_extension.ports["o2"])
    component.info["ridge_converter_end_um"] = [
        float(value) for value in ridge_converter.ports["o2"].center
    ]

    # Layer 10 stage 2 starts at the exact stage-1 port; use a separate
    # component so its 100 um length is independent of the final 10 um taper.
    # --- 5C continued. Layer 10: independent second stage ---
    middle_slab_taper = component.add_ref(
        design2_tapered_fiber_adiabatic(
            length=SLAB_STAGE2_LENGTH,
            width1=SLAB_EXTENSION_START_WIDTH,
            width2=SLAB_STAGE2_END_WIDTH,
            taper_length=SLAB_STAGE2_LENGTH,
            layer=slab_layer,
        )
    )
    middle_slab_taper.connect("o1", post_pnc_slab_taper.ports["o2"])

    # Layer 10 alone continues through the third stage. Its 200 nm port
    # connects to the original butterfly tether; Layer 15 ends 10 um earlier.
    # --- 5D. Layer 10: final slab taper, ending at the tether neck ---
    outer_slab_taper = component.add_ref(
        design2_tapered_fiber_adiabatic(
            length=SLAB_TAPER_LENGTH,
            width1=SLAB_STAGE2_END_WIDTH,
            width2=tether_connection_width,
            taper_length=SLAB_TAPER_LENGTH,
            layer=slab_layer,
        )
    )
    outer_slab_taper.connect("o1", middle_slab_taper.ports["o2"])

    # Continue Layer 2 from the 3 um starting point to the first-stage end.
    # Beyond it the shared window follows both remaining Layer-10 stages.
    # The second, laterally shifted ridge stage stays inside that opening.
    # --- 5E. Layer 2: transition windows and delayed slab copy ---
    release_tether = phononic_tether(
        width=tether_connection_width,
        segments=PHONONIC_SUPPORT_SEGMENTS,
        gaps=PHONONIC_SUPPORT_GAPS,
        overshoot=PHONONIC_SUPPORT_OVERSHOOT,
    )
    release_band_width = float(
        release_tether.ports["a2"].center[0] - release_tether.ports["a1"].center[0]
    )
    transition_x = np.linspace(
        0.0,
        RIDGE_TAPER_LENGTH,
        ADIABATIC_TAPER_POINTS,
    )
    transition_u = transition_x / RIDGE_TAPER_LENGTH
    transition_weight = 6 * transition_u**5 - 15 * transition_u**4 + 10 * transition_u**3
    release_x = PNC_TAPER_LENGTH + transition_x
    start_inner = PNC_TAPER_END_WIDTH / 2
    start_outer = start_inner + release_band_width
    release_target_half_y = SLAB_EXTENSION_START_WIDTH / 2 + (POSITIVE_WINDOW_MARGIN)
    release_inner = start_inner * (1 - transition_weight)
    release_outer = start_outer + (release_target_half_y - start_outer) * transition_weight
    upper_release = list(zip(release_x, release_inner, strict=True)) + list(
        zip(release_x[::-1], release_outer[::-1], strict=True)
    )
    lower_release = [(x_value, -y_value) for x_value, y_value in upper_release]
    component.add_polygon(upper_release, layer=FULL_ETCH_LAYER)
    component.add_polygon(lower_release, layer=FULL_ETCH_LAYER)

    # Copy the first Layer-10 section (3->1.65 um quintic taper) onto Layer 2,
    # but delay its beginning by 0.5 um so it does not reach into the PNC
    # boundary. The copied profile remains identical to Layer 10 thereafter.
    copy_x_local = np.linspace(
        LAYER2_LAYER10_COPY_DELAY,
        POST_PNC_TAPER_LENGTH,
        ADIABATIC_TAPER_POINTS,
    )
    copy_u = copy_x_local / POST_PNC_TAPER_LENGTH
    copy_weight = 6 * copy_u**5 - 15 * copy_u**4 + 10 * copy_u**3
    copy_width = (
        PNC_TAPER_END_WIDTH + (SLAB_EXTENSION_START_WIDTH - PNC_TAPER_END_WIDTH) * copy_weight
    )
    copy_x = PNC_TAPER_LENGTH + copy_x_local
    copy_upper = [
        (float(x_value), float(width_value / 2))
        for x_value, width_value in zip(copy_x, copy_width, strict=True)
    ]
    copy_lower = [
        (float(x_value), float(-width_value / 2))
        for x_value, width_value in zip(
            copy_x[::-1],
            copy_width[::-1],
            strict=True,
        )
    ]
    component.add_polygon(
        copy_upper + copy_lower,
        layer=FULL_ETCH_LAYER,
    )

    # --- 5F. Layers 11/14: PNC taper supports and full release opening ---
    if add_phononic_taper_support:
        component.add_ref(
            design2_phononic_taper_shield(
                length=PNC_TAPER_LENGTH,
                width1=SLAB_EXTENSION_END_WIDTH,
                width2=PNC_TAPER_END_WIDTH,
                cutoff_width=PNC_TAPER_SUPPORT_CUTOFF_WIDTH,
                pitch=REVIEW_SUPPORT_ANCHOR_PERIOD,
                tether_width=tether_connection_width,
                layer=FULL_ETCH_LAYER,
                pnc_layer=PNC_EBL_LAYER,
                layer14_fill_layer=Layers.via,
                end_extension=PNC_NARROW_END_EXTENSION,
            )
        )

    # A lateral support pair stabilizes the longer suspended transition. Each
    # tip overlaps the local slab edge slightly and extends away from both
    # sides of the waveguide.
    # --- 5G. Optional extra side-supports (disabled by current factory) ---
    support_u = SUPPORT_TIP_POSITION / PNC_TAPER_LENGTH
    support_weight = 6 * support_u**5 - 15 * support_u**4 + 10 * support_u**3
    local_slab_width = (
        SLAB_EXTENSION_END_WIDTH + (PNC_TAPER_END_WIDTH - SLAB_EXTENSION_END_WIDTH) * support_weight
    )
    support_size = (
        SUPPORT_TIP_WIDTH,
        SUPPORT_TIP_LENGTH + SUPPORT_TIP_OVERLAP,
    )
    support_center_y = local_slab_width / 2 + (SUPPORT_TIP_LENGTH - SUPPORT_TIP_OVERLAP) / 2
    if add_support_tips:
        for side in (-1, 1):
            support = component.add_ref(
                gf.components.rectangle(
                    size=support_size,
                    centered=True,
                    layer=slab_layer,
                )
            )
            support.move((SUPPORT_TIP_POSITION, side * support_center_y))

    # --- 5H. Ports and metadata used by coupler packaging ---
    component.add_port(
        "o1",
        center=(0.0, 0.0),
        width=SLAB_EXTENSION_END_WIDTH,
        orientation=180,
        layer=slab_layer,
    )
    component.add_port("o2", port=outer_slab_taper.ports["o2"])
    component.info["profile"] = "two_layer_quintic_adiabatic"
    component.info["pnc_taper_length"] = PNC_TAPER_LENGTH
    component.info["post_pnc_taper_length"] = POST_PNC_TAPER_LENGTH
    component.info["slab_extension_length"] = PNC_TAPER_LENGTH + POST_PNC_TAPER_LENGTH
    component.info["slab_taper_length"] = SLAB_TAPER_LENGTH
    component.info["layer10_stage_lengths_um"] = [
        POST_PNC_TAPER_LENGTH,
        SLAB_STAGE2_LENGTH,
        SLAB_TAPER_LENGTH,
    ]
    component.info["layer10_width_sequence_um"] = [
        PNC_TAPER_END_WIDTH,
        SLAB_EXTENSION_START_WIDTH,
        SLAB_STAGE2_END_WIDTH,
        tether_connection_width,
    ]
    component.info["layer15_stage_lengths_um"] = [ridge_taper_length, RIDGE_CONVERTER_LENGTH]
    component.info["layer15_width_sequence_um"] = [
        RIDGE_MEANDER_WIDTH,
        RIDGE_EXTENSION_START_WIDTH,
        RIDGE_CONVERTER_END_WIDTH,
    ]
    component.info["ridge_taper_length"] = ridge_taper_length
    component.info["support_tip_size"] = list(support_size)
    component.info["support_tip_position"] = SUPPORT_TIP_POSITION
    component.info["support_tips_added"] = add_support_tips
    component.info["phononic_taper_support_added"] = add_phononic_taper_support
    component.info["layer10_max_taper_width"] = PNC_TAPER_END_WIDTH
    component.info["layer2_end_distance"] = PNC_TAPER_LENGTH + RIDGE_TAPER_LENGTH
    component.info["layer2_end_y_size"] = 2 * release_target_half_y
    component.info["layer2_layer10_copy_delay"] = LAYER2_LAYER10_COPY_DELAY
    component.info["slab_widths_at_extension"] = [
        SLAB_EXTENSION_END_WIDTH,
        PNC_TAPER_END_WIDTH,
        SLAB_EXTENSION_START_WIDTH,
    ]
    component.info["ridge_widths_at_extension"] = [
        RIDGE_EXTENSION_START_WIDTH,
        RIDGE_EXTENSION_END_WIDTH,
    ]
    return component


# %% [06] COUPLER PACKAGING AND LONG CHANNELS
# =============================================================================
# Uses tapered_local.py for fiber guide, glue spots, ticks and channels.
def make_tapered_coupler(
    tether_connection_width: float,
    ridge_taper_length: float,
    add_support_tips: bool = True,
    add_phononic_taper_support: bool = False,
    alignment_tick_layer=NEGATIVE_RESIST_LAYER,
    slab_layer=NEGATIVE_RESIST_LAYER,
):
    """Return a coupler whose narrow taper end matches its tether."""
    if not 0 < tether_connection_width < S_BEND_WIDTH:
        raise ValueError(
            "tether_connection_width must be positive and narrower than " "S_BEND_WIDTH"
        )

    taper_spec = partial(
        design2_two_layer_tapered_fiber,
        tether_connection_width=tether_connection_width,
        ridge_taper_length=ridge_taper_length,
        add_support_tips=add_support_tips,
        add_phononic_taper_support=add_phononic_taper_support,
        slab_layer=slab_layer,
    )
    return partial(
        tapered_fiber_coupler,
        tapered_fiber_spec=partial(
            tapered_fiber_3d,
            taper_spec=taper_spec,
            add_doubleetch_strip=False,
            input_extension_length=0.0,
        ),
        output_bend=False,
        center_taper_length=0,
        center_straight_length=0,
        label="",
        waveguide_label_size=TETHER_LABEL_SIZE,
        doubleetch_label_size=TETHER_LABEL_SIZE,
        waveguide_label_layer=MARKER_LAYER,
        doubleetch_label_layer=MARKER_LAYER,
        alignment_tick_layer=alignment_tick_layer,
        add_alignment_ticks=True,
        alignment_tick_window_layer=POSITIVE_RESIST_LAYER,
        alignment_tick_window_width=ALIGNMENT_TIP_WINDOW_WIDTH,
        alignment_tick_window_length=ALIGNMENT_TIP_WINDOW_LENGTH,
        alignment_tick_window_inward_shift=(ALIGNMENT_TIP_WINDOW_INWARD_SHIFT),
        alignment_tick_x_offset=SLAB_EXTENSION_LENGTH,
        glue_layers=(POSITIVE_RESIST_LAYER, Layers.via),
        oxide_open_layer=POSITIVE_RESIST_LAYER,
        oxide_open_angle=MLA_CHANNEL_ANGLE,
        oxide_open_length=POSITIVE_CHANNEL_LENGTH,
        oxide_open_start_width=MLA_CHANNEL_START_WIDTH,
        oxide_open_start_at_taper_end=True,
        oxide_open_start_offset=MLA_CHANNEL_START_OFFSET,
        oxide_open_duplicate_layer=Layers.via,
        oxide_open_duplicate_angle=LAYER14_CHANNEL_ANGLE,
        oxide_open_duplicate_length=LAYER14_CHANNEL_LENGTH,
        oxide_open_duplicate_start_width=LAYER14_CHANNEL_START_WIDTH,
        oxide_open_duplicate_start_at_taper_end=False,
        oxide_open_duplicate_start_offset=(SLAB_EXTENSION_LENGTH + LAYER14_CHANNEL_START_OFFSET),
    )


# %% [07] EXPOSURE AND PROTECTION MASKS: LAYERS 2 / 14 / 16
# =============================================================================
# Flatten independent copies before editing. Do NOT recursively remove
# layers from shared/cached child cells: that can erase other devices.
def fill_layer14_taper_seams(
    layer14: gf.Component,
    connection_centers,
    connection_orientations,
) -> gf.Component:
    """Fill enclosed L14 holes only inside the two taper center regions.

    Detect actual hole polygons, rather than drawing a guessed rectangle.
    The allowed region follows each port's local outward direction, from
    its meander connection to its slab tip. Larger support-pad openings
    farther out are explicitly outside this region and remain unchanged.
    Input cells are only read; the returned polygons are independently owned.
    """
    layer_index = gf.get_layer(Layers.via)
    dbu = float(layer14.kcl.dbu)
    region = gf.kdb.Region(layer14.kdb_cell.begin_shapes_rec(layer_index)).merged()
    fills = gf.kdb.Region()
    filled_count = 0
    half_width = SLAB_EXTENSION_END_WIDTH / 2
    tolerance = dbu / 2
    frames = []
    for center, orientation in zip(connection_centers, connection_orientations, strict=True):
        angle = radians(float(orientation))
        outward = np.array([-cos(angle), -sin(angle)])
        normal = np.array([-outward[1], outward[0]])
        frames.append((np.asarray(center, dtype=float), outward, normal))

    for polygon in region.each():
        for hole_index in range(polygon.holes()):
            points = list(polygon.each_point_hole(hole_index))
            coords = np.asarray([(p.x * dbu, p.y * dbu) for p in points])
            for center, outward, normal in frames:
                relative = coords - center
                longitudinal = relative @ outward
                transverse = relative @ normal
                inside_taper = (
                    np.all(longitudinal >= -tolerance)
                    and np.all(longitudinal <= SLAB_EXTENSION_LENGTH + tolerance)
                    and np.all(np.abs(transverse) <= half_width + tolerance)
                )
                if inside_taper:
                    fills.insert(gf.kdb.Polygon(points))
                    filled_count += 1
                    break

    repaired = gf.Component()
    repaired.kdb_cell.shapes(layer_index).insert((region | fills).merged())
    repaired.info["taper_seam_holes_filled"] = filled_count
    repaired.info["taper_seam_fill_area_um2"] = float(fills.area() * dbu**2)
    return repaired


def add_shared_mla_window(
    device: gf.Component,
    window_source: gf.Component,
    margin: float = POSITIVE_WINDOW_MARGIN,
    end_extension: float = POSITIVE_WINDOW_END_EXTENSION,
) -> gf.Component:
    """Add the shared Layer-2 MLA release window around the optical taper."""
    if margin <= 0:
        raise ValueError("POSITIVE_WINDOW_MARGIN must be positive")
    if end_extension < margin:
        raise ValueError("POSITIVE_WINDOW_END_EXTENSION must be at least " "POSITIVE_WINDOW_MARGIN")

    negative_layer = tuple(NEGATIVE_RESIST_LAYER)
    positive_layer = tuple(POSITIVE_RESIST_LAYER)

    # Flatten only a temporary copy for the offset calculation.
    # --- 7A. Layer 2: expanded optical-taper exposure window ---
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

    # Start the expanded Layer-2 taper window at the END OF STAGE 1.
    # For Design 2 this is not the final ridge endpoint. There is no Layer 12.
    release_start_distance = PNC_TAPER_LENGTH + RIDGE_TAPER_LENGTH
    layer12_start_regions = gf.Component()
    clip_length = 6000.0
    clip_width = 6000.0
    for center, orientation in zip(
        window_source.info["s_bend_connection_centers"],
        window_source.info["s_bend_connection_orientations"],
        strict=True,
    ):
        angle = radians(orientation)
        outward = np.array([-cos(angle), -sin(angle)])
        release_start = np.asarray(center, dtype=float) + (release_start_distance * outward)
        clip_center = release_start + clip_length / 2 * outward
        clip = layer12_start_regions.add_ref(
            gf.components.rectangle(
                size=(clip_length, clip_width),
                centered=True,
                layer=positive_layer,
            )
        )
        clip.rotate(orientation + 180).move(tuple(clip_center))

    positive_window = gf.boolean(
        A=positive_window,
        B=layer12_start_regions,
        operation="and",
        layer=positive_layer,
        layer1=positive_layer,
        layer2=positive_layer,
    )

    # --- 7B. Layer 16: optical-ridge offset (before moving alignment tips) ---
    ridge_protection = gf.Component()
    ridge_protection.add_ref(device)
    ridge_protection.flatten()
    ridge_layer = tuple(RIDGE_LAYER)
    ridge_protection.remove_layers(
        [layer for layer in ridge_protection.layers if tuple(layer) != ridge_layer]
    )
    ridge_protection.offset(
        layer=ridge_layer,
        distance=RIDGE_PROTECTION_MARGIN,
        corner_mode=2,
    )
    ridge_protection.remap_layers(
        {ridge_layer: tuple(RIDGE_PROTECTION_LAYER)},
        recursive=False,
    )

    result = gf.Component()
    result.add_ref(device)
    result.add_ref(positive_window)
    result.add_ref(ridge_protection)

    # Layer 14 is the final MLA release exposure. The meander itself already
    # receives a SOLID Layer-14 opening from design2_phononic_meander_shield(), bounded
    # exactly by the OUTER edges of Layer 2. Here we keep the existing Layer-2
    # copy through the adjacent PNC/taper regions and clip it at the current
    # Layer-14 channel starts.
    # --- 7C. Layer 14: central opening and local start coordinates ---
    layer14_start_points = []
    for center, orientation in zip(
        window_source.info["s_bend_connection_centers"],
        window_source.info["s_bend_connection_orientations"],
        strict=True,
    ):
        angle = radians(orientation)
        outward = np.array([-cos(angle), -sin(angle)])
        layer14_start_points.append(
            np.asarray(center, dtype=float)
            + (SLAB_EXTENSION_LENGTH + LAYER14_CHANNEL_START_OFFSET) * outward
        )

    layer14_starts = np.asarray(layer14_start_points, dtype=float)
    layer14_xmin = float(layer14_starts[:, 0].min())
    layer14_xmax = float(layer14_starts[:, 0].max())

    # First copy the central PNC/taper bands. The COMPLETE cleaned
    # Layer-2 window is merged into Layer 14 later in this function.
    layer14_pattern_source = device.copy()
    layer14_pattern_source.flatten()
    layer14_pattern_source.remove_layers(
        [layer for layer in layer14_pattern_source.layers if tuple(layer) != positive_layer]
    )
    layer14_pattern_source.remap_layers(
        {positive_layer: tuple(Layers.via)},
        recursive=False,
    )

    layer14_clip = gf.Component()
    clip_height = 6000.0
    clip_center_y = float(layer14_starts[:, 1].mean())
    clip_ref = layer14_clip.add_ref(
        gf.components.rectangle(
            size=(layer14_xmax - layer14_xmin, clip_height),
            centered=True,
            layer=Layers.via,
        )
    )
    clip_ref.move(
        (
            (layer14_xmin + layer14_xmax) / 2,
            clip_center_y,
        )
    )
    layer14_central_pattern = gf.boolean(
        A=layer14_pattern_source,
        B=layer14_clip,
        operation="and",
        layer=Layers.via,
        layer1=Layers.via,
        layer2=Layers.via,
    )

    # Remove the companion-window blocks at BOTH taper tips.  Work in
    # device-local coordinates so every straight-length variant receives
    # exactly the same two-ended cleanup. Alignment tips are retained
    # here and moved from Layer 10 to Layer 15 at the end of make_device().
    # --- 7D. Both ends: remove only companion-window blocks ---
    layer14_tip_exclusions = gf.Component()
    tick_half_separation = 15.0 / 2 - ALIGNMENT_TIP_WINDOW_INWARD_SHIFT
    tick_window_x_offset = 0.25 + ALIGNMENT_TIP_WINDOW_WIDTH / 2
    tip_orientations = window_source.info["s_bend_connection_orientations"]
    for tip_start, tip_orientation in zip(layer14_start_points, tip_orientations, strict=True):
        tip_angle = radians(tip_orientation)
        outward = np.array([-cos(tip_angle), -sin(tip_angle)])
        normal = np.array([-outward[1], outward[0]])
        for transverse_side in (-1.0, 1.0):
            for longitudinal_side in (-1.0, 1.0):
                window_center = (
                    np.asarray(tip_start, dtype=float)
                    + longitudinal_side * tick_window_x_offset * outward
                    + transverse_side * tick_half_separation * normal
                )
                exclusion = layer14_tip_exclusions.add_ref(
                    gf.components.rectangle(
                        size=(ALIGNMENT_TIP_WINDOW_WIDTH, ALIGNMENT_TIP_WINDOW_LENGTH),
                        centered=True,
                        layer=Layers.via,
                    )
                )
                exclusion.rotate(tip_orientation + 180)
                exclusion.move(tuple(window_center))

    layer14_central_pattern = gf.boolean(
        A=layer14_central_pattern,
        B=layer14_tip_exclusions,
        operation="not",
        layer=Layers.via,
        layer1=Layers.via,
        layer2=Layers.via,
    )
    result.add_ref(layer14_central_pattern)

    # --- 7E. Layer 2: extend the outer window to meet its long channel ---
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
    # Final Layer-2 / Layer-14 strategy:
    # 1) Remove the companion-window blocks at both taper tips
    #    directly from Layer 2.
    # 2) Copy ALL of that cleaned Layer 2 onto Layer 14.
    # 3) Merge it with the existing Layer 14 geometry.
    #
    # The same local exclusion is also applied to the pre-existing Layer 14
    # before the union, so the old tip blocks cannot survive from another
    # Layer-14 construction path.
    # --- 7F. Clean Layer 2 and restore its continuous central openings ---
    layer2_tip_exclusions = layer14_tip_exclusions.copy()
    layer2_tip_exclusions.flatten()
    layer2_tip_exclusions.remap_layers(
        {tuple(Layers.via): positive_layer},
        recursive=False,
    )

    existing_layer2 = result.copy()
    existing_layer2.flatten()
    existing_layer2.remove_layers(
        [layer for layer in existing_layer2.layers if tuple(layer) != positive_layer]
    )
    cleaned_layer2 = gf.boolean(
        A=existing_layer2,
        B=layer2_tip_exclusions,
        operation="not",
        layer=positive_layer,
        layer1=positive_layer,
        layer2=positive_layer,
    )

    # Rebuild the central opening after companion-window removal. Keep the
    # central release opening as one smooth continuous bridge through this
    # local transition.  The bridge width is the intended taper opening
    # (stage-1 slab width + 2*MLA margin), so it does not recreate the
    # removed companion blocks.
    smooth_tip_bridge = gf.Component()
    bridge_length = 6.0
    bridge_width = SLAB_EXTENSION_START_WIDTH + 2 * margin
    for tip_start, tip_orientation in zip(layer14_start_points, tip_orientations, strict=True):
        bridge = smooth_tip_bridge.add_ref(
            gf.components.rectangle(
                size=(bridge_length, bridge_width),
                centered=True,
                layer=positive_layer,
            )
        )
        bridge.rotate(tip_orientation)
        bridge.move(tuple(np.asarray(tip_start, dtype=float)))

    cleaned_layer2 = gf.boolean(
        A=cleaned_layer2,
        B=smooth_tip_bridge,
        operation="or",
        layer=positive_layer,
        layer1=positive_layer,
        layer2=positive_layer,
    )

    # Replace Layer 2 with the cleaned, continuous version before creating
    # the Layer-14 copy. All other mask layers remain unchanged.
    result.flatten()
    result.remove_layers([positive_layer], recursive=False)
    result.add_ref(cleaned_layer2)

    # --- 7G. Layer 14: merge the ENTIRE cleaned Layer-2 mask ---
    layer2_full_copy = cleaned_layer2.copy()
    layer2_full_copy.flatten()
    layer2_full_copy.remap_layers(
        {positive_layer: tuple(Layers.via)},
        recursive=False,
    )

    existing_layer14 = result.copy()
    existing_layer14.flatten()
    existing_layer14.remove_layers(
        [layer for layer in existing_layer14.layers if tuple(layer) != tuple(Layers.via)]
    )
    existing_layer14 = gf.boolean(
        A=existing_layer14,
        B=layer14_tip_exclusions,
        operation="not",
        layer=Layers.via,
        layer1=Layers.via,
        layer2=Layers.via,
    )

    merged_layer14 = gf.boolean(
        A=existing_layer14,
        B=layer2_full_copy,
        operation="or",
        layer=Layers.via,
        layer1=Layers.via,
        layer2=Layers.via,
    )

    # Fill the tiny enclosed taper-junction holes on Layer 14 only.
    # Do this after the union so mask fragments cannot recreate the seams.
    # --- 7H. Layer 14 only: fill bounded taper-junction seams ---
    merged_layer14 = fill_layer14_taper_seams(
        merged_layer14,
        window_source.info["s_bend_connection_centers"],
        window_source.info["s_bend_connection_orientations"],
    )
    result.info["layer14_taper_holes_filled"] = int(merged_layer14.info["taper_seam_holes_filled"])
    result.info["layer14_taper_fill_area_um2"] = float(
        merged_layer14.info["taper_seam_fill_area_um2"]
    )
    result.info["tip_window_cleanup"] = "both_ends_with_continuous_central_opening"

    # Replace only Layer 14 with independently owned merged polygons.
    result.flatten()
    result.remove_layers([Layers.via], recursive=False)
    result.add_ref(merged_layer14)

    # The broad meander route is needed only to calculate ports and PNC
    # placement. Do not export that temporary geometry on any mask layer.
    result.flatten()
    result.remove_layers([MEANDER_CONSTRUCTION_LAYER], recursive=False)
    # --- 7I. Record final mask metadata; construction layer already removed ---
    result.info["mask_mode"] = "dual_resist_shared_etch"
    result.info["negative_resist_layer"] = negative_layer
    result.info["positive_resist_layer"] = positive_layer
    result.info["positive_window_margin"] = float(margin)
    result.info["positive_window_end_extension"] = float(end_extension)
    result.info["layer2_transition_endpoint_distance"] = release_start_distance
    result.info["layer2_release_target_y_size"] = SLAB_EXTENSION_START_WIDTH + 2 * margin
    result.info["layer2_adiabatic_taper_length"] = RIDGE_TAPER_LENGTH
    result.info["ridge_protection_layer"] = tuple(RIDGE_PROTECTION_LAYER)
    result.info["ridge_protection_margin"] = RIDGE_PROTECTION_MARGIN
    result.info["layer14_central_pattern_x_bounds"] = [
        layer14_xmin,
        layer14_xmax,
    ]
    result.flatten()
    return result


# %% [08] ALIGNMENT TICKS: MOVE LAYER 10 TO LAYER 15
# =============================================================================
# Two 0.5 x 10 um ticks per taper end; shapes/positions are unchanged.
# This late move prevents them from changing the Layer-16 protection mask.
def alignment_tip_region(tip_centers, orientations, dbu: float):
    """Exact 0.5 x 10 um alignment ticks: two at each taper end.

    Matches tapered_local.tapered_fiber_coupler's local tick definition.
    The returned region is in integer database units, with no size bias.
    """
    region = gf.kdb.Region()
    for tip, orientation in zip(tip_centers, orientations, strict=True):
        angle = radians(float(orientation) + 180.0)
        tangent = np.array([cos(angle), sin(angle)])
        normal = np.array([-sin(angle), cos(angle)])
        for side in (-1.0, 1.0):
            center = np.asarray(tip, dtype=float) + side * 7.5 * normal
            corners = [
                center + dx * tangent + dy * normal
                for dx, dy in ((-0.25, -5.0), (0.25, -5.0), (0.25, 5.0), (-0.25, 5.0))
            ]
            region.insert(
                gf.kdb.Polygon(
                    [
                        gf.kdb.Point(int(round(point[0] / dbu)), int(round(point[1] / dbu)))
                        for point in corners
                    ]
                )
            )
    return region.merged()


def move_alignment_tips_to_layer15(component, tip_centers, orientations):
    """Move only the four existing ticks after all exposure masks are built.

    Late reassignment deliberately leaves the L2/L14 windows and L16
    protection mask unchanged. The input must be independently flattened;
    no cached child cells or other Layer-10/Layer-15 polygons are edited.
    """
    component.flatten()
    expected = alignment_tip_region(tip_centers, orientations, float(component.kcl.dbu))
    expected_polygons = list(expected.each())
    source_shapes = component.kdb_cell.shapes(gf.get_layer(NEGATIVE_RESIST_LAYER))
    destination = component.kdb_cell.shapes(gf.get_layer(RIDGE_LAYER))
    selected = []
    selected_region = gf.kdb.Region()
    for shape in source_shapes.each():
        if not (shape.is_polygon() or shape.is_box()):
            continue
        polygon = gf.kdb.Polygon(shape.box) if shape.is_box() else shape.polygon
        for target in expected_polygons:
            if (
                polygon.bbox() == target.bbox()
                and (gf.kdb.Region(polygon) ^ gf.kdb.Region(target)).is_empty()
            ):
                selected.append(shape)
                selected_region.insert(polygon)
                break
    if len(selected) != 2 * len(tip_centers) or not (selected_region ^ expected).is_empty():
        raise ValueError("Expected exactly two complete Layer-10 ticks per taper end")
    # Copy each original shape exactly, then erase only that source shape.
    for shape in selected:
        destination.insert(shape)
        source_shapes.erase(shape)
    component.info["alignment_tip_layer"] = list(RIDGE_LAYER)
    component.info["alignment_tip_count"] = len(selected)


# %% [09] ONE COMPLETE DEVICE: make_device(meander_straight_length)
# =============================================================================
# This is the sole all-layer device factory. Every sweep entry calls it.
def make_device(meander_straight_length: float) -> gf.Component:
    """Build one complete, independent all-layer device in local coordinates.

    All THREE straight segments have length meander_straight_length.
    The two port X coordinates are -L/2 and +L/2. Both taper/coupler assemblies,
    butterfly tethers, support pads, alignment tips, exposure windows, and
    all device layers are constructed here through the same helper calls.
    No geometry is copied from another sweep device or cut at a global X.
    """
    meander_straight_length = float(meander_straight_length)
    if not np.isfinite(meander_straight_length) or meander_straight_length <= 0:
        raise ValueError("meander_straight_length must be finite and positive")
    if TAPER_ANGLE != 0.0:
        raise ValueError("Equal-straight layout requires TAPER_ANGLE = 0; rotate a complete device instead")
    # The native router sets each outer lead to (port_dx + middle_L) / 2.
    # Setting port_dx = L makes both outer leads exactly L as well.
    meander_half_length = meander_straight_length / 2.0
    gf.gpdk.PDK.activate()
    connection_width = REVIEW_TETHER_CONNECTION_WIDTH
    taper_angle = TAPER_ANGLE
    support_period = MEANDER_STRAIGHT_PNC_PERIOD
    # --- 9A. Configure butterfly tether and common meander geometry ---
    tether_spec = partial(
        TETHER_FACTORIES[TETHER_VARIANT],
        layer=NEGATIVE_RESIST_LAYER,
        pad_scale=TETHER_PAD_SCALE,
        connection_width=connection_width,
        connection_taper_length=TETHER_CONNECTION_BLEND_LENGTH,
    )
    common_device_kwargs = dict(
        chip_length=6000,
        tapered_fiber_angle=taper_angle,
        edge_coupler_spec=partial(
            edge_coupler,
            length=300,
            taper_length=100,
        ),
        include_edge_couplers=False,
        flatten=False,
        s_bend_width=S_BEND_WIDTH,
        s_bend_taper_length=S_BEND_TAPER_LENGTH,
        s_bend_half_length=meander_half_length,
        s_bend_compact_180_radius=S_BEND_COMPACT_180_RADIUS,
        s_bend_compact_euler_p=S_BEND_EULER_P,
        s_bend_compact_straight_length=meander_straight_length,
        s_bend_vertical_half_span=S_BEND_VERTICAL_HALF_SPAN,
        # Use the wide meander only as construction geometry. It is
        # removed before export and is not transferred to Layer 2.
        s_bend_layer=MEANDER_CONSTRUCTION_LAYER,
        ridge_s_bend_layer=RIDGE_LAYER,
        ridge_s_bend_width=RIDGE_EXTENSION_END_WIDTH,
    )
    # --- 9B. Draw BOTH ends, meander, PNCs, tethers, and support pads ---
    device = tapered_fiber_inout(
        **common_device_kwargs,
        s_bend_support_shield_spec=partial(
            design2_phononic_meander_shield,
            width=connection_width,
            pitch=support_period,
            bend_pitch=REVIEW_SUPPORT_BEND_PERIOD,
            slab_width=S_BEND_WIDTH,
            layer=FULL_ETCH_LAYER,
            pnc_layer=PNC_EBL_LAYER,
            layer14_fill_layer=Layers.via,
            end_extension=PNC_NARROW_END_EXTENSION,
        ),
        tapered_fiber_spec=make_tapered_coupler(
            connection_width,
            RIDGE_TAPER_LENGTH,
            add_support_tips=False,
            add_phononic_taper_support=True,
            slab_layer=NEGATIVE_RESIST_LAYER,
        ),
        extra_rectangle_size=SUPPORT_RECTANGLE_SIZE,
        extra_rectangle_gap=EXTRA_RECTANGLE_GAP,
        extra_rectangle_layer=POSITIVE_RESIST_LAYER,
        extra_rectangle_fill_layer=Layers.ring,
        extra_rectangle_margin=0.0,
        extra_rectangle_surround_slot_width=(SUPPORT_SURROUND_SLOT_WIDTH),
        tether_spec=tether_spec,
    )
    # --- 9C. Separate source for exposure-window construction ---
    window_source = tapered_fiber_inout(
        **common_device_kwargs,
        tapered_fiber_spec=make_tapered_coupler(
            connection_width,
            RIDGE_TAPER_LENGTH,
            add_support_tips=False,
            alignment_tick_layer=MARKER_LAYER,
            slab_layer=NEGATIVE_RESIST_LAYER,
        ),
    )

    # Snapshot placement information BEFORE assembling the independent masks.
    left_pnc_x = float(device.bbox(layer=gf.get_layer(PNC_EBL_LAYER)).left)
    connection_centers = [
        list(map(float, p)) for p in window_source.info["s_bend_connection_centers"]
    ]
    connection_orientations = list(map(float, window_source.info["s_bend_connection_orientations"]))
    # --- 9D. Complete Layer 2 / 14 / 16 masks on owned polygons ---
    component = add_shared_mla_window(device, window_source)
    component.flatten()

    # --- 9E. Length label and alignment-tip coordinates ---
    euler_length = float(
        gf.path.euler(
            radius=S_BEND_COMPACT_180_RADIUS, angle=180, p=S_BEND_EULER_P, use_eff=True
        ).length()
    )
    # Centerline between the two meander ports only; excludes end tapers/tethers.
    # Euler length is measured from the actual path, not approximated by pi*R.
    total_length = 3 * meander_straight_length + 2 * euler_length
    tips = []
    for center, orientation in zip(connection_centers, connection_orientations, strict=True):
        angle = radians(orientation)
        tips.append(
            np.asarray(center) + SLAB_EXTENSION_LENGTH * np.array([-cos(angle), -sin(angle)])
        )
    left_tip = min(tips, key=lambda point: float(point[0]))
    label = component.add_ref(
        gf.components.text(
            text=f"DESIGN 2_L{meander_straight_length:g}_ML{total_length:.1f}",
            size=DETAIL_LABEL_SIZE,
            layer=MARKER_LAYER,
        )
    )
    label.move(
        (
            float(left_tip[0]) + DETAIL_LABEL_TIP_X_OFFSET,
            float(left_tip[1]) + DETAIL_LABEL_TIP_Y_OFFSET,
        )
    )
    label_box = label.bbox()
    label_bounds = [
        float(label_box.left),
        float(label_box.bottom),
        float(label_box.right),
        float(label_box.top),
    ]
    component.flatten()
    # --- 9F. Move alignment ticks AFTER all exposure/protection masks ---
    move_alignment_tips_to_layer15(component, tips, connection_orientations)
    # --- 9G. Return complete independent geometry plus metadata ---
    component.info["builder"] = "make_device"
    component.info["design"] = DESIGN_ID
    component.info["layer10_final_taper_length_um"] = SLAB_TAPER_LENGTH
    component.info["ridge_stage1_length_um"] = RIDGE_TAPER_LENGTH
    component.info["ridge_stage2_length_um"] = RIDGE_CONVERTER_LENGTH
    component.info["ridge_width_sequence_um"] = [
        RIDGE_MEANDER_WIDTH,
        RIDGE_EXTENSION_START_WIDTH,
        RIDGE_CONVERTER_END_WIDTH,
    ]
    component.info["layer10_stage_lengths_um"] = [
        POST_PNC_TAPER_LENGTH,
        SLAB_STAGE2_LENGTH,
        SLAB_TAPER_LENGTH,
    ]
    component.info["layer10_width_sequence_um"] = [
        PNC_TAPER_END_WIDTH,
        SLAB_EXTENSION_START_WIDTH,
        SLAB_STAGE2_END_WIDTH,
        REVIEW_TETHER_CONNECTION_WIDTH,
    ]
    component.info["layer15_stage_lengths_um"] = [RIDGE_TAPER_LENGTH, RIDGE_CONVERTER_LENGTH]
    component.info["layer10_post_pnc_length_um"] = (
        POST_PNC_TAPER_LENGTH + SLAB_STAGE2_LENGTH + SLAB_TAPER_LENGTH
    )
    component.info["layer15_post_pnc_length_um"] = RIDGE_TAPER_LENGTH + RIDGE_CONVERTER_LENGTH
    component.info["ridge_tip_centerline_offset_um"] = RIDGE_CONVERTER_Y_SHIFT
    component.info["ridge_separation_definition"] = RIDGE_SEPARATION_DEFINITION
    converter_ends = []
    for center, orientation in zip(connection_centers, connection_orientations, strict=True):
        angle = radians(float(orientation) + 180.0)
        tangent = np.array([cos(angle), sin(angle)])
        normal = np.array([-sin(angle), cos(angle)])
        converter_ends.append(
            list(
                map(
                    float,
                    np.asarray(center)
                    + (PNC_TAPER_LENGTH + RIDGE_TAPER_LENGTH + RIDGE_CONVERTER_LENGTH) * tangent
                    + RIDGE_CONVERTER_Y_SHIFT * normal,
                )
            )
        )
    component.info["ridge_converter_end_centers_um"] = converter_ends
    component.info["meander_straight_length_um"] = meander_straight_length
    component.info["meander_straight_segment_lengths_um"] = [meander_straight_length] * 3
    component.info["meander_straight_total_um"] = 3 * meander_straight_length
    component.info["meander_half_port_separation_um"] = meander_half_length
    component.info["meander_180_bend_length_um"] = euler_length
    component.info["meander_length_formula"] = "3 * L + 2 * Euler_180_length"
    component.info["total_meander_length_um"] = total_length
    component.info["left_pnc_bend_x"] = left_pnc_x
    component.info["connection_centers"] = connection_centers
    component.info["connection_orientations"] = connection_orientations
    component.info["tip_centers"] = [list(map(float, p)) for p in tips]
    component.info["label_bbox_um"] = label_bounds
    component.info["mask_layers"] = [list(map(int, layer)) for layer in component.layers]
    return component


# %% [10] CHIP PLACEMENT AND GLOBAL MARKERS
# =============================================================================
# Only whole-device translations occur here. All layers move together.
def add_chip_outline_and_markers(component: gf.Component, device_count: int) -> None:
    """Add only the chip outline and global overlay markers to a layout.

    All devices are already placed. No device polygons are edited here.
    """
    # Chip outline and global overlay markers are not part of the length sweep.
    chip_shape = (
        CHIP_WIDTH,
        max(CHIP_MIN_HEIGHT, (device_count - 1) * DEVICE_Y_PITCH + CHIP_Y_ALLOWANCE),
    )
    component.add_ref(gf.components.rectangle(size=chip_shape, centered=True, layer=(0, 0)))
    for marker_chip_shape, arm_length in (
        (chip_shape, GLOBAL_MARKER_ARM_LENGTHS[0]),
        (
            (chip_shape[0] - INNER_MARKER_CHIP_WIDTH_REDUCTION, chip_shape[1]),
            GLOBAL_MARKER_ARM_LENGTHS[1],
        ),
    ):
        component.add_ref(
            calipermarker_array(
                chip_shape=marker_chip_shape,
                arm_length=arm_length,
                arm_width=LAYER10_MARKER_ARM_WIDTH,
                center_square_size=LAYER10_MARKER_CENTER_SQUARE_SIZE,
                layer1=NEGATIVE_RESIST_LAYER,
                layer2=None,
            )
        )
        component.add_ref(
            calipermarker_array(
                chip_shape=marker_chip_shape,
                arm_length=arm_length + MLA_MARKER_SIZE_INCREASE,
                arm_width=MLA_MARKER_ARM_WIDTH,
                center_square_size=MLA_MARKER_CENTER_SQUARE_SIZE,
                layer1=MARKER_LAYER,
                layer2=None,
            )
        )


def build_layout() -> gf.Component:
    """Call make_device for each length; align entire devices, then add chip marks."""
    validate_settings()
    gf.gpdk.PDK.activate()
    fin = gf.Component(TOP_CELL_NAME)
    lengths = tuple(float(length) for length in MEANDER_STRAIGHT_LENGTHS)
    devices = [make_device(length) for length in lengths]
    # Translate each COMPLETE device. Do not edit, clip, or remap any layer.
    # The reference is the leftmost Layer-11 PNC edge at the left Euler bend.
    # Device 1 stays in place; the longer devices move right to align here.
    x_shifts = [LEFT_BEND_ALIGNMENT_X - float(device.info["left_pnc_bend_x"]) for device in devices]
    y_positions = [BOTTOM_DEVICE_Y + i * DEVICE_Y_PITCH for i in range(len(devices))]
    for device, x_shift, y_position in zip(devices, x_shifts, y_positions, strict=True):
        fin.add_ref(device).move((x_shift, y_position))

    add_chip_outline_and_markers(fin, len(devices))
    fin.info["design"] = DESIGN_ID
    fin.info["device_count"] = len(devices)
    fin.info["device_cells"] = [device.name for device in devices]
    fin.info["meander_straight_lengths"] = list(lengths)
    fin.info["total_meander_lengths"] = [float(d.info["total_meander_length_um"]) for d in devices]
    fin.info["device_x_shifts"] = [float(value) for value in x_shifts]
    fin.info["left_bend_alignment_x"] = LEFT_BEND_ALIGNMENT_X
    fin.info["alignment_reference_layer"] = list(PNC_EBL_LAYER)
    fin.info["alignment_reference"] = "leftmost_pnc_edge_of_left_euler_bend"
    fin.info["device_y_positions"] = y_positions
    fin.info["device_y_pitch"] = DEVICE_Y_PITCH
    fin.info["bottom_device_y"] = BOTTOM_DEVICE_Y
    fin.info["placement_mode"] = "independent_all_layer_factory_left_bend_aligned"
    return fin


# %% [11] SETTINGS CHECKS AND GDS EXPORT
# =============================================================================
def validate_settings() -> None:
    """Catch empty sweeps and invalid common dimensions before writing a GDS.

    These are basic input checks, NOT fabrication DRC or optical validation.
    """
    if not MEANDER_STRAIGHT_LENGTHS:
        raise ValueError("MEANDER_STRAIGHT_LENGTHS must contain at least one device")
    for length in MEANDER_STRAIGHT_LENGTHS:
        if not np.isfinite(float(length)) or float(length) <= 0:
            raise ValueError("Meander straight lengths must be finite and positive")
    dimensions = {
        "MEANDER_STRAIGHT_PNC_PERIOD": MEANDER_STRAIGHT_PNC_PERIOD,
        "DEVICE_Y_PITCH": DEVICE_Y_PITCH,
        "S_BEND_COMPACT_180_RADIUS": S_BEND_COMPACT_180_RADIUS,
        "RIDGE_MEANDER_WIDTH": RIDGE_MEANDER_WIDTH,
        "RIDGE_TAPER_LENGTH": RIDGE_TAPER_LENGTH,
        "SLAB_TAPER_LENGTH": SLAB_TAPER_LENGTH,
        "REVIEW_TETHER_CONNECTION_WIDTH": REVIEW_TETHER_CONNECTION_WIDTH,
        "POSITIVE_WINDOW_MARGIN": POSITIVE_WINDOW_MARGIN,
        "CHIP_WIDTH": CHIP_WIDTH,
        "CHIP_MIN_HEIGHT": CHIP_MIN_HEIGHT,
    }
    dimensions.update(
        {
            "TAPER_STAGE2_LENGTH": TAPER_STAGE2_LENGTH,
            "SLAB_STAGE1_END_WIDTH": SLAB_STAGE1_END_WIDTH,
            "SLAB_STAGE2_END_WIDTH": SLAB_STAGE2_END_WIDTH,
            "RIDGE_CONVERTER_END_WIDTH": RIDGE_CONVERTER_END_WIDTH,
        }
    )
    for name, value in dimensions.items():
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if not 0 <= S_BEND_EULER_P <= 1:
        raise ValueError("S_BEND_EULER_P must be between 0 and 1")
    if POSITIVE_WINDOW_END_EXTENSION < POSITIVE_WINDOW_MARGIN:
        raise ValueError("POSITIVE_WINDOW_END_EXTENSION must be >= POSITIVE_WINDOW_MARGIN")
    if not all(np.isfinite(v) for v in (BOTTOM_DEVICE_Y, LEFT_BEND_ALIGNMENT_X)):
        raise ValueError("Device placement coordinates must be finite")


def export_gds(output_path: str | Path | None = None) -> Path:
    """Build and write the complete layout; return the actual saved path.

    Default: GDS / OUTPUT_GDS_FILENAME. Supplying a path does not change
    geometry. An existing file at that path is overwritten.
    """
    path = Path(output_path) if output_path is not None else GDS / OUTPUT_GDS_FILENAME
    if path.suffix.lower() != ".gds":
        raise ValueError("The output filename must end in .gds")
    component = build_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    options = gf.kdb.SaveLayoutOptions()
    options.gds2_max_vertex_count = GDS_MAX_VERTEX_COUNT
    component.write_gds(path, save_options=options, with_metadata=False)
    return path


if __name__ == "__main__":
    saved_path = export_gds()
    print(f"Wrote {saved_path}")
