"""Three-device two-EBL Euler-meander review.

Layer 10 is limited to taper geometry at or below 3 um width.
The broad 15 um meander body is omitted from the exported mask layers.
Literal GDS layer 11 carries each solid PNC cell for negative-resist EBL.
Layer 2 carries a continuous expanded positive-resist region for MLA exposure;
it has no PNC-shaped exclusions or openings.
The narrow PNC neck extends 0.5 um beyond both radial (local y) ends.
The positive release windows, long tapered slots, support surround slots,
glue holes, and former Layer-12 geometry are all assigned to Layer 2 for one
shared MLA exposure. No Boolean geometry merge is applied.
"""

from functools import partial
from math import cos, radians, sin
import os
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


# =============================================================================
# 1. MASK LAYERS AND FABRICATION TONES
# =============================================================================

NEGATIVE_RESIST_LAYER = Layers.wg       # GDS 10/0: negative-resist fine mask
# Temporary construction layer used to route and position the meander. It is
# removed recursively before the final GDS is written.
MEANDER_CONSTRUCTION_LAYER = (999, 0)
# This device intentionally repurposes literal GDS 2/0 as its second EBL,
# positive-resist full-etch mask around the PNC features.
FULL_ETCH_LAYER = Layers.globalmarker  # GDS 2/0: surrounding MLA open mask
PNC_EBL_LAYER = Layers.wg_doubleetch   # GDS 11/0: negative-resist PNC
PNC_NARROW_END_EXTENSION = 0.5
LAYER2_LAYER10_COPY_DELAY = 0.5
POSITIVE_RESIST_LAYER = FULL_ETCH_LAYER  # GDS 2/0: shared MLA open mask
MARKER_LAYER = FULL_ETCH_LAYER           # GDS 2/0: markers and labels
RIDGE_LAYER = Layers.mc                  # GDS 15/0: ridge taper
RIDGE_PROTECTION_LAYER = Layers.mc2      # GDS 16/0: ridge protection
RIDGE_PROTECTION_MARGIN = 0.2


# =============================================================================
# 2. MLA RELEASE-WINDOW AND LONG-CHANNEL SETTINGS
# =============================================================================

POSITIVE_WINDOW_MARGIN = 3.0
# Exact shared-MLA window endpoint beyond the physical taper end.
# Keep this at least as large as POSITIVE_WINDOW_MARGIN.
POSITIVE_WINDOW_END_EXTENSION = 15.0


# =============================================================================
# 3. OPTICAL TAPER AND BUTTERFLY TETHER SETTINGS
# =============================================================================

# Variant 2 is the library's split/forked butterfly-style tether.
TETHER_VARIANT = 2
REVIEW_TETHER_CONNECTION_WIDTH = 0.2
TAPER_ANGLE = 0.0
MEANDER_STRAIGHT_LENGTHS = (200.0, 800.0, 1400.0)
REVIEW_SUPPORT_ANCHOR_PERIOD = 15.0
REVIEW_SUPPORT_BEND_PERIOD = 3.0
PNC_TAPER_SUPPORT_CUTOFF_WIDTH = 3.0
# Phononic-crystal support copied from layout.acoustics.phononics. The thin
# beam follows the taper tip width and carries the library's six square masses.
PHONONIC_SUPPORT_SEGMENTS = (0.626, 0.600, 0.560, 0.540, 0.516, 0.504)
PHONONIC_SUPPORT_GAPS = 0.25
PHONONIC_SUPPORT_OVERSHOOT = PNC_NARROW_END_EXTENSION
RIDGE_TAPER_LENGTH = 40.0
SLAB_TAPER_LENGTH = 20.0
TETHER_CONNECTION_BLEND_LENGTH = 2.0
TETHER_PAD_SCALE = 1.0
# Independent full taper angles for the two long slot layers.
MLA_CHANNEL_ANGLE = 2.0
LAYER14_CHANNEL_ANGLE = 1.5
MLA_CHANNEL_START_WIDTH = 6.0
LAYER14_CHANNEL_START_WIDTH = 1.0
# Keep each Layer-14 tapered channel beginning exactly at the optical tip.
# Glue-circle positions are calculated independently and do not move.
LAYER14_CHANNEL_START_OFFSET = 0.0
# The shared-MLA and restored Layer-14 tapered slots are 2000 um long.
POSITIVE_CHANNEL_LENGTH = 2000.0
S_BEND_WIDTH = 15.0


# =============================================================================
# 4. EULER MEANDER GEOMETRY
# =============================================================================

S_BEND_TAPER_LENGTH = 10.0
S_BEND_HALF_LENGTH = 150.0
S_BEND_COMPACT_180_RADIUS = 100.0
# Euler fraction: 0=circular, 0.5=partial Euler, 1=full Euler.
S_BEND_EULER_P = 0.2
# Two opposite 180-degree bends shift the route by four center radii.
S_BEND_VERTICAL_HALF_SPAN = 2 * S_BEND_COMPACT_180_RADIUS
# The PNC section spans exactly three 15 um periods. At its final cell, the
# Layer-10 slab is 3 um wide; both post-PNC adiabatic tapers start there.
PNC_TAPER_LENGTH = 3 * REVIEW_SUPPORT_ANCHOR_PERIOD
POST_PNC_TAPER_LENGTH = RIDGE_TAPER_LENGTH
SLAB_EXTENSION_LENGTH = PNC_TAPER_LENGTH + POST_PNC_TAPER_LENGTH
SLAB_EXTENSION_START_WIDTH = 1.0
SLAB_EXTENSION_END_WIDTH = 15.0
PNC_TAPER_END_WIDTH = PNC_TAPER_SUPPORT_CUTOFF_WIDTH
RIDGE_EXTENSION_START_WIDTH = 0.1
RIDGE_MEANDER_WIDTH = 1.0
RIDGE_EXTENSION_END_WIDTH = RIDGE_MEANDER_WIDTH
SUPPORT_TIP_POSITION = 30.0
SUPPORT_TIP_WIDTH = 0.5
SUPPORT_TIP_LENGTH = 10.0
SUPPORT_TIP_OVERLAP = 0.2
# ACTIVE S-BEND SUPPORT GEOMETRY.
# Two shared-MLA side windows along each existing 0.5 x 10 um alignment tip.
ALIGNMENT_TIP_WINDOW_WIDTH = 2.0
ALIGNMENT_TIP_WINDOW_LENGTH = 8.0
ALIGNMENT_TIP_WINDOW_INWARD_SHIFT = 0.2
DEVICE_Y_PITCH = 700.0


# =============================================================================
# 5. DEVICE PLACEMENT AND 600 UM EBL SUBFIELDS
# =============================================================================

# Align the left meander connection for every swept straight length. Three
# 600 um EBL subfield columns can then share the same x boundaries.
LEFT_BEND_ALIGNMENT_X = -900.0
EBL_SUBFIELD_SIZE = 600.0
TETHER_LABEL_SIZE = 2.0
ADIABATIC_TAPER_POINTS = 201
# Detailed per-device label placed beside the left Layer-15 taper endpoint.
DETAIL_LABEL_TIP_X_OFFSET = 5.0
DETAIL_LABEL_TIP_Y_OFFSET = 10.0
DETAIL_LABEL_SIZE = 3.0
# Corner overlay markers: the MLA marker is larger than the Layer-10 marker.
LAYER10_MARKER_ARM_WIDTH = 2.0
LAYER10_MARKER_CENTER_SQUARE_SIZE = 3.0
MLA_MARKER_SIZE_INCREASE = 10.0
MLA_MARKER_ARM_WIDTH = 3.0
MLA_MARKER_CENTER_SQUARE_SIZE = 5.0

# The support island remains; only its border is exposed on the shared MLA mask.
SUPPORT_RECTANGLE_SIZE = (50.0, 40.0)
SUPPORT_SURROUND_SLOT_WIDTH = 5.0
EXTRA_RECTANGLE_GAP = 100.0
# Start the long shared-MLA tapered slot at the end of the short rectangular
# taper-window extension, not at the distant support rectangle.
MLA_CHANNEL_START_OFFSET = POSITIVE_WINDOW_END_EXTENSION


# =============================================================================
# 6. PHONONIC-CRYSTAL SUPPORT COMPONENTS
# =============================================================================

@gf.cell
def phononic_meander_support(
    width: float = REVIEW_TETHER_CONNECTION_WIDTH,
    segments: tuple[float, ...] = PHONONIC_SUPPORT_SEGMENTS,
    gaps: float = PHONONIC_SUPPORT_GAPS,
    overshoot: float = PHONONIC_SUPPORT_OVERSHOOT,
    layer=NEGATIVE_RESIST_LAYER,
) -> gf.Component:
    """Layer-remapped library PNC tether for the meander supports."""
    source = phononic_tether(
        width=width,
        segments=segments,
        gaps=gaps,
        overshoot=overshoot,
    ).copy()
    source.flatten()
    source.remap_layers(
        {Layers.wg_doubleetch: layer},
        recursive=False,
    )
    component = gf.Component()
    reference = component.add_ref(source)
    component.add_port("a1", port=reference.ports["a1"])
    component.add_port("a2", port=reference.ports["a2"])
    component.info["support_type"] = "phononic_crystal"
    return component


@gf.cell
def phononic_meander_shield(
    path,
    width: float = REVIEW_TETHER_CONNECTION_WIDTH,
    segments: tuple[float, ...] = PHONONIC_SUPPORT_SEGMENTS,
    gaps: float = PHONONIC_SUPPORT_GAPS,
    overshoot: float = PHONONIC_SUPPORT_OVERSHOOT,
    pitch: float = REVIEW_SUPPORT_ANCHOR_PERIOD,
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
            section_ranges.append(
                (section_start, index, bool(section_types[section_start]))
            )
            section_start = index
    section_ranges.append(
        (section_start, len(section_types), bool(section_types[section_start]))
    )

    tether = phononic_tether(
        width=width,
        segments=segments,
        gaps=gaps,
        overshoot=overshoot,
    )
    tether_length = float(
        tether.ports["a2"].center[0] - tether.ports["a1"].center[0]
    )
    cell_extent = max(float(max(segments)), width)
    cell_half_length = cell_extent / 2 + 1e-3

    # Make one continuous band around the complete meander. PNC cells are
    # subtracted only after they have been placed as complete library cells,
    # so no cell is clipped at a straight/Euler section boundary.
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
        segment_fraction = (
            path_position - cumulative_lengths[segment_index]
        ) / lengths[segment_index]
        center = (
            points[segment_index]
            + segment_fraction * vectors[segment_index]
        )
        tangent = (
            vectors[segment_index] / lengths[segment_index]
            if tangent_override is None
            else tangent_override / np.linalg.norm(tangent_override)
        )
        normal = np.array([-tangent[1], tangent[0]])
        tangent_angle = float(
            np.degrees(np.arctan2(tangent[1], tangent[0]))
        )
        for side in (-1, 1):
            anchor = center + side * normal * slab_width / 2
            full_cells.add_ref(tether).rotate(
                tangent_angle + side * 90
            ).move(tuple(anchor))
        support_positions.append(float(path_position))

    # A single shared PNC is placed exactly at every straight/bend boundary.
    # Keep the first interior PNC one complete local pitch away from it, rather
    # than only one cell length away, to avoid a dense three-PNC cluster.
    for start, stop, is_bend in section_ranges:
        section_start_length = cumulative_lengths[start]
        section_end_length = cumulative_lengths[stop]
        section_length = section_end_length - section_start_length
        requested_pitch = bend_pitch if is_bend else pitch
        start_clearance = (
            requested_pitch if start > 0 else cell_half_length
        )
        end_clearance = (
            requested_pitch
            if stop < len(lengths)
            else cell_half_length
        )
        usable_start = section_start_length + start_clearance
        usable_end = section_end_length - end_clearance
        if usable_end < usable_start:
            continue
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
        component.info["layer14_full_open_width"] = (
            slab_width + 2 * tether_length
        )
    component.info["pnc_narrow_end_extension"] = end_extension
    component.info["full_pnc_cells"] = True
    component.info["support_positions"] = support_positions
    component.info["connection_pnc_positions"] = connection_positions
    return component


@gf.cell
def phononic_taper_shield(
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
    end_extension: float = PNC_NARROW_END_EXTENSION,
    npoints: int = ADIABATIC_TAPER_POINTS,
) -> gf.Component:
    """Draw solid negative-resist PNC cells and surrounding taper MLA mask."""
    if not width2 <= cutoff_width < width1:
        raise ValueError("Require width2 <= cutoff_width < width1")
    if length <= 0 or pitch <= 0:
        raise ValueError("Taper length and PNC pitch must be positive")
    if end_extension < 0:
        raise ValueError("PNC narrow-end extension must be non-negative")

    # Solve the same quintic width profile used by tapered_fiber_adiabatic.
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
    tether_length = float(
        tether.ports["a2"].center[0] - tether.ports["a1"].center[0]
    )
    cell_extent = max(float(max(segments)), tether_width)
    cell_half_length = cell_extent / 2 + 1e-3

    outer_bands = gf.Component()
    inner_edge = widths / 2
    outer_edge = inner_edge + tether_length
    # Extend the PNC mask by half a cell at both longitudinal ends. The last
    # cell remains centered exactly at the 3 um taper endpoint but is no
    # longer cut in half by the end of the surrounding positive-resist band.
    band_end_extension = cell_half_length
    band_x = np.concatenate(
        ([-band_end_extension], x, [cutoff_x + band_end_extension])
    )
    band_inner_edge = np.concatenate(
        ([inner_edge[0]], inner_edge, [inner_edge[-1]])
    )
    band_outer_edge = np.concatenate(
        ([outer_edge[0]], outer_edge, [outer_edge[-1]])
    )
    upper_polygon = list(
        zip(band_x, band_inner_edge, strict=True)
    ) + list(
        zip(band_x[::-1], band_outer_edge[::-1], strict=True)
    )
    lower_polygon = [(xi, -yi) for xi, yi in upper_polygon]
    outer_bands.add_polygon(upper_polygon, layer=Layers.wg_doubleetch)
    outer_bands.add_polygon(lower_polygon, layer=Layers.wg_doubleetch)

    cuts = gf.Component()
    support_positions = np.arange(pitch, cutoff_x + 1e-9, pitch)
    for support_x in support_positions:
        support_u = support_x / length
        support_weight = (
            6 * support_u**5 - 15 * support_u**4 + 10 * support_u**3
        )
        support_width = width1 + (width2 - width1) * support_weight
        cuts.add_ref(tether).rotate(90).move(
            (support_x, support_width / 2)
        )
        cuts.add_ref(tether).rotate(-90).move(
            (support_x, -support_width / 2)
        )

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
    component.info["support_positions"] = [
        float(value) for value in support_positions
    ]
    component.info["pnc_ebl_layer"] = tuple(pnc_layer)
    component.info["pnc_mla_layer"] = tuple(layer)
    component.info["pnc_ebl_tone"] = "negative_solid_pnc"
    component.info["pnc_mla_tone"] = "positive_continuous_expanded_area"
    component.info["pnc_narrow_end_extension"] = end_extension
    return component


# =============================================================================
# 7. SUPPORT-ISLAND AND ADIABATIC-TAPER COMPONENTS
# =============================================================================

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


@gf.cell
def tapered_fiber_adiabatic(
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
    lower = [
        (float(xi), float(-wi / 2))
        for xi, wi in zip(x[::-1], widths[::-1])
    ]
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
def two_layer_tapered_fiber(
    tether_connection_width: float = 0.2,
    ridge_taper_length: float = RIDGE_TAPER_LENGTH,
    add_support_tips: bool = True,
    add_phononic_taper_support: bool = False,
    slab_layer=NEGATIVE_RESIST_LAYER,
) -> gf.Component:
    """Piecewise-adiabatic slab and ridge transition after the meander.

    Along +x, the 15->3 um PNC transition is excluded from Layer 10. Layer 10
    starts again at the 3 um endpoint, then tapers 3->1 um while Layer 15
    begins its aligned 1->0.1 um taper. Layer 10 finally continues from
    1 um to the tether width.
    """
    if not 0 < ridge_taper_length <= POST_PNC_TAPER_LENGTH:
        raise ValueError(
            "ridge_taper_length must satisfy 0 < length <= "
            "POST_PNC_TAPER_LENGTH"
        )
    component = gf.Component()

    # The first transition is intentionally not drawn on Layer 10 because its
    # width is greater than 3 um. Keep its longitudinal distance so all PNC,
    # ridge, and downstream taper coordinates remain unchanged.
    post_pnc_slab_taper = component.add_ref(
        tapered_fiber_adiabatic(
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
    ridge_pnc_section = component.add_ref(
        tapered_fiber_adiabatic(
            length=PNC_TAPER_LENGTH,
            width1=RIDGE_MEANDER_WIDTH,
            width2=RIDGE_MEANDER_WIDTH,
            taper_length=PNC_TAPER_LENGTH,
            layer=RIDGE_LAYER,
        )
    )

    ridge_extension = component.add_ref(
        tapered_fiber_adiabatic(
            length=ridge_taper_length,
            width1=RIDGE_EXTENSION_END_WIDTH,
            width2=RIDGE_EXTENSION_START_WIDTH,
            taper_length=ridge_taper_length,
            layer=RIDGE_LAYER,
        )
    )
    ridge_extension.connect("o1", ridge_pnc_section.ports["o2"])

    outer_slab_taper = component.add_ref(
        tapered_fiber_adiabatic(
            length=SLAB_TAPER_LENGTH,
            width1=SLAB_EXTENSION_START_WIDTH,
            width2=tether_connection_width,
            taper_length=SLAB_TAPER_LENGTH,
            layer=slab_layer,
        )
    )
    outer_slab_taper.connect("o1", post_pnc_slab_taper.ports["o2"])

    # Continue Layer 2 from the 3 um Layer-10 starting point to the Layer-15
    # endpoint with a quintic adiabatic profile. The two inner edges close to
    # the centerline while the outer edges reach a total y-size of 7 um.
    release_tether = phononic_tether(
        width=tether_connection_width,
        segments=PHONONIC_SUPPORT_SEGMENTS,
        gaps=PHONONIC_SUPPORT_GAPS,
        overshoot=PHONONIC_SUPPORT_OVERSHOOT,
    )
    release_band_width = float(
        release_tether.ports["a2"].center[0]
        - release_tether.ports["a1"].center[0]
    )
    transition_x = np.linspace(
        0.0,
        RIDGE_TAPER_LENGTH,
        ADIABATIC_TAPER_POINTS,
    )
    transition_u = transition_x / RIDGE_TAPER_LENGTH
    transition_weight = (
        6 * transition_u**5
        - 15 * transition_u**4
        + 10 * transition_u**3
    )
    release_x = PNC_TAPER_LENGTH + transition_x
    start_inner = PNC_TAPER_END_WIDTH / 2
    start_outer = start_inner + release_band_width
    release_target_half_y = SLAB_EXTENSION_START_WIDTH / 2 + (
        POSITIVE_WINDOW_MARGIN
    )
    release_inner = start_inner * (1 - transition_weight)
    release_outer = start_outer + (
        release_target_half_y - start_outer
    ) * transition_weight
    upper_release = list(zip(release_x, release_inner, strict=True)) + list(
        zip(release_x[::-1], release_outer[::-1], strict=True)
    )
    lower_release = [(x_value, -y_value) for x_value, y_value in upper_release]
    component.add_polygon(upper_release, layer=FULL_ETCH_LAYER)
    component.add_polygon(lower_release, layer=FULL_ETCH_LAYER)

    # Copy the first Layer-10 section (the 3->1 um quintic taper) onto Layer 2,
    # but delay its beginning by 0.5 um so it does not reach into the PNC
    # boundary. The copied profile remains identical to Layer 10 thereafter.
    copy_x_local = np.linspace(
        LAYER2_LAYER10_COPY_DELAY,
        POST_PNC_TAPER_LENGTH,
        ADIABATIC_TAPER_POINTS,
    )
    copy_u = copy_x_local / POST_PNC_TAPER_LENGTH
    copy_weight = 6 * copy_u**5 - 15 * copy_u**4 + 10 * copy_u**3
    copy_width = PNC_TAPER_END_WIDTH + (
        SLAB_EXTENSION_START_WIDTH - PNC_TAPER_END_WIDTH
    ) * copy_weight
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

    if add_phononic_taper_support:
        component.add_ref(
            phononic_taper_shield(
                length=PNC_TAPER_LENGTH,
                width1=SLAB_EXTENSION_END_WIDTH,
                width2=PNC_TAPER_END_WIDTH,
                cutoff_width=PNC_TAPER_SUPPORT_CUTOFF_WIDTH,
                pitch=REVIEW_SUPPORT_ANCHOR_PERIOD,
                tether_width=tether_connection_width,
                layer=FULL_ETCH_LAYER,
                pnc_layer=PNC_EBL_LAYER,
                end_extension=PNC_NARROW_END_EXTENSION,
            )
        )

    # A lateral support pair stabilizes the longer suspended transition. Each
    # tip overlaps the local slab edge slightly and extends away from both
    # sides of the waveguide.
    support_u = SUPPORT_TIP_POSITION / PNC_TAPER_LENGTH
    support_weight = 6 * support_u**5 - 15 * support_u**4 + 10 * support_u**3
    local_slab_width = SLAB_EXTENSION_END_WIDTH + (
        PNC_TAPER_END_WIDTH - SLAB_EXTENSION_END_WIDTH
    ) * support_weight
    support_size = (
        SUPPORT_TIP_WIDTH,
        SUPPORT_TIP_LENGTH + SUPPORT_TIP_OVERLAP,
    )
    support_center_y = (
        local_slab_width / 2
        + (SUPPORT_TIP_LENGTH - SUPPORT_TIP_OVERLAP) / 2
    )
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
    component.info["slab_extension_length"] = (
        PNC_TAPER_LENGTH + POST_PNC_TAPER_LENGTH
    )
    component.info["slab_taper_length"] = SLAB_TAPER_LENGTH
    component.info["ridge_taper_length"] = ridge_taper_length
    component.info["support_tip_size"] = list(support_size)
    component.info["support_tip_position"] = SUPPORT_TIP_POSITION
    component.info["support_tips_added"] = add_support_tips
    component.info["phononic_taper_support_added"] = (
        add_phononic_taper_support
    )
    component.info["layer10_max_taper_width"] = PNC_TAPER_END_WIDTH
    component.info["layer2_end_distance"] = (
        PNC_TAPER_LENGTH + RIDGE_TAPER_LENGTH
    )
    component.info["layer2_end_y_size"] = 2 * release_target_half_y
    component.info["layer2_layer10_copy_delay"] = (
        LAYER2_LAYER10_COPY_DELAY
    )
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


# =============================================================================
# 8. COUPLER ASSEMBLY AND SHARED MLA MASK
# =============================================================================

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
            "tether_connection_width must be positive and narrower than "
            "S_BEND_WIDTH"
        )

    taper_spec = partial(
        two_layer_tapered_fiber,
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
        alignment_tick_window_layer=POSITIVE_RESIST_LAYER,
        alignment_tick_window_width=ALIGNMENT_TIP_WINDOW_WIDTH,
        alignment_tick_window_length=ALIGNMENT_TIP_WINDOW_LENGTH,
        alignment_tick_window_inward_shift=(
            ALIGNMENT_TIP_WINDOW_INWARD_SHIFT
        ),
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
        oxide_open_duplicate_length=2000.0,
        oxide_open_duplicate_start_width=LAYER14_CHANNEL_START_WIDTH,
        oxide_open_duplicate_start_at_taper_end=False,
        oxide_open_duplicate_start_offset=(
            SLAB_EXTENSION_LENGTH + LAYER14_CHANNEL_START_OFFSET
        ),
    )


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

    # Clip the former Layer-12 shape at the Layer-15 endpoint before assigning
    # it to the shared Layer-2 MLA mask.
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
        release_start = np.asarray(center, dtype=float) + (
            release_start_distance * outward
        )
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

    ridge_protection = gf.Component()
    ridge_protection.add_ref(device)
    ridge_protection.flatten()
    ridge_layer = tuple(RIDGE_LAYER)
    ridge_protection.remove_layers(
        [
            layer
            for layer in ridge_protection.layers
            if tuple(layer) != ridge_layer
        ]
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
    # receives a SOLID Layer-14 opening from phononic_meander_shield(), bounded
    # exactly by the OUTER edges of Layer 2. Here we keep the existing Layer-2
    # copy through the adjacent PNC/taper regions and clip it at the current
    # Layer-14 channel starts.
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
            + (
                SLAB_EXTENSION_LENGTH
                + LAYER14_CHANNEL_START_OFFSET
            )
            * outward
        )

    layer14_starts = np.asarray(layer14_start_points, dtype=float)
    layer14_xmin = float(layer14_starts[:, 0].min())
    layer14_xmax = float(layer14_starts[:, 0].max())

    # Copy only the device's meander/PNC release geometry.  Deliberately do
    # not copy ``positive_window`` because that is the expanded taper-tip
    # region on Layer 2 and is not wanted on Layer 14.
    layer14_pattern_source = device.copy()
    layer14_pattern_source.flatten()
    layer14_pattern_source.remove_layers(
        [
            layer
            for layer in layer14_pattern_source.layers
            if tuple(layer) != positive_layer
        ]
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
    result.add_ref(layer14_central_pattern)

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
    # The broad meander route is needed only to calculate ports and PNC
    # placement. Do not export that temporary geometry on any mask layer.
    result.remove_layers([MEANDER_CONSTRUCTION_LAYER], recursive=True)
    result.info["mask_mode"] = "dual_resist_shared_etch"
    result.info["negative_resist_layer"] = negative_layer
    result.info["positive_resist_layer"] = positive_layer
    result.info["positive_window_margin"] = float(margin)
    result.info["positive_window_end_extension"] = float(end_extension)
    result.info["layer2_transition_endpoint_distance"] = (
        release_start_distance
    )
    result.info["layer2_release_target_y_size"] = (
        SLAB_EXTENSION_START_WIDTH + 2 * margin
    )
    result.info["layer2_adiabatic_taper_length"] = RIDGE_TAPER_LENGTH
    result.info["ridge_protection_layer"] = tuple(RIDGE_PROTECTION_LAYER)
    result.info["ridge_protection_margin"] = RIDGE_PROTECTION_MARGIN
    result.info["layer14_central_pattern_x_bounds"] = [
        layer14_xmin,
        layer14_xmax,
    ]
    return result


# =============================================================================
# 9. TOP-LEVEL DEVICE SWEEP AND PLACEMENT
# =============================================================================

def build_layout() -> gf.Component:
    """Build three review devices sweeping only meander straight length."""
    name = Path(__file__).stem
    top_cell_name = f"{name}_L10_L12"

    gf.gpdk.PDK.activate()
    fin = gf.Component(top_cell_name)
    sweep_parameters = MEANDER_STRAIGHT_LENGTHS
    euler_180_length = float(
        gf.path.euler(
            radius=S_BEND_COMPACT_180_RADIUS,
            angle=180,
            p=S_BEND_EULER_P,
            use_eff=True,
        ).length()
    )
    layer10_taper_length = POST_PNC_TAPER_LENGTH + SLAB_TAPER_LENGTH
    device_count = len(sweep_parameters)
    chip_height = max(10e3, (device_count - 1) * DEVICE_Y_PITCH + 1000)
    chip_shape = (5e3, chip_height)

    fin.add_ref(
        gf.components.rectangle(size=chip_shape, centered=True, layer=(0, 0))
    )
    for marker_chip_shape, arm_length in (
        ((chip_shape[0], chip_shape[1]), 50),
        ((chip_shape[0] - 1000, chip_shape[1]), 80),
    ):
        # Smaller first-layer reference marker.
        fin.add_ref(
            calipermarker_array(
                chip_shape=marker_chip_shape,
                arm_length=arm_length,
                arm_width=LAYER10_MARKER_ARM_WIDTH,
                center_square_size=LAYER10_MARKER_CENTER_SQUARE_SIZE,
                layer1=NEGATIVE_RESIST_LAYER,
                layer2=None,
            )
        )
        # Larger second-layer marker centered on the Layer-10 marker.
        fin.add_ref(
            calipermarker_array(
                chip_shape=marker_chip_shape,
                arm_length=arm_length + MLA_MARKER_SIZE_INCREASE,
                arm_width=MLA_MARKER_ARM_WIDTH,
                center_square_size=MLA_MARKER_CENTER_SQUARE_SIZE,
                layer1=MARKER_LAYER,
                layer2=None,
            )
        )

    y_positions = tuple(
        (index - (device_count - 1) / 2) * DEVICE_Y_PITCH
        for index in range(device_count)
    )

    device_x_shifts = []
    left_pnc_bend_x_before_alignment = []
    total_meander_lengths = []
    for y_position, meander_straight_length in zip(
        y_positions,
        sweep_parameters,
        strict=True,
    ):
        connection_width = REVIEW_TETHER_CONNECTION_WIDTH
        taper_angle = TAPER_ANGLE
        support_period = REVIEW_SUPPORT_ANCHOR_PERIOD
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
            s_bend_half_length=S_BEND_HALF_LENGTH,
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
        device = tapered_fiber_inout(
            **common_device_kwargs,
            s_bend_support_shield_spec=partial(
                phononic_meander_shield,
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
            extra_rectangle_surround_slot_width=(
                SUPPORT_SURROUND_SLOT_WIDTH
            ),
            tether_spec=tether_spec,
        )
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
        # The taper connection remains at the same x for every length, so it
        # is not a useful bend-alignment reference. Use the left edge of the
        # actual Layer-11 PNC region, which follows the left Euler bend.
        left_pnc_bend_x = float(
            device.bbox(layer=gf.get_layer(PNC_EBL_LAYER)).left
        )
        device_x_shift = LEFT_BEND_ALIGNMENT_X - left_pnc_bend_x
        placed_device = fin.add_ref(
            add_shared_mla_window(device, window_source)
        )
        placed_device.move((device_x_shift, y_position))
        device_x_shifts.append(float(device_x_shift))
        left_pnc_bend_x_before_alignment.append(left_pnc_bend_x)

        # This matches the compact-meander path construction in tapered_local:
        # two outer leads, one middle straight, and two 180-degree Euler bends.
        half_lead = (
            2 * S_BEND_HALF_LENGTH + meander_straight_length
        ) / 2
        total_meander_length = (
            2 * half_lead
            + meander_straight_length
            + 2 * euler_180_length
        )
        total_meander_lengths.append(float(total_meander_length))

        layer15_endpoints = []
        for center, orientation in zip(
            window_source.info["s_bend_connection_centers"],
            window_source.info["s_bend_connection_orientations"],
            strict=True,
        ):
            angle = radians(orientation)
            outward = np.array([-cos(angle), -sin(angle)])
            layer15_endpoints.append(
                np.asarray(center, dtype=float)
                + (PNC_TAPER_LENGTH + RIDGE_TAPER_LENGTH) * outward
            )
        left_layer15_endpoint = min(
            layer15_endpoints,
            key=lambda point: float(point[0]),
        )

        # Keep the written device label compact: total meander centerline
        # length and total Layer-10 taper length only.
        detail_text = (
            f"ML{total_meander_length:.1f}_"
            f"TL{layer10_taper_length:g}"
        )
        detail_label = fin.add_ref(
            gf.components.text(
                text=detail_text,
                size=DETAIL_LABEL_SIZE,
                layer=MARKER_LAYER,
            )
        )
        detail_label.move(
            (
                device_x_shift
                + float(left_layer15_endpoint[0])
                + DETAIL_LABEL_TIP_X_OFFSET,
                y_position
                + float(left_layer15_endpoint[1])
                + DETAIL_LABEL_TIP_Y_OFFSET,
            )
        )

    fin.info["device_count"] = device_count
    fin.info["left_bend_alignment_x"] = LEFT_BEND_ALIGNMENT_X
    fin.info["device_x_shifts"] = device_x_shifts
    fin.info["left_pnc_bend_x_before_alignment"] = (
        left_pnc_bend_x_before_alignment
    )
    fin.info["ebl_subfield_size"] = EBL_SUBFIELD_SIZE
    fin.info["ebl_subfield_columns"] = 3
    fin.info["total_meander_lengths"] = total_meander_lengths
    fin.info["layer10_taper_length"] = layer10_taper_length
    fin.info["euler_180_length"] = euler_180_length
    fin.info["tether_connection_widths"] = [
        REVIEW_TETHER_CONNECTION_WIDTH
    ]
    fin.info["tether_variant"] = TETHER_VARIANT
    fin.info["slab_taper_length"] = SLAB_TAPER_LENGTH
    fin.info["ridge_taper_length"] = RIDGE_TAPER_LENGTH
    fin.info["taper_angles"] = [TAPER_ANGLE]
    fin.info["s_bend_width"] = S_BEND_WIDTH
    fin.info["s_bend_half_length"] = S_BEND_HALF_LENGTH
    fin.info["taper_profile"] = "quintic_adiabatic"
    fin.info["adiabatic_taper_points"] = ADIABATIC_TAPER_POINTS
    fin.info["ridge_layer"] = tuple(RIDGE_LAYER)
    fin.info["ridge_protection_layer"] = tuple(RIDGE_PROTECTION_LAYER)
    fin.info["ridge_protection_margin"] = RIDGE_PROTECTION_MARGIN
    fin.info["slab_extension_length"] = SLAB_EXTENSION_LENGTH
    fin.info["pnc_taper_length"] = PNC_TAPER_LENGTH
    fin.info["post_pnc_taper_length"] = POST_PNC_TAPER_LENGTH
    fin.info["support_anchor_periods"] = [
        REVIEW_SUPPORT_ANCHOR_PERIOD
    ]
    fin.info["support_bend_period"] = REVIEW_SUPPORT_BEND_PERIOD
    fin.info["pnc_taper_support_cutoff_width"] = (
        PNC_TAPER_SUPPORT_CUTOFF_WIDTH
    )
    fin.info["meander_bend_center_radius"] = (
        S_BEND_COMPACT_180_RADIUS
    )
    fin.info["meander_euler_fraction"] = S_BEND_EULER_P
    fin.info["meander_turn_count"] = 2
    fin.info["meander_straight_lengths"] = list(
        MEANDER_STRAIGHT_LENGTHS
    )
    fin.info["slab_extension_widths"] = [
        SLAB_EXTENSION_END_WIDTH,
        PNC_TAPER_END_WIDTH,
        SLAB_EXTENSION_START_WIDTH,
    ]
    fin.info["ridge_extension_widths"] = [
        RIDGE_EXTENSION_START_WIDTH,
        RIDGE_EXTENSION_END_WIDTH,
    ]
    fin.info["alignment_tip_window_size"] = [
        ALIGNMENT_TIP_WINDOW_WIDTH,
        ALIGNMENT_TIP_WINDOW_LENGTH,
    ]
    fin.info["alignment_tip_slab_junction_offset"] = SLAB_EXTENSION_LENGTH
    fin.info["sweep_parameters"] = [
        [float(length)] for length in sweep_parameters
    ]
    fin.info["sweep_parameter_names"] = [
        "meander_straight_length_um",
    ]
    fin.info["mask_mode"] = "dual_resist_shared_etch"
    fin.info["interlayer_boolean_applied"] = False
    fin.info["negative_resist_layer"] = tuple(NEGATIVE_RESIST_LAYER)
    fin.info["full_etch_layer"] = tuple(FULL_ETCH_LAYER)
    fin.info["meander_layer"] = None
    fin.info["layer10_max_width"] = PNC_TAPER_END_WIDTH
    fin.info["pnc_ebl_layer"] = tuple(PNC_EBL_LAYER)
    fin.info["pnc_mla_layer"] = tuple(FULL_ETCH_LAYER)
    fin.info["pnc_ebl_tone"] = "negative_solid_pnc"
    fin.info["pnc_mla_tone"] = "positive_continuous_expanded_area"
    fin.info["pnc_narrow_end_extension"] = PNC_NARROW_END_EXTENSION
    fin.info["output_taper_layer"] = tuple(NEGATIVE_RESIST_LAYER)
    fin.info["positive_resist_layer"] = tuple(POSITIVE_RESIST_LAYER)
    fin.info["positive_window_margin"] = POSITIVE_WINDOW_MARGIN
    fin.info["positive_window_end_extension"] = (
        POSITIVE_WINDOW_END_EXTENSION
    )
    fin.info["marker_layer"] = tuple(MARKER_LAYER)
    fin.info["detail_label_tip_offset"] = [
        DETAIL_LABEL_TIP_X_OFFSET,
        DETAIL_LABEL_TIP_Y_OFFSET,
    ]
    fin.info["detail_label_size"] = DETAIL_LABEL_SIZE
    fin.info["tether_pad_scale"] = TETHER_PAD_SCALE
    fin.info["support_rectangle_size"] = list(SUPPORT_RECTANGLE_SIZE)
    fin.info["support_surround_slot_width"] = SUPPORT_SURROUND_SLOT_WIDTH
    fin.info["meander_support_type"] = (
        "periodic_phononic_crystal_outer_slab"
    )
    fin.info["phononic_connection_style"] = "solid_negative_resist_pnc"
    fin.info["phononic_support_segments"] = list(PHONONIC_SUPPORT_SEGMENTS)
    fin.info["phononic_support_gap"] = PHONONIC_SUPPORT_GAPS
    fin.info["phononic_support_neck_width"] = REVIEW_TETHER_CONNECTION_WIDTH
    fin.info["layer13_support_rectangles_only"] = True
    fin.info["mla_channel_angle"] = MLA_CHANNEL_ANGLE
    fin.info["layer14_channel_angle"] = LAYER14_CHANNEL_ANGLE
    fin.info["mla_channel_start_width"] = MLA_CHANNEL_START_WIDTH
    fin.info["layer14_channel_start_width"] = (
        LAYER14_CHANNEL_START_WIDTH
    )
    fin.info["layer14_channel_start_offset"] = (
        LAYER14_CHANNEL_START_OFFSET
    )
    fin.info["mla_channel_start_offset"] = MLA_CHANNEL_START_OFFSET
    fin.info["positive_channel_length"] = POSITIVE_CHANNEL_LENGTH
    fin.info["tether_rectangle_slot"] = False
    return fin


if __name__ == "__main__":
    GDS.mkdir(parents=True, exist_ok=True)
    fin = build_layout()
    gds_save_options = gf.kdb.SaveLayoutOptions()
    gds_save_options.gds2_max_vertex_count = 2000
    fin.write_gds(
        os.path.join(GDS, f"{fin.name}.gds"),
        save_options=gds_save_options,
        with_metadata=False,
    )
