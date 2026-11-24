"""Straight-coupler taper sweep, adiabatic vs linear tapers, derived from
gds/2611_straight_coupler_L15000_output_comparison.gds (+ .json).

Changes relative to that source layout
--------------------------------------
1. Device blocks (BLOCKS):
   - 500 um straights with tapers at both ends ("2T"), W = 1.60 um,
     L10 first-stage length T = 40/50/60/70/80 um; each row holds A, L and a
     copy of each.
   - 15000 um straights with an input taper + 2000 um plain output ("1T"),
     T = 60 um, W = 1.590/1.595/1.600/1.605/1.610 um (0.5 nm centreline offset
     where W/2 is a half nm, so the width is exact on the 1 nm grid).
2. Both taper profiles everywhere: "A" uses the source's adiabatic quintic smootherstep
   profile (zero width slope at both ends of every stage), "L" straight-sided
   linear stages. Same widths and lengths for both:
   - Layer 10: 15 um -> 0.8 um over T, then 0.8 um -> 0.200 um over 58 um.
   - Layer 15: W -> 0.200 um over T, ending at the L10 0.8 um point.
   - Layer-15 alignment tips at the 0.8 um point.
4. Device labels are drawn on Layer 15 (source: Layer 2).
5. Layer 10 now extends 2.0 um past the start of the Layer-2 window
   (source: 0.5 um, the JSON's layer2_start_offset_um).

Everything outward of the Layer-10 neck (butterfly tether, Layer-2/14 channels,
glue circles, Layer-13 pad) is copied unchanged from the source device and
shifted outward by the change in Layer-10 length (T + 58 - 80 um per end).

Only needs the bundled pure-Python gdsio_min.py (no klayout / gdsfactory).
Run from the Cading repository root:  python 2611_straight_coupler_taper_sweep.py
"""
from __future__ import annotations

# %% [00] Imports and paths
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from gdsio_min import Library, Cell, clip_halfplane, bbox  # noqa: E402

SOURCE_STEM = '2611_straight_coupler_L15000_output_comparison'
SOURCE_GDS = ROOT / 'gds' / (SOURCE_STEM + '.gds')
SOURCE_DEVICE = 'D1_SINGLE_L15000_W1620_2T'
SOURCE_LABEL_CELLS = ('LABEL_D1_SINGLE_L15000_W1620_2T',
                      'LABEL_D1_SINGLE_L15000_W1640_1T_OUT2000',
                      'LABEL_D1_SINGLE_L15000_W1680_2T')
STEM = '2611_straight_coupler_A-L_L300-2T-T40-80_L15000-1T-W1590-1610'
OUTPUT = ROOT / 'gds' / (STEM + '.gds')

# %% [01] USER PARAMETERS (um)
WG_WIDTH = 1.60                               # Layer 15 waveguide width (all devices)
L10_STAGE1_LENGTHS = (40.0, 50.0, 60.0, 70.0, 80.0)   # 15 -> 0.8 um stage; bottom to top in each block
PROFILES = {'A': 'smootherstep', 'L': 'linear'}   # label tag -> taper profile
# Each block gives one row per taper length T; members are placed left to right.
LONG_TAPER_T = 60.0                           # fixed L10 first-stage length for the long devices
LONG_WIDTHS = (1.590, 1.595, 1.600, 1.605, 1.610)   # long-device waveguide width sweep
# 'rows' lists (T, width) per row; 'label' says which parameter the label shows.
S_BEND_RADIUS = 150.0                         # short devices: the straight becomes a gentle S-bend:
S_BEND_ANGLE_DEG = 10.0                       # arc (R, this angle) - tilted straight - arc back
S_BEND_OFFSET = None                          # output-end y offset; None = follows from the angle
SHORT_COLUMN_PITCH = 4930.0                   # = 17 x 290 um field pitch: every column sits identically in the field grid
SHORT_LENGTH = 300.0                          # short devices: port-to-port x span (S-bend inside); with the
                                              # tapers every short device fits two 300 um fields (590 um)
BLOCKS = (dict(mode='2T', align='left', members=((SHORT_LENGTH, 'A'), (SHORT_LENGTH, 'L'),
                                                  (SHORT_LENGTH, 'A'), (SHORT_LENGTH, 'L')),
               rows=None, label='T', s_bend=True, max_fields=2, column_pitch=SHORT_COLUMN_PITCH),                          # rows: T sweep at WG_WIDTH
          dict(mode='1T', align='right', members=((15000.0, 'A'),),
               rows=tuple((LONG_TAPER_T, w) for w in LONG_WIDTHS), label='W', s_bend=False),
          dict(mode='1T', align='right', members=((15000.0, 'L'),),
               rows=tuple((LONG_TAPER_T, w) for w in LONG_WIDTHS), label='W', s_bend=False),
          dict(mode='1T', align='right', members=((15000.0, 'A'),), copy=2,        # duplicate long set
               rows=tuple((LONG_TAPER_T, w) for w in LONG_WIDTHS), label='W', s_bend=False),
          dict(mode='1T', align='right', members=((15000.0, 'L'),), copy=2,
               rows=tuple((LONG_TAPER_T, w) for w in LONG_WIDTHS), label='W', s_bend=False))
STRAIGHT_LENGTHS = tuple(sorted({L for b in BLOCKS for L, _ in b['members']}))
OUTPUT_EXTENSION_LENGTH = 2000.0              # 1T plain output beyond the output port
PORT_TO_L10_ENTRANCE = 45.0                   # plain Layer-15 straight between port and taper start

TAPER_PROFILE = 'smootherstep'                # set per device from PROFILES
TAPER_SAMPLE_STEP = 0.25                      # vertex spacing along curved tapers

L15_TIP_WIDTH = 0.200
L10_START_WIDTH = 15.0
L10_MID_WIDTH = 0.8                           # two stages: start -> 0.8 um -> tip
L10_STAGE2_LENGTH = 58.0                      # 0.8 -> 0.2 um
L10_MID_OFFSET = 0.0                          # L15 0.2 um end sits at the L10 0.8 um point (and the tips)
L10_TIP_WIDTH = 0.200
L14_LINEAR_START_BEFORE_MID = 3.153           # straight-sided L14 bridge starts this far before the 0.8 um point

# Per-device taper geometry, set by set_taper(T) before each device is built
L10_MID_POSITION = L10_TAPER_LENGTH = L15_TAPER_LENGTH = TICK_POSITION = L14_LINEAR_START = None


def set_taper(stage1_length):
    global L10_MID_POSITION, L10_TAPER_LENGTH, L15_TAPER_LENGTH, TICK_POSITION, L14_LINEAR_START
    L10_MID_POSITION = stage1_length
    L10_TAPER_LENGTH = stage1_length + L10_STAGE2_LENGTH
    L15_TAPER_LENGTH = stage1_length - L10_MID_OFFSET
    TICK_POSITION = L10_MID_POSITION          # alignment tips mark the L10 0.8 um point
    L14_LINEAR_START = L10_MID_POSITION - L14_LINEAR_START_BEFORE_MID

L2_MARGIN = 3.0                               # per side around Layer 10 (source: 2.5)
L2_START_OFFSET = 2.0                         # Layer 2 starts this far OUTWARD of the L10 entrance
L14_MARGIN = 1.0                              # per side around Layer 10
L14_START_OFFSET = 0.0

TICK_WIDTH = 0.5                              # Layer-15 alignment tips (ticks)
TICK_LENGTH = 10.0
TICK_GAP_TO_L10 = 2.0                         # same gap to the slab edge as the source (2.5 - 0.5)

BEAM_LAYER = (11, 0)
BEAM_WIDTH = 0.5
BEAM_OFFSETS_FROM_L2_START = (-0.5, 0.5, 1.5)  # same relation to the L2 window as the source
BEAM_SLAB_OVERLAP = 1.5                       # L11 beam reaches this far inside the L10 edge
BEAM_SPAN = 5.096 + 0.5                       # source anchor span + outer overlap

LABEL_LAYER = (15, 0)
LABEL_TIP_X_OFFSET = 2.0                      # from the alignment tips, toward the device
LABEL_Y_OFFSET = 15.0                         # source used 10; raised to clear the wider L11 beams (to 13.1 um)
LABEL_SOURCE_HEIGHT = 5.0                     # character height of the label font as copied from the source
LABEL_SOURCE_STROKE = 1.0                     # its stroke width
LABEL_HEIGHT = 2.0                            # drawn character height; outlines scale by LABEL_HEIGHT / 5
LABEL_LINE_WIDTH = 0.4                        # drawn stroke width (thinned only if scaling leaves it wider)
LABEL_GAP = 0.2                               # gap between characters (source 0.5 um at 5 um height)
LABEL_RASTER = 0.02                           # grid used to thin the strokes (20 nm)

SMALL_MARKER_SQUARE = 25.0                    # the smaller corner cross is replaced by this square (um side)
SMALL_MARKER_LAYER = (15, 0)
MARKER_X_INSET = 1000.0                       # every corner marker set moves this far toward the centre in x
ROW_PITCH = 300.0                             # device row spacing in y (glue circles 50 um apart)
GLUE_CIRCLE_DIAMETER = 250.0                  # the source's 300 um Layer-2/14 glue circles are redrawn at this size
SRC_GLUE_CIRCLE_R = 150.0
SRC_GLUE_CIRCLE_X = (8075.0, 8575.0)          # |x| of the source circle centres (on the waveguide axis)
FIRST_ROW_Y = -700.0
DEVICE_GAP = 500.0                            # minimum edge-to-edge gap between devices in one row
DEVICE_GAP_EBL_RANGE = 290.0                  # the 2nd device of a shared row may move right by up to this
                                              # (per row) so its tapers avoid EBL stitches; rows stay aligned
ROW_LEFT_EDGE_X = -6937.847                   # left edge of every left-aligned row (see ALIGN_BLOCK_STARTS)
ALIGN_BLOCK_STARTS = True                     # left-aligned blocks start exactly where the long devices start
ROW_RIGHT_EDGE_X = 11742.153                  # right edge (plain output end) of every right-aligned row
FRAME_X = (-7400.0, 12300.0)                  # same chip width as the source
FRAME_BOTTOM_MARGIN = 727.75                  # puts the corner markers at the same field phase as the rows
FRAME_TOP_MARGIN = 732.25                     # (marker centre 150 um inside the frame edge)

# EBL write-field placement (Layer 15). Fields are EBL_FIELD_SIZE square on a
# pitch of EBL_FIELD_SIZE - EBL_FIELD_OVERLAP, with a field corner at
# EBL_GRID_ORIGIN in GDS coordinates. The whole mask is shifted by
# (0..pitch) in x and y so Layer 15 needs the fewest fields; ties go to the
# placement whose field-overlap strips stay clear of tapers, tips and labels.
EBL_FIELD_SIZE = 300.0
EBL_FIELD_OVERLAP = 10.0
EBL_GRID_ORIGIN = (0.0, 0.0)
EBL_LAYER = (15, 0)
EBL_SEARCH_STEP = 0.5
EBL_CRITICAL_MARGIN = 2.0                     # keep overlap strips this far from critical zones
EBL_PIECE_STEP = 10.0                         # L15 shapes are sliced into pieces this long in x
EBL_Y_SEARCH_STEP = 1.0                       # y stitching is irrelevant: y only sets which fields are needed
EBL_Y_CANDIDATES = 6                          # y offsets (best by estimated count) searched fully in x
EBL_PRIORITY = 'no_taper_stitch'                       # 'fields': fewest fields first; 'no_taper_stitch': no stitch
                                              # through any L15 taper/alignment tip first, then fewest fields

# source geometry facts (measured from the source GDS)
SRC_HALF_L = 7500.0
SRC_L10_LENGTH = 80.0                          # source L10 entrance -> 0.2 um neck
DBU_PER_UM = 1000.0


# %% [02] Helpers
def um(v):
    return v * DBU_PER_UM


def profile(t):
    t = min(max(t, 0.0), 1.0)
    if TAPER_PROFILE == 'linear':
        return t
    return t * t * t * (t * (6 * t - 15) + 10)       # quintic smootherstep


def l10_half(s):
    s = min(max(s, 0.0), L10_TAPER_LENGTH)
    if s <= L10_MID_POSITION:
        return 0.5 * (L10_START_WIDTH + (L10_MID_WIDTH - L10_START_WIDTH) * profile(s / L10_MID_POSITION))
    t = (s - L10_MID_POSITION) / (L10_TAPER_LENGTH - L10_MID_POSITION)
    return 0.5 * (L10_MID_WIDTH + (L10_TIP_WIDTH - L10_MID_WIDTH) * profile(t))


def stations(s0, s1, breaks):
    """Edge vertices from s0 to s1: stage breaks always, plus dense samples for curved profiles."""
    pts = {s0, s1} | {b for b in breaks if s0 < b < s1}
    if TAPER_PROFILE != 'linear':
        n = max(1, int(math.ceil((s1 - s0) / TAPER_SAMPLE_STEP)))
        pts |= {s0 + (s1 - s0) * k / n for k in range(n + 1)}
    return sorted(pts)


def l10_stations(s0, s1):
    return stations(s0, s1, (0.0, L10_MID_POSITION, L10_TAPER_LENGTH))


def l15_half(s, w):
    s = min(max(s, 0.0), L15_TAPER_LENGTH)
    return 0.5 * (w + (L15_TIP_WIDTH - w) * profile(s / L15_TAPER_LENGTH))


def cut_far_piece(poly, side):
    """Clip a source polygon to the region outward of the source L10 neck.
    Returns ring (nm) ordered from the cut-top point, around the far side, to the cut-bottom point."""
    xc = um(-(SRC_HALF_L + 45.0 + SRC_L10_LENGTH)) if side < 0 else um(SRC_HALF_L + 45.0 + SRC_L10_LENGTH)
    ring = clip_halfplane(poly, xc, keep_left=(side < 0))
    on = [i for i, p in enumerate(ring) if abs(p[0] - xc) < 0.5]
    assert len(on) == 2, ('expected one crossing interval', len(on))
    i_hi, i_lo = sorted(on, key=lambda i: -ring[i][1])
    n = len(ring)
    # walk from i_hi away from i_lo
    step = 1 if (i_hi + 1) % n != i_lo else -1
    seq = [ring[i_hi]]
    i = i_hi
    while i != i_lo:
        i = (i + step) % n
        seq.append(ring[i])
    return xc, seq


def y_extent(poly, x):
    ys = []
    n = len(poly)
    for i in range(n):
        (x0, y0), (x1, y1) = poly[i], poly[(i + 1) % n]
        if (x0 < x) != (x1 < x):
            ys.append(y0 + (x - x0) / (x1 - x0) * (y1 - y0))
    return (min(ys), max(ys)) if ys else None


def mapper(side, half_L):
    """s (um, outward from the L10 entrance) -> x (nm) for this device end."""
    xe = side * (half_L + PORT_TO_L10_ENTRANCE)
    return lambda s: um(xe + side * s)


def shifted(seq, side, half_L):
    dx = um(side * (half_L - SRC_HALF_L + (L10_TAPER_LENGTH - SRC_L10_LENGTH)))
    return [(x + dx, y) for x, y in seq]


def merged_end_polygon(near_top, near_bot, far_seq):
    """near_top/near_bot: [(x,y)] from inner start to the cut; far_seq from cut-top to cut-bottom."""
    pts = near_top[:-1] + far_seq + list(reversed(near_bot[:-1]))
    return pts


def s_bend_offset(span, radius, angle_deg):
    th = math.radians(angle_deg)
    ls = (span - 2 * radius * math.sin(th)) / math.cos(th)
    return 2 * radius * (1 - math.cos(th)) + ls * math.sin(th)


def s_offset(L):
    return S_BEND_OFFSET if S_BEND_OFFSET is not None else s_bend_offset(L, S_BEND_RADIUS, S_BEND_ANGLE_DEG)


def s_bend(span, radius, offset):
    """Arc (radius, theta) - tilted straight - arc back over x-span `span` with y-offset `offset`.
    Returns theta (rad), tilted straight length, path length and a centreline [(x, y, heading)]."""
    def dy(th):
        ls = (span - 2 * radius * math.sin(th)) / math.cos(th)
        return 2 * radius * (1 - math.cos(th)) + ls * math.sin(th)
    lo, hi = 0.0, min(math.pi / 2 - 1e-9, math.asin(min(1.0, span / (2 * radius))))
    assert dy(hi) >= offset, 'offset too large for this span and radius'
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if dy(mid) < offset else (lo, mid)
    th = (lo + hi) / 2
    ls = (span - 2 * radius * math.sin(th)) / math.cos(th)
    x0 = -span / 2
    n_arc = max(8, int(math.ceil(radius * th / 0.5)))           # <= 0.5 um chords (0.2 nm sagitta)
    pts = [(x0 + radius * math.sin(th * k / n_arc), radius * (1 - math.cos(th * k / n_arc)), th * k / n_arc)
           for k in range(n_arc + 1)]
    x2, y2 = span / 2, offset
    pts += [(x2 - radius * math.sin(th * k / n_arc), offset - radius * (1 - math.cos(th * k / n_arc)), th * k / n_arc)
            for k in range(n_arc, -1, -1)]
    return th, ls, 2 * radius * th + ls, pts


def sample(s0, s1, n=2):
    return [s0 + (s1 - s0) * k / (n - 1) for k in range(n)]


# %% [03] Source data
def load_source():
    lib = Library.read(SOURCE_GDS)
    dev = lib.cells[SOURCE_DEVICE]
    far = {}   # (layer, side) -> list of (xc, seq) or plain polygons
    for l, d, p in dev.polygons:
        b = bbox(p)
        for side in (-1, 1):
            xc = um(side * (SRC_HALF_L + 45.0 + SRC_L10_LENGTH))
            fully = (b[2] <= xc) if side < 0 else (b[0] >= xc)
            crosses = (b[0] < xc < b[2])
            if l in (2, 10, 14) and crosses:
                far.setdefault((l, side), []).append(('cut',) + cut_far_piece(p, side))
            elif fully:
                far.setdefault((l, side), []).append(('whole', None, p))
    if GLUE_CIRCLE_DIAMETER != 2 * SRC_GLUE_CIRCLE_R:
        resize_glue_circles(far)
    glyphs = extract_glyphs(lib)
    frame = lib.cells['FRAME_AND_MARKERS']
    return lib, far, glyphs, frame


def resize_glue_circles(far):
    """Remove the source circle arcs from the Layer-2/14 channel outlines (each arc is replaced by
    the straight chord between its two channel junctions) and add separate circles of
    GLUE_CIRCLE_DIAMETER at the same centres."""
    r_old = um(SRC_GLUE_CIRCLE_R)
    for (layer, side), items in far.items():
        if layer not in (2, 14):
            continue
        centres = [(um(side * cx), 0) for cx in SRC_GLUE_CIRCLE_X]

        def on_circle(q):
            return any(abs(math.hypot(q[0] - c[0], q[1] - c[1]) - r_old) < 50 for c in centres)

        def strip(pts, ring):
            flags = [on_circle(q) for q in pts]
            n = len(pts)
            keep = []
            for k in range(n):
                if not flags[k]:
                    keep.append(pts[k])
                    continue
                prev_on = flags[k - 1] if (ring or k > 0) else False
                next_on = flags[(k + 1) % n] if (ring or k < n - 1) else False
                if not (prev_on and next_on):          # run ends (channel junctions) stay
                    keep.append(pts[k])
            return keep
        found = 0
        for n_item, (kind, xc, poly) in enumerate(items):
            before = len(poly)
            new = strip(poly, ring=(kind == 'whole'))
            found += before - len(new)
            items[n_item] = (kind, xc, new)
        assert found > 100, 'glue circle arcs not found on layer %d' % layer
        r = um(GLUE_CIRCLE_DIAMETER / 2)
        for c in centres:
            items.append(('whole', None, [(round(c[0] + r * math.cos(2 * math.pi * k / 720)),
                                           round(c[1] + r * math.sin(2 * math.pi * k / 720))) for k in range(720)]))


def extract_glyphs(lib):
    texts = {'LABEL_D1_SINGLE_L15000_W1620_2T': 'L15000_W1.62_2T',
             'LABEL_D1_SINGLE_L15000_W1640_1T_OUT2000': 'L15000_W1.64_1T',
             'LABEL_D1_SINGLE_L15000_W1680_2T': 'L15000_W1.68_2T'}
    glyphs = {}
    for cn, txt in texts.items():
        polys = sorted(lib.cells[cn].polygons, key=lambda t: bbox(t[2])[0])
        clusters = []
        for l, d, p in polys:
            b = bbox(p)
            if clusters and b[0] < clusters[-1]['x1'] + 100:   # < 0.1 um gap => same glyph
                clusters[-1]['polys'].append(p)
                clusters[-1]['x1'] = max(clusters[-1]['x1'], b[2])
            else:
                clusters.append({'x0': b[0], 'x1': b[2], 'polys': [p]})
        assert len(clusters) == len(txt), (cn, len(clusters), len(txt))
        for ch, cl in zip(txt, clusters):
            glyphs.setdefault(ch, ([[(x - cl['x0'], y) for x, y in p] for p in cl['polys']], cl['x1'] - cl['x0']))
    # '9' is the same font's '6' rotated by 180 degrees about the glyph box centre
    six, w6 = glyphs['6']
    ys = [y for p in six for _, y in p]
    h = max(ys) + min(ys)
    glyphs['9'] = ([[(w6 - x, h - y) for x, y in p] for p in six], w6)
    # 'A' does not occur in the source labels: stencil 'A' in the same style
    if 'A' not in glyphs:
        k = DBU_PER_UM
        glyphs['A'] = ([[(0, 0), (0, 4.5 * k), (0.5 * k, 5 * k), (3 * k, 5 * k), (3.5 * k, 4.5 * k),
                         (3.5 * k, 0), (2.5 * k, 0), (2.5 * k, 4 * k), (1 * k, 4 * k), (1 * k, 0)],
                        [(1.5 * k, 2 * k), (2.5 * k, 2 * k), (2.5 * k, 3 * k), (1.5 * k, 3 * k)]], 3.5 * k)
    # '3' does not occur in the source labels: chamfered like '2', 1 um strokes, 3.5 x 5 um
    if '3' not in glyphs:
        k = DBU_PER_UM
        glyphs['3'] = ([[(0, 0), (3 * k, 0), (3.5 * k, 0.5 * k), (3.5 * k, 4.5 * k), (3 * k, 5 * k), (0, 5 * k),
                         (0, 4 * k), (2.5 * k, 4 * k), (2.5 * k, 1 * k), (0, 1 * k)],
                        [(1 * k, 2 * k), (2.5 * k, 2 * k), (2.5 * k, 3 * k), (1 * k, 3 * k)]], 3.5 * k)
    # '7' does not occur in the source labels: drawn in the same style (1 um strokes, 3.5 x 5 um)
    if '7' not in glyphs:
        k = DBU_PER_UM
        glyphs['7'] = ([[(0, 4 * k), (0, 5 * k), (3.5 * k, 5 * k), (3.5 * k, 0),
                         (2.5 * k, 0), (2.5 * k, 4 * k)]], 3.5 * k)
    k = LABEL_HEIGHT / LABEL_SOURCE_HEIGHT
    glyphs = {ch: ([[(round(x * k), round(y * k)) for x, y in p] for p in polys], w * k)
              for ch, (polys, w) in glyphs.items()}
    excess = (LABEL_SOURCE_STROKE * k - LABEL_LINE_WIDTH) / 2
    assert excess > -1e-9, 'scaled strokes are already thinner than LABEL_LINE_WIDTH'
    if excess > 1e-9:
        glyphs = {ch: (thin_glyph(polys, excess), w) for ch, (polys, w) in glyphs.items()}
    return glyphs


def thin_glyph(polys, d):
    """Union of a glyph's polygons (nm), shrunk by d um on every side (Euclidean erosion on a
    LABEL_RASTER grid), traced back to simple polygons. Strokes 1 um -> 1 - 2d um wide."""
    import numpy as np
    from scipy import ndimage
    px = LABEL_RASTER
    xs = [x for p in polys for x, _ in p]
    ys = [y for p in polys for _, y in p]
    x0, y0 = min(xs) / DBU_PER_UM - px, min(ys) / DBU_PER_UM - px
    nx = int(math.ceil((max(xs) / DBU_PER_UM - x0) / px)) + 2
    ny = int(math.ceil((max(ys) / DBU_PER_UM - y0) / px)) + 2
    cx = x0 + (np.arange(nx) + 0.5) * px
    cy = y0 + (np.arange(ny) + 0.5) * px
    X, Y = np.meshgrid(cx, cy)
    img = np.zeros((ny, nx), bool)
    for p in polys:                                   # even-odd point-in-polygon per pixel centre
        pts = [(a / DBU_PER_UM, b / DBU_PER_UM) for a, b in p]
        inside = np.zeros_like(img)
        for (ax, ay), (bx, by) in zip(pts, pts[1:] + pts[:1]):
            if ay == by:
                continue
            cond = (ay > Y) != (by > Y)
            xint = ax + (Y - ay) * (bx - ax) / (by - ay)
            inside ^= cond & (X < xint)
        img |= inside
    dist = ndimage.distance_transform_edt(img) * px
    out_img = dist > d + px / 2
    # trace pixel-edge boundaries (interior on the left)
    edges = {}
    filled = np.pad(out_img, 1)
    for i, j in zip(*np.nonzero(filled)):
        # pixel (row i, col j) spans lattice corners (j, i)..(j+1, i+1); edges run CCW (interior on the left)
        if not filled[i - 1, j]:
            edges.setdefault((j, i), []).append((j + 1, i))              # bottom, +x
        if not filled[i + 1, j]:
            edges.setdefault((j + 1, i + 1), []).append((j, i + 1))      # top, -x
        if not filled[i, j - 1]:
            edges.setdefault((j, i + 1), []).append((j, i))              # left, -y
        if not filled[i, j + 1]:
            edges.setdefault((j + 1, i), []).append((j + 1, i + 1))      # right, +y
    loops = []
    while edges:
        start = next(iter(edges))
        loop, cur = [start], start
        prev_dir = None
        while True:
            nxts = edges[cur]
            if len(nxts) == 1:
                nxt = nxts.pop()
            else:                                            # pinch point: take the right-most turn,
                def turn(n):                                 # so diagonal neighbours stay separate
                    d1 = prev_dir
                    d2 = (n[0] - cur[0], n[1] - cur[1])
                    return d1[0] * d2[1] - d1[1] * d2[0]
                nxt = min(nxts, key=turn)
                nxts.remove(nxt)
            if not nxts:
                del edges[cur]
            prev_dir = (nxt[0] - cur[0], nxt[1] - cur[1])
            if nxt == start:
                break
            loop.append(nxt)
            cur = nxt
        loops.append(loop)
    result = []
    for loop in loops:
        pts = [(x0 + (c - 1) * px, y0 + (r - 1) * px) for c, r in loop]   # undo the 1-pixel pad
        simp = []
        n = len(pts)
        for k in range(n):                                   # drop collinear vertices
            a, b, c = pts[k - 1], pts[k], pts[(k + 1) % n]
            if abs((b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])) > 1e-12:
                simp.append(b)
        area = sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(simp, simp[1:] + simp[:1])) / 2
        assert area > 0, 'unexpected hole in a thinned glyph'
        result.append([(round(x * DBU_PER_UM), round(y * DBU_PER_UM)) for x, y in simp])
    return result


def label_polygons(text, glyphs, x0, y0):
    out, x = [], um(x0)
    for ch in text:
        polys, w = glyphs[ch]
        out += [[(x + px, um(y0) + py) for px, py in p] for p in polys]
        x += w + um(LABEL_GAP)
    return out, (x - um(LABEL_GAP) - um(x0)) / DBU_PER_UM


# %% [04] Device construction
def build_end(cell, side, half_L, width, far):
    X = mapper(side, half_L)
    S = L10_TAPER_LENGTH
    report = {}

    # Layer 10: linear taper merged with the copied neck + butterfly tether
    (_, xc, seq), = far[(10, side)]
    seq = shifted(seq, side, half_L)
    top = [(X(s), um(l10_half(s))) for s in l10_stations(0, S)]
    bot = [(X(s), -um(l10_half(s))) for s in l10_stations(0, S)]
    top[-1], bot[-1] = seq[0], seq[-1]
    cell.polygons.append((10, 0, merged_end_polygon(top, bot, seq)))

    # Layer 2: window = L10 + margin, starting L2_START_OFFSET outward of the L10 entrance
    s2 = L2_START_OFFSET
    (_, xc, seq), = [f for f in far[(2, side)] if f[0] == 'cut']
    seq = shifted(seq, side, half_L)
    top = [(X(s), um(l10_half(s) + L2_MARGIN)) for s in l10_stations(s2, S)]
    bot = [(X(s), -um(l10_half(s) + L2_MARGIN)) for s in l10_stations(s2, S)]
    # keep the margin along the copied 0.2 um neck until the copied L2 is wider on both sides
    h_neck = um(l10_half(S) + L2_MARGIN)
    s_join = None
    for k in range(1, 2001):
        e = y_extent(seq, X(S + 0.1 * k))
        if e and e[1] >= h_neck and e[0] <= -h_neck:
            s_join = 0.1 * k + 0.5
            break
    assert s_join is not None, 'copied Layer-2 channel never exceeds the new margin'
    top.append((X(S + s_join), h_neck))
    bot.append((X(S + s_join), -h_neck))
    cell.polygons.append((2, 0, top + list(reversed(bot))))
    cell.polygons.append((2, 0, seq))                     # copied channel (overlaps the extension)
    report['layer2_neck_extension_um'] = s_join
    report['layer2_start_width_um'] = 2 * (l10_half(s2) + L2_MARGIN)

    # Layer 14: L10 + margin up to L14_LINEAR_START, then straight-sided to the copied channel
    s14, sj = L14_START_OFFSET, L14_LINEAR_START
    (_, xc, seq), = [f for f in far[(14, side)] if f[0] == 'cut']
    seq = shifted(seq, side, half_L)
    top = [(X(s), um(l10_half(s) + L14_MARGIN)) for s in l10_stations(s14, sj)] + [seq[0]]
    bot = [(X(s), -um(l10_half(s) + L14_MARGIN)) for s in l10_stations(s14, sj)] + [seq[-1]]
    cell.polygons.append((14, 0, top[:-1] + seq + list(reversed(bot[:-1]))))
    report['layer14_start_width_um'] = 2 * (l10_half(s14) + L14_MARGIN)
    # safety: straight-sided L14 bridge must still cover L10 + margin
    for k in range(51):
        s = sj + (S - sj) * k / 50
        t = k / 50
        yt = (1 - t) * um(l10_half(sj) + L14_MARGIN) + t * seq[0][1]
        yb = (1 - t) * -um(l10_half(sj) + L14_MARGIN) + t * seq[-1][1]
        assert yt >= um(l10_half(s) + L14_MARGIN) - 2 and yb <= -um(l10_half(s) + L14_MARGIN) + 2

    # everything fully outward of the neck (glue circles, L13 pad, other L2/L14 pieces)
    for (l, sd), items in far.items():
        if sd != side:
            continue
        for kind, _, p in items:
            if kind == 'whole':
                cell.polygons.append((l, 0, shifted(p, side, half_L)))

    # Layer 11: local support beams at the L2 window entrance
    s2 = L2_START_OFFSET
    for off in BEAM_OFFSETS_FROM_L2_START:
        sc = s2 + off
        yin = l10_half(sc) - BEAM_SLAB_OVERLAP
        yout = yin + BEAM_SLAB_OVERLAP + BEAM_SPAN
        xa, xb = sorted((X(sc - BEAM_WIDTH / 2), X(sc + BEAM_WIDTH / 2)))
        for sg in (1, -1):
            ya, yb = sorted((um(sg * yin), um(sg * yout)))
            cell.polygons.append((BEAM_LAYER[0], BEAM_LAYER[1], [(xa, ya), (xb, ya), (xb, yb), (xa, yb)]))

    # Layer 15: alignment tips at the L10 1 um point
    st = TICK_POSITION
    yin = l10_half(st) + TICK_GAP_TO_L10
    xa, xb = sorted((X(st - TICK_WIDTH / 2), X(st + TICK_WIDTH / 2)))
    for sg in (1, -1):
        ya, yb = sorted((um(sg * yin), um(sg * (yin + TICK_LENGTH))))
        cell.polygons.append((15, 0, [(xa, ya), (xb, ya), (xb, yb), (xa, yb)]))
    report['tick_inner_y_um'] = yin
    return report


def build_device(name, L, width, mode, far, s_offset=0.0):
    if s_offset:
        return build_s_device(name, L, width, mode, far, s_offset)
    cell = Cell(name)
    h = L / 2
    xe = h + PORT_TO_L10_ENTRANCE
    xt = xe + L15_TAPER_LENGTH
    ss = stations(0.0, L15_TAPER_LENGTH, ())                       # outward from the entrance

    def lo(s):   # bottom/top edges on the 1 nm grid with the exact total width
        wn = int(round(2 * um(l15_half(s, width))))
        return -(wn // 2)

    def hi(s):
        wn = int(round(2 * um(l15_half(s, width))))
        return wn - wn // 2
    left = [(um(-(xe + s)), lo(s)) for s in reversed(ss)]     # tip -> entrance, bottom
    left_top = [(um(-(xe + s)), hi(s)) for s in ss]           # entrance -> tip, top
    if mode == '2T':
        right = ([(um(xe + s), lo(s)) for s in ss] +
                 [(um(xe + s), hi(s)) for s in reversed(ss)])
    else:
        xo = h + OUTPUT_EXTENSION_LENGTH
        right = [(um(xo), lo(0.0)), (um(xo), hi(0.0))]
    cell.polygons.append((15, 0, left + right + left_top))
    rep = build_end(cell, -1, h, width, far)
    if mode == '2T':
        build_end(cell, 1, h, width, far)
    return cell, rep


def build_s_device(name, L, width, mode, far, off):
    """Two-taper device whose straight is an S-bend: the output end assembly is raised by `off`."""
    assert mode == '2T'
    cell = Cell(name)
    h = L / 2
    xe = h + PORT_TO_L10_ENTRANCE
    th, ls, path, sb = s_bend(L, S_BEND_RADIUS, off)
    ss = stations(0.0, L15_TAPER_LENGTH, ())
    # centreline with half-widths: left taper tip -> entrance -> port, S, port -> entrance -> tip
    line = [(-(xe + s), 0.0, 0.0, l15_half(s, width)) for s in reversed(ss)]
    line += [(x, y, a, width / 2) for x, y, a in sb]
    line += [(xe + s, off, 0.0, l15_half(s, width)) for s in ss]
    top_e = [(um(x - hw * math.sin(a)), um(y + hw * math.cos(a))) for x, y, a, hw in line]
    bot_e = [(um(x + hw * math.sin(a)), um(y - hw * math.cos(a))) for x, y, a, hw in line]
    cell.polygons.append((15, 0, bot_e + list(reversed(top_e))))
    rep = build_end(cell, -1, h, width, far)
    n0 = len(cell.polygons)
    build_end(cell, 1, h, width, far)
    for k in range(n0, len(cell.polygons)):                     # output end assembly follows the S
        l, d, pts = cell.polygons[k]
        cell.polygons[k] = (l, d, [(x, y + um(off)) for x, y in pts])
    rep.update(s_bend_theta_deg=math.degrees(th), s_bend_tilted_straight_um=ls, s_bend_path_length_um=path)
    return cell, rep


# %% [05] Frame
def build_frame(src_frame, y0, y1):
    fx0, fx1 = FRAME_X
    cx, cy = (fx0 + fx1) / 2, (y0 + y1) / 2
    hw, hh = (fx1 - fx0) / 2, (y1 - y0) / 2
    src_hh = 2050.0
    c = Cell('FRAME_AND_MARKERS')
    c.polygons.append((0, 0, [(um(-hw), um(-hh)), (um(hw), um(-hh)), (um(hw), um(hh)), (um(-hw), um(hh))]))
    # group the source marker polygons into crosses around their centre squares
    centres = [tuple(sum(q[k] for q in p) / len(p) for k in (0, 1))
               for l, d, p in src_frame.polygons if l != 0 and len(p) == 8]
    groups = {cen: [] for cen in centres}
    for l, d, p in src_frame.polygons:
        if l == 0:
            continue
        mx, my = (sum(q[k] for q in p) / len(p) for k in (0, 1))
        cen = min(centres, key=lambda q: (q[0] - mx) ** 2 + (q[1] - my) ** 2)
        groups[cen].append((l, d, p))
    reach = {cen: max(max(abs(q[0] - cen[0]), abs(q[1] - cen[1])) for _, _, p in g for q in p)
             for cen, g in groups.items()}
    small = min(reach.values())
    n = 0
    for cen, g in groups.items():
        dy = um(hh - src_hh) * (1 if cen[1] > 0 else -1)
        mdx = -um(MARKER_X_INSET) * (1 if cen[0] > 0 else -1)
        if abs(reach[cen] - small) < 1 and SMALL_MARKER_SQUARE:
            h = um(SMALL_MARKER_SQUARE / 2)
            x0, y0_ = cen[0] + mdx, cen[1] + dy
            c.polygons.append((SMALL_MARKER_LAYER[0], SMALL_MARKER_LAYER[1],
                               [(x0 - h, y0_ - h), (x0 + h, y0_ - h), (x0 + h, y0_ + h), (x0 - h, y0_ + h)]))
            n += 1
            continue
        for l, d, p in g:
            c.polygons.append((l, d, [(x + mdx, y + dy) for x, y in p]))
            n += 1
    return c, (cx, cy), n


# %% [05b] EBL write-field placement
def lab_len_of(text, glyphs):
    return (sum(glyphs[c][1] for c in text) + um(LABEL_GAP) * (len(text) - 1)) / DBU_PER_UM


def _poly_y_range(poly, x0, x1):
    ys = [q[1] for q in poly if x0 <= q[0] <= x1]
    n = len(poly)
    for xc in (x0, x1):
        for i in range(n):
            (ax, ay), (bx, by) = poly[i], poly[(i + 1) % n]
            if (ax < xc) != (bx < xc):
                ys.append(ay + (xc - ax) / (bx - ax) * (by - ay))
    return (min(ys), max(ys)) if ys else None


def cell_pieces(cell, _cache={}):
    """Layer-15 shapes of a cell as small boxes (um): long shapes are sliced every EBL_PIECE_STEP in x."""
    if cell.name in _cache:
        return _cache[cell.name]
    out = []
    step = um(EBL_PIECE_STEP)
    for l, d, p in cell.polygons:
        if (l, d) != EBL_LAYER:
            continue
        b = bbox(p)
        if b[2] - b[0] <= step:
            out.append((b[0], b[2], b[1], b[3]))
            continue
        n = int(math.ceil((b[2] - b[0]) / step))
        for k in range(n):
            x0 = b[0] + (b[2] - b[0]) * k / n
            x1 = b[0] + (b[2] - b[0]) * (k + 1) / n
            yr = _poly_y_range(p, x0, x1)
            if yr:
                out.append((x0, x1, yr[0], yr[1]))
    out = [tuple(v / DBU_PER_UM for v in q) for q in out]
    _cache[cell.name] = out
    return out


def mask_pieces(lib, top):
    """All Layer-15 pieces of the (unshifted) mask: (x0, x1, y0, y1, owner)."""
    out = []
    for r in top.refs:
        ox, oy = r['xy'][0] / DBU_PER_UM, r['xy'][1] / DBU_PER_UM
        out += [(a + ox, b + ox, c + oy, d + oy, r['name']) for a, b, c, d in cell_pieces(lib.cells[r['name']])]
    return out


def piece_rows(y0, y1, py):
    """Field rows that write a piece: the row whose field holds it whole (preferring the row
    of its centre), otherwise every row it touches (a y stitch; harmless here)."""
    pitch = EBL_FIELD_SIZE - EBL_FIELD_OVERLAP
    jc = math.floor(((y0 + y1) / 2 - py) / pitch)
    for j in (jc, jc - 1, jc + 1):
        if py + pitch * j <= y0 and y1 <= py + pitch * j + EBL_FIELD_SIZE:
            return (j,)
    return tuple(range(math.floor((y0 - py) / pitch), math.floor((y1 - py) / pitch) + 1))


def row_intervals(pieces, py):
    rows = {}
    for x0, x1, y0, y1, *_ in pieces:
        for j in piece_rows(y0, y1, py):
            rows.setdefault(j, []).append((x0, x1))
    out = {}
    for j, xs in rows.items():
        xs.sort()
        merged = [list(xs[0])]
        for a, b in xs[1:]:
            if a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        out[j] = merged
    return out


def cover_x(intervals, px):
    """Minimal 1-D field cover (greedy is optimal for intervals on a fixed grid). Returns field indices."""
    pitch = EBL_FIELD_SIZE - EBL_FIELD_OVERLAP
    used, reach = [], -1e18
    for a, b in intervals:
        while reach < b:
            pnt = a if a > reach else reach
            i = math.floor((pnt - px) / pitch)
            used.append(i)
            reach = px + pitch * i + EBL_FIELD_SIZE
    return used


def used_fields(rowiv, px):
    return {(i, j) for j, iv in rowiv.items() for i in cover_x(iv, px)}


def zone_rows_for(critical, py):
    # a zone is written in the field row that holds the taper's centreline
    return [(piece_rows(yz - 1.0, yz + 1.0, py)[0], x0, x1, kind, case) for yz, x0, x1, kind, case in critical]


def stitch_hits(zone_rows, used, px, skip=frozenset()):
    """x stitches (overlap strip between two written fields of one row) through critical zones."""
    pitch = EBL_FIELD_SIZE - EBL_FIELD_OVERLAP
    hits = {'taper': 0, 'label': 0}
    clear = 1e9
    for j, x0, x1, kind, case in zone_rows:
        if case in skip:
            continue
        for i in range(math.floor((x0 - px) / pitch) - 1, math.floor((x1 - px) / pitch) + 2):
            if (i - 1, j) not in used or (i, j) not in used:
                continue                                     # no stitch: only one side is written
            e0 = px + pitch * i                              # stitch overlap strip [e0, e0 + overlap]
            gap = max(x0 - (e0 + EBL_FIELD_OVERLAP), e0 - x1)
            if gap < EBL_CRITICAL_MARGIN:
                hits[kind] += 1
            if kind == 'taper':
                clear = min(clear, gap)
    return hits, clear


RANK = {'fields': lambda r: (r['n'], r['taper'], r['label'], -r['clear']),
        'no_taper_stitch': lambda r: (r['taper'], r['n'], r['label'], -r['clear'])}


def one_field_row(pieces, py):
    """True if all pieces fit inside a single field row."""
    pitch = EBL_FIELD_SIZE - EBL_FIELD_OVERLAP
    y0 = min(q[2] for q in pieces)
    y1 = max(q[3] for q in pieces)
    j = math.floor((y0 - py) / pitch)
    return any(py + pitch * jj <= y0 and y1 <= py + pitch * jj + EBL_FIELD_SIZE for jj in (j - 1, j))


def place_for_ebl(lib, top, critical, skip=frozenset(), fixed=None, exclude=frozenset(), y_groups=(),
                  limits=()):
    """Shift of the whole mask (dx, dy in [0, pitch)) against the field grid. Pieces owned by
    `exclude` are ignored; each owner group in `y_groups` must fit inside one field row."""
    pitch = EBL_FIELD_SIZE - EBL_FIELD_OVERLAP
    gx, gy = EBL_GRID_ORIGIN
    allp = mask_pieces(lib, top)
    pieces = [q for q in allp if q[4] not in exclude]
    groups = [[q for q in allp if q[4] in g] for g in y_groups]
    lim = [([q[:4] for q in allp if q[4] in g], n) for g, n in limits]
    if fixed is not None:
        y_cands = [fixed[1]]
    else:   # rank y offsets by expected field count, then search x for the best few
        scored_y = []
        for k in range(int(round(pitch / EBL_Y_SEARCH_STEP))):
            dy = k * EBL_Y_SEARCH_STEP
            bad = sum(not one_field_row(g, gy - dy) for g in groups)
            rowiv = row_intervals(pieces, gy - dy)
            est = sum((b - a + pitch) / pitch for iv in rowiv.values() for a, b in iv)
            scored_y.append((bad, est, dy))
        scored_y.sort()
        y_cands = [dy for _, _, dy in scored_y[:EBL_Y_CANDIDATES]]
    results = []
    for dy in y_cands:
        py = gy - dy
        rowiv = row_intervals(pieces, py)
        zr = zone_rows_for(critical, py)
        xs = [fixed[0]] if fixed is not None else [k * EBL_SEARCH_STEP for k in range(int(round(pitch / EBL_SEARCH_STEP)))]
        lim_rows = [(row_intervals(ps, py), n) for ps, n in lim]
        for dx in xs:
            used = used_fields(rowiv, gx - dx)
            hits, clear = stitch_hits(zr, used, gx - dx, skip)
            over = sum(max(0, len(used_fields(riv, gx - dx)) - n) for riv, n in lim_rows)
            results.append(dict(over=over, n=len(used), taper=hits['taper'], label=hits['label'],
                                clear=clear, dx=dx, dy=dy))
    alternatives = {k: min(results, key=lambda r, f=f: (r['over'],) + f(r)) for k, f in RANK.items()}
    chosen = alternatives[EBL_PRIORITY]
    dx, dy = chosen['dx'], chosen['dy']
    fields = sorted(used_fields(row_intervals(allp, gy - dy), gx - dx))
    return dict(
        field_size_um=EBL_FIELD_SIZE, field_overlap_um=EBL_FIELD_OVERLAP, field_pitch_um=pitch,
        grid_origin_um=list(EBL_GRID_ORIGIN), layer=list(EBL_LAYER),
        mask_shift_um=[dx, dy], layer15_field_count=len(fields), priority=EBL_PRIORITY,
        devices_over_field_limit=chosen['over'],
        taper_zones_with_stitch=chosen['taper'], label_zones_with_stitch=chosen['label'],
        min_stitch_to_taper_zone_um=round(chosen['clear'], 3),
        alternatives={k: dict(fields=v['n'], taper_zones_with_stitch=v['taper'],
                              label_zones_with_stitch=v['label'], mask_shift_um=[v['dx'], v['dy']])
                      for k, v in alternatives.items()},
        fields_xy=[(gx + pitch * i, gy + pitch * j) for i, j in fields])


def device_fields(dev_pieces, px, py):
    return len(used_fields(row_intervals(dev_pieces, py), px))


def shift_devices(top, critical, devices, cases, delta):
    refs = {r['name']: r for r in top.refs}
    rec = {d['case_id']: d for d in devices}
    for c in cases:
        for r in (refs[c], refs['LABEL_' + c]):
            r['xy'] = (r['xy'][0] + um(delta), r['xy'][1])
        for n, z in enumerate(critical):
            if z[4] == c:
                critical[n] = (z[0], z[1] + delta, z[2] + delta, z[3], z[4])
        d = rec[c]
        d['origin_um'][0] = round(d['origin_um'][0] + delta, 3)
        d['label_origin_um'][0] = round(d['label_origin_um'][0] + delta, 3)
        d['global_x_extent_um'] = [round(v + delta, 3) for v in d['global_x_extent_um']]


def place_limited_block(lib, top, critical, devices, block_cases, lead_cases, limit, shift):
    """Slide a whole block of rows right by 0..pitch so that every lead device (first of its row)
    needs at most `limit` fields and has no x stitch through a taper."""
    pitch = EBL_FIELD_SIZE - EBL_FIELD_OVERLAP
    gx, gy = EBL_GRID_ORIGIN
    px, py = gx - shift[0], gy - shift[1]
    allp = mask_pieces(lib, top)
    own = {c: [q[:4] for q in allp if q[4] in (c, 'LABEL_' + c)] for c in lead_cases}
    zones = [z for z in critical if z[4] in lead_cases]
    best = None
    for k in range(int(round(pitch / EBL_SEARCH_STEP))):
        delta = k * EBL_SEARCH_STEP
        over, taper, clear = 0, 0, 1e9
        for c in lead_cases:
            ps = [(a + delta, b + delta, y0, y1) for a, b, y0, y1 in own[c]]
            used = used_fields(row_intervals(ps, py), px)
            over += max(0, len(used) - limit)
            zr = zone_rows_for([(z[0], z[1] + delta, z[2] + delta, z[3], z[4]) for z in zones if z[4] == c], py)
            h, cl = stitch_hits(zr, used, px)
            taper += h['taper']
            clear = min(clear, cl)
        key = (over, taper, -clear, delta)
        if best is None or key < best[0]:
            best = (key, delta)
    shift_devices(top, critical, devices, block_cases, best[1])
    return dict(block_shift_um=best[1], lead_devices_over_limit=best[0][0],
                lead_taper_stitches=best[0][1], lead_min_stitch_to_taper_um=round(-best[0][2], 3))


def adjust_movable_gaps(lib, top, critical, movable, devices, shift):
    """Slide each movable device (and the devices right of it in its row) right by
    0..DEVICE_GAP_EBL_RANGE so its taper zones avoid x stitches with the fewest fields."""
    dx, dy = shift
    gx, gy = EBL_GRID_ORIGIN
    px, py = gx - dx, gy - dy
    refs = {r['name']: r for r in top.refs}
    rec = {d['case_id']: d for d in devices}
    chosen = {}
    for case, info in movable.items():
        group = [c for c, inf in movable.items() if inf['row'] == info['row'] and inf['order'] >= info['order']]
        owners = set(group) | {'LABEL_' + c for c in group}
        pieces = [q for q in mask_pieces(lib, top) if abs((q[2] + q[3]) / 2 - info['row_y']) < 300.0]
        mine = [q for q in pieces if q[4] in owners]
        others = [q for q in pieces if q[4] not in owners]
        zones = [z for z in critical if z[4] == case]
        row_hi = max(refs[c]['xy'][0] / DBU_PER_UM + movable[c]['hi'] for c in group)
        best = None
        for k in range(int(round(DEVICE_GAP_EBL_RANGE / EBL_SEARCH_STEP))):
            delta = k * EBL_SEARCH_STEP
            if row_hi + delta > ROW_RIGHT_EDGE_X:
                break
            moved = [(a + delta, b + delta, c, d) for a, b, c, d, _ in mine]
            used = used_fields(row_intervals(others + moved, py), px)
            zr = zone_rows_for([(z[0], z[1] + delta, z[2] + delta, z[3], z[4]) for z in zones], py)
            hits, clear = stitch_hits(zr, used, px)
            over = 0
            if info.get('max_fields'):
                ps = [(a + delta, b + delta, c, d) for a, b, c, d, o in mine if o in (case, 'LABEL_' + case)]
                over = max(0, device_fields(ps, px, py) - info['max_fields'])
            key = (over, hits['taper'], len(used), hits['label'], -clear)
            if best is None or key < best[0]:
                best = (key, delta)
        delta = best[1]
        if delta:
            for c in group:
                for r in (refs[c], refs['LABEL_' + c]):
                    r['xy'] = (r['xy'][0] + um(delta), r['xy'][1])
                for n, z in enumerate(critical):
                    if z[4] == c:
                        critical[n] = (z[0], z[1] + delta, z[2] + delta, z[3], z[4])
                d = rec[c]
                d['origin_um'][0] = round(d['origin_um'][0] + delta, 3)
                d['label_origin_um'][0] = round(d['label_origin_um'][0] + delta, 3)
                d['global_x_extent_um'] = [round(v + delta, 3) for v in d['global_x_extent_um']]
        chosen[case] = round(DEVICE_GAP + delta, 3)
        rec[case]['gap_to_left_device_um'] = chosen[case]
    return chosen


def write_field_overlay(ebl):
    """Separate review file: used Layer-15 write fields as boxes on layer 99 (not for fabrication)."""
    lib = Library('LIB')
    c = Cell('EBL_L15_WRITE_FIELDS')
    s = EBL_FIELD_SIZE
    for x, y in ebl['fields_xy']:
        c.polygons.append((99, 0, [(um(x), um(y)), (um(x + s), um(y)), (um(x + s), um(y + s)), (um(x), um(y + s))]))
    lib.add(c)
    lib.write(OUTPUT.with_name(STEM + '_EBL_FIELDS_L15.gds'))


# %% [06] Assemble, write, report
def main():
    src, far, glyphs, src_frame = load_source()
    lib = Library('LIB')
    top = Cell('STRAIGHT_COUPLER_ADIABATIC_TAPER_SWEEP')
    devices = []
    row = 0
    critical = []                                  # (y, x0, x1, kind, case) Layer-15 taper/tip/label zones
    movable = {}                                   # 2nd device of left-aligned shared rows
    limited = {}                                   # blocks whose devices must fit `max_fields` fields
    global TAPER_PROFILE
    for block in BLOCKS:
        mode, align = block['mode'], block['align']
        if True:
          for T, w in (block['rows'] or tuple((t, WG_WIDTH) for t in L10_STAGE1_LENGTHS)):
              set_taper(T)
              key = f'T{int(T)}' if block['label'] == 'T' else f'W{w:.3f}'
              y = FIRST_ROW_Y + row * ROW_PITCH
              built = []
              seen = {}
              for L, tag in block['members']:
                  TAPER_PROFILE = PROFILES[tag]
                  seen[(L, tag)] = copy = seen.get((L, tag), 0) + block.get('copy', 1)
                  suffix = '' if copy == 1 else f'_{copy}'           # repeated members are copies
                  case = f'SC_L{int(L)}_{key.replace(".", "")}_{tag}_{mode}{suffix}'
                  off = s_offset(L) if block.get('s_bend') else 0.0
                  cell, rep = build_device(case, L, w, mode, far, off)
                  xs = [q[0] for _, _, pp in cell.polygons for q in pp]
                  built.append((L, case, cell, rep, min(xs) / DBU_PER_UM, max(xs) / DBU_PER_UM, tag, suffix))
              span = sum(b[5] - b[4] for b in built) + DEVICE_GAP * (len(built) - 1)
              x_left = ROW_LEFT_EDGE_X if align == 'left' else ROW_RIGHT_EDGE_X - span
              for k_dev, (L, case, cell, rep, lo, hi, tag, suffix) in enumerate(built):
                  off = s_offset(L) if block.get('s_bend') else 0.0
                  profile = PROFILES[tag]
                  if block.get('column_pitch'):                   # aligned columns: left edges on a pitch
                      x_left = ROW_LEFT_EDGE_X + k_dev * block['column_pitch']
                  if k_dev > 0 and align == 'left' and not block.get('column_pitch'):
                      movable[case] = dict(row_y=y, lo=lo, hi=hi, row=row, order=k_dev,
                                           max_fields=block.get('max_fields'))
                  if block.get('max_fields'):
                      limited.setdefault(id(block), dict(limit=block['max_fields'], cases=[], leads=[]))
                      limited[id(block)]['cases'].append(case)
                      if k_dev == 0 or block.get('column_pitch'):   # aligned columns: all devices
                          limited[id(block)]['leads'].append(case)
                  dev_x = round(x_left - lo, 3)
                  x_left += (hi - lo) + DEVICE_GAP
                  text = f'L{int(L)}_{key}_{tag}_{mode}{suffix}'
                  lab = Cell('LABEL_' + case)
                  tip_x = -(L / 2 + PORT_TO_L10_ENTRANCE + TICK_POSITION)
                  polys, lab_len = label_polygons(text, glyphs, 0.0, 0.0)
                  lab.polygons += [(LABEL_LAYER[0], LABEL_LAYER[1], p) for p in polys]
                  lib.add(lab)
                  lib.add(cell)
                  lx, ly = dev_x + tip_x + LABEL_TIP_X_OFFSET, y + LABEL_Y_OFFSET
                  xe = L / 2 + PORT_TO_L10_ENTRANCE
                  hi_s = TICK_POSITION + TICK_WIDTH / 2 + 1.0          # taper + alignment tips
                  zones = (('taper', -1.0, hi_s),
                           ('label', TICK_POSITION - LABEL_TIP_X_OFFSET - lab_len_of(text, glyphs) - 1.0,
                            TICK_POSITION - LABEL_TIP_X_OFFSET + 1.0))
                  for kind, a_s, b_s in zones:
                      critical.append((y, dev_x - xe - b_s, dev_x - xe - a_s, kind, case))
                      if mode == '2T':
                          critical.append((y + off, dev_x + xe + a_s, dev_x + xe + b_s, kind, case))
                  top.refs.append(dict(name=case, xy=(um(dev_x), um(y)), angle=0, mag=1.0, reflect=False))
                  top.refs.append(dict(name='LABEL_' + case, xy=(um(lx), um(ly)), angle=0, mag=1.0, reflect=False))
                  out_len = OUTPUT_EXTENSION_LENGTH if mode == '1T' else 0.0
                  devices.append(dict(
                      case_id=case, mode=mode, set=tag, taper_profile=profile, copy=int(suffix[1:] or 1),
                      label=text, label_layer=list(LABEL_LAYER),
                      width_um=w, l10_stage1_length_um=T, l10_taper_length_um=L10_TAPER_LENGTH,
                      l15_taper_length_um=L15_TAPER_LENGTH, alignment_tip_position_um=TICK_POSITION,
                      layer14_linear_start_distance_um=L14_LINEAR_START,
                      straight_length_um=L, output_extension_length_um=out_len,
                      path_length_excluding_tapers_um=round(rep.get('s_bend_path_length_um', L) + out_len, 4),
                      s_bend=(dict(radius_um=S_BEND_RADIUS, angle_deg=round(rep['s_bend_theta_deg'], 6),
                                   offset_um=round(off, 4)) if off else None),
                      output_end_y_offset_um=off,
                      origin_um=[dev_x, y], row_index=row, row_align=align,
                      global_x_extent_um=[round(dev_x + lo, 3), round(dev_x + hi, 3)], label_origin_um=[round(lx, 3), round(ly, 3)],
                      label_length_um=round(lab_len, 3),
                      input_port_local_um=[-L / 2, 0.0],
                      output_port_local_um=[L / 2 + out_len, off],
                      l10_entrance_local_x_um=-(L / 2 + PORT_TO_L10_ENTRANCE),
                      l15_tip_local_x_um=-(L / 2 + PORT_TO_L10_ENTRANCE + L15_TAPER_LENGTH),
                      alignment_tip_local_x_um=tip_x,
                      l10_tip_local_x_um=-(L / 2 + PORT_TO_L10_ENTRANCE + L10_TAPER_LENGTH),
                      **{k: round(v, 4) for k, v in rep.items()}))
              row += 1
    y0 = FIRST_ROW_Y - FRAME_BOTTOM_MARGIN
    y1 = FIRST_ROW_Y + (row - 1) * ROW_PITCH + FRAME_TOP_MARGIN
    frame, (cx, cy), nmark = build_frame(src_frame, y0, y1)
    lib.add(frame)
    top.refs.append(dict(name='FRAME_AND_MARKERS', xy=(um(cx), um(cy)), angle=0, mag=1.0, reflect=False))
    lib.add(top)

    # EBL field placement: global mask shift (ignoring the movable devices' zones), then
    # per-row gap for the movable devices, then report the final placement
    # every left-aligned block starts where the right-aligned (long) devices start
    starts = [d['global_x_extent_um'][0] for d in devices if d['row_align'] == 'right']
    if starts and ALIGN_BLOCK_STARTS:
        x_start = min(starts)
        for b in BLOCKS:
            if b['align'] != 'left':
                continue
            cases = [d['case_id'] for d in devices if d['row_align'] == 'left']
            delta = x_start - min(d['global_x_extent_um'][0] for d in devices if d['case_id'] in cases)
            shift_devices(top, critical, devices, cases, delta)
    lim_cases = {c for b in limited.values() for c in b['cases']}
    limits = [({c, 'LABEL_' + c}, b['limit']) for b in limited.values() for c in b['cases']]
    ebl0 = place_for_ebl(lib, top, critical, skip=set(movable),
                         y_groups=[{c, 'LABEL_' + c} for c in lim_cases], limits=limits)
    gaps = adjust_movable_gaps(lib, top, critical, movable, devices, ebl0['mask_shift_um'])
    ebl = place_for_ebl(lib, top, critical, fixed=ebl0['mask_shift_um'], limits=limits)
    ebl['alternatives_before_gap_adjustment'] = ebl0['alternatives']
    ebl.pop('alternatives')
    ebl['shared_row_gaps_um'] = gaps
    sx, sy = ebl['mask_shift_um']
    for r in top.refs:
        r['xy'] = (r['xy'][0] + um(sx), r['xy'][1] + um(sy))
    for d in devices:
        d['origin_um'] = [round(d['origin_um'][0] + sx, 3), round(d['origin_um'][1] + sy, 3)]
        d['label_origin_um'] = [round(d['label_origin_um'][0] + sx, 3), round(d['label_origin_um'][1] + sy, 3)]
        d['global_x_extent_um'] = [round(v + sx, 3) for v in d['global_x_extent_um']]
    y0, y1 = y0 + sy, y1 + sy
    write_field_overlay(ebl)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    lib.write(OUTPUT)

    params = dict(
        wg_width_um=WG_WIDTH, l10_stage1_lengths_um=list(L10_STAGE1_LENGTHS),
        long_device_taper_t_um=LONG_TAPER_T, long_device_widths_um=list(LONG_WIDTHS),
        profiles=PROFILES, blocks=[dict(mode=b['mode'], align=b['align'], label=b['label'],
                                         rows=[list(r) for r in b['rows']] if b['rows'] else 'T sweep',
                                         members=[[L, t] for L, t in b['members']]) for b in BLOCKS],
        layer15_taper=dict(start_width='waveguide width', end_width_um=L15_TIP_WIDTH,
                           length='stage-1 length - mid offset'),
        layer10_taper=dict(start_width_um=L10_START_WIDTH, mid_width_um=L10_MID_WIDTH,
                           stage1_length='swept', stage2_length_um=L10_STAGE2_LENGTH,
                           mid_offset_beyond_l15_tip_um=L10_MID_OFFSET, end_width_um=L10_TIP_WIDTH),
        port_to_l10_entrance_um=PORT_TO_L10_ENTRANCE,
        layer2_margin_um=L2_MARGIN, layer2_start_offset_um=L2_START_OFFSET,
        layer2_l10_overlap_um=L2_START_OFFSET,
        layer14_margin_um=L14_MARGIN, layer14_start_offset_um=L14_START_OFFSET,
        layer14_linear_start_before_mid_um=L14_LINEAR_START_BEFORE_MID,
        tick=dict(layer=[15, 0], position_from_l10_entrance_um=TICK_POSITION, width_um=TICK_WIDTH, length_um=TICK_LENGTH, gap_to_l10_edge_um=TICK_GAP_TO_L10),
        local_taper_beams=dict(layer=list(BEAM_LAYER), width_um=BEAM_WIDTH,
                               offsets_from_l2_start_um=list(BEAM_OFFSETS_FROM_L2_START),
                               slab_overlap_um=BEAM_SLAB_OVERLAP, span_um=BEAM_SPAN),
        label=dict(layer=list(LABEL_LAYER), tip_x_offset_um=LABEL_TIP_X_OFFSET, y_offset_um=LABEL_Y_OFFSET),
        far_assembly_shift_um='l10_taper_length - %g per end' % SRC_L10_LENGTH,
        row_pitch_um=ROW_PITCH,
        device_gap_um=DEVICE_GAP,
        row_left_edge_x_um=ROW_LEFT_EDGE_X + sx, row_right_edge_x_um=ROW_RIGHT_EDGE_X + sx,
        ebl_write_fields={k: v for k, v in ebl.items() if k != 'fields_xy'}, frame_bounds_um=[FRAME_X[0] + sx, y0, FRAME_X[1] + sx, y1],
        moved_marker_polygon_count=nmark)
    manifest = dict(
        source_gds=str(SOURCE_GDS.relative_to(ROOT)),
        source_gds_sha256=hashlib.sha256(SOURCE_GDS.read_bytes()).hexdigest(),
        source_device_copied_for_far_assemblies=SOURCE_DEVICE,
        gds=str(OUTPUT.relative_to(ROOT)), gds_sha256=hashlib.sha256(OUTPUT.read_bytes()).hexdigest(),
        top_cell=top.name, device_count=len(devices), parameters=params, devices=devices,
        optical_or_mechanical_simulation_performed=False)
    OUTPUT.with_suffix('.json').write_text(json.dumps(manifest, indent=2))
    with open(OUTPUT.with_suffix('.csv'), 'w', newline='') as f:
        cols = ['case_id', 'mode', 'width_um', 'straight_length_um', 'output_extension_length_um',
                'path_length_excluding_tapers_um', 'origin_um', 'label']
        wr = csv.DictWriter(f, fieldnames=cols)
        wr.writeheader()
        for d in devices:
            wr.writerow({k: d[k] for k in cols})
    print('wrote', OUTPUT, len(devices), 'devices')


# Layout variants written by one run (each overrides the parameters above)
VARIANTS = {
    '': {},                                            # 300 um rows, 250 um circles, long set duplicated
    '_pitch500': dict(ROW_PITCH=500.0, GLUE_CIRCLE_DIAMETER=300.0,   # 500 um rows, source 300 um
                      BLOCKS=tuple(b for b in BLOCKS if 'copy' not in b)),  # circles, no long copies
}


def run_variants():
    global STEM, OUTPUT
    base = dict(STEM=STEM, **{k: globals()[k] for v in VARIANTS.values() for k in v})
    for suffix, over in VARIANTS.items():
        globals().update(base)
        globals().update(over)
        STEM = base['STEM'] + suffix
        OUTPUT = ROOT / 'gds' / (STEM + '.gds')
        cell_pieces.__defaults__[0].clear()           # piece cache is per layout
        main()
    globals().update(base)


if __name__ == '__main__':
    run_variants()
