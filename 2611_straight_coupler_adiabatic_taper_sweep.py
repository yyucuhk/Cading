"""Straight-coupler sweep with adiabatic (quintic smootherstep) tapers, derived from
gds/2611_straight_coupler_L15000_output_comparison.gds (+ .json).

Changes relative to that source layout
--------------------------------------
1. Waveguide (Layer 15) width sweep: 1.58, 1.59, 1.60, 1.61, 1.62 um.
2. Straight length sweep: 500, 5000, 15000 um (both output styles kept:
   "2T" = taper at both ends, "1T" = input taper + 2000 um plain output).
3. Taper lengths/widths changed; each stage keeps the source's adiabatic
   quintic smootherstep profile (zero width slope at both ends of every stage;
   set TAPER_PROFILE = 'linear' for straight-sided tapers):
   - Layer 15: W -> 0.200 um over 60 um.
   - Layer 10: two linear stages over 120 um total (60 um under the Layer-15
     taper plus 60 um extra): 10 um -> 1 um at 62 um (2 um beyond the
     Layer-15 tip), then 1 um -> 0.200 um at 120 um.
4. Device labels are drawn on Layer 15 (source: Layer 2).
5. Layer 10 now extends 2.0 um past the start of the Layer-2 window
   (source: 0.5 um, the JSON's layer2_start_offset_um).

Everything outward of the Layer-10 neck (butterfly tether, Layer-2/14 channels,
glue circles, Layer-13 pad) is copied unchanged from the source device and
shifted outward by the extra Layer-10 length (40 um per end).

Only needs the bundled pure-Python gdsio_min.py (no klayout / gdsfactory).
Run from the Cading repository root:  python 2611_straight_coupler_adiabatic_taper_sweep.py
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
STEM = '2611_straight_coupler_W158-162_L500-15000_adiabatic_taper'
OUTPUT = ROOT / 'gds' / (STEM + '.gds')

# %% [01] USER PARAMETERS (um)
WG_WIDTHS = (1.58, 1.59, 1.60, 1.61, 1.62)   # Layer 15, bottom to top in each block
STRAIGHT_LENGTHS = (500.0, 5000.0, 15000.0)
DEVICE_SETS = ('2T', '1T')                    # 2T: two tapers; 1T: input taper + plain output
OUTPUT_EXTENSION_LENGTH = 2000.0              # 1T plain output beyond the output port
PORT_TO_L10_ENTRANCE = 45.0                   # plain Layer-15 straight between port and taper start

TAPER_PROFILE = 'smootherstep'                # 'smootherstep' (adiabatic, as source) or 'linear'
TAPER_SAMPLE_STEP = 0.25                      # vertex spacing along curved tapers

L15_TAPER_LENGTH = 60.0
L15_TIP_WIDTH = 0.200
L10_START_WIDTH = 10.0
L10_EXTRA_LENGTH = 60.0                       # beyond the Layer-15 tip
L10_TAPER_LENGTH = L15_TAPER_LENGTH + L10_EXTRA_LENGTH
L10_MID_WIDTH = 1.0                           # two linear stages: start -> 1 um -> tip
L10_MID_OFFSET = 2.0                          # 1 um point sits this far beyond the L15 tip
L10_MID_POSITION = L15_TAPER_LENGTH + L10_MID_OFFSET
L10_TIP_WIDTH = 0.200

L2_MARGIN = 2.5                               # per side around Layer 10
L2_START_OFFSET = 2.0                         # Layer 2 starts this far OUTWARD of the L10 entrance
L14_MARGIN = 1.0                              # per side around Layer 10
L14_START_OFFSET = 0.0
L14_LINEAR_START = 58.847                     # straight-sided L14 bridge from here to the neck

TICK_WIDTH = 0.5                              # Layer-15 alignment tips (ticks)
TICK_LENGTH = 10.0
TICK_GAP_TO_L10 = 2.0                         # same gap to the slab edge as the source (2.5 - 0.5)
TICK_POSITION = L10_MID_POSITION              # alignment tips mark the L10 1 um point

BEAM_LAYER = (11, 0)
BEAM_WIDTH = 0.5
BEAM_OFFSETS_FROM_L2_START = (-0.5, 0.5, 1.5)  # same relation to the L2 window as the source
BEAM_SLAB_OVERLAP = 0.5
BEAM_SPAN = 5.096 + 0.5                       # source anchor span + outer overlap

LABEL_LAYER = (15, 0)
LABEL_TIP_X_OFFSET = 2.0                      # from the alignment tips, toward the device
LABEL_Y_OFFSET = 12.0                         # source used 10; raised to clear the wider L11 beams
LABEL_GAP = 0.5

ROW_PITCH = 550.0
FIRST_ROW_Y = -700.0
ROW_GROUPS = ((500.0, 5000.0), (15000.0,))    # lengths sharing one row, placed left to right
DEVICE_GAP = 500.0                            # edge-to-edge gap between devices in one row
ROW_ALIGN = {'2T': 'left', '1T': 'right'}     # per device set
ROW_LEFT_EDGE_X = -6937.847                   # left edge of every left-aligned row
ROW_RIGHT_EDGE_X = 11742.153                  # right edge (plain output end) of every right-aligned row
FRAME_X = (-7400.0, 12300.0)                  # same chip width as the source
FRAME_BOTTOM_MARGIN = 700.0
FRAME_TOP_MARGIN = 650.0

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
    glyphs = extract_glyphs(lib)
    frame = lib.cells['FRAME_AND_MARKERS']
    return lib, far, glyphs, frame


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
    return glyphs


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
    top[-1], bot[-1] = seq[0], seq[-1]
    cell.polygons.append((2, 0, merged_end_polygon(top, bot, seq)))
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


def build_device(name, L, width, mode, far):
    cell = Cell(name)
    h = L / 2
    xe = h + PORT_TO_L10_ENTRANCE
    xt = xe + L15_TAPER_LENGTH
    ss = stations(0.0, L15_TAPER_LENGTH, ())                       # outward from the entrance
    left = [(um(-(xe + s)), -um(l15_half(s, width))) for s in reversed(ss)]   # tip -> entrance, bottom
    left_top = [(um(-(xe + s)), um(l15_half(s, width))) for s in ss]          # entrance -> tip, top
    if mode == '2T':
        right = ([(um(xe + s), -um(l15_half(s, width))) for s in ss] +
                 [(um(xe + s), um(l15_half(s, width))) for s in reversed(ss)])
    else:
        xo = h + OUTPUT_EXTENSION_LENGTH
        right = [(um(xo), -um(width / 2)), (um(xo), um(width / 2))]
    cell.polygons.append((15, 0, left + right + left_top))
    rep = build_end(cell, -1, h, width, far)
    if mode == '2T':
        build_end(cell, 1, h, width, far)
    return cell, rep


# %% [05] Frame
def build_frame(src_frame, y0, y1):
    fx0, fx1 = FRAME_X
    cx, cy = (fx0 + fx1) / 2, (y0 + y1) / 2
    hw, hh = (fx1 - fx0) / 2, (y1 - y0) / 2
    src_hh = 2050.0
    c = Cell('FRAME_AND_MARKERS')
    c.polygons.append((0, 0, [(um(-hw), um(-hh)), (um(hw), um(-hh)), (um(hw), um(hh)), (um(-hw), um(hh))]))
    n = 0
    for l, d, p in src_frame.polygons:
        if l == 0:
            continue
        ymid = sum(q[1] for q in p) / len(p)
        dy = um(hh - src_hh) * (1 if ymid > 0 else -1)
        c.polygons.append((l, d, [(x, y + dy) for x, y in p]))
        n += 1
    return c, (cx, cy), n


# %% [06] Assemble, write, report
def main():
    src, far, glyphs, src_frame = load_source()
    lib = Library('LIB')
    top = Cell('STRAIGHT_COUPLER_ADIABATIC_TAPER_SWEEP')
    devices = []
    row = 0
    assert sorted(L for g in ROW_GROUPS for L in g) == sorted(STRAIGHT_LENGTHS)
    for mode in DEVICE_SETS:
        for group in ROW_GROUPS:
          for w in WG_WIDTHS:
              y = FIRST_ROW_Y + row * ROW_PITCH
              built = []
              for L in group:
                  case = f'SC_L{int(L)}_W{int(round(w * 1000))}_{mode}'
                  cell, rep = build_device(case, L, w, mode, far)
                  xs = [q[0] for _, _, pp in cell.polygons for q in pp]
                  built.append((L, case, cell, rep, min(xs) / DBU_PER_UM, max(xs) / DBU_PER_UM))
              span = sum(b[5] - b[4] for b in built) + DEVICE_GAP * (len(built) - 1)
              x_left = ROW_LEFT_EDGE_X if ROW_ALIGN[mode] == 'left' else ROW_RIGHT_EDGE_X - span
              for L, case, cell, rep, lo, hi in built:
                  dev_x = round(x_left - lo, 3)
                  x_left += (hi - lo) + DEVICE_GAP
                  text = f'L{int(L)}_W{w:.2f}_{mode}'
                  lab = Cell('LABEL_' + case)
                  tip_x = -(L / 2 + PORT_TO_L10_ENTRANCE + TICK_POSITION)
                  polys, lab_len = label_polygons(text, glyphs, 0.0, 0.0)
                  lab.polygons += [(LABEL_LAYER[0], LABEL_LAYER[1], p) for p in polys]
                  lib.add(lab)
                  lib.add(cell)
                  lx, ly = dev_x + tip_x + LABEL_TIP_X_OFFSET, y + LABEL_Y_OFFSET
                  top.refs.append(dict(name=case, xy=(um(dev_x), um(y)), angle=0, mag=1.0, reflect=False))
                  top.refs.append(dict(name='LABEL_' + case, xy=(um(lx), um(ly)), angle=0, mag=1.0, reflect=False))
                  out_len = OUTPUT_EXTENSION_LENGTH if mode == '1T' else 0.0
                  devices.append(dict(
                      case_id=case, mode=mode, label=text, label_layer=list(LABEL_LAYER),
                      width_um=w, straight_length_um=L, output_extension_length_um=out_len,
                      path_length_excluding_tapers_um=L + out_len,
                      origin_um=[dev_x, y], row_index=row, row_align=ROW_ALIGN[mode],
                      global_x_extent_um=[round(dev_x + lo, 3), round(dev_x + hi, 3)], label_origin_um=[round(lx, 3), round(ly, 3)],
                      label_length_um=round(lab_len, 3),
                      input_port_local_um=[-L / 2, 0.0],
                      output_port_local_um=[L / 2 + out_len, 0.0],
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
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    lib.write(OUTPUT)

    params = dict(
        wg_widths_um=list(WG_WIDTHS), straight_lengths_um=list(STRAIGHT_LENGTHS), device_sets=list(DEVICE_SETS),
        taper_profile=TAPER_PROFILE,
        layer15_taper=dict(start_width='waveguide width', end_width_um=L15_TIP_WIDTH, length_um=L15_TAPER_LENGTH),
        layer10_taper=dict(start_width_um=L10_START_WIDTH, mid_width_um=L10_MID_WIDTH,
                           mid_position_um=L10_MID_POSITION, mid_offset_beyond_l15_tip_um=L10_MID_OFFSET,
                           end_width_um=L10_TIP_WIDTH, length_um=L10_TAPER_LENGTH,
                           extra_length_beyond_l15_tip_um=L10_EXTRA_LENGTH),
        port_to_l10_entrance_um=PORT_TO_L10_ENTRANCE,
        layer2_margin_um=L2_MARGIN, layer2_start_offset_um=L2_START_OFFSET,
        layer2_l10_overlap_um=L2_START_OFFSET,
        layer14_margin_um=L14_MARGIN, layer14_start_offset_um=L14_START_OFFSET,
        layer14_linear_start_distance_um=L14_LINEAR_START,
        tick=dict(layer=[15, 0], position_from_l10_entrance_um=TICK_POSITION, width_um=TICK_WIDTH, length_um=TICK_LENGTH, gap_to_l10_edge_um=TICK_GAP_TO_L10),
        local_taper_beams=dict(layer=list(BEAM_LAYER), width_um=BEAM_WIDTH,
                               offsets_from_l2_start_um=list(BEAM_OFFSETS_FROM_L2_START),
                               slab_overlap_um=BEAM_SLAB_OVERLAP, span_um=BEAM_SPAN),
        label=dict(layer=list(LABEL_LAYER), tip_x_offset_um=LABEL_TIP_X_OFFSET, y_offset_um=LABEL_Y_OFFSET),
        far_assembly_shift_um=L10_TAPER_LENGTH - SRC_L10_LENGTH,
        row_pitch_um=ROW_PITCH, row_groups_um=[list(g) for g in ROW_GROUPS], device_gap_um=DEVICE_GAP,
        row_align=ROW_ALIGN, row_left_edge_x_um=ROW_LEFT_EDGE_X, row_right_edge_x_um=ROW_RIGHT_EDGE_X, frame_bounds_um=[FRAME_X[0], y0, FRAME_X[1], y1],
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


if __name__ == '__main__':
    main()
