"""Design 1 L6000: four three-width sets (12 devices), all dimensions in um.

Two meander sets and two single-straight sets. New single straights are 6000 um
with two couplers, or one input coupler plus the retained 2000 um bare output.
Complete optical terminal assemblies are copied, never redrawn.
L2 starts 0.5 um outward from L10, at TOTAL width 15 um.
L14 starts with L10, at TOTAL width 8 um. Downstream margins are 2.5 um
per side on L2 and 1.0 um per side on L14. First-stage window envelopes blend smoothly over the original stage.
L16/L51 and the old PNC patterns are absent. L11 contains only local plain
beam supports at taper entrances; global marks are on L15 and labels on L2.
Requires the existing source modules and source width-sweep GDS/manifest.
"""
from __future__ import annotations

# %% [00] Imports and paths
import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
BASE_BUILDER = ROOT / '2611_design1_L6000_P10_width_sweep.py'
BASE_GDS = ROOT / 'gds/2611_design1_L6000_P10_width_sweep.gds'
STEM = '2611_design1_L6000_P10_straight_output_comparison'
OUTPUT = ROOT / 'gds' / (STEM + '.gds')
WORK = ROOT / 'validation_2611_straight_output'

# %% [01] USER PARAMETERS - EDIT HERE
STRAIGHT_LENGTH = 6000.0              # each original meander straight
WG_WIDTHS = (1.62, 1.64, 1.68)        # Layer 15; bottom to top in each set
OUTPUT_EXTENSION_LENGTH = 2000.0     # from the lower-right meander port
PNC_PITCH = 10.0                     # straight/bend; also added straight
INCLUDE_REFERENCE_DEVICES = True     # two-coupler set first, straight outputs above
COUPLER_END_WINDOWS_ONLY = True       # applies to BOTH device sets; crops L2 and L14
LAYER2_WINDOW_MARGIN = 2.5            # per-side taper-window expansion on Layer 2
LAYER14_WINDOW_MARGIN = 1.0           # independent Layer 14 expansion
RELEASE_WINDOW_START_OFFSET = 0.0     # legacy construction setting; see independent offsets below
LAYER2_START_OFFSET = 0.5             # preserve original L10-to-L2 delayed-copy interface
LAYER14_START_OFFSET = 0.0            # same entrance as L10
LAYER2_START_WIDTH = 15.0             # TOTAL entrance width, not extra margin
LAYER14_START_WIDTH = 8.0             # TOTAL entrance width, not extra margin
ROW_PITCH = 600.0
FIRST_ROW_Y = -700.0
LEFT_BEND_X = -900.0
LABEL_SIZE = 5.0
FRAME_MARGIN = 500.0
MIN_FRAME_HEIGHT = 4000.0
DBU = 0.001
FAB_LAYERS = (2, 10, 11, 13, 14, 15, 16)  # pre-export construction checks
REMOVED_OUTPUT_LAYERS = (16, 51)       # absent from every final output cell
CLEARED_SOURCE_PNC_LAYERS = (11,)     # remove all old PNCs before adding local beams
LOCAL_TAPER_BEAM_LAYER = (11, 0)      # same solid-support mask as original PNCs
LOCAL_TAPER_BEAM_WIDTH = 0.5          # longitudinal width of each transverse beam
LOCAL_TAPER_BEAM_PITCH = 1.0          # center-to-center, measured along the taper
LOCAL_TAPER_BEAM_COUNT = 3            # three stations; TWO beams per station
LOCAL_TAPER_BEAM_START = 0.0          # first center aligned with the L10 entrance
LOCAL_TAPER_ANCHOR_SPAN = 5.096       # original PNC anchor-to-anchor span
LOCAL_TAPER_ANCHOR_OVERLAP = 0.5      # overlap into slab and outer anchor
ALIGNMENT_MARKER_SOURCE_LAYER = (2, 0)
ALIGNMENT_MARKER_TARGET_LAYER = (15, 0)  # global alignment markers ONLY


def require(ok, message):
    if not ok:
        raise ValueError(message)


def check(checks, name, ok, **details):
    checks.append(dict(name=name, passed=bool(ok), **details))
    require(ok, name + ': ' + str(details))


def load_base():
    spec = importlib.util.spec_from_file_location('long_design1', BASE_BUILDER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.STRAIGHT_LENGTH = STRAIGHT_LENGTH
    module.STRAIGHT_PNC_PERIOD = module.BEND_PNC_PERIOD = PNC_PITCH
    return module


def box(bounds):
    import klayout.db as kdb
    return kdb.Region(kdb.Box(*(round(v / DBU) for v in bounds)))


def regions(layout, cell):
    import klayout.db as kdb
    result = {}
    for li in layout.layer_indexes():
        r = kdb.Region(cell.begin_shapes_rec(li)).merged()
        if not r.is_empty():
            info = layout.get_info(li)
            result[(info.layer, info.datatype)] = r
    return result


def read_gds(path):
    import klayout.db as kdb
    ly = kdb.Layout()
    ly.read(str(path))
    require(abs(ly.dbu - DBU) < 1e-12, 'Expected 1 nm database grid')
    require(len(ly.top_cells()) == 1, 'Expected exactly one top cell')
    return ly, ly.top_cell()


def insert_regions(ly, name, rr):
    c = ly.create_cell(name)
    for (n, d), r in rr.items():
        c.shapes(ly.layer(n, d)).insert(r)
    return c


def write_layout(ly, path):
    import klayout.db as kdb
    options = kdb.SaveLayoutOptions()
    options.gds2_max_vertex_count = 2000
    path.parent.mkdir(parents=True, exist_ok=True)
    ly.write(str(path), options)


def finalize_output_layers(layout, frame_cell):
    """Apply only the requested layer changes to the independently owned layout.

    The frame cell contains chip outline and alignment marks, never device
    windows or labels. Delete excluded layers across every cell, including
    unreferenced cells. Never remap all of Layer 2: it also carries release
    windows and device labels. Original source/baseline layouts are read-only.
    """
    import klayout.db as kdb
    require(not list(frame_cell.each_inst()), 'Expected an owned flat frame cell')
    source = layout.layer(*ALIGNMENT_MARKER_SOURCE_LAYER)
    target = layout.layer(*ALIGNMENT_MARKER_TARGET_LAYER)
    moved = kdb.Region(frame_cell.shapes(source)).merged()
    require(not moved.is_empty(), 'No source alignment markers were found')
    existing = kdb.Region(frame_cell.shapes(target)).merged()
    frame_cell.shapes(target).clear()
    frame_cell.shapes(target).insert((existing | moved).merged())
    frame_cell.shapes(source).clear()
    removed_shapes = {}
    for index in tuple(layout.layer_indexes()):
        info = layout.get_info(index)
        if info.layer in (*REMOVED_OUTPUT_LAYERS, *CLEARED_SOURCE_PNC_LAYERS):
            key = str(info.layer) + '/' + str(info.datatype)
            removed_shapes[key] = sum(cell.shapes(index).size() for cell in layout.each_cell())
            layout.delete_layer(index)
    return dict(removed_layers=list(REMOVED_OUTPUT_LAYERS),
                cleared_source_pnc_layers=list(CLEARED_SOURCE_PNC_LAYERS),
                removed_shape_counts=removed_shapes,
                marker_source_layer=list(ALIGNMENT_MARKER_SOURCE_LAYER),
                marker_target_layer=list(ALIGNMENT_MARKER_TARGET_LAYER),
                marker_cell=frame_cell.name,
                moved_marker_polygon_count=moved.count(),
                moved_marker_area_um2=moved.area() * layout.dbu ** 2,
                layer2_labels_preserved=True, release_masks_handled_separately=True)


def finalize_coupler_windows(layout, entries):
    """Apply the independent entrance widths, start positions and margins.

    L2 starts 0.5 um outward from L10, preserving the original delayed-copy
    interface. L14 starts with L10. Entrance widths are TOTAL widths, not
    extra clearances: L2=15 um and L14=8 um. Downstream margins: L2=2.5 um, L14=1 um.
    Only the first 40 um window envelope is replaced; downstream mask
    candidates still come from the original mask generator.
    """
    import klayout.db as kdb
    policies = {}
    require(COUPLER_END_WINDOWS_ONLY, 'This comparison requires terminal-only masks')
    for e in entries:
        cell = layout.cell(e['case_id'])
        rr = regions(layout,cell)
        ml,mc = read_gds(WORK/e['release_window_file'])
        candidates = regions(ml,mc)
        result = {2:kdb.Region(),14:kdb.Region()}
        ends=[]
        for side,port in (('input',e['input_port_local_um']),('output',e['output_port_local_um'])):
            px,cy = map(float,port)
            sign = -1 if px<0 else 1
            search = box((-20000,cy/8,px,20000) if sign<0 else (px,-20000,20000,cy/8))
            taper = (rr[(10,0)] & search).merged()
            if taper.is_empty():
                require(side=='output' and e['mode']=='straight_output','Unexpected absent taper')
                continue
            sx = (taper.bbox().right if sign<0 else taper.bbox().left)*DBU
            end = dict(side=side,center_y_um=cy,l10_start_local_x_um=sx)
            for layer,offset,start_width,margin in (
                (2,LAYER2_START_OFFSET,LAYER2_START_WIDTH,LAYER2_WINDOW_MARGIN),
                (14,LAYER14_START_OFFSET,LAYER14_START_WIDTH,LAYER14_WINDOW_MARGIN)):
                start = sx+sign*offset
                keep = box((-20000,cy/8,start,20000) if sign<0 else (start,-20000,20000,cy/8))
                retained = (rr[(layer,0)] & candidates[(layer,0)] & keep).merged()
                endpoint = sx+sign*40.0
                replace = box((min(sx,endpoint),cy-25,max(sx,endpoint),cy+25))
                shaped = first_stage_mask(sx,cy,sign,offset,start_width,1.0+2*margin)
                retained = ((retained-replace) | shaped).merged()
                require((retained-keep).is_empty(),'Mask extends before its entrance')
                result[layer] |= retained
                end[f'layer{layer}_start_local_x_um']=start
                end[f'layer{layer}_start_global_x_um']=start+e['origin_um'][0]
                end[f'layer{layer}_start_width_um']=start_width
                end[f'layer{layer}_offset_um']=offset
            ends.append(end)
        require(len(ends)==(1 if e['mode']=='straight_output' else 2),'Wrong number of couplers')
        for n,r in result.items():
            cell.shapes(layout.layer(n,0)).clear()
            cell.shapes(layout.layer(n,0)).insert(r.merged())
        policies[e['case_id']]=dict(ends=ends,layer2_margin_um=LAYER2_WINDOW_MARGIN,
            layer14_margin_um=LAYER14_WINDOW_MARGIN)
        e.update(final_pnc_present=False,final_release_window_count=len(ends),
                 layer2_margin_um=LAYER2_WINDOW_MARGIN,layer14_margin_um=LAYER14_WINDOW_MARGIN,
                 layer2_start_offset_um=LAYER2_START_OFFSET,layer14_start_offset_um=LAYER14_START_OFFSET,
                 layer2_start_width_um=LAYER2_START_WIDTH,layer14_start_width_um=LAYER14_START_WIDTH)
    return policies


def build_release_window_candidate(base, width):
    """Create isolated L2/L14 candidates using the existing mask construction.

    Other candidate layers, including labels and optical shapes, are discarded.
    Final intersection preserves the old channels, pads and terminal geometry.
    Each worker owns its configuration and uses a separate cached-cell namespace.
    """
    import gdsfactory as gf
    import klayout.db as kdb
    require(math.isfinite(RELEASE_WINDOW_START_OFFSET) and RELEASE_WINDOW_START_OFFSET >= 0,
            'Release start offset must be finite and non-negative')
    candidates = {}
    margins = {2: LAYER2_WINDOW_MARGIN, 14: LAYER14_WINDOW_MARGIN}
    for margin in sorted(set(margins.values())):
        require(math.isfinite(margin) and 0 < margin <= 3.0,
                'Release margins must be finite, positive and at most the original 3 um')
        mask_model = base.configured_design(width, {'POSITIVE_WINDOW_MARGIN': margin})
        mask_model.validate_settings()
        component = mask_model.make_device(STRAIGHT_LENGTH)
        for layer, requested in margins.items():
            if requested != margin:
                continue
            region = kdb.Region(component.kdb_cell.begin_shapes_rec(gf.get_layer((layer,0)))).merged()
            require(not region.is_empty(), f'No Layer{layer} mask candidate generated')
            candidates[(layer,0)] = region
    name = ('RELEASE_M2_' + str(LAYER2_WINDOW_MARGIN).replace('.', 'p')
            + '_M14_' + str(LAYER14_WINDOW_MARGIN).replace('.', 'p')
            + '_' + base.case_id(width))
    ly = kdb.Layout(); ly.dbu = DBU
    insert_regions(ly, name, candidates)
    write_layout(ly, WORK/(name+'.gds'))
    return name+'.gds'


def add_local_taper_beams(layout, entries, policies):
    """Add plain transverse beam pairs only at the optical taper entrances.

    Old PNCs have already been cleared. Layer 11 now contains these rectangles
    only: three stations at 0, 1 and 2 um outward from L10, on both sides.
    The 5.096 um anchor span and 0.5 um end overlap match the old PNC support.
    No optical polygon or release-mask boundary is edited here.
    """
    import klayout.db as kdb
    require(LOCAL_TAPER_BEAM_LAYER == (11, 0), 'Use the original solid support layer')
    require(LOCAL_TAPER_BEAM_COUNT >= 1 and isinstance(LOCAL_TAPER_BEAM_COUNT, int),
            'Beam count must be a positive integer')
    require(0 < LOCAL_TAPER_BEAM_WIDTH < LOCAL_TAPER_BEAM_PITCH,
            'Beam width must be positive and smaller than its pitch')
    require(LOCAL_TAPER_BEAM_START >= 0 and LOCAL_TAPER_ANCHOR_SPAN > 0 and
            LOCAL_TAPER_ANCHOR_OVERLAP > 0, 'Invalid support dimensions')
    last = LOCAL_TAPER_BEAM_START + (LOCAL_TAPER_BEAM_COUNT-1)*LOCAL_TAPER_BEAM_PITCH
    require(last+LOCAL_TAPER_BEAM_WIDTH/2 < 40, 'Support group must stay in the first taper stage')
    li = layout.layer(*LOCAL_TAPER_BEAM_LAYER)
    records = []
    for entry in entries:
        cell = layout.cell(entry['case_id'])
        require(cell.shapes(li).size()==0, 'Unexpected pre-existing PNC/support geometry')
        rr = regions(layout, cell)
        supports = kdb.Region()
        for end in policies[entry['case_id']]['ends']:
            sx, cy = float(end['l10_start_local_x_um']), float(end['center_y_um'])
            direction = -1.0 if end['side']=='input' else 1.0
            for station in range(LOCAL_TAPER_BEAM_COUNT):
                distance = LOCAL_TAPER_BEAM_START + station*LOCAL_TAPER_BEAM_PITCH
                x = sx + direction*distance
                x1, x2 = x-LOCAL_TAPER_BEAM_WIDTH/2, x+LOCAL_TAPER_BEAM_WIDTH/2
                slab = (rr[(10, 0)] & box((x1,cy-10,x2,cy+10))).merged()
                require(not slab.is_empty(), 'No L10 slab at the beam attachment')
                bb = slab.bbox()
                for side in (-1, 1):
                    edge = (bb.top if side > 0 else bb.bottom)*DBU
                    inner = edge-side*LOCAL_TAPER_ANCHOR_OVERLAP
                    outer = edge+side*(LOCAL_TAPER_ANCHOR_SPAN+LOCAL_TAPER_ANCHOR_OVERLAP)
                    bounds = (x1, min(inner,outer), x2, max(inner,outer))
                    beam = box(bounds)
                    require(not (beam & rr[(10,0)]).is_empty(), 'Beam is detached from L10')
                    require(not (beam - rr[(2,0)] - rr[(14,0)]).is_empty(),
                            'Beam does not reach the outer anchor beyond the release windows')
                    require((beam & rr[(15,0)]).is_empty(), 'Beam touches the optical ridge')
                    require((beam & supports).is_empty(), 'Adjacent beams intersect')
                    supports |= beam
                    records.append(dict(case_id=entry['case_id'],end=end['side'],side=side,
                        station=station,distance_from_l10_um=distance,center_local_x_um=x,
                        center_global_x_um=x+entry['origin_um'][0],bounds_local_um=list(bounds)))
        cell.shapes(li).insert(supports)
        expected = 2*LOCAL_TAPER_BEAM_COUNT*len(policies[entry['case_id']]['ends'])
        require(supports.count()==expected, 'Incorrect local beam population')
        entry['local_taper_beam_count'] = expected
        entry['local_taper_beam_layer'] = list(LOCAL_TAPER_BEAM_LAYER)
    expected_total = 2*LOCAL_TAPER_BEAM_COUNT*sum(len(p['ends']) for p in policies.values())
    require(len(records)==expected_total, 'Wrong total support count')
    return dict(layer=list(LOCAL_TAPER_BEAM_LAYER),width_um=LOCAL_TAPER_BEAM_WIDTH,
        pitch_um=LOCAL_TAPER_BEAM_PITCH,stations_per_end=LOCAL_TAPER_BEAM_COUNT,
        start_um=LOCAL_TAPER_BEAM_START,anchor_span_um=LOCAL_TAPER_ANCHOR_SPAN,
        anchor_overlap_um=LOCAL_TAPER_ANCHOR_OVERLAP,total_beams=len(records),beams=records)


# %% [02] Complete all-layer straight-output factory

def make_straight_output(reference, model, width, port, length):
    """Return independent all-layer regions; never mutate reference/shared cells."""
    import gdsfactory as gf
    import klayout.db as kdb
    import numpy as np
    px, py = map(float, port)
    require(px > 0 and py < 0, 'This variant uses the lower-right output only')
    require(length > PNC_PITCH and math.isfinite(length), 'Extension is too short')
    require(abs(length / PNC_PITCH - round(length / PNC_PITCH)) < 1e-8,
            'Use a full number of 10 um periods for this straight output')
    bounds = kdb.Box()
    for r in reference.values():
        bounds += r.bbox()
    # Exactly at the final meander port, below the middle arm. The right bend
    # occupies y >= 0; the lower output sits at y=-200. Nothing inside the
    # meander or input assembly is replaced. At y=-25 this also includes the
    # complete output glue pads (whose top is y=-50), not clipped pad caps.
    roi = box((px, bounds.bottom * DBU - 10,
               max(bounds.right * DBU, px + length) + 10, py / 8))
    end_x = px + length
    half = model.S_BEND_WIDTH / 2
    margin = model.RIDGE_PROTECTION_MARGIN
    tether = model.phononic_tether(
        width=model.REVIEW_TETHER_CONNECTION_WIDTH,
        segments=model.PHONONIC_SUPPORT_SEGMENTS,
        gaps=model.PHONONIC_SUPPORT_GAPS,
        overshoot=model.PHONONIC_SUPPORT_OVERSHOOT,
    )
    tether_length = float(tether.ports['a2'].center[0] - tether.ports['a1'].center[0])
    # Continue the same masks as the adjacent straight: L2 side bands,
    # L14 solid release opening, L15 optical ridge and L16 protection.
    extra = {
        (2, 0): box((px, py + half, end_x, py + half + tether_length)) |
                box((px, py - half - tether_length, end_x, py - half)),
        (14, 0): box((px, py - half - tether_length, end_x, py + half + tether_length)),
        (15, 0): box((px, py - width / 2, end_x, py + width / 2)),
        (16, 0): box((px, py - width / 2 - margin,
                      end_x + margin, py + width / 2 + margin)),
    }
    support = gf.Component()
    stations = np.arange(0.0, length - 1e-8, PNC_PITCH)
    # A full pair at x=0 replaces/rejoins the old port pair. The last pair
    # is 10 um before the flat end; no taper or 1 um terminal cluster is added.
    for x in stations:
        for side in (-1, 1):
            support.add_ref(tether).rotate(side * 90).move((px + float(x), py + side * half))
    extra[(11, 0)] = kdb.Region(support.kdb_cell.begin_shapes_rec(gf.get_layer((11, 0)))).merged()
    output = {layer: ((r - roi) | extra.get(layer, kdb.Region())).merged()
              for layer, r in reference.items()}
    return output, extra, roi, stations.tolist(), tether_length


# %% [03] Fresh-process case build and exact geometry regression

def build_case(width):
    import gdsfactory as gf
    import klayout.db as kdb
    import numpy as np
    gf.gpdk.PDK.activate()
    base = load_base()
    m = base.configured_design(width)
    m.validate_settings()
    # The approved GDS is the exact reference, not a new reconstruction.
    # Reading its independently owned polygons avoids export-grid differences
    # in unchanged corners and preserves every existing layer exactly.
    ref_id = base.case_id(width)
    baseline, _ = read_gds(BASE_GDS)
    saved = baseline.cell(ref_id)
    require(saved is not None, 'Reference GDS does not contain the requested width')
    manifest = json.loads(BASE_GDS.with_suffix('.json').read_text())
    require(manifest['gds_sha256'] == hashlib.sha256(BASE_GDS.read_bytes()).hexdigest(),
            'Approved GDS and manifest hashes disagree')
    record = next(v for v in manifest['devices'] if v['case_id'] == ref_id)
    require(record['width_um'] == width and record['straight_length_um'] == STRAIGHT_LENGTH,
            'Approved reference dimensions differ from requested dimensions')
    require(record['straight_pnc_period_um'] == PNC_PITCH,
            'Approved reference has a different straight PNC pitch')
    original = regions(baseline, saved)
    centers0 = record['connection_centers']
    tips0 = record['tip_centers']
    orientations = [math.degrees(math.atan2(c[1]-t[1], c[0]-t[0])) % 360
                    for c, t in zip(centers0, tips0, strict=True)]
    info = dict(connection_centers=centers0, tip_centers=tips0,
                connection_orientations=orientations,
                total_meander_length_um=record['total_meander_um'])
    centers = [list(map(float, p)) for p in info['connection_centers']]
    oi = int(np.argmax([p[0] for p in centers]))
    ii = 1 - oi
    port = centers[oi]
    require(np.allclose(port, [STRAIGHT_LENGTH / 2, -200], rtol=0, atol=.001),
            'Unexpected output port; inspect the new router before proceeding')
    changed, extra, roi, stations, tether_length = make_straight_output(
        original, m, width, port, OUTPUT_EXTENSION_LENGTH)
    checks = []
    ref_id = base.case_id(width)
    variant_id = ref_id + f'_OUT{OUTPUT_EXTENSION_LENGTH:g}'
    baseline, _ = read_gds(BASE_GDS)
    saved = baseline.cell(ref_id)
    require(saved is not None, 'Reference GDS does not contain the requested width')
    saved_regions = regions(baseline, saved)
    check(checks, 'all_reference_layers_present', set(original) == set(saved_regions))
    for layer, r in original.items():
        area = (r ^ saved_regions[layer]).area() * DBU ** 2
        check(checks, f'approved_reference_L{layer[0]}_exact', area == 0, xor_um2=area)
        area = ((changed[layer] ^ r) - roi).area() * DBU ** 2
        check(checks, f'input_meander_bends_L{layer[0]}_unchanged', area == 0, xor_um2=area)
        area = ((changed[layer] ^ extra.get(layer, kdb.Region())) & roi).area() * DBU ** 2
        check(checks, f'output_L{layer[0]}_only_intended_straight', area == 0, xor_um2=area)
    input_ticks = m.alignment_tip_region([info['tip_centers'][ii]],
                                         [info['connection_orientations'][ii]], DBU)
    all_ticks = m.alignment_tip_region(info['tip_centers'], info['connection_orientations'], DBU)
    optical = (changed[(15, 0)] - input_ticks).merged()
    check(checks, 'one_connected_optical_path', optical.count() == 1)
    check(checks, 'one_connected_protection_path', changed[(16, 0)].count() == 1)
    check(checks, 'only_two_input_ticks_retained', input_ticks.count() == 2 and
          (input_ticks - changed[(15, 0)]).is_empty() and
          ((all_ticks - input_ticks) & changed[(15, 0)]).is_empty())
    for mask in (14, 16):
        missing = (optical - changed[(mask, 0)]).area() * DBU ** 2
        check(checks, f'layer{mask}_covers_entire_optical_path', missing == 0, missing_um2=missing)
    px, py = port
    end_x = px + OUTPUT_EXTENSION_LENGTH
    for dx in [-.02, .02, 1, 10, 100, 500, 1000, 1500, OUTPUT_EXTENSION_LENGTH - .02]:
        sl = (optical & box((px + dx, py - 3, px + dx + DBU, py + 3))).merged()
        actual = sl.bbox().height() * DBU
        check(checks, f'output_width_at_{dx:g}', sl.count() == 1 and abs(actual - width) < .0011,
              width_um=actual)
    check(checks, 'exact_2000um_optical_extension',
          optical.bbox().right * DBU == end_x and
          (optical & box((end_x + DBU, py - 3, end_x + 10, py + 3))).is_empty(),
          port_x_um=px, end_x_um=end_x, length_um=OUTPUT_EXTENSION_LENGTH)
    # Inspect the EXPORTED-STAGE polygon population, not sampler metadata.
    pnc_local = (changed[(11, 0)] & box((px - 1, py - 30, end_x + 1, py + 30))).merged()
    centers_by_side = {1: [], -1: []}
    for poly in pnc_local.each():
        b = poly.bbox()
        side = 1 if b.center().y * DBU > py else -1
        centers_by_side[side].append((b.left + b.right) * DBU / 2 - px)
    for side, values in centers_by_side.items():
        values = np.sort(values)
        check(checks, f'extension_side{side}_exact_10um_PNC_grid',
              len(values) == len(stations) and np.allclose(values, stations, rtol=0, atol=.0011),
              count=len(values), first_um=float(values[0]), last_um=float(values[-1]),
              min_pitch_um=float(np.diff(values).min()), max_pitch_um=float(np.diff(values).max()))
    for layer in (10, 13, 51):
        check(checks, f'no_old_output_taper_pad_or_guide_L{layer}',
              (changed[(layer, 0)] & roi).is_empty())
    support_cut = changed[(14, 0)] & box((px + 100, py - 30, px + 100 + DBU, py + 30))
    check(checks, 'constant_15um_support_envelope',
          abs(support_cut.bbox().height() * DBU - (m.S_BEND_WIDTH + 2*tether_length)) < .0021)
    candidate_file=build_release_window_candidate(base,width)
    report = dict(width_um=width, reference_id=ref_id, variant_id=variant_id,
                  release_window_file=candidate_file,
                  layer2_margin_um=LAYER2_WINDOW_MARGIN,
                  layer14_margin_um=LAYER14_WINDOW_MARGIN,
                  release_start_offset_um=RELEASE_WINDOW_START_OFFSET,
                  input_port_um=centers[ii], output_port_um=port,
                  straight_output_endpoint_um=[end_x, py],
                  original_straight_lengths_um=[STRAIGHT_LENGTH]*3,
                  meander_length_um=float(info['total_meander_length_um']),
                  new_path_length_excluding_input_taper_um=float(info['total_meander_length_um'])+OUTPUT_EXTENSION_LENGTH,
                  output_extension_length_um=OUTPUT_EXTENSION_LENGTH,
                  output_pnc_pairs=len(stations), pnc_pitch_um=PNC_PITCH,
                  input_taper_total_um=float(m.RIDGE_TAPER_LENGTH+m.SLAB_TAPER_LENGTH),
                  input_pnc_transition_um=float(m.PNC_TAPER_LENGTH), input_pnc_cutoff_width_um=float(m.PNC_TAPER_END_WIDTH),
                  input_dense_pnc_periods=int(m.PNC_CONNECTION_PERIODS),
                  input_dense_pnc_pitch_um=float(m.PNC_CONNECTION_PITCH),
                  left_pnc_x_um=original[(11,0)].bbox().left*DBU,
                  changed_window_um=[roi.bbox().left*DBU, roi.bbox().bottom*DBU,
                                     roi.bbox().right*DBU, roi.bbox().top*DBU],
                  checks=checks, status='PASS')
    for name, rr in ((ref_id, original), (variant_id, changed)):
        ly=kdb.Layout(); ly.dbu=DBU
        insert_regions(ly, name, rr)
        write_layout(ly, WORK/(name+'.gds'))
        ey, ec=read_gds(WORK/(name+'.gds'))
        er=regions(ey,ec)
        for layer,r in rr.items():
            check(checks, f'{name}_export_L{layer[0]}',
                  (r ^ er[layer]).is_empty())
    (WORK/(ref_id+'.json')).write_text(json.dumps(report,indent=2))
    print(ref_id, 'PASS', len(checks), 'checks', flush=True)


# %% [04] Whole-device placement and exactly one label per device

def combine():
    """Place four sets of three widths; preserve the first six device origins."""
    import copy
    import gdsfactory as gf
    import klayout.db as kdb
    gf.gpdk.PDK.activate()
    base=load_base()
    ly=kdb.Layout();ly.dbu=DBU
    top=ly.create_cell('D1_L6000_STRAIGHT_OUT_COMPARISON')
    entries=[];checks=[];bounds=kdb.Box();texts=[]

    def place(ident,rr,entry):
        nonlocal bounds
        row=len(entries);y=FIRST_ROW_Y+row*ROW_PITCH
        x=entry['origin_um'][0]
        cell=insert_regions(ly,ident,rr)
        tr=kdb.Trans(round(x/DBU),round(y/DBU))
        top.insert(kdb.CellInstArray(cell.cell_index(),tr))
        entry.update(case_id=ident,origin_um=[tr.disp.x*DBU,tr.disp.y*DBU])
        for (n,d),r in rr.items():
            if n in FAB_LAYERS:bounds+=r.transformed(tr).bbox()
        text=gf.components.text(text=entry['label'],size=LABEL_SIZE,layer=(2,0))
        textreg=kdb.Region(text.kdb_cell.begin_shapes_rec(gf.get_layer((2,0)))).merged()
        labelcell=insert_regions(ly,'LABEL_'+ident,{(2,0):textreg})
        lt=kdb.Trans(round((LEFT_BEND_X+200)/DBU),round((y+250)/DBU))
        top.insert(kdb.CellInstArray(labelcell.cell_index(),lt))
        texts.append(textreg.transformed(lt));bounds+=texts[-1].bbox()
        entry['label_origin_um']=[lt.disp.x*DBU,lt.disp.y*DBU]
        entries.append(entry)

    # A/B: existing meander designs, retaining every original end and origin.
    modes=([False] if INCLUDE_REFERENCE_DEVICES else [])+[True]
    for one in modes:
        for width in WG_WIDTHS:
            meta=json.loads((WORK/(base.case_id(width)+'.json')).read_text())
            ident=meta['variant_id'] if one else meta['reference_id']
            src,sc=read_gds(WORK/(ident+'.gds'))
            rr=regions(src,sc)
            label=f'D1_W{width:.2f}_T80_P10'+(f'_OUT{OUTPUT_EXTENSION_LENGTH:g}' if one else '')
            place(ident,rr,dict(mode='straight_output' if one else 'original_two_tapers',
                topology='meander',set_name='meander_one_coupler' if one else 'meander_two_couplers',
                width_um=width,origin_um=[LEFT_BEND_X-meta['left_pnc_x_um'],0],label=label,
                straight_section_count=3,straight_length_um=STRAIGHT_LENGTH,
                bend_count=2,added_length_um=OUTPUT_EXTENSION_LENGTH if one else 0.,
                path_length_excluding_tapers_um=meta['meander_length_um']+(OUTPUT_EXTENSION_LENGTH if one else 0.),
                release_window_file=meta['release_window_file'],
                input_port_local_um=meta['input_port_um'],output_port_local_um=meta['output_port_um']))
    policies=finalize_coupler_windows(ly,entries)

    # C/D: one straight only; derive complete couplers from the updated A/B cells.
    for one in (False,True):
        for width in WG_WIDTHS:
            source_id=base.case_id(width)+(f'_OUT{OUTPUT_EXTENSION_LENGTH:g}' if one else '')
            require(ly.cell(source_id) is not None,'Enable reference devices for the straight two-coupler set')
            source_entry=next(e for e in entries if e['case_id']==source_id)
            rr=make_single_straight(regions(ly,ly.cell(source_id)),width,one)
            ident=f'D1_SINGLE_L{STRAIGHT_LENGTH:g}_W{round(width*1000)}'+(f'_1T_OUT{OUTPUT_EXTENSION_LENGTH:g}' if one else '_2T')
            entry=copy.deepcopy(source_entry)
            entry.update(mode='single_straight_one_coupler' if one else 'single_straight_two_couplers',
                topology='single_straight',set_name='single_straight_one_coupler' if one else 'single_straight_two_couplers',
                straight_section_count=1,bend_count=0,
                path_length_excluding_tapers_um=STRAIGHT_LENGTH+(OUTPUT_EXTENSION_LENGTH if one else 0.),
                input_port_local_um=[-STRAIGHT_LENGTH/2,0.],output_port_local_um=[STRAIGHT_LENGTH/2,0.],
                label=f'D1_S{STRAIGHT_LENGTH:g}_W{width:.2f}_T80'+(f'_1T_OUT{OUTPUT_EXTENSION_LENGTH:g}' if one else '_2T'),
                copied_from=source_id)
            place(ident,rr,entry)
            policies[ident]=copy.deepcopy(policies[source_id])
            for end in policies[ident]['ends']:end['center_y_um']=0.

    # Retain original outline width; extend height to contain the added rows.
    xmin=math.floor((bounds.left*DBU-FRAME_MARGIN)/100)*100
    xmax=math.ceil((bounds.right*DBU+FRAME_MARGIN)/100)*100
    ymin=math.floor((bounds.bottom*DBU-FRAME_MARGIN)/100)*100
    ymax=math.ceil((bounds.top*DBU+FRAME_MARGIN)/100)*100
    model=base.configured_design(WG_WIDTHS[0],dict(CHIP_WIDTH=xmax-xmin,CHIP_MIN_HEIGHT=max(MIN_FRAME_HEIGHT,ymax-ymin)))
    actual_height=max(MIN_FRAME_HEIGHT,ymax-ymin,(len(entries)-1)*model.DEVICE_Y_PITCH+model.CHIP_Y_ALLOWANCE)
    ymax=ymin+actual_height
    frame=gf.Component();model.add_chip_outline_and_markers(frame,len(entries));frame.flatten()
    fc=insert_regions(ly,'FRAME_AND_MARKERS',regions(frame.kcl.layout,frame.kdb_cell))
    ft=kdb.Trans(round((xmin+xmax)/2/DBU),round((ymin+ymax)/2/DBU))
    top.insert(kdb.CellInstArray(fc.cell_index(),ft))
    layer_policy=finalize_output_layers(ly,fc)
    beam_policy=add_local_taper_beams(ly,entries,policies)
    snapshots={e['case_id']:regions(ly,ly.cell(e['case_id'])) for e in entries}
    placed=[{n:r.transformed(kdb.Trans(round(e['origin_um'][0]/DBU),round(e['origin_um'][1]/DBU)))
             for (n,d),r in snapshots[e['case_id']].items()} for e in entries]
    for i,a in enumerate(placed):
        for j,b in enumerate(placed[i+1:],i+1):
            area=sum((a.get(n,kdb.Region())&b.get(n,kdb.Region())).area()*DBU**2 for n in (2,10,11,13,14,15))
            check(checks,f'rows_{i}_{j}_no_overlap',area==0,area_um2=area)
    for i,t in enumerate(texts):
        check(checks,f'label{i}_clear',all((t&r).is_empty() for p in placed for r in p.values()))
        check(checks,f'label{i}_unique',all((t&s).is_empty() for s in texts[i+1:]))
    fr=regions(ly,fc)
    for n in (10,15):
        marks=fr[(n,0)].transformed(ft)
        check(checks,f'L{n}_markers_clear',all((marks&r).is_empty() for p in placed for r in p.values()))
    write_layout(ly,OUTPUT)
    final,finaltop=read_gds(OUTPUT)
    insts={i.cell.name:i for i in finaltop.each_inst()}
    check(checks,'twelve_devices_and_labels',len(entries)==12 and len(insts)==25)
    check(checks,'excluded_layers_absent',not ({final.get_info(i).layer for i in final.layer_indexes()}&set(REMOVED_OUTPUT_LAYERS)))
    for e in entries:
        rr=regions(final,insts[e['case_id']].cell)
        check(checks,e['case_id']+'_saved_layers',set(rr)==set(snapshots[e['case_id']]))
        for key,r in snapshots[e['case_id']].items():
            check(checks,e['case_id']+f'_L{key[0]}_export_exact',(r^rr[key]).is_empty())
    workers=[json.loads((WORK/(base.case_id(w)+'.json')).read_text()) for w in WG_WIDTHS]
    report=dict(status='PASS',devices=entries,device_count=len(entries),set_count=4,
        width_sweep_um=list(WG_WIDTHS),meander_straight_length_um=STRAIGHT_LENGTH,
        single_straight_length_um=STRAIGHT_LENGTH,extension_length_um=OUTPUT_EXTENSION_LENGTH,
        row_pitch_um=ROW_PITCH,left_bend_x_um=LEFT_BEND_X,frame_bounds_um=[xmin,ymin,xmax,ymax],
        layer2_start_width_um=LAYER2_START_WIDTH,layer14_start_width_um=LAYER14_START_WIDTH,
        layer2_window_margin_um=LAYER2_WINDOW_MARGIN,layer14_window_margin_um=LAYER14_WINDOW_MARGIN,
        layer2_start_offset_um=LAYER2_START_OFFSET,layer14_start_offset_um=LAYER14_START_OFFSET,
        final_layer_policy=layer_policy,coupler_layer_policy=policies,local_taper_beams=beam_policy,checks=checks,workers=workers,
        check_count=len(checks)+sum(len(w['checks']) for w in workers),
        original_gds_sha256=hashlib.sha256(BASE_GDS.read_bytes()).hexdigest(),
        gds_sha256=hashlib.sha256(OUTPUT.read_bytes()).hexdigest(),
        optical_or_mechanical_simulation_performed=False)
    OUTPUT.with_suffix('.json').write_text(json.dumps(report,indent=2))
    with OUTPUT.with_suffix('.csv').open('w',newline='') as fp:
        fields=['case_id','set_name','width_um','origin_um','straight_section_count','bend_count','added_length_um','path_length_excluding_tapers_um','label']
        writer=csv.DictWriter(fp,fieldnames=fields);writer.writeheader()
        writer.writerows({k:e[k] for k in fields} for e in entries)
    print('PASS',report['check_count'],'construction/export checks;',OUTPUT,flush=True)


# %% [05] Views rendered directly from the exported GDS

def render(ly,top,result):
    """Render actual exported polygons using the independent reader."""
    from validate_2611_four_sets import render as render_saved
    render_saved(ROOT)


def first_stage_mask(sx, cy, outward, start_offset, start_width, end_width, length=40.0):
    """Smooth full-width release envelope, with exact widths at both ends."""
    import klayout.db as kdb
    import numpy as np
    require(0 <= start_offset < length, 'Invalid first-stage start offset')
    require(start_width >= end_width > 0, 'Mask must narrow monotonically')
    s = np.linspace(start_offset, length, 401)
    t = (s - start_offset) / (length - start_offset)
    f = t**3 * (10 - 15*t + 6*t*t)
    w = start_width + (end_width-start_width)*f
    upper = [(sx+outward*float(a), cy+float(b)/2) for a,b in zip(s,w)]
    lower = [(sx+outward*float(a), cy-float(b)/2) for a,b in zip(s[::-1],w[::-1])]
    return kdb.Region(kdb.Polygon([kdb.Point(round(x/DBU),round(y/DBU)) for x,y in upper+lower]))

def make_single_straight(source_regions, width, single_coupler):
    """Copy complete approved terminal assemblies and replace only routing.

    New ports are (-3000,0) and (+3000,0). Move the complete upper-left end
    down 200 um and lower-right end up 200 um; no optical taper is redrawn.
    The one-coupler version retains the 2000 um taper-free output extension.
    """
    import klayout.db as kdb
    left = box((-20000,25,-STRAIGHT_LENGTH/2,20000))
    right = box((STRAIGHT_LENGTH/2,-20000,20000,-25))
    down = kdb.Trans(0,round(-200/DBU))
    up = kdb.Trans(0,round(200/DBU))
    out={}
    for key,r in source_regions.items():
        if key[0] in (*REMOVED_OUTPUT_LAYERS, *CLEARED_SOURCE_PNC_LAYERS):
            continue
        out[key]=((r&left).transformed(down) | (r&right).transformed(up)).merged()
    ridge=box((-STRAIGHT_LENGTH/2,-width/2,STRAIGHT_LENGTH/2,width/2))
    out[(15,0)]=(out.get((15,0),kdb.Region())|ridge).merged()
    return {key:r for key,r in out.items() if not r.is_empty()}

# %% [06] Entry point: one isolated process per width

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker-width',type=float,help=argparse.SUPPRESS)
    args=parser.parse_args();WORK.mkdir(parents=True,exist_ok=True)
    if args.worker_width is not None:
        build_case(args.worker_width);return
    require(len(WG_WIDTHS)>0 and len(set(WG_WIDTHS))==len(WG_WIDTHS),'Widths must be distinct')
    require(all(math.isfinite(x) and x>0 for x in (*WG_WIDTHS,STRAIGHT_LENGTH,OUTPUT_EXTENSION_LENGTH,PNC_PITCH,ROW_PITCH)),
            'All dimensions must be finite and positive')
    for width in WG_WIDTHS:
        run=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker-width',str(width)],
                           cwd=ROOT,text=True,capture_output=True,timeout=300)
        (WORK/(f'W{width:g}.log')).write_text(run.stdout+'\n'+run.stderr)
        require(run.returncode==0,f'Width {width} failed:\n'+run.stdout[-1000:]+run.stderr[-6000:])
        print(run.stdout,flush=True)
    combine()
    from validate_2611_four_sets import validate, render as render_saved
    validate(ROOT)
    render_saved(ROOT)


if __name__=='__main__':
    main()
