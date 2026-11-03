"""Build two 3x3 CAD sweeps in one GDS, without sharing cached device cells.

Run from Cading: python 2611_build_combined_sweep.py
Edit only Section 01 for the sweep. Source geometry remains in the two
numbered Design-1/Design-2 scripts. Each case evaluates a parameter-only
copy of its source in a FRESH Python process, before defaults/derived values
are bound. No source file is rewritten, and no device is clipped or spliced.
"""
from __future__ import annotations

# %% [00] Imports and paths
import argparse
import ast
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import types

ROOT = Path(__file__).resolve().parent

# %% [01] USER PARAMETERS (all dimensions in um)
STRAIGHT_LENGTH = 1400.0  # EACH of three straight segments
BEND_PNC_PERIODS = (2.0, 3.0, 5.0)
POST_PNC_TAPER_LENGTHS = (80.0, 100.0, 120.0)  # total Layer-10 length; excludes PNC
DESIGN_X_OFFSET = 5000.0  # Design 2 relative to Design 1
ROW_PITCH = 600.0
FIRST_ROW_Y = -700.0
LEFT_BEND_X = -900.0  # leftmost L11 PNC edge in the Design-1 column
# Preserve the existing relative stage lengths, NOT absolute stage lengths.
# Change these explicit ratios to choose a different taper-length allocation.
STAGE_RATIOS = {1: (1.0, 1.0), 2: (2.0, 10.0, 1.0)}
SOURCES = {1: '2611_tapered_yue_dual_resist.py', 2: '2611_tapered_yue_design2.py'}
OUTPUT = ROOT / 'gds/2611_design1_design2_L1400_18device_sweep.gds'
MANIFEST = OUTPUT.with_suffix('.json')
CSV_PATH = OUTPUT.with_suffix('.csv')
WORK = ROOT / 'validation_2611_combined_sweep'
LABEL_SIZE = 5.0
# A single encompassing outline/marker set, not 18 overlapping chip frames.
FRAME_MARGIN = 500.0
MIN_FRAME_HEIGHT = 10000.0
DBU = 0.001
FAB_LAYERS = (2, 10, 11, 13, 14, 15, 16)


def require(condition, message):
    if not condition:
        raise ValueError(message)


# %% [02] Cartesian-product cases and safe parameter loading

def cases():
    require(STRAIGHT_LENGTH > 0 and ROW_PITCH > 0 and DESIGN_X_OFFSET > 0,
            'Lengths and spacing must be positive')
    rows = [(p, t) for p in BEND_PNC_PERIODS for t in POST_PNC_TAPER_LENGTHS]
    result = []
    for design in (1, 2):
        for row, (period, total) in enumerate(rows):
            ratio = STAGE_RATIOS[design]
            require(period > 0 and total > 0 and min(ratio) > 0, 'Invalid sweep setting')
            stages = [total * r / sum(ratio) for r in ratio]
            result.append(dict(design=design, row=row, bend_period_um=period,
                               taper_total_um=total, stage_lengths_um=stages,
                               straight_um=STRAIGHT_LENGTH,
                               case_id=f'D{design}_{row+1:02d}_B{period:g}_T{total:g}'))
    return result


def configured_module(case, extra=None):
    """Replace only named TOP-LEVEL parameter assignments in a memory copy.

    All geometry functions and expressions are left untouched. This is not
    runtime mutation of a loaded module: defaults and decorators see the new
    values on their very first evaluation. The caller isolates each case in
    a new process so GDSFactory caches cannot cross-contaminate the sweep.
    """
    design = case['design']
    path = ROOT / SOURCES[design]
    values = dict(MEANDER_STRAIGHT_LENGTHS=(case['straight_um'],),
                  REVIEW_SUPPORT_BEND_PERIOD=case['bend_period_um'])
    stages = case['stage_lengths_um']
    if design == 1:
        values.update(RIDGE_TAPER_LENGTH=stages[0], SLAB_TAPER_LENGTH=stages[1])
    else:
        values.update(TAPER_STAGE1_LENGTH=stages[0], TAPER_STAGE2_LENGTH=stages[1],
                      TAPER_STAGE3_LENGTH=stages[2])
    values.update(extra or {})
    tree = ast.parse(path.read_text(), filename=str(path))
    found = {key: 0 for key in values}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name) and target.id in values:
                node.value = ast.parse(repr(values[target.id]), mode='eval').body
                found[target.id] += 1
    require(all(v == 1 for v in found.values()), f'Parameter assignment changed: {found}')
    ast.fix_missing_locations(tree)
    module = types.ModuleType('sweep_' + case['case_id'])
    module.__file__ = str(path)
    sys.modules[module.__name__] = module
    exec(compile(tree, str(path), 'exec'), module.__dict__)
    return module


# %% [03] Build and measure one COMPLETE all-layer device (child process)

def worker(case):
    import gdsfactory as gf
    import klayout.db as kdb
    import numpy as np
    m = configured_module(case)
    m.validate_settings()
    capture = {}
    shield_name = 'phononic_meander_shield' if case['design'] == 1 else 'design2_phononic_meander_shield'
    original = getattr(m, shield_name)

    def measured_shield(path, **kwargs):
        shield = original(path=path, **kwargs)
        capture['points'] = np.asarray(path.points, dtype=float)
        capture['positions'] = list(shield.info['support_positions'])
        capture['bend_pitch'] = float(kwargs['bend_pitch'])
        return shield

    measured_shield.__name__ = original.__name__
    setattr(m, shield_name, measured_shield)
    d = m.make_device(case['straight_um'])
    d.flatten()
    d.name = case['case_id']
    dbu = float(d.kcl.dbu)
    require(abs(dbu - DBU) < 1e-12, 'Unexpected database grid')
    regions = {tuple(layer): kdb.Region(d.kdb_cell.begin_shapes_rec(gf.get_layer(layer))).merged()
               for layer in d.layers}
    required = {(n, 0) for n in (*FAB_LAYERS, 51)}
    require(set(regions) == required, f'Unexpected layers: {set(regions)}')
    ticks = m.alignment_tip_region(d.info['tip_centers'], d.info['connection_orientations'], dbu)
    optical = (regions[(15, 0)] - ticks).merged()
    require(optical.count() == 1, 'Optical ridge is disconnected')
    require(ticks.count() == 4 and (ticks - regions[(15, 0)]).is_empty(), 'Missing alignment tips')
    require((ticks & regions[(10, 0)]).is_empty(), 'Alignment tips remain on L10')
    require(regions[(16, 0)].count() == 1, 'Protection is disconnected')
    for layer in (14, 16):
        require((optical - regions[(layer, 0)]).area() * dbu**2 < 0.005,
                f'Layer {layer} does not cover the optical ridge')

    # Measure all five actual path sections, not metadata alone.
    seg = np.linalg.norm(np.diff(capture['points'], axis=0), axis=1)
    flags = seg <= max(1.0, 4.0 * float(np.median(seg)))
    edges = np.r_[0, np.nonzero(flags[1:] != flags[:-1])[0] + 1, len(seg)]
    cumulative = np.r_[0., np.cumsum(seg)]
    sections = [(int(a), int(b), bool(flags[a])) for a, b in zip(edges[:-1], edges[1:])]
    require([b for _, _, b in sections] == [False, True, False, True, False], 'Wrong path topology')
    actual_lengths = [float(seg[a:b].sum()) for a, b, _ in sections]
    require(np.allclose(actual_lengths[::2], case['straight_um'], atol=1e-6), 'Unequal straights')
    require(abs(sum(actual_lengths) - float(d.info['total_meander_length_um'])) < .003,
            'Incorrect total-length label')
    bend_measurements = []
    for a, b, is_bend in sections:
        if is_bend:
            positions = sorted(p for p in capture['positions']
                               if cumulative[a] + 1e-6 < p < cumulative[b] - 1e-6)
            diffs = np.diff(positions)
            require(len(diffs) > 1 and np.min(diffs) >= case['bend_period_um'] - 1e-6,
                    'Incorrect bend pitch')
            require(np.max(diffs) < case['bend_period_um'] * 1.03, 'Bend pitch not applied')
            bend_measurements.append(dict(pairs=len(positions), fitted_pitch_um=float(np.mean(diffs))))

    # Check the physical taper profile and stage boundaries in local coordinates.
    taper_name = 'two_layer_tapered_fiber' if case['design'] == 1 else 'design2_two_layer_tapered_fiber'
    taper = getattr(m, taper_name)(add_support_tips=False, add_phononic_taper_support=True)
    taper_regions = {n: kdb.Region(taper.kdb_cell.begin_shapes_rec(gf.get_layer((n, 0)))).merged()
                     for n in (10, 15)}
    stages = case['stage_lengths_um']
    widths10 = [3., 1., .2] if case['design'] == 1 else [3., 1.65, .98, .2]
    widths15 = [1., .1] if case['design'] == 1 else [1., .5, .2]
    lengths15 = stages[:1] if case['design'] == 1 else stages[:2]
    samples = []
    for n, lengths, widths in ((10, stages, widths10), (15, lengths15, widths15)):
        xs = float(m.PNC_TAPER_LENGTH) + np.r_[0., np.cumsum(lengths)]
        require(abs(taper_regions[n].bbox().right * dbu - xs[-1]) < .0031,
                f'Layer {n} endpoint wrong')
        require(taper_regions[n].count() == 1, f'Layer {n} taper junction disconnected')
        for j, (x, width) in enumerate(zip(xs, widths)):
            ix = round((x + (.002 if j == 0 else -.002)) / dbu)
            rr = (taper_regions[n] & kdb.Region(kdb.Box(ix, -5000, ix + 1, 5000))).merged()
            require(rr.count() == 1, 'Invalid taper cross-section')
            box = rr.bbox()
            measured = box.height() * dbu
            require(abs(measured - width) <= .004, f'Layer {n} width {measured} != {width}')
            samples.append(dict(layer=n, distance_um=float(x-m.PNC_TAPER_LENGTH),
                                measured_width_um=measured, target_width_um=width))
    if case['design'] == 2:
        require(m.RIDGE_CONVERTER_Y_SHIFT == .840, 'Ridge offset changed')
        xend = (m.PNC_TAPER_LENGTH + sum(lengths15))
        ix = round((xend-.002)/dbu)
        rr = taper_regions[15] & kdb.Region(kdb.Box(ix, -5000, ix+1, 5000))
        require(abs(rr.bbox().center().y*dbu - .840) < .002, 'Wrong terminal ridge shift')

    out = WORK / 'cases'
    out.mkdir(parents=True, exist_ok=True)
    options = kdb.SaveLayoutOptions()
    options.gds2_max_vertex_count = 2000
    d.write_gds(out / (case['case_id']+'.gds'), save_options=options, with_metadata=False)
    data = dict(case, source=SOURCES[case['design']],
                source_sha256=hashlib.sha256((ROOT/SOURCES[case['design']]).read_bytes()).hexdigest(),
                left_pnc_x_um=regions[(11,0)].bbox().left*dbu,
                measured_straights_um=actual_lengths[::2], measured_bends_um=actual_lengths[1::2],
                meander_length_um=sum(actual_lengths), bend_pnc_measurements=bend_measurements,
                straight_pnc_period_um=float(m.REVIEW_SUPPORT_ANCHOR_PERIOD),
                layer10_widths_um=widths10, layer15_widths_um=widths15,
                layer15_lengths_um=lengths15, taper_samples=samples,
                tip_centers=d.info['tip_centers'], orientations=d.info['connection_orientations'],
                status='PASS')
    (out/(case['case_id']+'.json')).write_text(json.dumps(data, indent=2))
    print(case['case_id'], 'PASS', flush=True)


# %% [04] Combine owned cells by translation; preserve every layer

def read_gds(path):
    import klayout.db as kdb
    ly = kdb.Layout()
    ly.read(str(path))
    require(abs(ly.dbu-DBU) < 1e-12, 'Grid mismatch during GDS merge')
    return ly, ly.top_cell()


def copy_flat(dst_layout, src_layout, src_cell, name):
    import klayout.db as kdb
    cell = dst_layout.create_cell(name)
    for idx in src_layout.layer_indexes():
        info = src_layout.get_info(idx)
        region = kdb.Region(src_cell.begin_shapes_rec(idx))
        if not region.is_empty():
            cell.shapes(dst_layout.layer(info)).insert(region)
    return cell


def cell_regions(ly, cell):
    import klayout.db as kdb
    return {(ly.get_info(i).layer, ly.get_info(i).datatype):
            kdb.Region(cell.begin_shapes_rec(i)).merged() for i in ly.layer_indexes()
            if not kdb.Region(cell.begin_shapes_rec(i)).is_empty()}


def combine(entries):
    import gdsfactory as gf
    import klayout.db as kdb
    gf.gpdk.PDK.activate()
    ly = kdb.Layout()
    ly.dbu = DBU
    top = ly.create_cell('D1_D2_L1400_18_SWEEP')
    placed = []
    locals_by_case = {}
    fabrication_bbox = kdb.Box()
    labels = gf.Component('SWEEP_CASE_LABELS')
    for case in entries:
        ident = case['case_id']
        meta = json.loads((WORK/'cases'/f'{ident}.json').read_text())
        source, src_cell = read_gds(WORK/'cases'/f'{ident}.gds')
        cell = copy_flat(ly, source, src_cell, ident)
        src = cell_regions(source, src_cell)
        copied = cell_regions(ly, cell)
        require(set(src) == set(copied) and all((src[l]^copied[l]).is_empty() for l in src),
                f'{ident}: lost geometry while merging')
        locals_by_case[ident] = copied
        x = LEFT_BEND_X + (case['design']-1)*DESIGN_X_OFFSET - meta['left_pnc_x_um']
        y = FIRST_ROW_Y + case['row']*ROW_PITCH
        trans = kdb.Trans(round(x/DBU), round(y/DBU))
        top.insert(kdb.CellInstArray(cell.cell_index(), trans))
        meta['origin_um'] = [trans.disp.x*DBU, trans.disp.y*DBU]
        meta['left_bend_reference_x_um'] = meta['left_pnc_x_um'] + meta['origin_um'][0]
        physical = {n: copied[(n,0)].transformed(trans) for n in FAB_LAYERS}
        for r in physical.values():
            fabrication_bbox += r.bbox()
        placed.append((meta, physical))
        # Keep the existing length label; add explicit sweep ID above each row.
        labels.add_ref(gf.components.text(
            text=f'{ident}_L{STRAIGHT_LENGTH:g}', size=LABEL_SIZE, layer=(2,0)
        )).move((LEFT_BEND_X+(case['design']-1)*DESIGN_X_OFFSET+200, y+250))

    # Only Layer 11 may change when only the bend period is varied.
    for design in (1,2):
        for total in POST_PNC_TAPER_LENGTHS:
            family=[e for e in entries if e['design']==design and e['taper_total_um']==total]
            ref=locals_by_case[family[0]['case_id']]
            counts=[]
            for e in family:
                other=locals_by_case[e['case_id']]
                require(all((ref[l]^other[l]).is_empty() for l in ref if l!=(11,0)),
                        'Bend-period sweep changed a non-PNC layer')
                counts.append(other[(11,0)].count())
            require(counts[0]>counts[1]>counts[2], 'PNC populations do not follow period sweep')

    overlaps=[]
    for i,(a, ar) in enumerate(placed):
        for b,br in placed[i+1:]:
            for n in FAB_LAYERS:
                if ar[n].bbox().overlaps(br[n].bbox()):
                    area=(ar[n]&br[n]).area()*DBU**2
                    if area > .001:
                        overlaps.append(dict(a=a['case_id'],b=b['case_id'],layer=n,area_um2=area))
    require(not overlaps, f'Fabrication-mask collision: {overlaps}')
    for i in range(len(entries)//2):
        a=placed[i][0]; b=placed[i+len(entries)//2][0]
        require(abs(b['origin_um'][0]-a['origin_um'][0]-DESIGN_X_OFFSET)<.0011,
                'Design column displacement is not 5000 um')
        require(b['origin_um'][1]==a['origin_um'][1], 'Design columns have unequal row coordinates')

    # Enclose fabrication masks, NOT long nonfabrication fiber-guide polygons.
    xmin=math.floor((fabrication_bbox.left*DBU-FRAME_MARGIN)/100)*100
    xmax=math.ceil((fabrication_bbox.right*DBU+FRAME_MARGIN)/100)*100
    ymin=min(-MIN_FRAME_HEIGHT/2, math.floor((fabrication_bbox.bottom*DBU-FRAME_MARGIN)/100)*100)
    ymax=max(MIN_FRAME_HEIGHT/2, math.ceil((fabrication_bbox.top*DBU+FRAME_MARGIN)/100)*100)
    m=configured_module(entries[0], dict(CHIP_WIDTH=xmax-xmin, CHIP_MIN_HEIGHT=ymax-ymin))
    frame=gf.Component('SWEEP_FRAME_SOURCE')
    m.add_chip_outline_and_markers(frame, len(entries)//2)
    frame.add_ref(labels).move((-(xmax+xmin)/2, -(ymax+ymin)/2))
    frame.flatten()
    frame.write_gds(WORK/'frame.gds', with_metadata=False)
    f_ly,f_cell=read_gds(WORK/'frame.gds')
    f=copy_flat(ly,f_ly,f_cell,'SWEEP_FRAME_AND_LABELS')
    top.insert(kdb.CellInstArray(f.cell_index(),kdb.Trans(round((xmax+xmin)/2/DBU),round((ymax+ymin)/2/DBU))))
    OUTPUT.parent.mkdir(exist_ok=True)
    options=kdb.SaveLayoutOptions();options.gds2_max_vertex_count=2000
    ly.write(str(OUTPUT),options)
    exported, ex_top=read_gds(OUTPUT)
    actual=[i for i in ex_top.each_inst() if i.cell.name.startswith(('D1_','D2_'))]
    require(len(actual)==len(entries)==18,'Incorrect exported device count')
    for inst in actual:
        a=locals_by_case[inst.cell.name]; b=cell_regions(exported,inst.cell)
        require(set(a)==set(b) and all((a[l]^b[l]).area()*DBU**2 < .005 for l in a),
                f'{inst.cell.name}: GDS export altered a layer')
    manifest=dict(status='PASS', devices=[a for a,_ in placed], device_count=len(entries),
                  row_pitch_um=ROW_PITCH, design_x_offset_um=DESIGN_X_OFFSET,
                  taper_allocation='proportional stage scaling', stage_ratios=STAGE_RATIOS,
                  taper_length_definition='full Layer-10 taper after the unchanged 45 um PNC',
                  frame_bounds_um=[xmin,ymin,xmax,ymax], fabrication_mask_overlaps=overlaps,
                  notes=['Bend period is the existing nominal parameter. Interior spacings are fitted to the arc; measured values are listed per case.',
                         'Layer 51 is a nonfabrication fiber guide and may cross columns. No fabrication layer is clipped.'])
    MANIFEST.write_text(json.dumps(manifest,indent=2))
    fields=['case_id','design','row','bend_period_um','taper_total_um','stage_lengths_um',
            'layer15_lengths_um','straight_um','meander_length_um','origin_um','left_bend_reference_x_um']
    with CSV_PATH.open('w',newline='') as fp:
        writer=csv.DictWriter(fp,fieldnames=fields);writer.writeheader()
        for meta,_ in placed:
            writer.writerow({k:meta[k] for k in fields})
    render(exported, ex_top, manifest)
    return manifest


# %% [05] Preview of the ACTUAL exported GDS

def render(ly, top, manifest):
    import klayout.db as kdb
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    fig,ax=plt.subplots(figsize=(13,10))
    for n in (14,11,15):
        lines=[]
        for poly in kdb.Region(top.begin_shapes_rec(ly.layer(n,0))).merged().each():
            pts=[(p.x*DBU,p.y*DBU) for p in poly.each_point_hull()]
            if pts: lines.append(pts+[pts[0]])
        ax.add_collection(LineCollection(lines,linewidths=.35,label=f'Layer {n}'))
    for e in manifest['devices']:
        x=e['left_bend_reference_x_um'];y=e['origin_um'][1]
        ax.text(x+200,y+255,e['case_id'],fontsize=7)
    ax.set(xlim=(-3200,8200),ylim=(-1200,4550),xlabel='X (um)',ylabel='Y (um)',
           title='Design 1 (left) / Design 2 (right): 9 cases each, 5000 um column offset')
    ax.set_aspect('equal');fig.tight_layout();fig.savefig(WORK/'combined_sweep.png',dpi=200)
    plt.close(fig)


# %% [06] Run sequential isolated workers, then combine and validate

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker',help=argparse.SUPPRESS)
    args=parser.parse_args()
    if args.worker:
        worker(json.loads(args.worker));return
    WORK.mkdir(exist_ok=True)
    entries=cases()
    for case in entries:
        result=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker',json.dumps(case)],
                              cwd=ROOT, text=True,capture_output=True,timeout=300)
        (WORK/(case['case_id']+'.log')).write_text(result.stdout+'\n'+result.stderr)
        if result.returncode:
            raise RuntimeError(f"{case['case_id']} failed; see {WORK/case['case_id']}.log\n{result.stderr[-5000:]}")
        print(result.stdout.strip(),flush=True)
    report=combine(entries)
    print(f"PASS: {report['device_count']} devices. Wrote {OUTPUT}",flush=True)


if __name__=='__main__':
    main()
