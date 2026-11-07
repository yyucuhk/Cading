"""Separate Design-1 width sweep; run this file from the Cading root.

All dimensions are um. Each of THREE meander straights is 6000 um.
Layer 15 optical ridge widths: 1.62, 1.64, 1.68. Both meander straight
and bend PNC pitches are 10 um. Existing Design-1 tapers are preserved:
L10 5 -> 1 -> 0.2 over 40/40; L15 W -> 0.1 over 40 after the shortened PNC.
The PNC ends at the 5 um width point of the original 45 um reference profile.
This builder does not rewrite the source or the other GDS outputs.

Bend pitch is exact CENTERLINE ARC-LENGTH pitch between interior PNC pairs.
The noninteger arc-length remainder is split equally at the junctions;
shared straight/bend boundary pairs are retained. This is not edge clearance
or constant spacing along the inner/outer offset curves.

Each width runs in a fresh process, then complete all-layer devices are
translated into one column. Checks and previews use the actual polygons.
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
SOURCE = ROOT / '2611_tapered_yue_dual_resist.py'

# %% [01] USER PARAMETERS - EDIT HERE (um)
STRAIGHT_LENGTH = 6000.0  # EACH of three straights, not total meander length
WG_WIDTHS = (1.62, 1.64, 1.68)  # Layer 15 ridge; ordered bottom to top
STRAIGHT_PNC_PERIOD = 10.0
BEND_PNC_PERIOD = 10.0
ROW_PITCH = 600.0
FIRST_ROW_Y = -700.0
LEFT_BEND_X = -900.0
LABEL_SIZE = 5.0
FRAME_MARGIN = 500.0
MIN_FRAME_HEIGHT = 4000.0
STEM = '2611_design1_L6000_P10_width_sweep'
OUTPUT = ROOT / 'gds' / (STEM + '.gds')
WORK = ROOT / 'validation_2611_L6000_width_sweep'
DBU = 0.001
FAB_LAYERS = (2, 10, 11, 13, 14, 15, 16)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def case_id(width):
    return f'D1_L{STRAIGHT_LENGTH:g}_W{round(width * 1000):d}_P{BEND_PNC_PERIOD:g}'


def check(report, name, ok, **details):
    report.append(dict(name=name, passed=bool(ok), **details))
    require(ok, name + ': ' + str(details))


# %% [02] Load an isolated Design-1 parameterization, before defaults bind

def exact_bend_positions(start, end, pitch):
    import numpy as np
    count = int(np.floor((end-start)/pitch + 1e-9)) + 1
    return (start+end)/2 + (np.arange(count)-(count-1)/2)*pitch


def configured_design(width, extra=None):
    values = dict(MEANDER_STRAIGHT_LENGTHS=(STRAIGHT_LENGTH,),
                  RIDGE_MEANDER_WIDTH=float(width),
                  MEANDER_STRAIGHT_PNC_PERIOD=STRAIGHT_PNC_PERIOD,
                  REVIEW_SUPPORT_BEND_PERIOD=BEND_PNC_PERIOD)
    values.update(extra or {})
    tree = ast.parse(SOURCE.read_text(), filename=str(SOURCE))
    found = {key: 0 for key in values}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name) and target.id in values:
                node.value = ast.parse(repr(values[target.id]), mode='eval').body
                found[target.id] += 1
    require(all(n == 1 for n in found.values()), f'Source parameters changed: {found}')
    label_count = sampler_count = 0
    for function in tree.body:
        if not isinstance(function, ast.FunctionDef):
            continue
        if function.name == 'make_device':
            for node in ast.walk(function):
                if isinstance(node, ast.Call) and ast.unparse(node.func) == 'gf.components.text':
                    kw = next(k for k in node.keywords if k.arg == 'layer')
                    require(ast.unparse(kw.value) == 'MARKER_LAYER', 'Unexpected label layer')
                    kw.value = ast.parse('(998, 0)', mode='eval').body
                    label_count += 1
        if function.name == 'phononic_meander_shield':
            for node in ast.walk(function):
                if (isinstance(node, ast.Assign) and len(node.targets) == 1
                    and isinstance(node.targets[0], ast.Name)
                    and node.targets[0].id == 'path_positions'
                    and isinstance(node.value, ast.IfExp)
                    and 'np.linspace' in ast.unparse(node.value)):
                    node.value = ast.parse(
                        '_exact_bend_positions(usable_start, usable_end, requested_pitch)',
                        mode='eval').body
                    sampler_count += 1
    require(label_count == sampler_count == 1, 'Unexpected label or bend sampler structure')
    ast.fix_missing_locations(tree)
    module = types.ModuleType('long_' + case_id(width))
    module.__file__ = str(SOURCE)
    module._exact_bend_positions = exact_bend_positions
    sys.modules[module.__name__] = module
    exec(compile(tree, str(SOURCE), 'exec'), module.__dict__)
    return module


# %% [03] Polygon utilities and measured cross-sections

def regions(layout, cell):
    import klayout.db as kdb
    result = {}
    for idx in layout.layer_indexes():
        r = kdb.Region(cell.begin_shapes_rec(idx)).merged()
        if not r.is_empty():
            info = layout.get_info(idx)
            result[(info.layer, info.datatype)] = r
    return result


def box(coords):
    import klayout.db as kdb
    return kdb.Region(kdb.Box(*(round(v / DBU) for v in coords)))


def contains_point(region, point):
    x, y = point
    return not (region & box((x, y, x+DBU, y+DBU))).is_empty()


def read_gds(path):
    import klayout.db as kdb
    ly = kdb.Layout()
    ly.read(str(path))
    require(abs(ly.dbu - DBU) < 1e-12, 'Unexpected database grid')
    require(len(ly.top_cells()) == 1, 'Expected one top cell')
    return ly, ly.top_cell()


def copy_flat(dst, src, cell, name):
    out = dst.create_cell(name)
    for (n, d), r in regions(src, cell).items():
        out.shapes(dst.layer(n, d)).insert(r)
    return out


# %% [04] Build one complete all-layer device and verify every section

def worker(width):
    import gdsfactory as gf
    import klayout.db as kdb
    import numpy as np
    m = configured_design(width)
    m.validate_settings()
    capture = {}
    original = m.phononic_meander_shield

    def measured_shield(path, **kwargs):
        shield = original(path=path, **kwargs)
        capture['points'] = np.asarray(path.points, dtype=float)
        capture['positions'] = np.asarray(shield.info['support_positions'], dtype=float)
        return shield

    measured_shield.__name__ = original.__name__
    m.phononic_meander_shield = measured_shield
    device = m.make_device(STRAIGHT_LENGTH)
    device.flatten()
    require(not kdb.Region(device.kdb_cell.begin_shapes_rec(gf.get_layer((998,0)))).is_empty(),
            'Missing isolated original label')
    device.remove_layers([(998,0)], recursive=False)
    device.name = case_id(width)
    checks = []
    r = regions(device.kcl.layout, device.kdb_cell)
    check(checks, 'all_device_layers_present', set(r) == {(n,0) for n in (*FAB_LAYERS,51)})
    ticks = m.alignment_tip_region(device.info['tip_centers'],
                                  device.info['connection_orientations'], DBU)
    optical = (r[(15,0)] - ticks).merged()
    check(checks, 'ridge_is_one_connected_body', optical.count() == 1)
    check(checks, 'protection_connected', r[(16,0)].count() == 1)
    check(checks, 'four_alignment_ticks_only_on_layer15', ticks.count() == 4 and
          (ticks-r[(15,0)]).is_empty() and (ticks&r[(10,0)]).is_empty())
    for mask in (14,16):
        missing = (optical-r[(mask,0)]).area()*DBU**2
        check(checks, f'layer{mask}_covers_ridge', missing < .005, missing_um2=missing)

    pts = capture['points']
    vectors = np.diff(pts, axis=0)
    lengths = np.linalg.norm(vectors, axis=1)
    cumulative = np.r_[0., np.cumsum(lengths)]
    flags = lengths <= max(1., 4*np.median(lengths))
    cuts = np.r_[0, np.flatnonzero(flags[1:] != flags[:-1])+1, len(lengths)]
    sections = [(int(a),int(b),bool(flags[a])) for a,b in zip(cuts[:-1],cuts[1:])]
    check(checks, 'three_straights_two_bends', [v for _,_,v in sections] == [False,True,False,True,False])
    positions = np.sort(capture['positions'])
    check(checks, 'no_duplicate_PNC_pairs', np.min(np.diff(positions)) > .001)
    section_data = []

    def path_point(s):
        i = min(len(lengths)-1, max(0, int(np.searchsorted(cumulative,s,side='right')-1)))
        t = vectors[i]/lengths[i]
        return pts[i]+(s-cumulative[i])*t, t

    for index,(a,b,bend) in enumerate(sections):
        length = float(lengths[a:b].sum())
        pitch = BEND_PNC_PERIOD if bend else STRAIGHT_PNC_PERIOD
        interior = positions[(positions > cumulative[a]+1e-6) & (positions < cumulative[b]-1e-6)]
        check(checks, f'section{index}_exact_pitch', len(interior)>1 and
              np.allclose(np.diff(interior),pitch,rtol=0,atol=1e-7),
              pitch_um=pitch, interior_pairs=len(interior))
        clearances = [float(interior[0]-cumulative[a]), float(cumulative[b]-interior[-1])]
        if bend:
            check(checks, f'bend{index}_symmetric_junction_clearance',
                  abs(clearances[0]-clearances[1]) < 1e-7 and min(clearances) >= pitch-1e-7,
                  junction_clearances_um=clearances)
        else:
            check(checks, f'straight{index}_6000um', abs(length-STRAIGHT_LENGTH)<1e-7, measured_um=length)
            xmin,xmax = sorted((float(pts[a,0]),float(pts[b,0])))
            cy = float(pts[a,1])
            local = (r[(11,0)] & box((xmin+1,cy-30,xmax-1,cy+30))).merged()
            centers = {1:[], -1:[]}
            for poly in local.each():
                bb = poly.bbox()
                sign = 1 if bb.center().y*DBU > cy else -1
                centers[sign].append((bb.left+bb.right)*DBU/2)
            expected = np.arange(xmin+pitch,xmax-1,pitch)
            for sign,values in centers.items():
                values = np.asarray(sorted(values))
                check(checks, f'straight{index}_side{sign}_exportable_PNC_centers',
                      len(values)==len(expected) and np.allclose(values,expected,rtol=0,atol=.0011),
                      count=len(values), expected_count=len(expected))
            sample = (optical & box(((xmin+xmax)/2,cy-3,(xmin+xmax)/2+DBU,cy+3))).merged()
            got = sample.bbox().height()*DBU
            check(checks, f'straight{index}_measured_ridge_width', sample.count()==1 and abs(got-width)<.0021,
                  measured_width_um=got, requested_width_um=width)
        # Confirm physical PNC necks exist at every measured path position,
        # on both sides, rather than trusting a parameter or count alone.
        all_present = True
        for s in interior:
            center,t = path_point(float(s))
            normal = np.array([-t[1],t[0]])
            for sign in (-1,1):
                all_present &= contains_point(r[(11,0)], center+sign*(m.S_BEND_WIDTH/2+.10)*normal)
        check(checks, f'section{index}_all_PNC_necks_present', all_present)
        if bend:
            correct_edges = True
            for s in np.linspace(cumulative[a]+10,cumulative[b]-10,11):
                center,t = path_point(float(s)); normal = np.array([-t[1],t[0]])
                for sign in (-1,1):
                    correct_edges &= contains_point(optical, center+sign*(width/2-.008)*normal)
                    correct_edges &= not contains_point(optical, center+sign*(width/2+.008)*normal)
            check(checks, f'bend{index}_ridge_width_at_11_sections', correct_edges)
        section_data.append(dict(kind='bend' if bend else 'straight', length_um=length,
                                 start_um=pts[a].tolist(),end_um=pts[b].tolist(),
                                 interior_PNC_pairs=len(interior),pitch_um=pitch,
                                 junction_clearances_um=clearances))
    total = float(lengths.sum())
    check(checks, 'correct_centerline_length', abs(total-device.info['total_meander_length_um']) < .003,
          measured_total_um=total)
    # Entire end assemblies, including all tethers/pads/windows, must match
    # under 180-degree rotation. There are no device text polygons here.
    roi = box((STRAIGHT_LENGTH/2+.001,-1000,STRAIGHT_LENGTH/2+11000,-25))
    for layer,reg in r.items():
        area = ((reg^reg.transformed(kdb.Trans(kdb.Trans.R180))) & roi).area()*DBU**2
        check(checks, f'layer{layer[0]}_complete_end_symmetry', area < .005, xor_um2=area)
    for side in (-1,1):
        rr = r if side==1 else {k:v.transformed(kdb.Trans(kdb.Trans.R180)) for k,v in r.items()}
        pnc_end = float(m.PNC_TAPER_LENGTH)
        first_end = pnc_end + float(m.RIDGE_TAPER_LENGTH)
        last_end = first_end + float(m.SLAB_TAPER_LENGTH)
        for layer,x,w in ((10,pnc_end+.002,float(m.PNC_TAPER_END_WIDTH)),
                          (10,first_end-.002,1.),(10,last_end-.002,.2),
                          (15,pnc_end+.002,width),(15,first_end-.002,.1)):
            x += STRAIGHT_LENGTH/2
            half_clip = max(2.0, w/2 + .02)
            section = (rr[(layer,0)] & box((x,-200-half_clip,x+DBU,-200+half_clip))).merged()
            measured = section.bbox().height()*DBU
            check(checks, f'end{side}_L{layer}_x{x:g}_taper_width', section.count()==1 and abs(measured-w)<.0041,
                  measured_width_um=measured, expected_width_um=w)
    check(checks, 'taper_PNC_stops_at_5um_and_matches_straight',
          0 < m.PNC_TAPER_LENGTH < m.PNC_TAPER_TRANSITION_LENGTH and
          m.PNC_TAPER_SUPPORT_CUTOFF_WIDTH == 5.0 and
          m.REVIEW_SUPPORT_ANCHOR_PERIOD == STRAIGHT_PNC_PERIOD)
    check(checks, 'Design1_taper_lengths_preserved', m.RIDGE_TAPER_LENGTH==40 and m.SLAB_TAPER_LENGTH==40)
    options = kdb.SaveLayoutOptions(); options.gds2_max_vertex_count=2000
    WORK.mkdir(exist_ok=True)
    device.write_gds(WORK/(device.name+'.gds'), save_options=options, with_metadata=False)
    exly,excell = read_gds(WORK/(device.name+'.gds'))
    exr = regions(exly,excell)
    for layer in r:
        area = (r[layer]^exr[layer]).area()*DBU**2
        check(checks, f'export_layer{layer[0]}_unchanged', area < .005, xor_um2=area)
    data = dict(case_id=device.name, width_um=width, straight_length_um=STRAIGHT_LENGTH,
                total_meander_um=total, section_measurements=section_data,
                left_pnc_x_um=r[(11,0)].bbox().left*DBU,
                tip_centers=device.info['tip_centers'],
                connection_centers=device.info['connection_centers'],
                straight_pnc_period_um=float(m.MEANDER_STRAIGHT_PNC_PERIOD),
                taper_pnc_period_um=float(m.REVIEW_SUPPORT_ANCHOR_PERIOD),
                pnc_taper_length_um=float(m.PNC_TAPER_LENGTH),
                pnc_cutoff_width_um=float(m.PNC_TAPER_SUPPORT_CUTOFF_WIDTH),
                pnc_reference_profile_length_um=float(m.PNC_TAPER_TRANSITION_LENGTH),
                connection_pnc_pitch_um=float(m.PNC_CONNECTION_PITCH),
                connection_pnc_periods=int(m.PNC_CONNECTION_PERIODS),
                connection_pnc_span_um=float(m.PNC_CONNECTION_PITCH*m.PNC_CONNECTION_PERIODS),
                connection_pnc_region="narrow_pnc_end",
                connection_pnc_start_um=float(m.PNC_TAPER_LENGTH-m.PNC_CONNECTION_PITCH*m.PNC_CONNECTION_PERIODS),
                connection_pnc_end_um=float(m.PNC_TAPER_LENGTH),
                post_pnc_taper_total_um=80., label=f'D1_W{width:.2f}_T80_P{BEND_PNC_PERIOD:g}',
                checks=checks,status='PASS')
    (WORK/(device.name+'.json')).write_text(json.dumps(data,indent=2))
    print(device.name, 'PASS', len(checks), 'checks', flush=True)


# %% [05] Complete-device placement, one label each, and enclosing frame

def combine():
    import gdsfactory as gf
    import klayout.db as kdb
    gf.gpdk.PDK.activate()
    ly=kdb.Layout();ly.dbu=DBU
    top=ly.create_cell('D1_L6000_P10_WIDTH_SWEEP')
    records=[]; snapshots={}; placed=[]; text_regions=[]; bounds=kdb.Box(); checks=[]
    for row,width in enumerate(WG_WIDTHS):
        ident=case_id(width)
        meta=json.loads((WORK/(ident+'.json')).read_text())
        old,src=read_gds(WORK/(ident+'.gds'))
        cell=copy_flat(ly,old,src,ident); r=regions(ly,cell)
        snapshots[ident]=r
        x=LEFT_BEND_X-meta['left_pnc_x_um'];y=FIRST_ROW_Y+row*ROW_PITCH
        trans=kdb.Trans(round(x/DBU),round(y/DBU))
        top.insert(kdb.CellInstArray(cell.cell_index(),trans))
        meta['origin_um']=[trans.disp.x*DBU,trans.disp.y*DBU]
        physical={n:r[(n,0)].transformed(trans) for n in FAB_LAYERS}
        placed.append(physical)
        for region in physical.values(): bounds+=region.bbox()
        text=gf.components.text(text=meta['label'],size=LABEL_SIZE,layer=(2,0))
        reg=kdb.Region(text.kdb_cell.begin_shapes_rec(gf.get_layer((2,0))))
        lab=ly.create_cell('LABEL_'+ident);lab.shapes(ly.layer(2,0)).insert(reg)
        lt=kdb.Trans(round((LEFT_BEND_X+200)/DBU),round((y+250)/DBU))
        top.insert(kdb.CellInstArray(lab.cell_index(),lt))
        text_regions.append(reg.transformed(lt));meta['label_origin_um']=[lt.disp.x*DBU,lt.disp.y*DBU]
        bounds+=text_regions[-1].bbox();records.append(meta)
    for ident,r in snapshots.items():
        reference=next(iter(snapshots.values()))
        for layer in r:
            if layer not in ((15,0),(16,0)):
                area=(r[layer]^reference[layer]).area()*DBU**2
                check(checks, ident+f'_width_sweep_preserves_layer{layer[0]}',area==0,xor_um2=area)
    for i,a in enumerate(placed):
        for j,b in enumerate(placed[i+1:],i+1):
            overlap=sum((a[n]&b[n]).area()*DBU**2 for n in FAB_LAYERS)
            check(checks,f'rows{i}_{j}_no_mask_overlaps',overlap==0,overlap_um2=overlap)
    for i,text in enumerate(text_regions):
        check(checks,f'label{i}_clear_of_all_devices',all((text&r).is_empty() for p in placed for r in p.values()))
    xmin=math.floor((bounds.left*DBU-FRAME_MARGIN)/100)*100
    xmax=math.ceil((bounds.right*DBU+FRAME_MARGIN)/100)*100
    ymin=math.floor((bounds.bottom*DBU-FRAME_MARGIN)/100)*100
    ymax=math.ceil((bounds.top*DBU+FRAME_MARGIN)/100)*100
    if ymax-ymin < MIN_FRAME_HEIGHT:
        middle=(ymin+ymax)/2;ymin=middle-MIN_FRAME_HEIGHT/2;ymax=middle+MIN_FRAME_HEIGHT/2
    m=configured_design(WG_WIDTHS[0],dict(CHIP_WIDTH=xmax-xmin,CHIP_MIN_HEIGHT=ymax-ymin))
    frame=gf.Component('LONG_MEANDER_FRAME')
    m.add_chip_outline_and_markers(frame,len(records));frame.flatten()
    frame.write_gds(WORK/'frame.gds',with_metadata=False)
    fl,fc=read_gds(WORK/'frame.gds');out=copy_flat(ly,fl,fc,'FRAME_AND_MARKERS')
    ft=kdb.Trans(round((xmin+xmax)/2/DBU),round((ymin+ymax)/2/DBU))
    top.insert(kdb.CellInstArray(out.cell_index(),ft))
    fr=regions(ly,out)
    for n in (2,10):
        markers=fr[(n,0)].transformed(ft)
        check(checks,f'global_L{n}_markers_clear',all((markers&p[n]).is_empty() for p in placed))
    OUTPUT.parent.mkdir(exist_ok=True)
    options=kdb.SaveLayoutOptions();options.gds2_max_vertex_count=2000
    ly.write(str(OUTPUT),options)
    final,final_top=read_gds(OUTPUT)
    insts={i.cell.name:i for i in final_top.each_inst() if i.cell.name in snapshots}
    labels=[i for i in final_top.each_inst() if i.cell.name.startswith('LABEL_')]
    check(checks,'three_complete_devices_and_three_labels',len(insts)==len(labels)==len(WG_WIDTHS)==3)
    for e in records:
        inst=insts[e['case_id']]; rr=regions(final,inst.cell)
        for layer,r in snapshots[e['case_id']].items():
            area=(r^rr[layer]).area()*DBU**2
            check(checks,e['case_id']+f'_final_export_layer{layer[0]}',area < .005,xor_um2=area)
        left=(rr[(11,0)].bbox().left+inst.trans.disp.x)*DBU
        check(checks,e['case_id']+'_left_bend_aligned',abs(left-LEFT_BEND_X)<1e-8,x_um=left)
        check(checks,e['case_id']+'_actual_y_origin',abs(inst.trans.disp.y*DBU-e['origin_um'][1])<1e-8)
    result=dict(status='PASS',source_file=SOURCE.name,
                source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                devices=records,frame_bounds_um=[xmin,ymin,xmax,ymax],row_pitch_um=ROW_PITCH,
                straight_pnc_pitch_um=STRAIGHT_PNC_PERIOD,bend_pnc_pitch_um=BEND_PNC_PERIOD,
                pitch_definition='center-to-center along meander centerline; symmetric bend end remainders',
                checks=checks,check_count=len(checks)+sum(len(e['checks']) for e in records),
                gds_sha256=hashlib.sha256(OUTPUT.read_bytes()).hexdigest())
    OUTPUT.with_suffix('.json').write_text(json.dumps(result,indent=2))
    with OUTPUT.with_suffix('.csv').open('w',newline='') as fp:
        fields=['case_id','width_um','straight_length_um','total_meander_um','post_pnc_taper_total_um','origin_um','label']
        w=csv.DictWriter(fp,fieldnames=fields);w.writeheader()
        w.writerows({k:e[k] for k in fields} for e in records)
    render(final,final_top,result)
    print('PASS',result['check_count'],'checks;',OUTPUT,flush=True)


# %% [06] Actual exported previews (not illustrative drawings)

def render(ly,top,result):
    import klayout.db as kdb
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    def view(cell,name,bounds,title,layers=(14,11,15)):
        fig,ax=plt.subplots(figsize=(7,8) if name=='bend_detail' else (15,5))
        for n in layers:
            reg=(kdb.Region(cell.begin_shapes_rec(ly.layer(n,0))) & box(bounds)).merged()
            lines=[]
            for p in reg.each():
                pts=[(v.x*DBU,v.y*DBU) for v in p.each_point_hull()]
                if pts:lines.append(pts+[pts[0]])
            flat = [point for line in lines for point in line + [(float('nan'),float('nan'))]]
            if flat:
                ax.plot([v[0] for v in flat],[v[1] for v in flat],linewidth=.6,label=f'Layer {n}')
        ax.set(xlim=(bounds[0],bounds[2]),ylim=(bounds[1],bounds[3]),
               xlabel='X (um)',ylabel='Y (um)',title=title)
        if name == 'overview':
            for entry in result['devices']:
                x,y=entry['label_origin_um']
                ax.text(x,y+20,entry['label'],fontsize=8)
        ax.set_aspect('equal');ax.legend(loc='upper right');fig.tight_layout()
        fig.savefig(WORK/(name+'.png'),dpi=200);plt.close(fig)

    f=result['frame_bounds_um']
    view(top,'overview',(f[0],FIRST_ROW_Y-450,f[2],FIRST_ROW_Y+(len(WG_WIDTHS)-1)*ROW_PITCH+450),
         'Design 1: three 6000 um straights; widths 1.62 / 1.64 / 1.68 um (bottom to top)')
    first=next(i.cell for i in top.each_inst() if i.cell.name==case_id(WG_WIDTHS[0]))
    view(first,'straight_detail',(-50,-220,50,-180),'Exported GDS: 10 um straight PNC pitch, 1.62 um ridge')
    view(first,'bend_detail',(STRAIGHT_LENGTH/2-20,-25,STRAIGHT_LENGTH/2+170,225),
         'Exported GDS: 10 um centerline arc-length PNC pitch in bend')
    view(first,'right_end',(STRAIGHT_LENGTH/2+20,-235,STRAIGHT_LENGTH/2+300,-165),
         'Design 1: complete right taper, alignment tips and tether',layers=(2,10,11,14,15,16))


# %% [07] Run each width in a fresh interpreter, then merge all layers

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker-width',type=float,help=argparse.SUPPRESS)
    args=parser.parse_args()
    if args.worker_width is not None:
        worker(args.worker_width);return
    require(len(WG_WIDTHS)==3 and len(set(WG_WIDTHS))==3,'Specify three distinct widths')
    require(all(math.isfinite(v) and v>0 for v in (*WG_WIDTHS,STRAIGHT_LENGTH,STRAIGHT_PNC_PERIOD,BEND_PNC_PERIOD,ROW_PITCH)),
            'Dimensions must be finite and positive')
    WORK.mkdir(exist_ok=True)
    for width in WG_WIDTHS:
        run=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker-width',str(width)],
                           cwd=ROOT,text=True,capture_output=True,timeout=300)
        (WORK/(case_id(width)+'.log')).write_text(run.stdout+'\n'+run.stderr)
        require(run.returncode==0,f'Width {width} failed:\n{run.stderr[-5000:]}')
        print(run.stdout,flush=True)
    combine()


if __name__=='__main__':
    main()
