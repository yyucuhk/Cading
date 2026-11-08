"""Design 1: original L6000 devices plus a one-sided 2000 um straight-output set.

Run from the Cading root. Existing source/GDS files are never overwritten.
The lower-right output is replaced; the upper-left input remains unchanged.
A straight output continues the existing constant-width Layer-15 ridge and
15 um support envelope, with 10 um PNC pitch and matching release/protection
masks. It has NO 15-to-5 transition, optical taper, dense end cluster,
butterfly tether, output support island, alignment ticks, or fiber guide.

Only owned, flattened polygons beyond the output port are replaced. The
replacement window is below the middle arm, so neither Euler bend is cut.
Every polygon outside this terminal window is compared against the reference.
All dimensions are micrometres. Optical/fabrication performance is not simulated.
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
INCLUDE_REFERENCE_DEVICES = True     # originals first, new variants above
ROW_PITCH = 600.0
FIRST_ROW_Y = -700.0
LEFT_BEND_X = -900.0
LABEL_SIZE = 5.0
FRAME_MARGIN = 500.0
MIN_FRAME_HEIGHT = 4000.0
DBU = 0.001
FAB_LAYERS = (2, 10, 11, 13, 14, 15, 16)


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
    report = dict(width_um=width, reference_id=ref_id, variant_id=variant_id,
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
    import gdsfactory as gf
    import klayout.db as kdb
    gf.gpdk.PDK.activate()
    base=load_base()
    ly=kdb.Layout();ly.dbu=DBU
    top=ly.create_cell('D1_L6000_STRAIGHT_OUT_COMPARISON')
    entries=[]; placed=[]; texts=[]; snapshots={}; checks=[]; bounds=kdb.Box()
    modes=([False] if INCLUDE_REFERENCE_DEVICES else [])+[True]
    for straight_output in modes:
        for width in WG_WIDTHS:
            meta=json.loads((WORK/(base.case_id(width)+'.json')).read_text())
            ident=meta['variant_id'] if straight_output else meta['reference_id']
            src,sc=read_gds(WORK/(ident+'.gds'))
            rr=regions(src,sc)
            cell=insert_regions(ly,ident,rr);snapshots[ident]=rr
            row=len(entries); y=FIRST_ROW_Y+row*ROW_PITCH
            x=LEFT_BEND_X-meta['left_pnc_x_um']
            tr=kdb.Trans(round(x/DBU),round(y/DBU))
            top.insert(kdb.CellInstArray(cell.cell_index(),tr))
            physical={n:rr[(n,0)].transformed(tr) for n in FAB_LAYERS}
            placed.append(physical)
            for r in physical.values():bounds+=r.bbox()
            label=f'D1_W{width:.2f}_T80_P10'+(f'_OUT{OUTPUT_EXTENSION_LENGTH:g}' if straight_output else '')
            text=gf.components.text(text=label,size=LABEL_SIZE,layer=(2,0))
            textreg=kdb.Region(text.kdb_cell.begin_shapes_rec(gf.get_layer((2,0)))).merged()
            labelcell=insert_regions(ly,'LABEL_'+ident,{(2,0):textreg})
            lt=kdb.Trans(round((LEFT_BEND_X+200)/DBU),round((y+250)/DBU))
            top.insert(kdb.CellInstArray(labelcell.cell_index(),lt))
            texts.append(textreg.transformed(lt));bounds+=texts[-1].bbox()
            entry=dict(case_id=ident,mode='straight_output' if straight_output else 'original_two_tapers',
                       width_um=width,origin_um=[tr.disp.x*DBU,tr.disp.y*DBU], label=label,
                       label_origin_um=[lt.disp.x*DBU,lt.disp.y*DBU],
                       added_length_um=OUTPUT_EXTENSION_LENGTH if straight_output else 0.,
                       path_length_excluding_tapers_um=meta['meander_length_um']+(OUTPUT_EXTENSION_LENGTH if straight_output else 0.),
                       input_port_local_um=meta['input_port_um'],output_port_local_um=meta['output_port_um'])
            entries.append(entry)
    for i,a in enumerate(placed):
        for j,b in enumerate(placed[i+1:],i+1):
            area=sum((a[n]&b[n]).area()*DBU**2 for n in FAB_LAYERS)
            check(checks,f'rows{i}_{j}_no_fabrication_overlap',area==0,area_um2=area)
    for i,t in enumerate(texts):
        check(checks,f'label{i}_clear_of_devices',all((t&r).is_empty() for p in placed for r in p.values()))
        check(checks,f'label{i}_clear_of_other_labels',all((t&s).is_empty() for s in texts[i+1:]))
    xmin=math.floor((bounds.left*DBU-FRAME_MARGIN)/100)*100
    xmax=math.ceil((bounds.right*DBU+FRAME_MARGIN)/100)*100
    ymin=math.floor((bounds.bottom*DBU-FRAME_MARGIN)/100)*100
    ymax=math.ceil((bounds.top*DBU+FRAME_MARGIN)/100)*100
    if ymax-ymin < MIN_FRAME_HEIGHT:
        middle=(ymin+ymax)/2;ymin=middle-MIN_FRAME_HEIGHT/2;ymax=middle+MIN_FRAME_HEIGHT/2
    m=base.configured_design(WG_WIDTHS[0],dict(CHIP_WIDTH=xmax-xmin,CHIP_MIN_HEIGHT=ymax-ymin))
    frame=gf.Component();m.add_chip_outline_and_markers(frame,len(entries));frame.flatten()
    fr=regions(frame.kcl.layout,frame.kdb_cell)
    fc=insert_regions(ly,'FRAME_AND_MARKERS',fr)
    ft=kdb.Trans(round((xmin+xmax)/2/DBU),round((ymin+ymax)/2/DBU))
    top.insert(kdb.CellInstArray(fc.cell_index(),ft))
    for n in (2,10):
        markers=fr[(n,0)].transformed(ft)
        check(checks,f'global_L{n}_markers_clear',all((markers&p[n]).is_empty() for p in placed))
    write_layout(ly,OUTPUT)
    final,top=read_gds(OUTPUT)
    insts={i.cell.name:i for i in top.each_inst()}
    check(checks,'correct_device_and_label_count',
          len(insts)==2*len(entries)+1 and sum(n.startswith('LABEL_') for n in insts)==len(entries))
    for i,e in enumerate(entries):
        inst=insts[e['case_id']]; rr=regions(final,inst.cell)
        for layer,r in snapshots[e['case_id']].items():
            check(checks,e['case_id']+f'_saved_L{layer[0]}_exact',(r^rr[layer]).is_empty())
        check(checks,e['case_id']+'_actual_x_alignment',
              abs((rr[(11,0)].bbox().left+inst.trans.disp.x)*DBU-LEFT_BEND_X)<1e-8)
        check(checks,e['case_id']+'_actual_y_position',
              abs(inst.trans.disp.y*DBU-(FIRST_ROW_Y+i*ROW_PITCH))<1e-8)
    workers=[json.loads((WORK/(base.case_id(w)+'.json')).read_text()) for w in WG_WIDTHS]
    result=dict(status='PASS',devices=entries,device_count=len(entries),
                reference_count=len(WG_WIDTHS) if INCLUDE_REFERENCE_DEVICES else 0,
                variant_count=len(WG_WIDTHS), width_sweep_um=list(WG_WIDTHS),
                meander_straight_length_um=STRAIGHT_LENGTH,extension_length_um=OUTPUT_EXTENSION_LENGTH,
                pnc_pitch_um=PNC_PITCH,row_pitch_um=ROW_PITCH,left_bend_x_um=LEFT_BEND_X,
                side='lower-right output; upper-left input unchanged',
                extension_layers=[2,11,14,15,16],frame_bounds_um=[xmin,ymin,xmax,ymax],
                original_gds_sha256=hashlib.sha256(BASE_GDS.read_bytes()).hexdigest(),
                gds_sha256=hashlib.sha256(OUTPUT.read_bytes()).hexdigest(),
                checks=checks,workers=workers,
                check_count=len(checks)+sum(len(w['checks']) for w in workers),
                optical_or_mechanical_simulation_performed=False)
    OUTPUT.with_suffix('.json').write_text(json.dumps(result,indent=2))
    with OUTPUT.with_suffix('.csv').open('w',newline='') as fp:
        fields=['case_id','mode','width_um','origin_um','added_length_um','path_length_excluding_tapers_um','label']
        writer=csv.DictWriter(fp,fieldnames=fields);writer.writeheader()
        writer.writerows({k:e[k] for k in fields} for e in entries)
    render(final,top,result)
    print('PASS',result['check_count'],'checks;',OUTPUT,flush=True)


# %% [05] Views rendered directly from the exported GDS

def render(ly,top,result):
    import klayout.db as kdb
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    def view(cell,name,bounds,title,layers=(14,11,15)):
        fig,ax=plt.subplots(figsize=(16,7) if name=='overview' else (14,5))
        for n in layers:
            rr=(kdb.Region(cell.begin_shapes_rec(ly.layer(n,0)))&box(bounds)).merged()
            coords=[]
            for p in rr.each():
                points=[(q.x*DBU,q.y*DBU) for q in p.each_point_hull()]
                if points:coords+=points+[points[0],(float('nan'),float('nan'))]
            if coords:ax.plot([p[0] for p in coords],[p[1] for p in coords],linewidth=.7,label=f'Layer {n}')
        ax.set(xlim=(bounds[0],bounds[2]),ylim=(bounds[1],bounds[3]),xlabel='X (um)',ylabel='Y (um)',title=title)
        if name=='overview':
            for e in result['devices']:
                x,y=e['label_origin_um'];ax.text(x,y+20,e['label'],fontsize=8)
        if name!='output_connection':ax.set_aspect('equal')
        ax.legend(loc='upper right');fig.tight_layout();fig.savefig(WORK/(name+'.png'),dpi=200);plt.close(fig)
    f=result['frame_bounds_um']
    view(top,'overview',(f[0],FIRST_ROW_Y-400,f[2],FIRST_ROW_Y+(len(result['devices'])-1)*ROW_PITCH+400),
         'Design 1: original devices below; three straight-output variants above')
    e=next(e for e in result['devices'] if e['mode']=='straight_output')
    c=ly.cell(e['case_id']);px,py=e['output_port_local_um']
    view(c,'output_connection',(px-30,py-23,px+180,py+23),
         'No output taper: constant-width ridge and 10 um PNC supports (axes expanded)',layers=(2,10,11,14,15,16))
    view(c,'full_output',(px-20,py-40,px+OUTPUT_EXTENSION_LENGTH+30,py+40),
         'New straight output: 2000 um from the original meander port',layers=(14,11,15))
    ix,iy=e['input_port_local_um']
    view(c,'preserved_input',(ix-310,iy-35,ix+20,iy+35),
         'Input preserved: 15-to-5 um PNC transition, dense end group, original optical taper',layers=(2,10,11,14,15,16))


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


if __name__=='__main__':
    main()
