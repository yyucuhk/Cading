"""Audit exported two-coupler L11 removal and terminal-only L14 windows.

Reads GDSII records and uses Shapely; does not import the builder or KLayout.
The baseline is the preceding six-device comparison, before this change.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from shapely.geometry import box, LineString
from validate_2611_straight_output_raw import read, union

STEM = '2611_design1_L6000_P10_straight_output_comparison'

def validate(root: Path, baseline: Path) -> dict:
    current = root/'gds'/(STEM+'.gds')
    old, ot = read(baseline)
    new, nt = read(current)
    checks, measurements = [], []
    def check(name, ok, **detail):
        checks.append(dict(name=name, passed=bool(ok), **detail))
        if not ok:
            raise AssertionError((name, detail))
    def exact(name, a, b):
        area = a.symmetric_difference(b).area
        check(name, area < 1e-6, xor_um2=area)
    check('same_top_cell', nt == ot)
    check('same_complete_cell_inventory', set(new) == set(old))
    refs = new[nt]['refs']
    devices = [(n,xy) for n,xy in refs if n.startswith('D1_L')]
    labels = [(n,xy) for n,xy in refs if n.startswith('LABEL_')]
    check('six_devices_and_six_labels', len(devices)==6 and len(labels)==6 and len(refs)==13)
    check('same_600um_y_pitch', sorted(xy[1] for _,xy in devices)==[-700000,-100000,500000,1100000,1700000,2300000])
    edited = {f'D1_L6000_W{w}_P10' for w in (1620,1640,1680)}
    for name, cell in new.items():
        check(name+'_same_references_and_transforms', sorted(cell['refs'])==sorted(old[name]['refs']))
        check(name+'_no_L16_L51', not ({16,51} & set(cell['polys'])))
        expected = set(old[name]['polys']) - ({11} if name in edited else set())
        check(name+'_exact_layer_inventory', set(cell['polys'])==expected)
        for layer in sorted(expected):
            if name in edited and layer==14:
                continue
            p0=sorted(tuple(p) for p in old[name]['polys'].get(layer,[]))
            p1=sorted(tuple(p) for p in cell['polys'].get(layer,[]))
            if p0==p1:
                check(name+f'_L{layer}_raw_polygons_identical', True)
            else:
                exact(name+f'_L{layer}_unchanged',union(old[name],layer),union(cell,layer))
    frame = new['FRAME_AND_MARKERS']
    check('markers_still_on_L15_not_L2', bool(frame['polys'].get(15)) and not frame['polys'].get(2))
    check('all_labels_preserved_on_L2', all(set(new[n]['polys'])=={2} for n,_ in labels))
    placed = dict(devices)
    for name in sorted(edited):
        cell = new[name]
        before = union(old[name],14)
        after = union(cell,14)
        l10 = union(cell,10)
        left = l10.intersection(box(-20000,25,-3000,20000))
        right = l10.intersection(box(3000,-20000,20000,-25))
        check(name+'_both_L10_tapers_present',not left.is_empty and not right.is_empty)
        lx, rx = left.bounds[2], right.bounds[0]
        check(name+'_actual_L10_entrances',abs(lx+3045)<1e-9 and abs(rx-3045)<1e-9,
              input_local_x_um=lx,output_local_x_um=rx)
        keep_left = box(-20000,25,lx,20000)
        keep_right = box(rx,-20000,20000,-25)
        keep = keep_left.union(keep_right)
        exact(name+'_L14_equals_old_end_windows_exactly',after,before.intersection(keep))
        check(name+'_L14_no_new_area',after.difference(before).area<1e-6)
        check(name+'_L14_no_geometry_outside_end_windows',after.difference(keep).area<1e-6)
        check(name+'_L14_no_meander_or_PNC_transition',after.intersection(box(-3045,-1000,3045,1000)).area<1e-6)
        bend_rois=box(-3200,-225,-3000,25).union(box(3000,-25,3200,225))
        check(name+'_L14_no_left_or_right_bend_fragments',after.intersection(bend_rois).area<1e-6)
        check(name+'_L11_absent',not cell['polys'].get(11))
        check(name+'_two_complete_terminal_regions',after.geom_type=='MultiPolygon' and len(after.geoms)==2)
        for title,window,start,y,sign in [('input',keep_left,lx,200,-1),('output',keep_right,rx,-200,1)]:
            kept=after.intersection(window)
            edge=kept.bounds[2] if sign<0 else kept.bounds[0]
            check(name+'_'+title+'_L14_start_equals_L10_start',abs(edge-start)<1e-9)
            width=l10.intersection(LineString([(start+sign*.001,y-10),(start+sign*.001,y+10)])).length
            check(name+'_'+title+'_unchanged_5um_L10_entrance',abs(width-5)<.002,measured_width_um=width)
            exact(name+'_'+title+'_all_downstream_L14_preserved',kept,before.intersection(window))
        # L15 is unchanged; inspect its three straight cross-sections as well.
        ridge=union(cell,15)
        width=int(name.split('_W')[1].split('_')[0])/1000
        for y in (-200,0,200):
            measured=ridge.intersection(LineString([(0,y-3),(0,y+3)])).length
            check(name+f'_straight_ridge_width_y{y}',abs(measured-width)<1e-8,measured_um=measured)
        geoms=list(ridge.geoms) if hasattr(ridge,'geoms') else [ridge]
        check(name+'_one_optical_body_four_existing_ticks',len(geoms)==5)
        dx,dy=placed[name]
        measurements.append(dict(case_id=name,removed_layer11=True,
            input_L14_start_global_x_um=lx+dx*.001,
            output_L14_start_global_x_um=rx+dx*.001,
            origin_um=[dx*.001,dy*.001],
            old_L14_area_um2=before.area,new_L14_area_um2=after.area,
            removed_L14_area_um2=before.area-after.area))
    result=dict(status='PASS',check_count=len(checks),checks=checks,measurements=measurements,
        scope='Only three two-coupler devices in the six-device L6000 comparison',
        straight_output_variants_unchanged=True,all_other_layers_and_labels_unchanged=True,
        original_source_GDS_unchanged=True,
        baseline_sha256=hashlib.sha256(baseline.read_bytes()).hexdigest(),
        gds_sha256=hashlib.sha256(current.read_bytes()).hexdigest(),
        method='Independent raw GDSII records and Shapely; no generator or KLayout',
        optical_or_mechanical_simulation_performed=False)
    (root/'gds/2611_two_coupler_windows_validation.json').write_text(json.dumps(result,indent=2))
    print('PASS',len(checks),'independent two-coupler window checks')
    return result

def render(root: Path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    c,t=read(root/'gds'/(STEM+'.gds'))
    n='D1_L6000_W1620_P10'
    out=root/'validation_2611_two_coupler_windows';out.mkdir(exist_ok=True)
    views=[('two_coupler_overview',(-5200,-400,5200,400),'Two-coupler device: Layer 11 removed; Layer 14 only at the two ends'),
           ('input_window',(-3200,170,-2990,230),'Input: Layer 14 starts at the Layer 10 entrance'),
           ('output_window',(2990,-230,3200,-170),'Output: Layer 14 starts at the Layer 10 entrance')]
    for filename,bounds,title in views:
        fig,ax=plt.subplots(figsize=(14,5))
        for layer in (2,10,14,15):
            r=union(c[n],layer).intersection(box(*bounds))
            gs=list(r.geoms) if hasattr(r,'geoms') else [r]
            xs,ys=[],[]
            for p in gs:
                if p.is_empty or p.geom_type!='Polygon':continue
                for ring in [p.exterior,*p.interiors]:
                    x,y=ring.xy;xs+=list(x)+[float('nan')];ys+=list(y)+[float('nan')]
            if xs:ax.plot(xs,ys,label=f'Layer {layer}',linewidth=1)
        if filename!='two_coupler_overview':
            ax.axvline(-3045 if filename=='input_window' else 3045,linestyle='--',linewidth=.7,label='L10 / L14 entrance')
        ax.set(xlim=(bounds[0],bounds[2]),ylim=(bounds[1],bounds[3]),xlabel='Local X (um)',ylabel='Local Y (um)',title=title)
        ax.legend(loc='upper right');fig.tight_layout();fig.savefig(out/(filename+'.png'),dpi=180);plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--render',action='store_true')
    a=p.parse_args();validate(a.root,a.baseline)
    if a.render:render(a.root)
