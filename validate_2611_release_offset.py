"""Measure both 1 um release margins and the 10 um longitudinal mask setback.

Reads actual GDSII polygons through the independent raw parser, not the CAD
builder. Optional --baseline additionally verifies exact preservation of all
non-mask layers, labels, markers, transforms, and the remaining Layer14 shape.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from shapely.geometry import box, Polygon, LineString
from shapely.affinity import scale
from validate_2611_straight_output_raw import read, union

STEM = '2611_design1_L6000_P10_straight_output_comparison'
OFFSET = 10.0
MARGIN = 1.0


def validate(root: Path, baseline: Path | None = None) -> dict:
    path = root/'gds'/(STEM+'.gds')
    new, top = read(path)
    old, oldtop = read(baseline) if baseline else (None, None)
    checks, measured, cache = [], [], {}

    def check(name, ok, **details):
        checks.append(dict(name=name, passed=bool(ok), **details))
        if not ok:
            raise AssertionError((name, details))

    def exact(name, a, b):
        area = a.symmetric_difference(b).area
        check(name, area < 1e-6, xor_um2=area)

    refs = new[top]['refs']
    devices = {n:p for n,p in refs if n.startswith('D1_L')}
    labels = {n:p for n,p in refs if n.startswith('LABEL_')}
    check('six_devices_six_labels_one_frame', len(devices)==6 and len(labels)==6 and len(refs)==13)
    check('preserved_600um_vertical_pitch', sorted(p[1] for p in devices.values())==[-700000,-100000,500000,1100000,1700000,2300000])
    check('preserved_horizontal_positions', all(p[0]==2242153 for p in devices.values()))
    if old:
        check('same_top_and_cell_inventory', top==oldtop and set(new)==set(old))
    for name, cell in new.items():
        check(name+'_excluded_layers_absent', not ({11,16,51}&set(cell['polys'])))
        if old:
            check(name+'_transforms_unchanged', sorted(cell['refs'])==sorted(old[name]['refs']))
            for layer in set(cell['polys'])|set(old[name]['polys']):
                if name in devices and layer in (2,14):
                    continue
                a = sorted(tuple(x) for x in cell['polys'].get(layer,[]))
                b = sorted(tuple(x) for x in old[name]['polys'].get(layer,[]))
                check(name+f'_L{layer}_raw_polygon_coordinates_unchanged', a==b)
    check('labels_only_on_layer2', all(set(new[n]['polys'])=={2} for n in labels))
    check('global_marks_only_moved_to_layer15', bool(new['FRAME_AND_MARKERS']['polys'].get(15)) and not new['FRAME_AND_MARKERS']['polys'].get(2))

    for name, translation in devices.items():
        cell = new[name]
        variant = name.endswith('_OUT2000')
        width = int(name.split('_W')[1].split('_')[0])/1000
        rs = {n:union(cell,n) for n in (2,10,14,15)}
        cache[name] = rs
        check(name+'_correct_remaining_device_layers', set(cell['polys'])=={2,10,13,14,15})
        ends = [('input',-1,200)] + ([] if variant else [('output',1,-200)])
        keep = Polygon()
        end_boxes = []
        for side, sign, y in ends:
            key = name+'_'+side
            search = box(-20000,25,-3000,20000) if sign<0 else box(3000,-20000,20000,-25)
            slab = rs[10].intersection(search)
            check(key+'_optical_taper_present', not slab.is_empty)
            sx = slab.bounds[2] if sign<0 else slab.bounds[0]
            check(key+'_L10_entrance_unchanged', abs(sx-sign*3045)<1e-9, x_um=sx)
            start = sx + sign*OFFSET
            end = box(-20000,25,start,20000) if sign<0 else box(start,-20000,20000,-25)
            end_boxes.append(end)
            keep = keep.union(end)
            margins = {2:[],14:[]}
            for layer in (2,14):
                region = rs[layer].intersection(end)
                check(key+f'_L{layer}_terminal_window_present', not region.is_empty)
                edge = region.bounds[2] if sign<0 else region.bounds[0]
                check(key+f'_L{layer}_exact_10um_start', abs(edge-start)<1e-9,
                      l10_entrance_um=sx, mask_start_um=edge, setback_um=abs(edge-sx))
                lead = box(min(sx,start)+.001,y-20,max(sx,start)-.001,y+20)
                check(key+f'_L{layer}_no_mask_in_first10um', rs[layer].intersection(lead).area<1e-9)
                for distance in (10.1,15,20,30,35,39):
                    x = sx+sign*distance
                    cut = rs[layer].intersection(LineString([(x,y-20),(x,y+20)]))
                    u = distance/40
                    f = 6*u**5-15*u**4+10*u**3
                    expected = 7.596*(1-f)+1.5*f
                    check(key+f'_L{layer}_smooth_first_stage_{distance}',
                          cut.geom_type=='LineString' and abs(cut.length/2-expected)<.0021,
                          measured_halfwidth_um=cut.length/2, expected_halfwidth_um=expected)
                for distance in (50,60,70,75):
                    x = sx+sign*distance
                    line = LineString([(x,y-20),(x,y+20)])
                    mask = rs[layer].intersection(line)
                    optical = rs[10].intersection(line)
                    check(key+f'_L{layer}_one_window_cross_section_{distance}', mask.geom_type=='LineString')
                    lower = optical.bounds[1]-mask.bounds[1]
                    upper = mask.bounds[3]-optical.bounds[3]
                    margins[layer].extend((lower,upper))
                    check(key+f'_L{layer}_1um_on_each_side_{distance}',
                          .9989<=lower<=1.0021 and .9989<=upper<=1.0021,
                          lower_um=lower, upper_um=upper)
                    check(key+f'_L{layer}_covers_taper_after_setback_{distance}', optical.difference(mask).length<1e-9)
                if old:
                    far = box(-20000,25,sx-95.01,20000) if sign<0 else box(sx+95.01,-20000,20000,-25)
                    exact(key+f'_L{layer}_far_channels_pads_unchanged',
                          rs[layer].intersection(far),union(old[name],layer).intersection(far))
            for distance in (0.1,5,9.9):
                x = sx+sign*distance
                check(key+f'_waveguide_not_shortened_{distance}',
                      rs[10].intersection(LineString([(x,y-5),(x,y+5)])).length>0)
            measured.append(dict(case_id=name,side=side,l10_start_local_x_um=sx,
                mask_start_local_x_um=start,mask_start_global_x_um=start+translation[0]*.001,
                offset_um=OFFSET,layer2_margin_min_um=min(margins[2]),layer2_margin_max_um=max(margins[2]),
                layer14_margin_min_um=min(margins[14]),layer14_margin_max_um=max(margins[14])))
        for layer in (2,14):
            check(name+f'_L{layer}_no_mask_elsewhere', rs[layer].difference(keep).area<1e-6)
            if old:
                before = union(old[name],layer)
                check(name+f'_L{layer}_only_reduced',rs[layer].difference(before).area<1e-6)
                if layer==14:
                    exact(name+'_L14_only10um_setback_no_other_change',rs[14],before.intersection(keep))
                else:
                    check(name+'_L2_expansion_reduced_not_just_cropped',
                          before.intersection(keep).area-rs[2].area>100,
                          removed_after_crop_um2=before.intersection(keep).area-rs[2].area)
        for y in (-200,0,200):
            length = rs[15].intersection(LineString([(0,y-4),(0,y+4)])).length
            check(name+f'_ridge_width_{y}',abs(length-width)<1e-8,width_um=length)
        pieces = list(rs[15].geoms) if hasattr(rs[15],'geoms') else [rs[15]]
        check(name+'_connected_ridge_plus_original_ticks',len(pieces)==(3 if variant else 5))
        if variant:
            roi = box(3000,-400,20000,-25)
            exact(name+'_straight_output2000um_unchanged',rs[15].intersection(roi),box(3000,-200-width/2,5000,-200+width/2))
            for layer in (2,10,14):
                check(name+f'_bare_straight_output_no_L{layer}',rs[layer].intersection(roi).area<1e-6)
        else:
            for layer in (2,14):
                left = scale(rs[layer].intersection(end_boxes[0]),xfact=-1,yfact=-1,origin=(0,0))
                right = rs[layer].intersection(end_boxes[1])
                error = left.hausdorff_distance(right)
                check(name+f'_L{layer}_ends_symmetric_within_grid',error<=.00101,max_edge_difference_um=error)
    for w in (1620,1640,1680):
        a=f'D1_L6000_W{w}_P10';b=a+'_OUT2000';roi=box(-20000,25,-3045,20000)
        for layer in (2,10,14,15):
            exact(a+f'_same_input_across_sets_L{layer}',cache[a][layer].intersection(roi),cache[b][layer].intersection(roi))
    check('all_nine_retained_windows_measured',len(measured)==9)
    report = dict(status='PASS',check_count=len(checks),checks=checks,measurements=measured,
        device_count=6,retained_coupler_windows=9,layer2_margin_um=MARGIN,layer14_margin_um=MARGIN,
        release_start_offset_um=OFFSET,offset_reference='Outward along taper from actual Layer10 entrance',
        removed_layers=[11,16,51],baseline_comparison_performed=bool(old),
        baseline_sha256=hashlib.sha256(baseline.read_bytes()).hexdigest() if baseline else None,
        gds_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        method='Independent raw GDSII polygons and Shapely; no CAD builder or KLayout import',
        optical_or_mechanical_simulation_performed=False)
    (root/'gds/2611_release_offset_validation.json').write_text(json.dumps(report,indent=2))
    print('PASS',len(checks),'independent mask-margin and setback checks')
    return report


def render(root: Path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cells,_=read(root/'gds'/(STEM+'.gds'))
    dest=root/'validation_2611_release_offset';dest.mkdir(exist_ok=True)
    views=[('two_coupler_input','D1_L6000_W1620_P10',(-3140,190,-3035,210),-3045,200,-1),
           ('two_coupler_output','D1_L6000_W1620_P10',(3035,-210,3140,-190),3045,-200,1),
           ('single_coupler_input','D1_L6000_W1620_P10_OUT2000',(-3140,190,-3035,210),-3045,200,-1)]
    for name,cell,bounds,sx,y,sign in views:
        fig,ax=plt.subplots(figsize=(13,5))
        for layer,style in ((2,'--'),(10,'-'),(14,':'),(15,'-')):
            region=union(cells[cell],layer).intersection(box(*bounds))
            pieces=list(region.geoms) if hasattr(region,'geoms') else [region]
            xs=[];ys=[]
            for p in pieces:
                if p.is_empty or p.geom_type!='Polygon':continue
                for ring in (p.exterior,*p.interiors):
                    x,z=ring.xy;xs+=list(x)+[float('nan')];ys+=list(z)+[float('nan')]
            if xs:ax.plot(xs,ys,style,label=f'Layer {layer}',linewidth=1.2)
        ax.axvline(sx,linestyle='--',linewidth=.6)
        ax.axvline(sx+sign*10,linestyle=':',linewidth=.6)
        ax.annotate('10 um mask setback',xy=(sx+sign*5,y+8),xytext=(sx+sign*45,y+8),
                    arrowprops={'arrowstyle':'->'},fontsize=10)
        ax.set(xlim=(bounds[0],bounds[2]),ylim=(bounds[1],bounds[3]),xlabel='Local X (um)',ylabel='Local Y (um)',
               title='Both release masks: 1 um expansion; start 10 um after Layer 10 begins')
        ax.legend(loc='lower right');fig.tight_layout();fig.savefig(dest/(name+'.png'),dpi=180);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--baseline',type=Path)
    p.add_argument('--render',action='store_true')
    a=p.parse_args();validate(a.root,a.baseline)
    if a.render:render(a.root)
