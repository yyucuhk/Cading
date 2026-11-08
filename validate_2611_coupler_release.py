"""Independent raw-GDS/Shapely audit of the final six-device release masks.

No layout generator or KLayout is imported. --baseline adds full before/after
polygon preservation checks. The layer-14 expansion applies to the existing
smooth opening: its broad entrance joins the inherited transition band, and
its narrow taper expansion is 1 um per side. Channels and pads are unchanged.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from shapely.geometry import box, LineString, Polygon
from shapely.affinity import scale
from validate_2611_straight_output_raw import read, union

STEM='2611_design1_L6000_P10_straight_output_comparison'

def validate(root: Path, baseline: Path | None = None) -> dict:
    path=root/'gds'/(STEM+'.gds')
    new,top=read(path)
    old,ot=read(baseline) if baseline else (None,None)
    checks=[];measured=[]
    def check(name,ok,**kw):
        checks.append(dict(name=name,passed=bool(ok),**kw))
        if not ok: raise AssertionError((name,kw))
    def exact(name,a,b,tol=1e-6):
        ar=a.symmetric_difference(b).area
        check(name,ar<tol,xor_um2=ar)
    refs=new[top]['refs']
    devs={n:p for n,p in refs if n.startswith('D1_L')}
    labels={n:p for n,p in refs if n.startswith('LABEL_')}
    check('six_devices_six_labels_one_frame',len(devs)==6 and len(labels)==6 and len(refs)==13)
    check('same_y_pitch',sorted(p[1] for p in devs.values())==[-700000,-100000,500000,1100000,1700000,2300000])
    check('same_x_origins',all(p[0]==2242153 for p in devs.values()))
    if old:
        check('same_top_cell',top==ot)
        check('same_cell_inventory',set(new)==set(old))
    for name,c in new.items():
        check(name+'_no11_16_51',not ({11,16,51}&set(c['polys'])))
        if old:
            check(name+'_unchanged_transforms',sorted(c['refs'])==sorted(old[name]['refs']))
            for layer in set(c['polys'])|set(old[name]['polys']):
                if name in devs and layer in (2,11,14): continue
                if sorted(tuple(x) for x in c['polys'].get(layer,[]))==sorted(tuple(x) for x in old[name]['polys'].get(layer,[])):
                    check(name+f'_L{layer}_raw_coordinates_identical',True)
                else:
                    exact(name+f'_L{layer}_unchanged',union(c,layer),union(old[name],layer))
    check('labels_remain_on_L2',all(set(new[n]['polys'])=={2} for n in labels))
    check('global_marks_remain_L15',bool(new['FRAME_AND_MARKERS']['polys'].get(15)) and not new['FRAME_AND_MARKERS']['polys'].get(2))
    cache={}
    for name in devs:
        c=new[name];variant=name.endswith('_OUT2000')
        width=int(name.split('_W')[1].split('_')[0])/1000
        rs={n:union(c,n) for n in (2,10,14,15)};cache[name]=rs
        check(name+'_expected_device_layers',set(c['polys'])=={2,10,13,14,15})
        sides=[('input',-1,200)] + ([] if variant else [('output',1,-200)])
        keep=Polygon();end_regions=[]
        for side,sign,y in sides:
            search=box(-20000,25,-3000,20000) if sign<0 else box(3000,-20000,20000,-25)
            slab=rs[10].intersection(search)
            check(name+'_'+side+'_L10_present',not slab.is_empty)
            sx=slab.bounds[2] if sign<0 else slab.bounds[0]
            check(name+'_'+side+'_L10_start_unchanged',abs(sx-sign*3045)<1e-9,start_local_um=sx)
            end=box(-20000,25,sx,20000) if sign<0 else box(sx,-20000,20000,-25)
            keep=keep.union(end);end_regions.append(end)
            for layer in (2,14):
                r=rs[layer].intersection(end)
                check(name+'_'+side+f'_L{layer}_nonempty',not r.is_empty)
                edge=r.bounds[2] if sign<0 else r.bounds[0]
                check(name+'_'+side+f'_L{layer}_starts_at_L10',abs(edge-sx)<1e-9)
            # Source transition side band has half-width 2.5+5.096 at entry.
            # Its downstream target is now 0.5+1, not 0.5+3. The original
            # smootherstep and 40 um first-stage length must be retained.
            first_errors=[]
            for distance in (1,5,10,20,30,35,39):
                x=sx+sign*distance;line=LineString([(x,y-20),(x,y+20)])
                half=rs[14].intersection(line).length/2
                u=distance/40;f=6*u**5-15*u**4+10*u**3
                expected=7.596*(1-f)+1.5*f
                err=abs(half-expected);first_errors.append(err)
                check(name+'_'+side+f'_smooth_L14_first_stage_{distance}',err<.0021,
                      measured_halfwidth_um=half,expected_halfwidth_um=expected)
            # Measure each side separately, rather than only total width.
            margins=[]
            for distance in (50,60,70,75):
                x=sx+sign*distance;line=LineString([(x,y-20),(x,y+20)])
                mask=rs[14].intersection(line);sl=rs[10].intersection(line)
                check(name+'_'+side+f'_single_L14_cross_section_{distance}',mask.geom_type=='LineString')
                lower=sl.bounds[1]-mask.bounds[1];upper=mask.bounds[3]-sl.bounds[3]
                margins.extend([lower,upper])
                check(name+'_'+side+f'_one_um_margin_both_sides_{distance}',
                      .9989<=lower<=1.0021 and .9989<=upper<=1.0021,
                      lower_um=lower,upper_um=upper)
                check(name+'_'+side+f'_taper_covered_{distance}',sl.difference(mask).length<1e-9)
            if old:
                far=box(-20000,25,sx-95.01,20000) if sign<0 else box(sx+95.01,-20000,20000,-25)
                exact(name+'_'+side+'_downstream_channels_pads_preserved',
                      rs[14].intersection(far),union(old[name],14).intersection(far))
            dx,dy=devs[name]
            measured.append(dict(case_id=name,side=side,start_global_x_um=sx+dx*.001,
                                 layer14_margin_min_um=min(margins),layer14_margin_max_um=max(margins),
                                 max_first_stage_halfwidth_error_um=max(first_errors)))
        for layer in (2,14):
            check(name+f'_L{layer}_nothing_outside_coupler_windows',rs[layer].difference(keep).area<1e-6)
            if old:
                before=union(old[name],layer)
                check(name+f'_L{layer}_only_removed_area',rs[layer].difference(before).area<1e-6)
                if layer==2:exact(name+'_L2_exact_crop_no_margin_change',rs[2],before.intersection(keep))
        r14=rs[14]
        pieces=list(r14.geoms) if hasattr(r14,'geoms') else [r14]
        check(name+'_one_L14_region_per_retained_coupler',len(pieces)==len(sides))
        ridge=rs[15];pieces=list(ridge.geoms) if hasattr(ridge,'geoms') else [ridge]
        check(name+'_one_optical_body_and_existing_ticks',len(pieces)==(3 if variant else 5))
        for y in (-200,0,200):
            w=ridge.intersection(LineString([(0,y-4),(0,y+4)])).length
            check(name+f'_ridge_width_y{y}',abs(w-width)<1e-8,width_um=w)
        if variant:
            roi=box(3000,-400,20000,-25)
            exact(name+'_unchanged_2000um_straight_ridge',ridge.intersection(roi),box(3000,-200-width/2,5000,-200+width/2))
            for n in (2,10,14):check(name+f'_straight_output_no_L{n}',rs[n].intersection(roi).area<1e-6)
        else:
            left=rs[14].intersection(end_regions[0]);right=rs[14].intersection(end_regions[1])
            exact(name+'_L14_two_ends_rotationally_symmetric',scale(left,xfact=-1,yfact=-1,origin=(0,0)),right)
    for w in (1620,1640,1680):
        n=f'D1_L6000_W{w}_P10';v=n+'_OUT2000';input_roi=box(-20000,25,-3045,20000)
        for layer in (2,10,14,15):exact(n+f'_both_sets_same_input_L{layer}',cache[n][layer].intersection(input_roi),cache[v][layer].intersection(input_roi))
    report=dict(status='PASS',check_count=len(checks),checks=checks,measurements=measured,
                scope='Both sets in six-device L6000 comparison',device_count=6,retained_coupler_windows=9,
                removed_layers=[11,16,51],cropped_release_layers=[2,14],layer14_margin_um=1.0,
                layer2_inside_retained_windows_unchanged=True,
                baseline_comparison_performed=bool(old),
                baseline_sha256=hashlib.sha256(baseline.read_bytes()).hexdigest() if baseline else None,
                gds_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                method='Raw GDSII parser and independent Shapely geometry; no generator or KLayout',
                optical_or_mechanical_simulation_performed=False)
    (root/'gds/2611_coupler_release_validation.json').write_text(json.dumps(report,indent=2))
    print('PASS',len(checks),'independent all-coupler release checks')
    return report

def render(root: Path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cs,_=read(root/'gds'/(STEM+'.gds'))
    out=root/'validation_2611_coupler_release';out.mkdir(exist_ok=True)
    views=[('two_coupler_input','D1_L6000_W1620_P10',(-3170,185,-3025,215),'Two-coupler input: retained L2 and independently narrowed L14'),
           ('two_coupler_output','D1_L6000_W1620_P10',(3025,-215,3170,-185),'Two-coupler output: 1 um Layer 14 expansion'),
           ('single_coupler_input','D1_L6000_W1620_P10_OUT2000',(-3170,185,-3025,215),'Straight-output variant: same input window policy'),
           ('straight_output','D1_L6000_W1620_P10_OUT2000',(2980,-208,5020,-192),'Straight output preserved: no PNC or release mask')]
    for fn,name,bounds,title in views:
        fig,ax=plt.subplots(figsize=(14,5))
        for n in (2,10,14,15):
            r=union(cs[name],n).intersection(box(*bounds));pieces=list(r.geoms) if hasattr(r,'geoms') else [r]
            xs=[];ys=[]
            for p in pieces:
                if p.is_empty or p.geom_type!='Polygon':continue
                for ring in [p.exterior,*p.interiors]:
                    x,y=ring.xy;xs+=list(x)+[float('nan')];ys+=list(y)+[float('nan')]
            if xs:ax.plot(xs,ys,label=f'Layer {n}',linewidth=1)
        ax.set(xlim=(bounds[0],bounds[2]),ylim=(bounds[1],bounds[3]),xlabel='Local X (um)',ylabel='Local Y (um)',title=title)
        ax.legend(loc='upper right');fig.tight_layout();fig.savefig(out/(fn+'.png'),dpi=180);plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--baseline',type=Path)
    p.add_argument('--render',action='store_true')
    a=p.parse_args();validate(a.root,a.baseline)
    if a.render:render(a.root)
