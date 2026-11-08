"""Independent raw-GDS/Shapely checks for original Design-1/2 L14 windows.

Requires --baseline with the three previous GDS files. Does not import the
new window generator or use GDSFactory/KLayout. Geometry checks only.
"""
from __future__ import annotations
import argparse,json,math,hashlib
from pathlib import Path
from shapely.geometry import Polygon,box,LineString
from shapely.affinity import affine_transform
from shapely.ops import unary_union
from validate_2611_straight_output_raw import read,union

FILES=('2611_tapered_yue_dual_resist_L10_L12.gds','2611_tapered_yue_design2.gds',
       '2611_design1_design2_L1400_18device_sweep.gds')

def validate(root,baseline,render=False):
    checks=[];ends=[];previews={}
    def ck(name,ok,**detail):
        checks.append(dict(name=name,passed=bool(ok),**detail))
        if not ok: raise AssertionError((name,detail))
    def same(name,a,b,tol=1e-7):
        area=a.symmetric_difference(b).area;ck(name,area<tol,xor_um2=area)
    def edges(shape,x):
        cut=shape.intersection(LineString([(x,-175),(x,175)]))
        if cut.is_empty or cut.geom_type!='LineString':raise AssertionError((x,cut.geom_type))
        return cut.bounds[1],cut.bounds[3]
    manifest=json.loads((root/'gds'/FILES[2]).with_suffix('.json').read_text())
    cases={e['case_id']:e for e in manifest['devices']}
    for fn in FILES:
        before,bt=read(baseline/fn);after,at=read(root/'gds'/fn)
        ck(fn+'_same_top_and_cells',bt==at and set(before)==set(after))
        for name in before:
            a,b=before[name],after[name];key=fn+'/'+name
            ck(key+'_references',a['refs']==b['refs'])
            for layer in set(a['polys'])|set(b['polys']):
                if layer==14:continue
                ck(key+f'_L{layer}_raw_polygons_identical',a['polys'].get(layer,[])==b['polys'].get(layer,[]))
            if not all(l in a['polys'] for l in (10,11,14)):
                same(key+'_nondevice_L14',union(a,14),union(b,14));continue
            if fn==FILES[0]:design=1;stages=[40.,40.]
            elif fn==FILES[1]:design=2;stages=[20.,100.,10.]
            else:design=cases[name]['design'];stages=cases[name]['stage_lengths_um']
            first,total=stages[0],sum(stages);old=union(a,14);new=union(b,14);slab=union(a,10)
            rois=[]
            for side,sign,cy in [('input',-1,200.),('output',1,-200.)]:
                sx=sign*745.;tr=[sign,0,0,sign,-sign*sx,-sign*cy]
                back=[sign,0,0,sign,sx,cy];roi=box(0,-175,20000,175)
                rois.append(affine_transform(roi,back))
                ow=affine_transform(old,tr).intersection(roi);nw=affine_transform(new,tr).intersection(roi)
                sw=affine_transform(slab,tr).intersection(box(0,-175,20000,175))
                k=key+'/'+side
                ck(k+'_connected_no_holes',nw.geom_type=='Polygon' and len(nw.interiors)==0)
                ck(k+'_PNC_end_position',abs(nw.bounds[0])<1e-8)
                lo0,hi0=edges(ow,0);lo,hi=edges(nw,0)
                ck(k+'_PNC_endpoint_width_preserved',abs(lo-lo0)<1e-8 and abs(hi-hi0)<1e-8)
                exp=sw.intersection(box(0,-30,total,30)).buffer(1,join_style=2)
                lo1,hi1=edges(exp,first)
                for frac in (.05,.125,.25,.5,.75,.95):
                    x=first*frac;f=frac**3*(10-15*frac+6*frac**2);l,h=edges(nw,x)
                    err=max(abs(l-(lo0+(lo1-lo0)*f)),abs(h-(hi0+(hi1-hi0)*f)))
                    ck(k+f'_first_stage_blend_{frac}',err<.0016,max_edge_error_um=err)
                narrow=box(first+.01,-30,total-.01,30)
                a1=nw.intersection(narrow);b1=exp.intersection(narrow)
                delta=a1.symmetric_difference(b1)
                ck(k+'_one_um_offset_window',a1.boundary.hausdorff_distance(b1.boundary)<.0016 and
                   delta.difference(b1.boundary.buffer(.0016)).area<1e-7,
                   max_edge_error_um=a1.boundary.hausdorff_distance(b1.boundary))
                centers=sorted({round(x,6) for x,y in ow.exterior.coords if abs(y-150)<1e-8})
                ncenters=sorted({round(x,6) for x,y in nw.exterior.coords if abs(y-150)<1e-8})
                ck(k+'_two_300um_circles',len(centers)==2 and centers==ncenters and
                   abs(nw.bounds[1]+150)<1e-8 and abs(nw.bounds[3]-150)<1e-8)
                arcs=box(0,100,20000,175).union(box(0,-175,20000,-100))
                same(k+'_circle_external_arcs_unchanged',ow.intersection(arcs),nw.intersection(arcs))
                far=ow.bounds[2];flo,fhi=edges(ow,far);slo,shi=edges(exp,total)
                ck(k+'_outer_endpoint_preserved',abs(nw.bounds[2]-far)<1e-8)
                for j in range(20):
                    x=total+(far-total)*j/19
                    if any(abs(x-cx)<151 for cx in centers):continue
                    l,h=edges(nw,x);f=(x-total)/(far-total)
                    err=max(abs(l-(slo+(flo-slo)*f)),abs(h-(shi+(fhi-shi)*f)))
                    ck(k+f'_linear_channel_{j}',err<.0016,max_edge_error_um=err)
                channel=Polygon([(total,slo),(far,flo),(far,fhi),(total,shi)])
                circles=[Polygon([(round((cx+150*math.cos(math.radians(i*2.5)))*1000)/1000,
                                   round(150*math.sin(math.radians(i*2.5))*1000)/1000)
                                  for i in range(144)]) for cx in centers]
                expected=unary_union([channel,*circles]);actual=nw.intersection(box(total,-175,20000,175))
                delta=actual.symmetric_difference(expected)
                ck(k+'_no_middle_frame_or_extra_opening',delta.difference(expected.boundary.buffer(.0016)).area<1e-7 and
                   actual.boundary.hausdorff_distance(expected.boundary)<.0016,
                   xor_um2=delta.area, max_edge_error_um=actual.boundary.hausdorff_distance(expected.boundary))
                ends.append(dict(file=fn,cell=name,design=design,side=side,
                    pnc_end_width_um=hi0-lo0,pnc_length_um=45.,
                    margin_um=1.,linear_start_from_l10_um=total,
                    linear_start_width_um=shi-slo,circle_centers_um=centers))
                if fn in FILES[:2] and side=='output':previews[design]=(ow,nw,sw,stages)
            changed=unary_union(rois)
            same(key+'_entire_PNC_meander_and_bends_unchanged',old.difference(changed),new.difference(changed))
            for sign,cy in [(-1,200.),(1,-200.)]:
                collar=box(sign*745.-.101,cy-30,sign*745.+.101,cy+30)
                same(key+f'_endpoint_PNC_collar_{sign}',old.intersection(collar),new.intersection(collar))
    ck('twenty_devices_forty_taper_ends',len(ends)==40,taper_ends=len(ends))
    if render:
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        out=root/'validation_original_design_windows';out.mkdir(exist_ok=True)
        for d,(old,new,slab,stages) in previews.items():
            for full in (False,True):
                fig,ax=plt.subplots(figsize=(12,4.5))
                ax.plot(*old.exterior.xy,'--',label='Previous Layer 14',linewidth=1.3)
                ax.plot(*new.exterior.xy,label='Updated Layer 14',linewidth=1.7)
                if not full:
                    for geom in ([slab] if slab.geom_type=='Polygon' else slab.geoms):
                        ax.plot(*geom.exterior.xy,label='Layer 10',linewidth=.8)
                ax.axvline(0,linestyle=':',linewidth=.9)
                ax.set(xlim=(-5,new.bounds[2]+20) if full else (-5,sum(stages)+170),
                    ylim=(-170,170) if full else (-28,28),xlabel='Outward distance from post-PNC Layer 10 start (um)',
                    ylabel='Transverse coordinate (um)',title=f'Design {d}: Layer 14 only')
                ax.legend(loc='upper left');fig.tight_layout();fig.savefig(out/f'design{d}_{"full" if full else "taper"}.png',dpi=170);plt.close(fig)
    report=dict(status='PASS',check_count=len(checks),checks=checks,windows=ends,
        scope='Original Design 1, original Design 2, and combined 18-device sweep only',
        changed_layers=[14],all_non_layer14_polygons_unchanged=True,pnc_region_unchanged=True,
        optical_or_mechanical_simulation_performed=False,
        files={fn:hashlib.sha256((root/'gds'/fn).read_bytes()).hexdigest() for fn in FILES})
    (root/'gds/2611_original_design_windows_validation.json').write_text(json.dumps(report,indent=2))
    print('PASS',len(checks),'checks;',len(ends),'taper ends')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--baseline',type=Path,required=True);p.add_argument('--render',action='store_true')
    a=p.parse_args();validate(a.root,a.baseline,a.render)
