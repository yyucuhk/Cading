"""Independent raw-GDS/Shapely audit of the Layer-14-only linear window."""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
from shapely.geometry import Polygon, box, LineString
from shapely.affinity import affine_transform
from shapely.ops import unary_union
from validate_2611_straight_output_raw import read, union

STEM='2611_design1_L6000_P10_straight_output_comparison'
JOIN=58.847


def validate(root: Path, baseline: Path):
    path=root/'gds'/(STEM+'.gds')
    cells,top=read(path);old,oldtop=read(baseline)
    meta=json.loads(path.with_suffix('.json').read_text())
    checks=[];records=[]
    def check(n,ok,**details):
        checks.append(dict(name=n,passed=bool(ok),**details))
        if not ok:raise AssertionError((n,details))
    check('linear_revision_metadata',meta.get('linear_layer14_windows',{}).get('enabled',False))
    check('all_cells_and_top_preserved',set(cells)==set(old) and top==oldtop)
    for n,c in cells.items():
        check(n+'_references_preserved',c['refs']==old[n]['refs'])
        check(n+'_layer_set_preserved',set(c['polys'])==set(old[n]['polys']))
        for l,polys in c['polys'].items():
            if l==14:continue
            check(n+f'_L{l}_raw_coordinates_exact',sorted(map(tuple,polys))==sorted(map(tuple,old[n]['polys'][l])))
        check(n+'_no_L16_L51',not ({16,51}&set(c['polys'])))
    refs=dict(cells[top]['refs']);devices={n:p for n,p in refs.items() if n.startswith('D1_')}
    check('12_devices_and_12_labels',len(devices)==12 and sum(n.startswith('LABEL_') for n in refs)==12)
    for name,pos in devices.items():
        single=name.startswith('D1_SINGLE_');one='_1T_' in name if single else name.endswith('_OUT2000')
        ends=[('input',-1,0 if single else 200)]+([] if one else [('output',1,0 if single else -200)])
        g=union(cells[name],14);b=union(old[name],14)
        for end,sign,cy in ends:
            sx=sign*3045.;matrix=[sign,0,0,sign,-sign*sx,-sign*cy]
            now=affine_transform(g,matrix).intersection(box(0,-175,20000,175))
            before=affine_transform(b,matrix).intersection(box(0,-175,20000,175))
            key=name+'_'+end
            check(key+'_single_connected_window',now.geom_type=='Polygon' and len(now.interiors)==0)
            cut=before.intersection(LineString([(JOIN,-175),(JOIN,175)]))
            lo,hi=round(cut.bounds[1],3),round(cut.bounds[3],3)
            far=before.bounds[2]
            tail=before.intersection(LineString([(far,-175),(far,175)]))
            elo,ehi=tail.bounds[1],tail.bounds[3]
            channel=Polygon([(JOIN,lo),(far,elo),(far,ehi),(JOIN,hi)])
            circles=unary_union([Polygon([(round(cx+150*math.cos(math.radians(i*2.5)),3),
                                          round(150*math.sin(math.radians(i*2.5)),3)) for i in range(144)]) for cx in (530.,1030.)])
            expected=before.intersection(box(0,-175,JOIN,175)).union(channel).union(circles)
            area=now.symmetric_difference(expected).area
            distance=now.boundary.hausdorff_distance(expected.boundary)
            # Boolean intersections are snapped to the 1 nm GDS grid. A
            # sub-grid edge displacement can accumulate area over millimetres;
            # bound both pointwise error and area using the same grid scale.
            grid_tolerance = .0011
            difference = now.symmetric_difference(expected)
            outside_grid_band = difference.difference(
                expected.boundary.buffer(grid_tolerance)).area
            area_bound = expected.boundary.length * grid_tolerance
            check(key+'_only_prefix_linear_channel_two_circles',
                  area <= area_bound and distance < grid_tolerance and outside_grid_band < 1e-8,
                  xor_um2=area, area_bound_um2=area_bound,
                  boundary_distance_um=distance, outside_grid_band_um2=outside_grid_band)
            upstream=box(0,-175,JOIN-.5,175)
            check(key+'_upstream_exact',now.intersection(upstream).symmetric_difference(before.intersection(upstream)).area<1e-7)
            actual=now.intersection(LineString([(JOIN,-175),(JOIN,175)]))
            check(key+'_no_width_step_at_join',abs(actual.bounds[1]-lo)<.0011 and abs(actual.bounds[3]-hi)<.0011,
                  preserved_width_um=hi-lo)
            check(key+'_outer_endpoint_unchanged',abs(now.bounds[2]-far)<1e-9)
            for q in (JOIN+.1,60,65,75,80,90,95,100,150,175,200,250,300,750,800,1250,1500,2000,far):
                fraction=(q-JOIN)/(far-JOIN)
                lower=lo+(elo-lo)*fraction;upper=hi+(ehi-hi)*fraction
                cut=now.intersection(LineString([(q,-175),(q,175)]))
                check(key+f'_linear_edges_{q:g}',cut.geom_type=='LineString' and
                      abs(cut.bounds[1]-lower)<.0015 and abs(cut.bounds[3]-upper)<.0015,
                      lower_um=cut.bounds[1],upper_um=cut.bounds[3])
            frame=box(175,-25,235,25)
            excess=now.intersection(frame).difference(channel.buffer(.0011)).area
            check(key+'_no_intermediate_frame',excess<1e-8,extra_area_um2=excess)
            for cx in (530.,1030.):
                section=now.intersection(LineString([(cx,-175),(cx,175)]))
                check(key+f'_circle_{cx:g}_diameter_300',abs(section.length-300)<1e-9)
                for angle in (30,60,90,120,150,210,240,270,300,330):
                    x=round(cx+150*math.cos(math.radians(angle)),3);y=round(150*math.sin(math.radians(angle)),3)
                    from shapely.geometry import Point
                    check(key+f'_circle_{cx:g}_arc_{angle}',now.boundary.distance(Point(x,y))<.0011)
            check(key+'_circle_extents_unchanged',now.bounds[1]==before.bounds[1] and now.bounds[3]==before.bounds[3])
            records.append(dict(case_id=name,end=end,start_global_x_um=pos[0]*.001+sx+sign*JOIN,
                                start_width_um=hi-lo,circle_diameter_um=300.,
                                circle_centers_global_um=[[pos[0]*.001+sx+sign*q,pos[1]*.001+cy] for q in (530.,1030.)]))
    check('18_windows_36_circles',len(records)==18)
    report=dict(status='PASS',check_count=len(checks),checks=checks,windows=records,
                scope='Layer 14 only, all 18 existing coupler ends in the 12-device comparison',
                unchanged_layers='Every non-L14 raw polygon and every cell transform',
                linear_start_from_l10_um=JOIN,circle_diameter_um=300.,
                method='Independent raw GDSII parser and Shapely; no generator or KLayout import',
                baseline_sha256=hashlib.sha256(baseline.read_bytes()).hexdigest(),
                gds_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                optical_or_mechanical_simulation_performed=False,
                grid_note='Independent unsnapped geometry compared to 1 nm GDS polygons; max edge error and the entire difference region are bounded within 1.1 nm.')
    (root/'gds/2611_linear_window_validation.json').write_text(json.dumps(report,indent=2))
    print('PASS',len(checks),'independent linear-window checks')
    return report


def render(root: Path, baseline: Path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cells,top=read(root/'gds'/(STEM+'.gds'));old,_=read(baseline)
    name='D1_L6000_W1620_P10';pos=dict(cells[top]['refs'])[name]
    out=root/'validation_2611_linear_window';out.mkdir(exist_ok=True)
    for side,sign,cy in [('output',1,-200),('input',-1,200)]:
        matrix=[1,0,0,1,pos[0]*.001,-cy]
        roi=box(3000,-375,5200,-25) if sign>0 else box(-5200,25,-3000,375)
        oldg=affine_transform(union(old[name],14).intersection(roi),matrix)
        newg=affine_transform(union(cells[name],14).intersection(roi),matrix)
        for zoom in (False,True):
            fig,ax=plt.subplots(figsize=(13,3.6 if not zoom else 4.5))
            for geom,style,label in [(oldg,'--','Previous Layer 14'),(newg,'-','Updated Layer 14')]:
                parts=list(geom.geoms) if hasattr(geom,'geoms') else [geom]
                for i,p in enumerate(parts):
                    x,y=p.exterior.xy;ax.plot(x,y,ls=style,lw=1.4,label=label if i==0 else None)
            start=pos[0]*.001+sign*(3045+JOIN)
            ax.axvline(start,ls=':',lw=.9)
            if zoom:
                ax.set_xlim(*sorted((start-sign*50,start+sign*210)));ax.set_ylim(-30,30)
            else:
                ax.set_xlim(newg.bounds[0]-10,newg.bounds[2]+10);ax.set_ylim(-165,165);ax.set_aspect('equal',adjustable='box')
            ax.set_xlabel('Global x (um)');ax.set_ylabel('Offset from taper centerline (um)')
            ax.legend(loc='best');ax.grid(alpha=.15);fig.tight_layout()
            fig.savefig(out/(side+('_junction' if zoom else '_window')+'.png'),dpi=170);plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--baseline',type=Path,required=True);p.add_argument('--render',action='store_true')
    a=p.parse_args();validate(a.root,a.baseline)
    if a.render:render(a.root,a.baseline)
