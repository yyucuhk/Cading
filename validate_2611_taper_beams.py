"""Independent GDSII/Shapely checks for plain local taper supports.

No import of the CAD generator, KLayout, or its support-placement helper.
The fixed expected dimensions are the requested release parameters.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from shapely.geometry import Polygon, box
from shapely.affinity import translate, scale
from shapely.ops import unary_union
from validate_2611_straight_output_raw import read, union

STEM = '2611_design1_L6000_P10_straight_output_comparison'
DBU = 0.001


def validate(root: Path, baseline: Path | None = None) -> dict:
    path = root/'gds'/(STEM+'.gds')
    cells, top = read(path)
    refs = dict(cells[top]['refs'])
    devices = {name: pos for name,pos in refs.items() if name.startswith('D1_')}
    checks = []
    measurements = []
    cache = {}
    def check(name, ok, **details):
        checks.append(dict(name=name,passed=bool(ok),**details))
        if not ok:
            raise AssertionError((name,details))
    def geom(name,layer):
        key = (name,layer)
        if key not in cache:
            cache[key] = union(cells[name],layer)
        return cache[key]
    check('twelve_devices',len(devices)==12)
    check('no_L16_or_L51_anywhere',all(not ({16,51}&set(c['polys'])) for c in cells.values()))
    check('L11_only_in_device_cells',all(not c['polys'].get(11) for n,c in cells.items() if n not in devices))
    total_beams = 0
    total_ends = 0
    for name,pos in devices.items():
        single = name.startswith('D1_SINGLE_')
        one = '_1T_' in name if single else name.endswith('_OUT2000')
        ends = [('input',-1,0.0 if single else 200.0)]
        if not one:
            ends.append(('output',1,0.0 if single else -200.0))
        support = geom(name,11)
        parts = list(support.geoms) if hasattr(support,'geoms') else [support]
        check(name+'_beam_count',len(parts)==6*len(ends),count=len(parts))
        total_beams += len(parts)
        allowed = Polygon()
        for end,direction,cy in ends:
            total_ends += 1
            sx = direction*3045.0
            roi = box(min(sx-0.251,sx+direction*2.251),cy-9,
                      max(sx+0.251,sx+direction*2.251),cy+9)
            allowed = allowed.union(roi)
            group = [p for p in parts if p.intersects(roi)]
            check(name+'_'+end+'_six_plain_rectangles',len(group)==6)
            for side in (-1,1):
                row = sorted((p for p in group if (p.centroid.y-cy)*side>0),
                             key=lambda p:(p.centroid.x-sx)*direction)
                check(name+'_'+end+f'_side{side}_three_beams',len(row)==3)
                centers = [(p.bounds[0]+p.bounds[2])/2 for p in row]
                check(name+'_'+end+f'_side{side}_pitch',
                      all(abs(abs(centers[i+1]-centers[i])-1.0)<1e-8 for i in range(2)))
                for station,p in enumerate(row):
                    key = name+'_'+end+f'_side{side}_beam{station}'
                    x1,y1,x2,y2 = p.bounds
                    distance = ((x1+x2)/2-sx)*direction
                    check(key+'_center',abs(distance-station)<1e-8,distance_um=distance)
                    check(key+'_width',abs(x2-x1-0.5)<1e-8,width_um=x2-x1)
                    check(key+'_rectangle_no_PNC_masses',p.symmetric_difference(box(*p.bounds)).area<1e-9)
                    check(key+'_original_radial_span',abs(y2-y1-6.096)<.0011,length_um=y2-y1)
                    # Original 5->1 quintic slab law; max width occurs at the
                    # first point of the 0.5 um-wide station, away from the tip.
                    near = max(0.0,station-0.25)
                    t=near/40.0
                    smooth=6*t**5-15*t**4+10*t**3
                    halfwidth=(5.0-4.0*smooth)/2
                    near_y=min(abs(y1-cy),abs(y2-cy))
                    far_y=max(abs(y1-cy),abs(y2-cy))
                    check(key+'_follows_slab_edge',abs(near_y-(halfwidth-.5))<.0021 and
                          abs(far_y-(halfwidth+5.596))<.0021)
                    contact=p.intersection(geom(name,10)).area
                    check(key+'_positive_slab_overlap',contact>.12,overlap_area_um2=contact)
                    check(key+'_clear_of_ridge',p.intersection(geom(name,15)).area<1e-9)
                    # Half-micron strip at its outer end must be outside BOTH
                    # exposure windows, not simply at the smaller L14 boundary.
                    anchor=box(x1,y2-.5,x2,y2) if side>0 else box(x1,y1,x2,y1+.5)
                    release=geom(name,2).union(geom(name,14))
                    check(key+'_outer_anchor_beyond_both_masks',anchor.intersection(release).area<1e-9)
                    measurements.append(dict(case_id=name,end=end,side=side,station=station,
                        width_um=x2-x1,length_um=y2-y1,distance_from_l10_um=distance,
                        center_global_x_um=(x1+x2)/2+pos[0]*DBU,bounds_local_um=list(p.bounds)))
            slabs=geom(name,10).intersection(box(min(sx,sx+direction*4),cy-10,max(sx,sx+direction*4),cy+10))
            connected=slabs.union(unary_union(group))
            check(name+'_'+end+'_one_connected_slab_plus_beams',connected.geom_type=='Polygon')
        check(name+'_no_L11_elsewhere',support.difference(allowed).area<1e-9)
        if one:
            cy=0 if single else -200
            check(name+'_bare_output_has_no_beams',support.intersection(box(3000,cy-20,5100,cy+20)).is_empty)
        if not one:
            left=support.intersection(box(-4000,-400,-3000,400))
            right=support.intersection(box(3000,-400,4000,400))
            error=scale(left,xfact=-1,yfact=-1,origin=(0,0)).symmetric_difference(right).area
            check(name+'_opposite_ends_symmetric',error<1e-8,xor_um2=error)
        for i,p in enumerate(parts):
            check(name+f'_beam{i}_separated',all(p.distance(q)>.49 for q in parts[i+1:]))
    check('18_supported_couplers',total_ends==18,count=total_ends)
    check('108_plain_beams',total_beams==108,count=total_beams)
    if baseline:
        before,oldtop=read(baseline)
        check('all_cells_preserved',set(before)==set(cells))
        for name,c in before.items():
            check(name+'_references_unchanged',c['refs']==cells[name]['refs'])
            for layer,polygons in c['polys'].items():
                a=sorted(tuple(p) for p in polygons)
                b=sorted(tuple(p) for p in cells[name]['polys'].get(layer,[]))
                check(name+f'_L{layer}_raw_polygons_unchanged',a==b)
            added=set(cells[name]['polys'])-set(c['polys'])
            check(name+'_only_L11_added',added==({11} if name in devices else set()))
    report=dict(status='PASS',check_count=len(checks),checks=checks,beams=measurements,
                total_beams=total_beams,couplers=total_ends,support_layer=11,
                width_um=.5,pitch_um=1.,station_offsets_um=[0.,1.,2.],
                anchor_to_anchor_span_um=5.096,anchor_overlap_um=.5,
                layer14_margin_um=1.,layer2_margin_um=2.5,
                baseline_compared=baseline is not None,
                gds_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                method='Raw GDSII parser and independent Shapely geometry; no CAD generator',
                optical_or_mechanical_simulation_performed=False)
    (root/'gds/2611_taper_beams_validation.json').write_text(json.dumps(report,indent=2))
    print('PASS',len(checks),'independent local-support checks')
    return report


def render(root: Path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    cells,top=read(root/'gds'/(STEM+'.gds'))
    refs=dict(cells[top]['refs'])
    name='D1_L6000_W1620_P10'
    xshift=refs[name][0]*DBU
    dest=root/'validation_2611_taper_beams'
    dest.mkdir(exist_ok=True)
    for side,direction,cy in [('input',-1,200.),('output',1,-200.)]:
        sx=direction*3045
        roi=box(min(sx-1,sx+direction*4),cy-9,max(sx+1,sx+direction*4),cy+9)
        fig,ax=plt.subplots(figsize=(10,6))
        for layer,style in [(2,'--'),(14,'-'),(10,'-'),(15,':'),(11,'-')]:
            r=translate(union(cells[name],layer).intersection(roi),xoff=xshift,yoff=-cy)
            geoms=list(r.geoms) if hasattr(r,'geoms') else [r]
            if layer==11:
                patches=[list(p.exterior.coords) for p in geoms if p.geom_type=='Polygon']
                ax.add_collection(PolyCollection(patches,alpha=.35,label='L11: plain 0.5 um beams'))
                continue
            first=True
            for p in geoms:
                if p.geom_type!='Polygon':continue
                xx,yy=p.exterior.xy
                ax.plot(xx,yy,linestyle=style,linewidth=1.3,label=f'Layer {layer}' if first else None)
                first=False
        for i in range(3):
            x=sx+xshift+direction*i
            ax.text(x,8.5,f'{i} um',ha='center',fontsize=10)
        ax.set(xlabel='Global X (um)',ylabel='Y relative to waveguide center (um)',
               title=f'{side.title()} taper: 3 beam pairs, 0.5 um width, 1 um pitch')
        ax.legend(loc='center left',bbox_to_anchor=(1.01,.5));ax.grid(alpha=.2)
        fig.tight_layout();fig.savefig(dest/f'{side}_beams.png',dpi=170);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--baseline',type=Path)
    parser.add_argument('--render',action='store_true')
    args=parser.parse_args()
    validate(args.root,args.baseline)
    if args.render:render(args.root)
