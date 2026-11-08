"""Independent raw-GDS/Shapely checks for the twelve-device L6000 layout."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from shapely.geometry import box, LineString, Polygon
from shapely.affinity import translate, scale
from shapely.ops import unary_union
from validate_2611_straight_output_raw import read, union

STEM='2611_design1_L6000_P10_straight_output_comparison'
DBU=.001


def validate(root: Path, baseline: Path | None = None) -> dict:
    path=root/'gds'/(STEM+'.gds')
    cells,top=read(path)
    old,oldtop=read(baseline) if baseline else (None,None)
    meta=json.loads(path.with_suffix('.json').read_text())
    refs=dict(cells[top]['refs']);checks=[];measured=[];cache={}
    devices={n:p for n,p in refs.items() if n.startswith('D1_')}
    labels={n:p for n,p in refs.items() if n.startswith('LABEL_')}
    def check(n,ok,**detail):
        checks.append(dict(name=n,passed=bool(ok),**detail))
        if not ok:raise AssertionError((n,detail))
    def exact(n,a,b,tol=1e-6):
        area=a.symmetric_difference(b).area
        check(n,area<tol,xor_um2=area)
    def geom(n,l):
        key=(n,l)
        if key not in cache:cache[key]=union(cells[n],l)
        return cache[key]
    check('12_devices_12_labels_one_frame',len(devices)==12 and len(labels)==12 and len(refs)==25)
    check('600um_row_pitch',sorted(p[1] for p in devices.values())==[-700000+i*600000 for i in range(12)])
    check('original_x_alignment_preserved',all(p[0]==2242153 for p in devices.values()))
    check('metadata_requested_dimensions',meta['layer2_start_width_um']==15 and meta['layer14_start_width_um']==8 and
          meta['layer2_window_margin_um']==2.5 and meta['layer14_window_margin_um']==1.0 and
          meta['layer2_start_offset_um']==.5 and meta['layer14_start_offset_um']==0)
    for n,c in cells.items():
        check(n+'_removed_layers_absent',not ({16,51}&set(c['polys'])))
    check('labels_on_layer2',all(set(cells[n]['polys'])=={2} for n in labels))
    check('global_markers_on15_not2',bool(cells['FRAME_AND_MARKERS']['polys'].get(15)) and not cells['FRAME_AND_MARKERS']['polys'].get(2))
    if old:
        for name,pos in old[oldtop]['refs']:
            check(name+'_original_placement',refs[name]==pos)
            for layer in old[name]['polys']:
                a=sorted(tuple(x) for x in old[name]['polys'][layer])
                b=sorted(tuple(x) for x in cells[name]['polys'].get(layer,[]))
                check(name+f'_L{layer}_original_raw_coordinates',a==b)
    if old:
        check('all_cell_names_preserved',set(cells)==set(old))
        check('all_cell_references_preserved',all(cells[n]['refs']==old[n]['refs'] for n in cells))
    for name,pos in devices.items():
        single=name.startswith('D1_SINGLE_')
        one=('_1T_' in name) if single else name.endswith('_OUT2000')
        width=int(name.split('_W')[1].split('_')[0])/1000
        check(name+'_requested_layers',set(cells[name]['polys'])=={2,10,11,13,14,15})
        ends=[('input',-1,0 if single else 200)]+([] if one else [('output',1,0 if single else -200)])
        keep={2:Polygon(),14:Polygon()}
        for side,sign,cy in ends:
            key=name+'_'+side;sx=sign*3045
            expected_widths={2:15.,14:8.};offsets={2:.5,14:0.};margins={2:2.5,14:1.0}
            search=box(-20000,-1000 if single else 25,-3000,1000) if sign<0 else box(3000,-1000,20000,1000 if single else -25)
            slab=geom(name,10).intersection(search)
            actual_sx=slab.bounds[2] if sign<0 else slab.bounds[0]
            check(key+'_L10_entrance',abs(actual_sx-sx)<1e-8,x_um=actual_sx)
            for layer in (2,14):
                offset=offsets[layer];start=sx+sign*offset
                end=box(-20000,-1000 if single else 25,start,1000) if sign<0 else box(start,-1000,20000,1000 if single else -25)
                keep[layer]=keep[layer].union(end)
                mask=geom(name,layer).intersection(end)
                edge=mask.bounds[2] if sign<0 else mask.bounds[0]
                check(key+f'_L{layer}_start',abs(edge-start)<1e-8,local_x_um=edge,global_x_um=edge+pos[0]*DBU)
                pre=box(min(sx-sign,start)+.001,cy-20,max(sx-sign,start)-.001,cy+20)
                check(key+f'_L{layer}_no_early_mask',geom(name,layer).intersection(pre).area<1e-9)
                for dist in (offset+.001,offset+.1,1.,5.,10.,20.,30.,39.9):
                    x=sx+sign*dist
                    cut=geom(name,layer).intersection(LineString([(x,cy-20),(x,cy+20)]))
                    t=(dist-offset)/(40-offset);f=t**3*(10-15*t+6*t*t)
                    expected=expected_widths[layer]+(1+2*margins[layer]-expected_widths[layer])*f
                    check(key+f'_L{layer}_profile_{dist:g}',cut.geom_type=='LineString' and abs(cut.length-expected)<.0025,
                          width_um=cut.length,expected_um=expected)
                entrance=geom(name,layer).intersection(LineString([(start,cy-20),(start,cy+20)]))
                check(key+f'_L{layer}_exact_entrance_width',abs(entrance.length-expected_widths[layer])<1e-8,width_um=entrance.length)
                for dist in (50,60,70,75):
                    x=sx+sign*dist;line=LineString([(x,cy-20),(x,cy+20)])
                    optical=geom(name,10).intersection(line);m=geom(name,layer).intersection(line)
                    lower=optical.bounds[1]-m.bounds[1];upper=m.bounds[3]-optical.bounds[3]
                    check(key+f'_L{layer}_margin_{dist}',abs(lower-margins[layer])<.0025 and abs(upper-margins[layer])<.0025,
                          lower_um=lower,upper_um=upper)
                measured.append(dict(case_id=name,side=side,layer=layer,start_global_x_um=edge+pos[0]*DBU,
                                     start_width_um=entrance.length,offset_um=offset))
                if old and not single:
                    far=box(-20000,25,sx-96,20000) if sign<0 else box(sx+96,-20000,20000,-25)
                    exact(key+f'_L{layer}_far_channel_pad_preserved',geom(name,layer).intersection(far),union(old[name],layer).intersection(far))
        for layer in (2,14):
            exact(name+f'_L{layer}_terminal_only',geom(name,layer).difference(keep[layer]),Polygon())
        r=geom(name,15);parts=list(r.geoms) if hasattr(r,'geoms') else [r]
        check(name+'_connected_path_and_ticks',len(parts)==(3 if one else 5))
        if single:
            source='D1_L6000_W'+str(round(width*1000))+'_P10'+('_OUT2000' if one else '')
            for layer in (2,10,11,13,14,15):
                left=translate(geom(source,layer).intersection(box(-20000,25,-3000,20000)),yoff=-200)
                right=translate(geom(source,layer).intersection(box(3000,-20000,20000,-25)),yoff=200)
                expected=left.union(right)
                if layer==15:expected=expected.union(box(-3000,-width/2,3000,width/2))
                exact(name+f'_L{layer}_exact_ends_and_single_route',geom(name,layer),expected)
            middle=box(-3000,-100,3000,100)
            exact(name+'_6000um_rectangle_no_bends',r.intersection(middle),box(-3000,-width/2,3000,width/2))
            for layer in (2,10,13,14):
                check(name+f'_no_L{layer}_on_center_straight',geom(name,layer).intersection(box(-2999.99,-180,2999.99,180)).area<1e-9)
            if one:
                exact(name+'_2000um_extension',r.intersection(box(3000,-100,20000,100)),box(3000,-width/2,5000,width/2))
        for x in (-2500,0,2500):
            for cy in ([0] if single else [-200,0,200]):
                cut=r.intersection(LineString([(x,cy-4),(x,cy+4)]))
                check(name+f'_ridge_width_{x}_{cy}',abs(cut.length-width)<1e-8,width_um=cut.length)
        if not one:
            for layer in (2,14):
                left=geom(name,layer).intersection(keep[layer].intersection(box(-20000,-20000,-3000,20000)))
                right=geom(name,layer).intersection(box(3000,-20000,20000,20000))
                err=scale(left,xfact=-1,yfact=-1,origin=(0,0)).hausdorff_distance(right)
                check(name+f'_L{layer}_symmetric_on_grid',err<.0011,edge_error_um=err)
    check('18_coupler_windows_36_mask_entrances',len(measured)==36)
    placed={n:unary_union([translate(geom(n,l),xoff=p[0]*DBU,yoff=p[1]*DBU) for l in (2,10,11,13,14,15)]) for n,p in devices.items()}
    pairs=list(placed.items())
    for i,(n,a) in enumerate(pairs):
        for m,b in pairs[i+1:]:check(n+'_'+m+'_no_overlap',a.intersection(b).area<1e-6)
    for n,p in labels.items():
        label=translate(union(cells[n],2),xoff=p[0]*DBU,yoff=p[1]*DBU)
        check(n+'_clear',all(label.intersection(r).area<1e-6 for r in placed.values()))
    for e in meta['devices']:
        single=e['topology']=='single_straight';one=e['added_length_um']>0
        check(e['case_id']+'_path_metadata',e['straight_section_count']==(1 if single else 3) and e['bend_count']==(0 if single else 2))
        if single:check(e['case_id']+'_single_path_length',e['path_length_excluding_tapers_um']==(8000 if one else 6000))
    from validate_2611_taper_beams import validate as validate_beams
    beam_report=validate_beams(root,baseline)
    checks.extend(beam_report['checks'])
    report=dict(status='PASS',check_count=len(checks),checks=checks,measured_entrances=measured,
                device_count=12,coupler_count=18,removed_layers=[16,51],
                layer2_start_width_um=15.,layer14_start_width_um=8.,downstream_margin_per_side_um={'2':2.5,'14':1.0},
                layer2_outward_offset_um=.5,layer14_outward_offset_um=0.,
                baseline_compared=bool(old),gds_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                method='Independent raw GDSII parser and Shapely; no CAD generator or KLayout',
                optical_or_mechanical_simulation_performed=False)
    (root/'gds/2611_four_sets_validation.json').write_text(json.dumps(report,indent=2))
    print('PASS',len(checks),'independent checks')
    return report


def render(root: Path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cells,top=read(root/'gds'/(STEM+'.gds'))
    dest=root/'validation_2611_four_sets';dest.mkdir(exist_ok=True)
    def draw(ax,r,label,style='-',lw=1):
        geoms=list(r.geoms) if hasattr(r,'geoms') else [r];xs=[];ys=[]
        for p in geoms:
            if p.geom_type!='Polygon':continue
            for ring in (p.exterior,*p.interiors):
                x,y=ring.xy;xs.extend(list(x)+[float('nan')]);ys.extend(list(y)+[float('nan')])
        if xs:ax.plot(xs,ys,style,linewidth=lw,label=label)
    name='D1_L6000_W1620_P10'
    for side,sign,cy in [('input',-1,200),('output',1,-200)]:
        sx=sign*3045
        fig,ax=plt.subplots(figsize=(13,5))
        for layer,style in ((2,'--'),(14,'-'),(10,':'),(15,'-')):
            draw(ax,union(cells[name],layer).intersection(box(min(sx-5,sx+sign*82),cy-10,max(sx+5,sx+sign*82),cy+10)),f'Layer {layer}',style)
        ax.axvline(sx,linestyle=':',linewidth=.8);ax.axvline(sx+sign*.5,linestyle='--',linewidth=.8)
        ax.set(xlabel='Local X (um)',ylabel='Local Y (um)',title=f'{side.title()}: L2 15 um inlet, 0.5 um offset; L14 8 um inlet; downstream margins L2=2.5 um / L14=1 um')
        ax.legend(loc='lower right');fig.tight_layout();fig.savefig(dest/f'{side}_masks.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(16,10))
    meta=json.loads((root/'gds'/(STEM+'.json')).read_text())
    refs=dict(cells[top]['refs'])
    for layer in (14,10,15):
        shapes=[]
        for e in meta['devices']:
            x,y=refs[e['case_id']];shapes.append(translate(union(cells[e['case_id']],layer),xoff=x*DBU,yoff=y*DBU))
        draw(ax,unary_union(shapes),f'Layer {layer}',lw=.6)
    for e in meta['devices']:
        x,y=e['label_origin_um'];ax.text(x,y+15,e['label'],fontsize=7)
    ax.set(xlabel='X (um)',ylabel='Y (um)',title='12 devices: meander two/one coupler, then single straight two/one coupler')
    ax.set_aspect('equal');ax.legend();fig.tight_layout();fig.savefig(dest/'overview.png',dpi=190);plt.close(fig)
    for name in ('D1_SINGLE_L6000_W1620_2T','D1_SINGLE_L6000_W1620_1T_OUT2000'):
        fig,ax=plt.subplots(figsize=(16,3))
        for layer in (2,14,10,15):draw(ax,union(cells[name],layer),f'Layer {layer}',lw=.8)
        ax.set(xlabel='Local X (um)',ylabel='Local Y (um)',title=name)
        ax.legend();fig.tight_layout();fig.savefig(dest/(name+'.png'),dpi=170);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--baseline',type=Path);p.add_argument('--render',action='store_true');a=p.parse_args()
    validate(a.root,a.baseline)
    if a.render:render(a.root)
