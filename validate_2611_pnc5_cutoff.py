"""Validate all four exported layouts after truncating PNC at 5 um width.

Run with 'capture' before changing geometry, then without arguments after
regenerating all four GDSs. Baseline comparisons are performed only when
baseline files exist. Checks are geometric, not optical or foundry DRC.
"""
from pathlib import Path
import hashlib
import json
import shutil
import sys
import klayout.db as kdb
import numpy as np

ROOT = Path(__file__).resolve().parent
WORK = ROOT / 'validation_2611_pnc5_cutoff'
BASE = WORK / 'baseline'
END = 31.379
OLD_END = 45.0
DENSE = np.arange(6, dtype=float) + END - 5.0
FILES = [
    ('design1','2611_tapered_yue_dual_resist_L10_L12.gds',1400.,5.,1,1),
    ('design2','2611_tapered_yue_design2.gds',1400.,5.,1,2),
    ('combined','2611_design1_design2_L1400_18device_sweep.gds',1400.,5.,18,None),
    ('long_design1','2611_design1_L6000_P10_width_sweep.gds',6000.,10.,3,1),
]
REPORT = {'checks': [], 'failures': [], 'measurements': [], 'outputs': {}}


def check(name, ok, **details):
    REPORT['checks'].append(dict(name=name,passed=bool(ok),**details))
    if not ok:
        REPORT['failures'].append(name)
        print('FAIL',name,details,flush=True)


def read(path):
    ly=kdb.Layout();ly.read(str(path))
    if len(ly.top_cells())!=1: raise ValueError('Expected one top cell: '+str(path))
    return ly,ly.top_cell()


def region(ly,cell,n):
    i=ly.find_layer(n,0)
    return kdb.Region() if i is None else kdb.Region(cell.begin_shapes_rec(i)).merged()


def regs(ly,cell):
    return {n:region(ly,cell,n) for n in (2,10,11,13,14,15,16,51)}


def box(coords,dbu):
    return kdb.Region(kdb.Box(*(round(float(v)/dbu) for v in coords)))


def devices(ly,top):
    return sorted([i for i in top.each_inst() if not region(ly,i.cell,51).is_empty()
                   and not region(ly,i.cell,15).is_empty()],key=lambda i:i.trans.disp.y)


def pnc_polygons(r,L,pitch,dbu):
    roi=box((L/2-3*pitch-.4,-230,L/2+OLD_END+.4,-170),dbu)
    result={1:[],-1:[]}
    for p in (r&roi).merged().each():
        b=p.bbox();side=1 if (b.bottom+b.top)*dbu/2 > -200 else -1
        result[side].append(((b.left+b.right)*dbu/2-L/2,p))
    return {side:sorted(v,key=lambda z:z[0]) for side,v in result.items()}


def capture():
    BASE.mkdir(parents=True,exist_ok=True)
    for _,fn,*_ in FILES: shutil.copy2(ROOT/'gds'/fn,BASE/fn)
    print('Captured four GDS baselines.',flush=True)


def render(ly,inst,tag,L):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    r=regs(ly,inst.cell);dbu=ly.dbu
    roi=box((L/2-10,-216,L/2+110,-184),dbu)
    fig,ax=plt.subplots(figsize=(12,4.8))
    for n,alpha in ((14,.12),(2,.18),(10,.55),(11,.85),(15,.9)):
        vertices=[]
        for p in (r[n]&roi).each():
            vertices.append([(q.x*dbu-L/2,q.y*dbu+200) for q in p.each_point_hull()])
        if vertices:
            # Default Matplotlib color cycle; no rendered geometry is guessed.
            c=next(ax._get_lines._cycler_items)['color'] if False else None
            col=PolyCollection(vertices,alpha=alpha,label=f'Layer {n}')
            ax.add_collection(col)
    ax.axvline(END,linestyle='--',linewidth=1)
    ax.axvline(OLD_END,linestyle=':',linewidth=1)
    ax.text(END+1,14.7,'New 5 um endpoint',fontsize=9)
    ax.text(OLD_END+1,12.7,'Previous 3 um endpoint',fontsize=8)
    ax.set(xlim=(-10,110),ylim=(-16,17),xlabel='Distance from meander connection (um)',
           ylabel='Transverse coordinate (um)',title=f'{tag}: actual exported taper and masks')
    ax.legend(loc='lower right',ncol=5,fontsize=8)
    fig.tight_layout();fig.savefig(WORK/(tag+'_5um_endpoint.png'),dpi=180);plt.close(fig)


def run():
    WORK.mkdir(exist_ok=True)
    checked=0
    for tag,fn,L,pitch,count,design_default in FILES:
        ly,top=read(ROOT/'gds'/fn);dbu=ly.dbu
        insts=devices(ly,top)
        check(tag+'_device_count',len(insts)==count)
        check(tag+'_1nm_grid',abs(dbu-.001)<1e-12)
        oldly=oldtop=None;oldinsts=[]
        if (BASE/fn).exists():
            oldly,oldtop=read(BASE/fn);oldinsts=devices(oldly,oldtop)
            check(tag+'_baseline_count',len(oldinsts)==len(insts))
        meta_path=(ROOT/'gds'/fn).with_suffix('.json')
        metadata={}
        if tag in ('combined','long_design1'):
            metadata={e['case_id']:e for e in json.loads(meta_path.read_text())['devices']}
        for i,inst in enumerate(insts):
            checked+=1;name=f'{tag}_{i+1}'
            r=regs(ly,inst.cell)
            check(name+'_all_eight_layers',all(not v.is_empty() for v in r.values()))
            design=design_default or int(inst.cell.name[1])
            e=metadata.get(inst.cell.name,{})
            stages=e.get('stage_lengths_um', [40.,40.] if design==1 else [20.,100.,10.])
            wg=e.get('width_um',1.)
            widths10=[5.,1.,.2] if design==1 else [5.,1.65,.98,.2]
            widths15=[wg,.1] if design==1 else [1.,.5,.2]
            lens15=stages[:1] if design==1 else stages[:2]
            check(name+'_aligned_left_bend',abs((r[11].bbox().left+inst.trans.disp.x)*dbu+900)<.0011)
            check(name+'_expected_y',abs(inst.trans.disp.y*dbu-(-700+600*i))<.0011)
            if oldly is not None and i<len(oldinsts):
                old=regs(oldly,oldinsts[i].cell)
                check(name+'_placement_unchanged',inst.trans==oldinsts[i].trans)
                right=box((L/2-.5,-500,20000,-25),dbu)
                allowed=right|right.transformed(kdb.Trans(kdb.Trans.R180))
                # Only the two end assemblies are permitted to change.
                for n in r:
                    area=((r[n]^old[n])-allowed).area()*dbu**2
                    check(name+f'_L{n}_outside_ends_unchanged',area<.000001,xor_um2=area)
            else: old=None
            for side in ('right','left'):
                rr=r if side=='right' else {n:v.transformed(kdb.Trans(kdb.Trans.R180)) for n,v in r.items()}
                oo=old if side=='right' or old is None else {n:v.transformed(kdb.Trans(kdb.Trans.R180)) for n,v in old.items()}
                prefix=name+'_'+side
                pncs=pnc_polygons(rr[11],L,pitch,dbu)
                regular=np.arange(0,END-5-1e-9,pitch)
                expected=np.r_[-np.arange(3,0,-1)*pitch,regular,DENSE]
                for sign,items in pncs.items():
                    xs=np.array([v[0] for v in items])
                    check(prefix+f'_PNC_side{sign}_positions',len(xs)==len(expected) and
                          np.allclose(xs,expected,atol=.0011,rtol=0),centers_um=xs.tolist())
                    terminal=[p for x,p in items if abs(x-END)<.0011]
                    check(prefix+f'_last_PNC_side{sign}_exists',len(terminal)==1)
                    if oo is not None and len(terminal)==1:
                        previous=pnc_polygons(oo[11],L,pitch,dbu)[sign]
                        last=[p for x,p in previous if abs(x-OLD_END)<.0011]
                        if len(last)==1:
                            # Moving the previous 3 um-wide endpoint anchors out
                            # by 1 um on each side gives the exact 5 um junction.
                            tr=kdb.Trans(round((END-OLD_END)/dbu),round(sign*1./dbu))
                            delta=kdb.Region(terminal[0])^kdb.Region(last[0]).transformed(tr)
                            check(prefix+f'_end_PNC_side{sign}_5um_anchor',delta.is_empty(),xor_um2=delta.area()*dbu**2)
                gap=rr[11]&box((L/2+END+.4,-230,L/2+OLD_END+.4,-170),dbu)
                check(prefix+'_no_PNC_in_removed_narrow_region',gap.is_empty())
                # Independent widths from exported taper polygons, not metadata.
                for n,lengths,widths in ((10,stages,widths10),(15,lens15,widths15)):
                    edges=END+np.r_[0.,np.cumsum(lengths)]
                    half=3. if n==10 else 2.
                    for j,(x,w) in enumerate(zip(edges,widths)):
                        x+=L/2+(.002 if j==0 else -.002)
                        cr=(rr[n]&box((x,-200-half,x+dbu,-200+half),dbu)).merged()
                        got=cr.bbox().height()*dbu if cr.count()==1 else None
                        check(prefix+f'_L{n}_stage{j}_width',got is not None and abs(got-w)<.0041,
                              measured_um=got,target_um=w)
                slab_roi=box((L/2+END+.001,-203,L/2+END+sum(stages)-.001,-197),dbu)
                slab=(rr[10]&slab_roi).merged()
                check(prefix+'_slab_taper_connected',slab.count()==1)
                check(prefix+'_L14_covers_new_slab_taper',(slab-rr[14]).area()*dbu**2<.005)
                before=box((L/2+.5,-202.49,L/2+END-.002,-197.51),dbu)
                check(prefix+'_L10_starts_only_at_new_endpoint',(rr[10]&before).is_empty())
                join=box((L/2+END-.15,-202.49,L/2+END+.15,-197.51),dbu)
                check(prefix+'_L14_full_5um_junction',(join-rr[14]).is_empty())
                check(prefix+'_local_L16_covers_ridge',((rr[15]&box((L/2-1,-202,L/2+END+1,-198),dbu))-rr[16]).area()*dbu**2<.005)
                REPORT['measurements'].append(dict(device=name,end=side,cutoff_um=END,width_um=5.,
                    dense_stations_um=DENSE.tolist(),stage_lengths_um=stages))
        if tag in ('design1','design2') and insts: render(ly,insts[0],tag,L)
        REPORT['outputs'][fn]={'sha256':hashlib.sha256((ROOT/'gds'/fn).read_bytes()).hexdigest(),
                              'devices':len(insts),'cutoff_width_um':5.,'pnc_length_um':END}
    REPORT.update(status='PASS' if not REPORT['failures'] else 'FAIL',
                  check_count=len(REPORT['checks']),devices_checked=checked,taper_ends_checked=2*checked,
                  cutoff_width_um=5.,pnc_length_um=END,previous_length_um=OLD_END,
                  shortening_um=OLD_END-END,dense_pitch_um=1.,dense_periods=5,
                  baseline_comparison_performed=all((BASE/f).exists() for _,f,*_ in FILES))
    (ROOT/'gds/2611_pnc5_cutoff_validation.json').write_text(json.dumps(REPORT,indent=2))
    (ROOT/'gds/2611_taper_pitch_validation.json').write_text(json.dumps(REPORT,indent=2))
    (WORK/'validation_report.json').write_text(json.dumps(REPORT,indent=2))
    print(REPORT['status'],REPORT['check_count'],'checks;',len(REPORT['failures']),'failures',flush=True)
    if REPORT['failures']: raise SystemExit(1)


if __name__=='__main__':
    capture() if len(sys.argv)>1 and sys.argv[1]=='capture' else run()
