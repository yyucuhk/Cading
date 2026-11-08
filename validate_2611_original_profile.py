"""Validate full original-profile 15->5 um PNC and downstream tapers.

Re-read saved polygons; never repair them. Optional original3/ baselines
refer to commit 6d4fcb56, before truncation and derivative matching.
"""
from pathlib import Path
import ast
import hashlib
import json
import numpy as np

ROOT=Path(__file__).resolve().parent
WORK=ROOT/'validation_2611_original_profile'
FILES=[('design1','2611_tapered_yue_dual_resist_L10_L12.gds',1400.,1),
       ('design2','2611_tapered_yue_design2.gds',1400.,1),
       ('combined','2611_design1_design2_L1400_18device_sweep.gds',1400.,18),
       ('long','2611_design1_L6000_P10_width_sweep.gds',6000.,3)]

def smooth_width(x,length,w0,w1):
    u=np.asarray(x)/length
    return w0+(w1-w0)*(10*u**3-15*u**4+6*u**5)

def run():
    import klayout.db as kdb
    WORK.mkdir(exist_ok=True)
    report=dict(status='RUNNING',checks=[],failures=[],outputs={},profiles=[])
    def check(name,ok,**data):
        report['checks'].append(dict(name=name,passed=bool(ok),**data))
        if not ok:report['failures'].append(name)
    coeff=np.array([0.,0.,0.,10.,-15.,6.])
    for order,expected in [(0,[0.,1.]),(1,[0.,0.]),(2,[0.,0.])]:
        values=np.polynomial.polynomial.polyval([0.,1.],np.polynomial.polynomial.polyder(coeff,order))
        check(f'original_quintic_endpoint_derivative_{order}',np.allclose(values,expected,atol=1e-12))
    check('original_quintic_monotonic',np.min(30*np.linspace(0,1,1001)**2*(1-np.linspace(0,1,1001))**2)>=0)
    for design,fn in [(1,'2611_tapered_yue_dual_resist.py'),(2,'2611_tapered_yue_design2.py')]:
        tree=ast.parse((ROOT/fn).read_text())
        fs={n.name for n in tree.body if isinstance(n,ast.FunctionDef)}
        check(f'design{design}_no_truncation_or_matched_Hermite_helpers',not fs.intersection(
            {'pnc_transition_cutoff_length','pnc_endpoint_width_derivatives','quintic_taper_coefficients','quintic_taper_widths'}))
        assignments={n.targets[0].id:ast.unparse(n.value) for n in tree.body
                     if isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name)}
        check(f'design{design}_full_45um_PNC',assignments.get('PNC_TAPER_LENGTH')=='PNC_TAPER_TRANSITION_LENGTH'
              and assignments.get('PNC_TAPER_TRANSITION_LENGTH')=='45.0')
        check(f'design{design}_5um_endpoint',assignments.get('PNC_TAPER_SUPPORT_CUTOFF_WIDTH')=='5.0')
    audit=WORK/'source_scope_audit.json'
    if audit.exists():
        for fn,d in json.loads(audit.read_text()).items():
            check(fn+'_original_drawing_functions_preserved',d['all_original_function_bodies_preserved'])
    def read(path):
        ly=kdb.Layout();ly.read(str(path));return ly,ly.top_cell()
    def region(ly,c,n):
        i=ly.find_layer(n,0)
        return kdb.Region() if i is None else kdb.Region(c.begin_shapes_rec(i)).merged()
    def insts(ly,c):
        return [i for i in c.each_inst() if not region(ly,i.cell,51).is_empty() and not region(ly,i.cell,15).is_empty()]
    previews=[]
    for tag,filename,L,number in FILES:
        path=ROOT/'gds'/filename;ly,top=read(path);dbu=ly.dbu
        check(tag+'_single_top',len(ly.top_cells())==1)
        check(tag+'_1nm_grid',abs(dbu-.001)<1e-12)
        def box(a):return kdb.Region(kdb.Box(*(round(v/dbu) for v in a)))
        devices=insts(ly,top);check(tag+'_device_count',len(devices)==number)
        origpath=WORK/'original3'/filename
        origly=orig=None
        if origpath.exists():
            origly,orig=read(origpath)
            for n in (0,13,15,16,51):
                delta=(region(ly,top,n)^region(origly,orig,n)).area()*dbu**2
                check(tag+f'_L{n}_exactly_matches_original_layout',delta==0.,xor_um2=delta)
        latestpath=WORK/'baseline'/filename
        latest=[]
        if latestpath.exists():
            latestly,latesttop=read(latestpath);latest=insts(latestly,latesttop)
        mpath=path.with_suffix('.json')
        entries=json.loads(mpath.read_text()).get('devices',[]) if mpath.exists() else []
        rows=[]
        for j,inst in enumerate(devices):
            name=f'{tag}_{j+1}';regs={n:region(ly,inst.cell,n) for n in (2,10,11,13,14,15,16,51)}
            check(name+'_all_device_layers',all(not r.is_empty() for r in regs.values()))
            entry=next((e for e in entries if e.get('case_id')==inst.cell.name),{})
            design=entry.get('design',2 if tag=='design2' or inst.cell.name.startswith('D2_') else 1)
            lengths=entry.get('stage_lengths_um',[40.,40.] if design==1 else [20.,100.,10.])
            if tag=='long':lengths=[40.,40.]
            length=lengths[0];w1=1. if design==1 else 1.65
            if latest:check(name+'_origin_unchanged_from_latest',any(i.trans==inst.trans for i in latest))
            original_inst=next((i for i in insts(origly,orig) if i.trans==inst.trans),None) if origly else None
            if origly:
                check(name+'_origin_matches_original',original_inst is not None)
                if original_inst:
                    allowed=box((L/2-.32,-240,L/2+45+length+.02,-160))
                    allowed+=allowed.transformed(kdb.Trans(kdb.Trans.R180))
                    pnc_allowed=box((L/2-.32,-240,L/2+45.32,-160))
                    pnc_allowed+=pnc_allowed.transformed(kdb.Trans(kdb.Trans.R180))
                    for n in (2,10,11,14):
                        delta=regs[n]^region(origly,original_inst.cell,n)
                        outside=(delta-(pnc_allowed if n==11 else allowed)).area()*dbu**2
                        check(name+f'_L{n}_changes_confined_to_width_transitions',outside<1e-6,outside_um2=outside)
            for end in ('right','left'):
                rr=regs if end=='right' else {n:r.transformed(kdb.Trans(kdb.Trans.R180)) for n,r in regs.items()}
                x0=L/2+45
                samples=[]
                for off in np.linspace(.01,length-.01,61):
                    ix=round((x0+off)/dbu)
                    cut=(rr[10]&kdb.Region(kdb.Box(ix,round(-203/dbu),ix+1,round(-197/dbu)))).merged()
                    got=cut.bbox().height()*dbu if cut.count()==1 else -100
                    x=(ix+.5)*dbu-x0;expected=float(smooth_width(x,length,5.,w1))
                    samples.append([x,got,expected])
                err=max(abs(a[1]-a[2]) for a in samples)
                check(name+'_'+end+'_first_stage_original_profile',err<.0041,max_width_error_um=err,sample_count=len(samples))
                check(name+'_'+end+'_monotonic_first_stage',all(b[1]<=a[1]+.0021 for a,b in zip(samples,samples[1:])))
                ix=round((x0-.001)/dbu)
                r=(rr[14]&kdb.Region(kdb.Box(ix,round(-230/dbu),ix+1,round(-170/dbu)))).merged()
                end_full=r.bbox().height()*dbu
                errs=[]
                for x in np.linspace(.02,44.98,61):
                    ix=round((L/2+x)/dbu)
                    r=(rr[14]&kdb.Region(kdb.Box(ix,round(-230/dbu),ix+1,round(-170/dbu)))).merged()
                    observed=r.bbox().height()*dbu-end_full+5 if r.count()==1 else -100
                    expected=float(smooth_width((ix+.5)*dbu-L/2,45.,15.,5.))
                    errs.append(abs(observed-expected))
                check(name+'_'+end+'_PNC_envelope_full_15to5_original_profile',max(errs)<.0061,max_width_error_um=max(errs))
                pitch=10. if tag=='long' else 5.
                expected=np.r_[np.arange(0,40-1e-9,pitch),np.arange(40.,46.)]
                roi=box((L/2-.4,-230,L/2+45.4,-170));stations={1:[],-1:[]}
                for p in (rr[11]&roi).merged().each():
                    b=p.bbox();side=1 if (b.bottom+b.top)*dbu/2>-200 else -1
                    stations[side].append((b.left+b.right)*dbu/2-L/2)
                for side,got in stations.items():
                    got=sorted(got)
                    check(name+f'_{end}_PNC_stations_side{side}',len(got)==len(expected) and np.allclose(got,expected,atol=.0011,rtol=0),measured_um=got)
                report['profiles'].append(dict(device=name,end=end,design=design,pnc_length_um=45.,first_stage_length_um=length,widths_um=[15.,5.,w1],max_first_stage_error_um=err))
            rows.append(inst.trans.disp.y*dbu)
            if tag in ('design1','design2'):previews.append((tag,regs,L/2,length,dbu))
        if number>1:
            check(tag+'_600um_y_pitch',np.allclose(np.diff(sorted(rows)),600,atol=1e-9,rtol=0))
            check(tag+'_single_label_per_device',sum(i.cell.name.startswith('LABEL_') for i in top.each_inst())==number)
        report['outputs'][filename]=dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest(),devices=len(devices))
    report.update(status='FAIL' if report['failures'] else 'PASS',check_count=len(report['checks']),
                  geometry_profile='original quintic smootherstep',pnc_length_um=45.,pnc_widths_um=[15.,5.],
                  junction_width_slope=0.,junction_width_curvature_per_um=0.,dense_station_positions_um=list(range(40,46)),
                  original_reference_commit='6d4fcb56e649e80e0270c3af37be16a102799c00',
                  optical_adiabaticity_simulated=False)
    (WORK/'validation_report.json').write_text(json.dumps(report,indent=2))
    (ROOT/'gds/2611_original_profile_validation.json').write_text(json.dumps(report,indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for tag,regs,xport,length,dbu in previews:
        fig,ax=plt.subplots(figsize=(12,4.5))
        roi=box((xport-5,-230,xport+45+length+5,-170))
        def curves(r):
            X=[];Y=[]
            for p in r.each():
                pts=list(p.each_point_hull());pts.append(pts[0])
                X.extend([q.x*dbu-xport for q in pts]+[np.nan]);Y.extend([q.y*dbu+200 for q in pts]+[np.nan])
            return X,Y
        for n in (14,11,10,15,16):
            ax.plot(*curves((regs[n]&roi).merged()),label=f'Layer {n}',linewidth=.7 if n!=10 else 1.2)
        x=np.linspace(0,45,201);w=smooth_width(x,45,15,5)
        ax.plot(x,w/2,'--',linewidth=1.,label='15 to 5 um envelope')
        ax.plot(x,-w/2,'--',linewidth=1.)
        ax.axvline(45,linestyle=':',linewidth=.7)
        ax.annotate('5 um junction\nzero slope and curvature',(45,2.5),xytext=(50,12),arrowprops={'arrowstyle':'->'})
        ax.set(xlim=(-5,45+length+5),ylim=(-15,15),xlabel='Distance from PNC entrance (um)',ylabel='Transverse position (um)',title=tag.title()+': original full-length taper profile (expanded transverse scale)')
        ax.legend(ncol=3,fontsize=8,loc='lower left');fig.tight_layout();fig.savefig(WORK/(tag+'_original_profile.png'),dpi=180);plt.close(fig)
    print(report['status'],report['check_count'],'checks',report['failures'],flush=True)
    if report['failures']:raise SystemExit(1)

if __name__=='__main__':run()
