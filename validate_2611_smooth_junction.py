"""Check C2 design equations and re-read all four production GDSs; no repairs.

Requires baseline copies in validation_2611_smooth_junction/baseline/ for
preservation checks. Geometry verification is not optical simulation.
"""
from pathlib import Path
import ast
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parent
WORK = ROOT / 'validation_2611_smooth_junction'
FILES = [
    ('design1', '2611_tapered_yue_dual_resist_L10_L12.gds', 1400., 1),
    ('design2', '2611_tapered_yue_design2.gds', 1400., 1),
    ('combined', '2611_design1_design2_L1400_18device_sweep.gds', 1400., 18),
    ('long', '2611_design1_L6000_P10_width_sweep.gds', 6000., 3),
]


def reference_endpoint():
    lo, hi = 0., 45.
    for _ in range(80):
        x = (lo + hi)/2; t=x/45
        if 15-12*(6*t**5-15*t**4+10*t**3)>5: lo=x
        else: hi=x
    x=(lo+hi)/2;t=x/45
    return x, -12/45*(30*t**4-60*t**3+30*t**2), -12/45**2*(120*t**3-180*t**2+60*t)


def independent_coefficients(length, target, slope, curvature):
    # Solve all six endpoint conditions, independently of production formulas.
    A=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,2,0,0,0],
                [1,1,1,1,1,1],[0,1,2,3,4,5],[0,0,2,6,12,20]],float)
    return np.linalg.solve(A, [5,slope*length,curvature*length**2,target,0,0])


def run():
    import klayout.db as kdb
    endpoint,slope,curvature=reference_endpoint()
    report=dict(status='RUNNING',checks=[],failures=[],outputs={},analytic_profiles=[])
    def check(name,ok,**data):
        report['checks'].append(dict(name=name,passed=bool(ok),**data))
        if not ok: report['failures'].append(name)
    helper_names={'quintic_taper_coefficients','quintic_taper_widths','pnc_endpoint_width_derivatives'}
    for design,fn in [(1,'2611_tapered_yue_dual_resist.py'),(2,'2611_tapered_yue_design2.py')]:
        tree=ast.parse((ROOT/fn).read_text())
        module=ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in helper_names],type_ignores=[])
        ns=dict(np=np,PNC_TAPER_LENGTH=endpoint,PNC_TAPER_TRANSITION_LENGTH=45.,
                PNC_TAPER_PROFILE_END_WIDTH=3.,SLAB_EXTENSION_END_WIDTH=15.)
        exec(compile(module,fn,'exec'),ns)
        got=ns['pnc_endpoint_width_derivatives']()
        check(f'design{design}_incoming_derivatives',np.allclose(got,[slope,curvature],rtol=0,atol=1e-13))
        lengths=[40.,50.,60.] if design==1 else [20.,80*2/13,100*2/13,120*2/13]
        target=1. if design==1 else 1.65
        for length in lengths:
            name=f'design{design}_L{length:.9g}'
            coeff=ns['quintic_taper_coefficients'](length,5.,target,*got)
            expected=independent_coefficients(length,target,slope,curvature)
            check(name+'_independent_six_constraint_solution',np.allclose(coeff,expected,rtol=0,atol=1e-10))
            bounds=[]
            for order,values in [(0,[5.,target]),(1,[slope,0.]),(2,[curvature,0.])]:
                measured=np.polynomial.polynomial.polyval([0.,1.],np.polynomial.polynomial.polyder(coeff,order))/length**order
                check(name+f'_C{order}_entrance_and_exit',np.allclose(measured,values,rtol=0,atol=1e-10))
                bounds.append(measured.tolist())
            d1=np.polynomial.polynomial.polyder(coeff)
            roots=np.polynomial.polynomial.polyroots(np.polynomial.polynomial.polyder(d1))
            coords=[0.,1.]+[r.real for r in roots if abs(r.imag)<1e-8 and 0<r.real<1]
            rates=np.polynomial.polynomial.polyval(coords,d1)/length
            check(name+'_monotonic_without_undershoot',float(max(rates))<=1e-10)
            report['analytic_profiles'].append(dict(design=design,length_um=length,
                input_width_um=5.,output_width_um=target,start_slope=slope,start_curvature_per_um=curvature,
                max_abs_width_slope=float(max(abs(rates))),endpoint_values=bounds))
        try: ns['quintic_taper_coefficients'](40.,5.,5.,slope,curvature)
        except ValueError: rejected=True
        else: rejected=False
        check(f'design{design}_invalid_profile_rejected',rejected)
    def read(path):
        ly=kdb.Layout();ly.read(str(path));return ly,ly.top_cell()
    def region(ly,cell,n):
        i=ly.find_layer(n,0)
        return kdb.Region() if i is None else kdb.Region(cell.begin_shapes_rec(i)).merged()
    def instances(ly,top):
        return [i for i in top.each_inst() if not region(ly,i.cell,51).is_empty() and not region(ly,i.cell,15).is_empty()]
    preview_data=[]
    for tag,filename,L,count in FILES:
        path=ROOT/'gds'/filename;old_path=WORK/'baseline'/filename
        ly,top=read(path);dbu=ly.dbu
        check(tag+'_one_top_cell',len(ly.top_cells())==1)
        check(tag+'_database_grid',abs(dbu-.001)<1e-12)
        check(tag+'_baseline_present',old_path.exists())
        if not old_path.exists():continue
        old_ly,old_top=read(old_path);ins=instances(ly,top);oldins=instances(old_ly,old_top)
        check(tag+'_device_count',len(ins)==count and len(oldins)==count)
        def box(a): return kdb.Region(kdb.Box(*(round(v/dbu) for v in a)))
        for n in (0,11,13,15,16,51):
            delta=(region(ly,top,n)^region(old_ly,old_top,n)).area()*dbu**2
            check(tag+f'_entire_L{n}_unchanged',delta==0.,xor_um2=delta)
        manifest_path=path.with_suffix('.json')
        entries=json.loads(manifest_path.read_text()).get('devices',[]) if manifest_path.exists() else []
        rows=[]
        for index,inst in enumerate(ins):
            old=next((i for i in oldins if i.trans==inst.trans),None)
            name=f'{tag}_{index+1}'
            check(name+'_placement_unchanged',old is not None)
            if old is None:continue
            regs={n:region(ly,inst.cell,n) for n in (2,10,11,13,14,15,16,51)}
            previous={n:region(old_ly,old.cell,n) for n in regs}
            entry=next((e for e in entries if e.get('case_id')==inst.cell.name),{})
            design=entry.get('design',2 if tag=='design2' or inst.cell.name.startswith('D2_') else 1)
            stages=entry.get('stage_lengths_um',[40.,40.] if design==1 else [20.,100.,10.])
            length=stages[0];target=1. if design==1 else 1.65
            coeff=independent_coefficients(length,target,slope,curvature)
            local_start=L/2+round(endpoint/dbu)*dbu
            allowed=box((local_start-.01,-232,local_start+length+.01,-168))
            allowed+=allowed.transformed(kdb.Trans(kdb.Trans.R180))
            for n in (2,10,14):
                diff=regs[n]^previous[n]
                outside=(diff-allowed).area()*dbu**2
                check(name+f'_L{n}_changes_only_first_stage',outside<1e-6,outside_um2=outside,
                      total_changed_um2=diff.area()*dbu**2)
            for end in ('right','left'):
                rr=regs if end=='right' else {n:r.transformed(kdb.Trans(kdb.Trans.R180)) for n,r in regs.items()}
                offsets=np.linspace(.01,length-.01,41)
                samples=[]
                for offset in offsets:
                    ix=round((local_start+offset)/dbu)
                    cut=(rr[10]&kdb.Region(kdb.Box(ix,round(-203/dbu),ix+1,round(-197/dbu)))).merged()
                    measured=cut.bbox().height()*dbu if cut.count()==1 else -100.
                    at=(ix+.5)*dbu-local_start
                    expected=float(np.polynomial.polynomial.polyval(at/length,coeff))
                    samples.append([float(at),measured,expected])
                errors=[abs(a[1]-a[2]) for a in samples]
                check(name+'_'+end+'_L10_matches_C2_profile',max(errors)<.0041,
                      max_full_width_error_um=float(max(errors)),sample_count=len(samples))
                check(name+'_'+end+'_no_width_step_or_reversal',all(b[1]<=a[1]+.0021 for a,b in zip(samples,samples[1:])))
                vals=[]
                for off in (.01,1.01):
                    ix=round((local_start+off)/dbu)
                    cut=(rr[10]&kdb.Region(kdb.Box(ix,round(-203/dbu),ix+1,round(-197/dbu)))).merged()
                    vals.append(((ix+.5)*dbu-local_start,cut.bbox().height()*dbu))
                observed=(vals[1][1]-vals[0][1])/(vals[1][0]-vals[0][0])
                predicted=np.diff(np.polynomial.polynomial.polyval(np.array([v[0] for v in vals])/length,coeff))[0]/(vals[1][0]-vals[0][0])
                check(name+'_'+end+'_exported_entrance_slope',abs(observed-predicted)<.0041,
                      measured_first_micron_slope=float(observed),expected=float(predicted))
                roi=box((local_start-.5,-215,local_start+length+.5,-185))
                missing=((rr[10]&roi)-rr[14]).area()*dbu**2
                check(name+'_'+end+'_L14_covers_modified_slab',missing<.005,missing_um2=missing)
            parts=sorted(regs[15].each(),key=lambda p:p.area(),reverse=True)
            optical=kdb.Region(parts[0])
            check(name+'_one_optical_body_and_four_ticks',len(parts)==5)
            check(name+'_connected_protection',regs[16].count()==1)
            for n in (14,16):
                check(name+f'_L{n}_covers_optical',(optical-regs[n]).area()*dbu**2<.005)
            for n in (10,11,13,14,15,16,51):
                roi=box((L/2+.001,-1000,L/2+11000,-25))
                delta=((regs[n]^regs[n].transformed(kdb.Trans(kdb.Trans.R180)))&roi).area()*dbu**2
                check(name+f'_L{n}_both_ends_symmetric',delta<.005,xor_um2=delta)
            rows.append(inst.trans.disp.y*dbu)
            if tag in ('design1','design2'):
                preview_data.append((tag,regs,previous,local_start,length,dbu))
        rows.sort()
        if len(rows)>1:
            check(tag+'_600um_pitch',np.allclose(np.diff(rows),600,rtol=0,atol=1e-9))
            check(tag+'_one_label_per_device',sum(i.cell.name.startswith('LABEL_') for i in top.each_inst())==count)
        report['outputs'][filename]=dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest(),devices=len(ins))
    report.update(status='FAIL' if report['failures'] else 'PASS',check_count=len(report['checks']),
                  cutoff_width_um=5.,pnc_endpoint_um=endpoint,
                  old_width_slope_jump=abs(slope),old_width_curvature_jump_per_um=abs(curvature),
                  new_analytic_continuity='C2',optical_adiabaticity_simulated=False,
                  note='C2 applies to the analytical envelope; exported polygons are sampled on a 1 nm grid.')
    WORK.mkdir(exist_ok=True)
    (WORK/'validation_report.json').write_text(json.dumps(report,indent=2))
    (ROOT/'gds/2611_smooth_junction_validation.json').write_text(json.dumps(report,indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for tag,regs,previous,x0,length,dbu in preview_data:
        fig,ax=plt.subplots(figsize=(12,4.5))
        def curves(reg):
            X=[];Y=[]
            for poly in reg.each():
                pts=list(poly.each_point_hull())
                if not pts:continue
                pts.append(pts[0]);X.extend([p.x*dbu-x0 for p in pts]+[np.nan]);Y.extend([p.y*dbu+200 for p in pts]+[np.nan])
            return X,Y
        for n in (14,11,10,15,16):
            ax.plot(*curves(regs[n]),label=f'Layer {n}',linewidth=1.3 if n==10 else .7)
        ax.plot(*curves(previous[10]),linestyle='--',linewidth=1,label='Previous Layer 10')
        ax.axvline(0,linestyle=':',linewidth=.7)
        ax.annotate('5 um junction\nmatched slope + curvature',(0,2.5),xytext=(3,10),arrowprops={'arrowstyle':'->'})
        ax.set(xlim=(-8,length+6),ylim=(-13,13),xlabel='Distance from 5 um junction (um)',ylabel='Transverse position (um)',
               title=f'{tag.title()}: exported C2-matched first taper (vertical scale expanded)')
        ax.legend(ncol=3,fontsize=8,loc='lower left');fig.tight_layout()
        fig.savefig(WORK/(tag+'_smooth_junction.png'),dpi=180);plt.close(fig)
    print(report['status'],report['check_count'],'checks')
    if report['failures']:
        print('\n'.join(report['failures']));raise SystemExit(1)


if __name__=='__main__':run()
