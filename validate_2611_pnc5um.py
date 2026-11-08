"""Read and validate the actual saved 5 um-cutoff GDS files; never repair them."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent
WORK = ROOT / 'validation_2611_pnc5um'
FILES = [
    ('design1', '2611_tapered_yue_dual_resist_L10_L12.gds', 1400., 5., 1),
    ('design2', '2611_tapered_yue_design2.gds', 1400., 5., 1),
    ('combined', '2611_design1_design2_L1400_18device_sweep.gds', 1400., 5., 18),
    ('long', '2611_design1_L6000_P10_width_sweep.gds', 6000., 10., 3),
]


def cutoff():
    # Full original profile now ends at 5 um; no interior-width root.
    return 45.0


def run():
    import klayout.db as kdb
    import numpy as np
    report = dict(status='RUNNING', checks=[], failures=[], outputs={}, baseline_files=[])
    endpoint = cutoff()
    def check(name, ok, **data):
        report['checks'].append(dict(name=name, passed=bool(ok), **data))
        if not ok: report['failures'].append(name)
    def read(path):
        ly = kdb.Layout(); ly.read(str(path))
        if len(ly.top_cells()) != 1: raise ValueError('Expected one top cell: '+str(path))
        return ly, ly.top_cell()
    def region(ly, cell, n):
        idx = ly.find_layer(n, 0)
        return kdb.Region() if idx is None else kdb.Region(cell.begin_shapes_rec(idx)).merged()
    def device_instances(ly, top):
        return [i for i in top.each_inst() if not region(ly,i.cell,51).is_empty()
                and not region(ly,i.cell,15).is_empty()]
    for tag,filename,L,pitch,number in FILES:
        path = ROOT/'gds'/filename
        if not path.exists():
            check(tag+'_output_exists', False); continue
        ly,top=read(path);dbu=ly.dbu
        def box(coords):
            return kdb.Region(kdb.Box(*(round(v/dbu) for v in coords)))
        insts=device_instances(ly,top)
        check(tag+'_device_count',len(insts)==number,actual=len(insts))
        check(tag+'_grid',abs(dbu-.001)<1e-12,dbu_um=dbu)
        old_path=WORK/'baseline'/filename
        before=[];old_ly=None
        if old_path.exists():
            old_ly,old_top=read(old_path);before=device_instances(old_ly,old_top)
            report['baseline_files'].append(filename)
        rows=[];lefts=[]
        for j,inst in enumerate(insts):
            name=tag+'_'+str(j+1)
            regs={n:region(ly,inst.cell,n) for n in (2,10,11,13,14,15,16,51)}
            check(name+'_all_layers',all(not r.is_empty() for r in regs.values()))
            parts=sorted(regs[15].each(),key=lambda p:p.area(),reverse=True)
            check(name+'_one_optical_body_four_ticks',len(parts)==5)
            if not parts: continue
            optical=kdb.Region(parts[0]); ticks=kdb.Region()
            for p in parts[1:]: ticks.insert(p)
            check(name+'_four_complete_ticks',len(parts)==5 and
                  all(abs(p.area()*dbu**2-5)<.002 for p in parts[1:]))
            check(name+'_ticks_not_on_L10',(ticks & regs[10]).is_empty())
            check(name+'_L16_connected',regs[16].count()==1)
            for n in (14,16):
                missing=(optical-regs[n]).area()*dbu**2
                check(name+f'_L{n}_covers_optical',missing<.005,missing_um2=missing)
            meta_path=path.with_suffix('.json');entries=[]
            if meta_path.exists():entries=json.loads(meta_path.read_text()).get('devices',[])
            entry=next((e for e in entries if e.get('case_id')==inst.cell.name),{})
            design=entry.get('design',2 if tag=='design2' or inst.cell.name.startswith('D2_') else 1)
            stages=entry.get('stage_lengths_um',[40.,40.] if design==1 else [20.,100.,10.])
            if tag=='long':stages=[40.,40.]
            first_stage=stages[0]
            expected=np.r_[np.arange(-3*pitch,-1e-6,pitch),
                           np.arange(0,endpoint-5-1e-9,pitch),
                           endpoint-np.arange(5,-1,-1)]
            for end in ('right','left'):
                rr={n:r if end=='right' else r.transformed(kdb.Trans(kdb.Trans.R180)) for n,r in regs.items()}
                roi=box((L/2-3*pitch-.4,-230,L/2+45.4,-170))
                centers={1:[],-1:[]}
                for p in (rr[11]&roi).merged().each():
                    b=p.bbox();side=1 if (b.bottom+b.top)*dbu/2>-200 else -1
                    centers[side].append((b.left+b.right)*dbu/2-L/2)
                for side,values in centers.items():
                    got=np.asarray(sorted(values))
                    check(name+f'_{end}_side{side}_PNC_stations',
                          len(got)==len(expected) and np.allclose(got,expected,atol=.0011,rtol=0),
                          measured_um=got.tolist(),expected_um=expected.tolist())
                    check(name+f'_{end}_side{side}_none_beyond_5um_cutoff',
                          bool(len(got)) and max(got)<=endpoint+.0011)
                ix=round((L/2+endpoint+.002)/dbu)
                slab=(rr[10] & kdb.Region(kdb.Box(ix,round(-202.6/dbu),ix+1,round(-197.4/dbu)))).merged()
                w=slab.bbox().height()*dbu if not slab.is_empty() else 0.
                check(name+f'_{end}_L10_starts_at_5um',slab.count()==1 and abs(w-5)<.0041,measured_width_um=w)
                tr=ticks if end=='right' else ticks.transformed(kdb.Trans(kdb.Trans.R180))
                near=tr & box((L/2+endpoint+first_stage-1,-218,L/2+endpoint+first_stage+1,-182))
                check(name+f'_{end}_ticks_follow_new_junction',near.merged().count()==2)
                junction=box((L/2+endpoint-.5,-215,L/2+endpoint+first_stage+.5,-185))
                missed=(rr[10]&junction)-rr[14]
                check(name+f'_{end}_L14_covers_L10_transition',missed.area()*dbu**2<.005,missing_um2=missed.area()*dbu**2)
            for n in (10,11,13,14,15,16,51):
                roi=box((L/2+.001,-1000,L/2+11000,-25))
                area=((regs[n]^regs[n].transformed(kdb.Trans(kdb.Trans.R180)))&roi).area()*dbu**2
                check(name+f'_L{n}_end_symmetry',area<.005,xor_um2=area)
            old=next((i for i in before if i.trans==inst.trans),None)
            if old_ly is not None:
                check(name+'_baseline_placement_matches',old is not None)
            if old is not None:
                core=box((-L/2+.5,-240,L/2-.5,240))
                core+=box((L/2-.3,-25,L/2+220,225))
                core+=box((-L/2-220,-225,-L/2+.3,25))
                for n in (2,10,11,14,15,16):
                    old_region=region(old_ly,old.cell,n)
                    change=(regs[n]^old_region)&core
                    raw_delta=change.area()*dbu**2
                    delta=raw_delta
                    if tag=='design2' and n==2:
                        from validate_2611_pnc5um_label_shift import exact_standalone_text
                        old_text=exact_standalone_text(L,45.,first_stage,dbu)
                        new_text=exact_standalone_text(L,endpoint,first_stage,dbu)
                        # Compare exact glyph regions, not a rectangular exclusion.
                        # Restoring BOTH known labels must recover identical core masks.
                        delta=(((regs[n]|old_text)^(old_region|new_text))&core).area()*dbu**2
                        outside=(change-(old_text|new_text)).area()*dbu**2
                        check(name+'_L2_difference_is_exact_label_translation',outside<.005,
                              outside_glyph_xor_um2=outside,raw_xor_um2=raw_delta)
                    check(name+f'_L{n}_meander_preserved',delta<.005,
                          xor_um2=delta,raw_xor_um2=raw_delta)
            rows.append(inst.trans.disp.y*dbu)
            lefts.append((regs[11].bbox().left+inst.trans.disp.x)*dbu)
        check(tag+'_left_bends_aligned',all(abs(x+900)<.0011 for x in lefts))
        if number>1:
            ys=sorted(rows)
            check(tag+'_600um_pitch',all(abs(b-a-600)<.0011 for a,b in zip(ys,ys[1:])))
            check(tag+'_one_label_per_device',sum(i.cell.name.startswith('LABEL_') for i in top.each_inst())==number)
        report['outputs'][filename]=dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest(),devices=len(insts))
    report.update(status='FAIL' if report['failures'] else 'PASS',
                  cutoff_width_um=5.,active_pnc_length_um=endpoint,
                  dense_start_um=endpoint-5,dense_end_um=endpoint,
                  check_count=len(report['checks']))
    WORK.mkdir(exist_ok=True)
    (WORK/'validation_report.json').write_text(json.dumps(report,indent=2))
    (ROOT/'gds/2611_pnc5um_validation.json').write_text(json.dumps(report,indent=2))
    print(report['status'],report['check_count'],'checks')
    if report['failures']:
        print('\n'.join(report['failures']));raise SystemExit(1)

if __name__=='__main__':run()
