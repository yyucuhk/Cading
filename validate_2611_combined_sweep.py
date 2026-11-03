"""Independent checks on the combined, re-read GDS (does not regenerate it).
Run after 2611_build_combined_sweep.py. Geometry checks are not foundry DRC.
"""
from pathlib import Path
import json
import math
import gdsfactory as gf
import klayout.db as kdb

ROOT=Path(__file__).resolve().parent
STEM=ROOT/'gds/2611_design1_design2_L1400_18device_sweep'
OUT=ROOT/'gds/2611_combined_sweep_validation.json'
REPORT={'checks':[], 'failures':[]}


def check(name, ok, **details):
    REPORT['checks'].append(dict(name=name, passed=bool(ok), **details))
    if not ok: REPORT['failures'].append(name)


def run():
    gf.gpdk.PDK.activate()
    manifest=json.loads(STEM.with_suffix('.json').read_text())
    ly=kdb.Layout();ly.read(str(STEM.with_suffix('.gds')));dbu=ly.dbu
    top=ly.top_cell()
    instances={i.cell.name:i for i in top.each_inst() if i.cell.name.startswith(('D1_','D2_'))}
    check('18_complete_device_instances',len(instances)==18)
    check('manifest_matches_actual_cells',set(instances)=={d['case_id'] for d in manifest['devices']})
    def region(cell,n):
        idx=ly.find_layer(n,0)
        return kdb.Region() if idx is None else kdb.Region(cell.begin_shapes_rec(idx)).merged()
    def box(coords):
        return kdb.Region(kdb.Box(*(round(v/dbu) for v in coords)))
    rows={1:[],2:[]}
    for e in manifest['devices']:
        name=e['case_id'];inst=instances[name];cell=inst.cell
        regs={n:region(cell,n) for n in (2,10,11,13,14,15,16,51)}
        check(name+'_all_layers_present',all(not r.is_empty() for r in regs.values()))
        actual=(inst.trans.disp.x*dbu,inst.trans.disp.y*dbu)
        check(name+'_placement',all(abs(a-b)<1e-8 for a,b in zip(actual,e['origin_um'])))
        left=(regs[11].bbox().left+inst.trans.disp.x)*dbu
        check(name+'_left_bend_aligned',abs(left-(-900+(e['design']-1)*5000))<1e-8,x_um=left)
        rows[e['design']].append((e['row'],*actual))
        ticks=kdb.Region()
        for (tx,ty),ori in zip(e['tip_centers'],e['orientations']):
            a=math.radians(ori+180);t=(math.cos(a),math.sin(a));n=(-t[1],t[0])
            for sign in (-1,1):
                cx=tx+sign*7.5*n[0];cy=ty+sign*7.5*n[1]
                ticks.insert(kdb.Polygon([kdb.Point(round((cx+dx*t[0]+dy*n[0])/dbu),
                                                        round((cy+dx*t[1]+dy*n[1])/dbu))
                    for dx,dy in ((-.25,-5),(.25,-5),(.25,5),(-.25,5))]))
        optical=(regs[15]-ticks).merged()
        check(name+'_connected_optical_path',optical.count()==1)
        check(name+'_four_tips_on_L15_only',ticks.count()==4 and
              (ticks-regs[15]).is_empty() and (ticks&regs[10]).is_empty())
        check(name+'_connected_L16',regs[16].count()==1)
        for mask in (14,16):
            area=(optical-regs[mask]).area()*dbu**2
            check(name+f'_L{mask}_covers_ridge',area<.005,uncovered_um2=area)

        # Full end-assembly symmetry: exclude only the exact length-text shapes.
        left_tip=min(e['tip_centers'],key=lambda p:p[0])
        label=(f"ML{e['meander_length_um']:.1f}_TL{e['taper_total_um']:g}" if e['design']==1 else
               f"DESIGN 2_L{e['straight_um']:g}_ML{e['meander_length_um']:.1f}")
        tc=gf.components.text(text=label,size=3.,layer=(2,0))
        text=kdb.Region(tc.kdb_cell.begin_shapes_rec(gf.get_layer((2,0)))).transformed(
            kdb.Trans(round((left_tip[0]+5.)/dbu),round((left_tip[1]+10.)/dbu)))
        text+=text.transformed(kdb.Trans(kdb.Trans.R180))
        roi=box((e['straight_um']/2+.001,-1000,10000,-25))
        for n,r in regs.items():
            delta=(r^r.transformed(kdb.Trans(kdb.Trans.R180)))&roi
            if n==2: delta-=text
            area=delta.area()*dbu**2
            check(name+f'_L{n}_both_complete_ends_symmetric',area<.005,xor_um2=area)

        # Measure all junction widths in BOTH complete device ends in the GDS.
        # The transverse +/-2 um clip avoids the separate alignment ticks.
        for side in ('right','left'):
            for n,lengths,widths in ((10,e['stage_lengths_um'],e['layer10_widths_um']),
                                     (15,e['layer15_lengths_um'],e['layer15_widths_um'])):
                r=regs[n] if side=='right' else regs[n].transformed(kdb.Trans(kdb.Trans.R180))
                x=e['straight_um']/2+45.
                coords=[x]
                for length in lengths: coords.append(coords[-1]+length)
                for j,(x,width) in enumerate(zip(coords,widths)):
                    ix=round((x+(.002 if j==0 else -.002))/dbu)
                    cr=(r & kdb.Region(kdb.Box(ix,round(-202/dbu),ix+1,round(-198/dbu)))).merged()
                    measured=cr.bbox().height()*dbu if cr.count()==1 else None
                    check(name+f'_{side}_L{n}_junction{j}',
                          measured is not None and abs(measured-width)<.0041,
                          expected_width_um=width,measured_width_um=measured)
    for design, positions in rows.items():
        positions.sort()
        check(f'design{design}_nine_rows',len(positions)==9)
        check(f'design{design}_600um_y_pitch',all(abs(b[2]-a[2]-600)<1e-8
              for a,b in zip(positions,positions[1:])))
    check('matching_rows_5000um_apart',all(abs(b[1]-a[1]-5000)<1e-8 and b[2]==a[2]
           for a,b in zip(rows[1],rows[2])))
    REPORT['status']='PASS' if not REPORT['failures'] else 'FAIL'
    REPORT['check_count']=len(REPORT['checks'])
    REPORT['gds_file']=STEM.with_suffix('.gds').name
    OUT.write_text(json.dumps(REPORT,indent=2))
    print(REPORT['status'],REPORT['check_count'],'checks',REPORT['failures'],flush=True)
    if REPORT['failures']: raise SystemExit(1)


if __name__=='__main__':
    run()
