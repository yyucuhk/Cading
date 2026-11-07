"""Measure straight-to-taper PNC continuity in all four exported GDSs.

Run 'python validate_2611_taper_pitch.py capture' BEFORE editing to retain a
baseline, then run without arguments AFTER regenerating the four layouts.
Only Layer 11 inside the first 5 um of each existing PNC taper may change.
Five 1 um intervals are checked, along with the unchanged surrounding grid. This does not
replace fabrication DRC or optical simulation.
"""
from pathlib import Path
import hashlib
import json
import shutil
import sys

import klayout.db as kdb
import numpy as np

ROOT = Path(__file__).resolve().parent
WORK = ROOT / 'validation_2611_taper_pitch'
DENSE_PITCH = 1.0
DENSE_PERIODS = 5
DENSE_SPAN = DENSE_PITCH * DENSE_PERIODS
BASE = WORK / 'baseline'
FILES = [
    ('design1', '2611_tapered_yue_dual_resist_L10_L12.gds', 1400.0, 5.0, 1),
    ('design2', '2611_tapered_yue_design2.gds', 1400.0, 5.0, 1),
    ('combined', '2611_design1_design2_L1400_18device_sweep.gds', 1400.0, 5.0, 18),
    ('long_design1', '2611_design1_L6000_P10_width_sweep.gds', 6000.0, 10.0, 3),
]
REPORT = {'checks': [], 'failures': [], 'measurements': [], 'outputs': {}}


def check(name, ok, **details):
    REPORT['checks'].append(dict(name=name, passed=bool(ok), **details))
    if not ok:
        REPORT['failures'].append(name)
    print(name, bool(ok), flush=True)


def read(path):
    layout = kdb.Layout()
    layout.read(str(path))
    if len(layout.top_cells()) != 1:
        raise ValueError(f'{path}: expected one top cell')
    return layout, layout.top_cell()


def region(layout, cell, layer):
    idx = layout.find_layer(*layer)
    return kdb.Region() if idx is None else kdb.Region(cell.begin_shapes_rec(idx)).merged()


def all_regions(layout, cell):
    result = {}
    for idx in layout.layer_indexes():
        info = layout.get_info(idx)
        reg = kdb.Region(cell.begin_shapes_rec(idx)).merged()
        if not reg.is_empty():
            result[(info.layer, info.datatype)] = reg
    return result


def box(bounds, dbu):
    return kdb.Region(kdb.Box(*(round(v / dbu) for v in bounds)))


def devices(layout, top):
    return [i for i in top.each_inst()
            if not region(layout, i.cell, (51, 0)).is_empty()
            and not region(layout, i.cell, (15, 0)).is_empty()]


def pnc_centers(reg, L, pitch, dbu):
    # Includes three unchanged straight stations and the complete taper grid.
    roi = box((L/2-3*pitch-.4, -230, L/2+45.4, -170), dbu)
    sides = {1: [], -1: []}
    for poly in (reg & roi).merged().each():
        b = poly.bbox()
        side = 1 if (b.bottom+b.top)*dbu/2 > -200 else -1
        sides[side].append((b.left+b.right)*dbu/2 - L/2)
    return {k: np.asarray(sorted(v)) for k, v in sides.items()}


def capture():
    BASE.mkdir(parents=True, exist_ok=True)
    for _, filename, _, _, _ in FILES:
        shutil.copy2(ROOT/'gds'/filename, BASE/filename)
    print('Captured four approved GDS baselines.', flush=True)


def run():
    WORK.mkdir(exist_ok=True)
    count = 0
    for tag, filename, L, pitch, expected_count in FILES:
        layout, top = read(ROOT/'gds'/filename)
        dbu = layout.dbu
        insts = devices(layout, top)
        check(tag+'_device_count', len(insts) == expected_count, count=len(insts))
        check(tag+'_database_grid', abs(dbu-.001) < 1e-12)
        regular = np.arange(0, 45+1e-9, pitch)
        dense = DENSE_PITCH * np.arange(DENSE_PERIODS + 1)
        expected = np.r_[np.arange(-3*pitch, -.0001, pitch), dense,
                         regular[regular > DENSE_SPAN + 1e-9]]
        allowed_global = kdb.Region()
        local_allowed = box((L/2-.4, -230, L/2+DENSE_SPAN+.4, -170), dbu)
        local_allowed += local_allowed.transformed(kdb.Trans(kdb.Trans.R180))
        old_layout = old_top = None
        if (BASE/filename).exists():
            old_layout, old_top = read(BASE/filename)
        before_by_origin = {} if old_layout is None else {
            (i.trans.disp.x, i.trans.disp.y): i for i in devices(old_layout, old_top)}
        for di, inst in enumerate(insts):
            count += 1
            ident = f'{tag}_{di+1}'
            reg = region(layout, inst.cell, (11,0))
            allowed_global += local_allowed.transformed(inst.trans)
            origin = (inst.trans.disp.x, inst.trans.disp.y)
            old_inst = before_by_origin.get(origin)
            if old_layout is not None:
                check(ident+'_same_placement', old_inst is not None)
            for end in ('right', 'left'):
                rr = reg if end == 'right' else reg.transformed(kdb.Trans(kdb.Trans.R180))
                values = pnc_centers(rr, L, pitch, dbu)
                for side, centers in values.items():
                    spacing = np.diff(centers)
                    ok = len(centers) == len(expected) and np.allclose(centers, expected, atol=.0011, rtol=0)
                    check(f'{ident}_{end}_{side}_actual_PNC_grid', ok,
                          measured_centers_from_port_um=centers.tolist(),
                          expected_pitch_um=pitch)
                    check(f'{ident}_{end}_{side}_correct_piecewise_intervals',
                          len(spacing) == len(expected)-1 and np.allclose(spacing, np.diff(expected), atol=.0011, rtol=0),
                          measured_min_um=float(min(spacing)) if len(spacing) else None,
                          measured_max_um=float(max(spacing)) if len(spacing) else None)
                    dense_centers = centers[(centers >= -.0011) & (centers <= DENSE_SPAN+.0011)]
                    check(f'{ident}_{end}_{side}_exactly_five_1um_periods',
                          len(dense_centers) == 6 and np.allclose(dense_centers, dense, atol=.0011, rtol=0),
                          stations_um=dense_centers.tolist())
                    check(f'{ident}_{end}_{side}_one_PNC_at_interface',
                          sum(abs(centers) < .0011) == 1)
                previous = None
                if old_inst is not None:
                    old_r = region(old_layout, old_inst.cell, (11,0))
                    if end == 'left':
                        old_r = old_r.transformed(kdb.Trans(kdb.Trans.R180))
                    previous = {str(k): v.tolist() for k,v in pnc_centers(old_r,L,pitch,dbu).items()}
                REPORT['measurements'].append(dict(layout=tag,cell=inst.cell.name,end=end,
                    pitch_um=pitch,taper_length_um=45.0,
                    dense_pitch_um=DENSE_PITCH,dense_periods=DENSE_PERIODS,dense_span_um=DENSE_SPAN,
                    before_centers_from_port_um=previous,
                    after_centers_from_port_um={str(k):v.tolist() for k,v in values.items()}))
        labels = [i for i in top.each_inst() if i.cell.name.startswith('LABEL_')]
        if expected_count > 1:
            check(tag+'_one_label_per_device',len(labels)==expected_count)
        if old_layout is not None:
            old = all_regions(old_layout, old_top)
            new = all_regions(layout, top)
            xor = {f'{n}/{d}': float((old.get((n,d),kdb.Region()) ^ new.get((n,d),kdb.Region())).area())*dbu**2
                   for n,d in sorted(set(old)|set(new)) if (n,d)!=(11,0)}
            check(tag+'_every_non_PNC_layer_exactly_unchanged', all(v==0 for v in xor.values()),xor_um2=xor)
            delta = (old[(11,0)] ^ new[(11,0)])
            outside = (delta-allowed_global).area()*dbu**2
            check(tag+'_PNC_changes_only_in_first_5um',outside==0,
                  outside_dense_region_xor_um2=float(outside),total_PNC_change_um2=float(delta.area()*dbu**2))
            check(tag+'_no_existing_PNC_removed', (old[(11,0)]-new[(11,0)]).is_empty())
            check(tag+'_physical_PNC_update_present',not delta.is_empty())
        REPORT['outputs'][filename] = {'sha256':hashlib.sha256((ROOT/'gds'/filename).read_bytes()).hexdigest(),
                                      'devices':len(insts),'pitch_um':pitch}
    REPORT.update(status='PASS' if not REPORT['failures'] else 'FAIL',
                  check_count=len(REPORT['checks']),devices_checked=count,
                  taper_ends_checked=2*count,
                  dense_pitch_um=DENSE_PITCH,dense_periods=DENSE_PERIODS,dense_span_um=DENSE_SPAN,
                  baseline_comparison_performed=all((BASE/f).exists() for _,f,_,_,_ in FILES))
    (ROOT/'gds/2611_taper_pitch_validation.json').write_text(json.dumps(REPORT,indent=2))
    (WORK/'validation_report.json').write_text(json.dumps(REPORT,indent=2))
    print(REPORT['status'],REPORT['check_count'],'checks',REPORT['failures'],flush=True)
    if REPORT['failures']:
        raise SystemExit(1)
    render()


def render():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    name='2611_design1_design2_L1400_18device_sweep.gds'
    for design in (1,2):
        for state, path in [('before',BASE/name),('after',ROOT/'gds'/name)]:
            if not path.exists():
                continue
            ly,top=read(path);dbu=ly.dbu
            inst=next(i for i in devices(ly,top) if i.cell.name.startswith(f'D{design}_01'))
            x0=700+inst.trans.disp.x*dbu
            y0=-200+inst.trans.disp.y*dbu
            bounds=(x0-21,y0-17,x0+51,y0+17)
            roi=box(bounds,dbu)
            fig,ax=plt.subplots(figsize=(12,5))
            for n in (14,11,15):
                r=region(ly,inst.cell,(n,0)).transformed(inst.trans)&roi
                lines=[]
                for p in r.each():
                    pts=[(q.x*dbu,q.y*dbu-y0) for q in p.each_point_hull()]
                    if pts:lines.append(pts+[pts[0]])
                if lines:
                    collection=LineCollection(lines,linewidths=.8,label=f'Layer {n}')
                    # Use the default matplotlib cycle, without a custom palette.
                    collection.set_color(ax._get_lines.get_next_color())
                    ax.add_collection(collection)
            ax.axvline(x0,linestyle=':',linewidth=1)
            ax.text(x0+.8,15,'Straight / taper junction',fontsize=9)
            ax.set(xlim=bounds[::2],ylim=(-17,17),xlabel='Global X (um)',
                   ylabel='Y relative to ridge (um)',
                   title=f'Design {design}: {state} - '+('5 um regular pitch' if state=='before' else 'five 1 um periods at connection; rest unchanged'))
            ax.set_aspect('equal');ax.legend(loc='lower left');fig.tight_layout()
            fig.savefig(WORK/f'design{design}_{state}.png',dpi=180);plt.close(fig)


if __name__ == '__main__':
    if len(sys.argv)>1 and sys.argv[1]=='capture':
        capture()
    else:
        run()
