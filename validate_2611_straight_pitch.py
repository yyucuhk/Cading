"""Measure the 5 um straight PNC pitch in the actual exported 18-device GDS.

Run after 2611_build_combined_sweep.py. Optional baseline files in
validation_2611_straight_pitch enable geometry-preservation comparisons.
These geometric checks are not full fabrication DRC or optical simulation.
"""
from pathlib import Path
import hashlib
import json
import klayout.db as kdb
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'validation_2611_straight_pitch'
STEM = ROOT / 'gds/2611_design1_design2_L1400_18device_sweep'
REPORT = {'checks': [], 'failures': []}


def check(name, passed, **details):
    REPORT['checks'].append(dict(name=name, passed=bool(passed), **details))
    if not passed:
        REPORT['failures'].append(name)


def read(path):
    layout = kdb.Layout()
    layout.read(str(path))
    return layout, layout.top_cell()


def regions(layout, cell):
    result = {}
    for idx in layout.layer_indexes():
        reg = kdb.Region(cell.begin_shapes_rec(idx)).merged()
        if not reg.is_empty():
            info = layout.get_info(idx)
            result[(info.layer, info.datatype)] = reg
    return result


def box(coords, dbu):
    return kdb.Region(kdb.Box(*(round(v / dbu) for v in coords)))


def all_layer_differences(old, new, dbu, exclude=()):
    return {f'{n}/{d}': float((old.get((n,d), kdb.Region()) ^
                             new.get((n,d), kdb.Region())).area()) * dbu**2
            for n,d in sorted(set(old) | set(new)) if (n,d) not in exclude}


def run():
    OUT.mkdir(exist_ok=True)
    manifest = json.loads(STEM.with_suffix('.json').read_text())
    layout, top = read(STEM.with_suffix('.gds'))
    dbu = layout.dbu
    instances = {i.cell.name: i for i in top.each_inst() if i.cell.name.startswith(('D1_', 'D2_'))}
    check('exactly_18_complete_devices', len(instances) == 18)
    old_layout = None
    old_cells = {}
    old_entries = {}
    if (OUT / 'before.gds').exists():
        old_layout, old_top = read(OUT / 'before.gds')
        old_cells = {i.cell.name: i.cell for i in old_top.each_inst()}
        old_entries = {e['case_id']: e for e in json.loads((OUT / 'before.json').read_text())['devices']}
        diff = all_layer_differences(regions(old_layout, old_top), regions(layout, top), dbu,
                                     exclude=((11,0),))
        check('entire_layout_all_non_PNC_layers_exactly_unchanged',
              all(v == 0 for v in diff.values()), xor_um2=diff)

    for entry in manifest['devices']:
        ident = entry['case_id']
        inst = instances[ident]
        regs = regions(layout, inst.cell)
        check(ident + '_independent_pitch_parameters',
              entry['straight_pnc_period_um'] == 5.0 and entry['taper_pnc_period_um'] == 5.0
              and abs(entry['pnc_taper_length_um']-45.0)<1e-9 and entry['pnc_cutoff_width_um']==5.0)
        check(ident + '_three_straight_sections', len(entry['straight_pnc_measurements']) == 3)
        allowed = kdb.Region()
        for index, section in enumerate(entry['straight_pnc_measurements']):
            a, b = section['start_um'], section['end_um']
            check(ident + f'_straight{index}_length', abs(a[1]-b[1]) < 1e-6 and
                  abs(abs(a[0]-b[0])-1400.0) < 1e-6)
            xmin, xmax = sorted((a[0], b[0]))
            cy = (a[1]+b[1])/2
            allowed += box((xmin-.001, cy-30, xmax+.001, cy+30), dbu)
            local = regs[(11,0)] & box((xmin+1, cy-30, xmax-1, cy+30), dbu)
            side_centers = {1: [], -1: []}
            for poly in local.merged().each():
                bounds = poly.bbox()
                y = (bounds.bottom + bounds.top)*dbu/2
                side_centers[1 if y > cy else -1].append((bounds.left + bounds.right)*dbu/2)
            expected = np.arange(xmin+5, xmax-1, 5.0)
            for side, values in side_centers.items():
                values = np.asarray(sorted(values))
                good_count = len(values) == len(expected) == 279
                spacings = np.diff(values)
                check(ident + f'_straight{index}_side{side}_actual_5um_grid',
                      good_count and np.allclose(values, expected, atol=.0011, rtol=0)
                      and np.allclose(spacings, 5.0, atol=.0011, rtol=0),
                      pncs=int(len(values)),
                      measured_min_pitch_um=float(min(spacings)) if len(spacings) else None,
                      measured_max_pitch_um=float(max(spacings)) if len(spacings) else None)
        if old_layout is not None:
            baseline = regions(old_layout, old_cells[ident])
            change = (baseline[(11,0)] ^ regs[(11,0)])
            outside = (change - allowed).area()*dbu**2
            check(ident + '_PNC_changes_confined_to_straights', outside == 0,
                  outside_straights_xor_um2=float(outside), total_PNC_change_um2=float(change.area()*dbu**2))
            previous = old_entries[ident]
            keys = ('stage_lengths_um', 'layer10_widths_um', 'layer15_widths_um',
                    'layer15_lengths_um', 'bend_pnc_measurements', 'tip_centers',
                    'straight_um', 'origin_um', 'label', 'label_origin_um')
            check(ident + '_bends_tapers_placement_and_labels_preserved',
                  all(entry[k] == previous[k] for k in keys))

    # The independent standalone files must agree with the same pitch change.
    for filename in ('2611_tapered_yue_dual_resist_L10_L12.gds', '2611_tapered_yue_design2.gds'):
        old_path = OUT / filename
        if old_path.exists():
            old, old_top = read(old_path)
            new, new_top = read(ROOT / 'gds' / filename)
            diff = all_layer_differences(regions(old, old_top), regions(new, new_top), dbu,
                                         exclude=((11,0),))
            check(filename + '_all_non_PNC_layers_unchanged', all(v == 0 for v in diff.values()), xor_um2=diff)

    # Retain the single-label layout and current device count.
    labels = [i for i in top.each_inst() if i.cell.name.startswith('LABEL_')]
    check('exactly_18_labels', len(labels) == 18)
    check('same_bend_and_taper_sweep',
          {(e['bend_period_um'], e['taper_total_um']) for e in manifest['devices']} ==
          {(p,t) for p in (2.,3.,5.) for t in (80.,100.,120.)})
    check('uniform_600um_vertical_pitch', np.allclose(np.diff(sorted(
          i.trans.disp.y*dbu for i in instances.values())), 600., atol=1e-8, rtol=0))
    check('all_left_bends_at_minus900', all(abs(
          regions(layout, i.cell)[(11,0)].bbox().left*dbu+i.trans.disp.x*dbu+900) < 1e-8
          for i in instances.values()))
    check('no_temporary_layers', all(layout.get_info(i).layer not in (998,999)
          or kdb.Region(top.begin_shapes_rec(i)).is_empty() for i in layout.layer_indexes()))
    REPORT.update(status='FAIL' if REPORT['failures'] else 'PASS', check_count=len(REPORT['checks']),
                  baseline_comparison_performed=old_layout is not None,
                  straight_pitch_definition='center_to_center', straight_pitch_um=5.0,
                  device_count=len(instances), straight_sections_measured=3*len(instances),
                  gds_sha256=hashlib.sha256(STEM.with_suffix('.gds').read_bytes()).hexdigest())
    output = ROOT / 'gds/2611_straight_pitch_validation.json'
    output.write_text(json.dumps(REPORT, indent=2))
    (OUT / 'validation_report.json').write_text(json.dumps(REPORT, indent=2))
    print(REPORT['status'], REPORT['check_count'], 'checks', REPORT['failures'], flush=True)
    if REPORT['failures']:
        raise SystemExit(1)
    render_detail(layout, instances[manifest['devices'][0]['case_id']].cell,
                  manifest['devices'][0]['straight_pnc_measurements'][0], dbu)


def render_detail(layout, cell, section, dbu):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    cy = section['start_um'][1]
    xmin = min(section['start_um'][0], section['end_um'][0]) + 20
    roi = box((xmin, cy-18, xmin+40, cy+18), dbu)
    fig, ax = plt.subplots(figsize=(11, 6))
    for layer in (14,11,15):
        lines = []
        reg = (kdb.Region(cell.begin_shapes_rec(layout.layer(layer,0))) & roi).merged()
        for poly in reg.each():
            pts = [(p.x*dbu, p.y*dbu-cy) for p in poly.each_point_hull()]
            if pts:
                lines.append(pts+[pts[0]])
        ax.add_collection(LineCollection(lines, linewidths=1.0, label=f'Layer {layer}'))
    ax.annotate('', xy=(xmin+10,15.5), xytext=(xmin+5,15.5),
                arrowprops={'arrowstyle':'<->'})
    ax.text(xmin+7.5,16.0,'5 um pitch',ha='center')
    ax.set(xlim=(xmin,xmin+40), ylim=(-18,18), xlabel='Local X (um)',
           ylabel='Distance from waveguide center (um)',
           title='Exported GDS: exact 5 um straight PNC pitch (center to center)')
    ax.set_aspect('equal')
    ax.legend(loc='lower right')
    fig.tight_layout()
    fig.savefig(OUT/'straight_pnc_5um.png',dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    run()
