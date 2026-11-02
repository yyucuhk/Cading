"""Verify that only the alignment-tip layer changes in the exported GDS.

Run after generation:
    python validate_2611_tip_layers.py BEFORE.gds AFTER.gds
The baseline must be the approved pre-transfer GDS. Every device is compared
in local coordinates and its actual placement is checked separately.
"""
from pathlib import Path
import json
import sys
import klayout.db as kdb
import numpy as np
from validate_2611_sweep import file_snapshot, compare, box_region

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'validation_2611_tips'
OUT.mkdir(exist_ok=True)
REPORT = {'checks': [], 'failures': []}


def check(name, passed, **details):
    REPORT['checks'].append({'name': name, 'passed': bool(passed), **details})
    if not passed:
        REPORT['failures'].append(name)
    (OUT / 'validation_report.json').write_text(json.dumps(REPORT, indent=2))
    print(json.dumps(REPORT['checks'][-1]), flush=True)


def read(path):
    layout = kdb.Layout()
    layout.read(str(path))
    layer = layout.find_layer(51, 0)
    instances = [inst for inst in layout.top_cell().each_inst()
                 if layer is not None and
                 not kdb.Region(inst.cell.begin_shapes_rec(layer)).is_empty()]
    instances.sort(key=lambda inst: inst.trans.disp.y)
    return layout, instances


def render(regions, dbu):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for side, bounds in [('right', (220, -216, 250, -184)),
                         ('left', (-250, 184, -220, 216))]:
        fig, ax = plt.subplots(figsize=(7, 7))
        clip = box_region(bounds, dbu)
        for layer in (2, 10, 15):
            points = []
            for polygon in (regions[(layer, 0)] & clip).merged().each():
                loops = [list(polygon.each_point_hull())]
                loops.extend(list(polygon.each_point_hole(i)) for i in range(polygon.holes()))
                for loop in loops:
                    vertices = [(p.x*dbu, p.y*dbu) for p in loop]
                    if vertices:
                        points.extend(vertices + [vertices[0], (np.nan, np.nan)])
            if points:
                points = np.asarray(points)
                ax.plot(points[:, 0], points[:, 1], linewidth=1.2,
                        label=f'Layer {layer}')
        x0, y0, x1, y1 = bounds
        ax.set(xlim=(x0, x1), ylim=(y0, y1), xlabel='Device-local x (um)',
               ylabel='Device-local y (um)',
               title=f'{side.title()} taper end: both alignment tips on Layer 15')
        ax.set_aspect('equal')
        ax.legend()
        fig.tight_layout()
        fig.savefig(OUT / f'{side}_tips_layer15.png', dpi=180)
        plt.close(fig)


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    old, before = read(Path(sys.argv[1]))
    new, after = read(Path(sys.argv[2]))
    dbu = new.dbu
    check('same_database_grid', old.dbu == dbu, dbu_um=dbu)
    check('exactly_three_devices', len(before) == len(after) == 3)
    if len(before) != 3 or len(after) != 3:
        raise AssertionError('Expected exactly three complete devices')
    # Independent expected locations from the original coupler tick definition.
    expected = kdb.Region()
    for tx, ty in [(-235.0, 200.0), (235.0, -200.0)]:
        for dy in (-7.5, 7.5):
            expected += box_region((tx-0.25, ty+dy-5, tx+0.25, ty+dy+5), dbu)
    check('four_original_ticks_per_device', expected.count() == 4,
          count=expected.count(), area_um2=float(expected.area()*dbu**2))
    placements = []
    for i, (old_inst, new_inst) in enumerate(zip(before, after, strict=True), 1):
        a = file_snapshot(old, old_inst.cell)
        b = file_snapshot(new, new_inst.cell)
        check(f'device_{i}_original_tips_identified',
              (a[(10,0)] & expected).area() == expected.area() and
              (a[(15,0)] & expected).is_empty())
        removed = a[(10,0)] - b[(10,0)]
        added = b[(15,0)] - a[(15,0)]
        check(f'device_{i}_only_four_tips_removed_from_layer10',
              (removed ^ expected).is_empty() and (b[(10,0)] - a[(10,0)]).is_empty(),
              removed_area_um2=float(removed.area()*dbu**2),
              tips_remaining_on_layer10_um2=float((b[(10,0)] & expected).area()*dbu**2))
        check(f'device_{i}_same_four_tips_added_to_layer15',
              (added ^ expected).is_empty() and (a[(15,0)] - b[(15,0)]).is_empty(),
              added_area_um2=float(added.area()*dbu**2), count=len(list(added.merged().each())))
        check(f'device_{i}_tip_shapes_and_positions_exactly_preserved',
              (removed ^ added).is_empty(), xor_um2=float((removed ^ added).area()*dbu**2))
        unaffected = {layer: float((a.get(layer, kdb.Region()) ^ b.get(layer, kdb.Region())).area()*dbu**2)
                      for layer in set(a) | set(b) if layer not in ((10,0), (15,0))}
        check(f'device_{i}_every_other_layer_exactly_unchanged',
              all(v == 0 for v in unaffected.values()),
              xor_um2={f'{l[0]}/{l[1]}': v for l,v in sorted(unaffected.items())})
        combined_diff = (a[(10,0)] + a[(15,0)]) ^ (b[(10,0)] + b[(15,0)])
        check(f'device_{i}_combined_layer10_15_geometry_unchanged', combined_diff.is_empty())
        check(f'device_{i}_placement_unchanged', old_inst.trans == new_inst.trans,
              origin_um=[new_inst.trans.disp.x*dbu, new_inst.trans.disp.y*dbu])
        placements.append([new_inst.trans.disp.x*dbu, new_inst.trans.disp.y*dbu])
        ridge = (b[(15,0)] - expected).merged()
        check(f'device_{i}_original_ridge_still_connected', ridge.count() == 1,
              ridge_components=ridge.count(), isolated_alignment_tips=4)
        check(f'device_{i}_layer16_still_connected', b[(16,0)].merged().count() == 1)
        check(f'device_{i}_left_bend_remains_aligned',
              abs(b[(11,0)].bbox().left*dbu + new_inst.trans.disp.x*dbu + 900) < 1e-9)
        if i == 1:
            render(b, dbu)
    check('600um_y_pitch_preserved',
          np.allclose(np.diff(np.asarray(placements)[:,1]),600,rtol=0,atol=1e-9),
          positions_um=placements)
    # Global chip/marker geometry outside the three device cells is unchanged.
    def global_regions(layout, devices):
        regions = file_snapshot(layout, layout.top_cell())
        for instance in devices:
            for layer, region in file_snapshot(layout, instance.cell).items():
                regions[layer] -= region.transformed(instance.trans)
        return regions
    delta = compare(global_regions(old, before), global_regions(new, after), dbu)
    check('chip_outline_and_global_markers_unchanged', all(v == 0 for v in delta.values()), xor_um2=delta)
    REPORT['baseline_gds'] = str(sys.argv[1])
    REPORT['output_gds'] = str(sys.argv[2])
    REPORT['status'] = 'FAIL' if REPORT['failures'] else 'PASS'
    (OUT / 'validation_report.json').write_text(json.dumps(REPORT, indent=2))
    if REPORT['failures']:
        raise AssertionError(', '.join(REPORT['failures']))
    print('TIP_LAYER_TRANSFER_VALIDATION_PASS', flush=True)


if __name__ == '__main__':
    main()
