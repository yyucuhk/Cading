"""Check placement-only 2611 alignment against a prior exported GDS.

Usage: python validate_2611_alignment.py path/to/before.gds
Compares every device layer in local coordinates and every placed device
against an exact translation of its prior instance. Global markers stay fixed.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import klayout.db as kdb
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'validation_2611_alignment'
OUT.mkdir(exist_ok=True)
REPORT = {'checks': [], 'failures': []}


def check(name, passed, **details):
    item = {'name': name, 'passed': bool(passed), **details}
    REPORT['checks'].append(item)
    if not passed:
        REPORT['failures'].append(name)
    (OUT / 'validation_report.json').write_text(json.dumps(REPORT, indent=2))
    print(json.dumps(item), flush=True)


def load(path):
    layout = kdb.Layout()
    layout.read(str(path))
    return layout


def regions(layout, cell):
    result = {}
    for index in layout.layer_indices():
        region = kdb.Region(cell.begin_shapes_rec(index)).merged()
        if not region.is_empty():
            info = layout.get_info(index)
            result[(info.layer, info.datatype)] = region
    return result


def instances(layout):
    index = layout.find_layer(51, 0)
    if index is None:
        return []
    result = [inst for inst in layout.top_cell().each_inst()
              if not kdb.Region(inst.cell.begin_shapes_rec(index)).is_empty()]
    return sorted(result, key=lambda inst: inst.trans.disp.y)


def differences(a, b, dbu):
    return {f'{layer[0]}/{layer[1]}': float(
            (a.get(layer, kdb.Region()) ^ b.get(layer, kdb.Region())).area() * dbu**2)
            for layer in sorted(set(a) | set(b))}


def transformed(values, transform):
    return {layer: region.transformed(transform) for layer, region in values.items()}


def global_only(layout, devices):
    result = regions(layout, layout.top_cell())
    for instance in devices:
        for layer, region in regions(layout, instance.cell).items():
            result[layer] = result.get(layer, kdb.Region()) - region.transformed(instance.trans)
    return result


def preview(layout, devices):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(15, 9))
    dbu = layout.dbu
    clip = kdb.Region(kdb.Box(-1020000, -1050000, 1000000, 850000))
    for layer in (2, 10, 11, 13, 14, 15, 16):
        reg = kdb.Region()
        for instance in devices:
            local = regions(layout, instance.cell).get((layer, 0), kdb.Region())
            reg += local.transformed(instance.trans)
        reg = (reg & clip).merged()
        points = []
        for polygon in reg.each():
            loops = [list(polygon.each_point_hull())]
            loops.extend(list(polygon.each_point_hole(i)) for i in range(polygon.holes()))
            for loop in loops:
                pp = [(p.x * dbu, p.y * dbu) for p in loop]
                if pp:
                    points.extend(pp + [pp[0], (np.nan, np.nan)])
        if points:
            pp = np.asarray(points)
            ax.plot(pp[:, 0], pp[:, 1], linewidth=0.65, label=f'Layer {layer}')
    ax.axvline(-900, linestyle='--', linewidth=1.2, label='Left PNC bend edge: x = -900 um')
    for number, (length, inst) in enumerate(zip((200, 800, 1400), devices), 1):
        ax.text(-870, inst.trans.disp.y * dbu + 270,
                f'Device {number} | straight {length} um', fontsize=10)
    ax.set(xlim=(-1020, 1000), ylim=(-1050, 850), xlabel='Chip x (um)', ylabel='Chip y (um)',
           title='Left-bend aligned devices | 600 um vertical pitch\nAll layers translated together; long outer channels continue outside this view')
    ax.set_aspect('equal', adjustable='box')
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1))
    fig.tight_layout()
    fig.savefig(OUT / 'aligned_devices.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline', type=Path)
    args = parser.parse_args()
    before = load(args.baseline)
    after = load(ROOT / 'gds/2611_tapered_yue_dual_resist_L10_L12.gds')
    check('same_database_grid', before.dbu == after.dbu, dbu_um=after.dbu)
    old, new = instances(before), instances(after)
    check('exactly_three_devices_before_and_after', len(old) == len(new) == 3,
          before_count=len(old), after_count=len(new))
    if len(old) != 3 or len(new) != 3:
        raise AssertionError('Expected three complete device instances')
    dbu = after.dbu
    expected_dx = (0.0, 300.0, 600.0)
    origins, left_edges = [], []
    for number, (a, b, dx) in enumerate(zip(old, new, expected_dx), 1):
        ra, rb = regions(before, a.cell), regions(after, b.cell)
        delta = differences(ra, rb, dbu)
        check(f'device_{number}_every_layer_exactly_unchanged',
              set(ra) == set(rb) and all(value == 0 for value in delta.values()), xor_um2=delta)
        displacement = [(b.trans.disp.x-a.trans.disp.x)*dbu,
                        (b.trans.disp.y-a.trans.disp.y)*dbu]
        check(f'device_{number}_rigid_x_translation_only',
              b.trans.angle == a.trans.angle == 0 and
              not b.trans.is_mirror() and not a.trans.is_mirror() and
              np.allclose(displacement, (dx, 0), atol=1e-9, rtol=0),
              translation_um=displacement)
        moved = transformed(ra, a.trans)
        moved = transformed(moved, kdb.Trans(int(round(dx/dbu)), 0))
        placed_delta = differences(moved, transformed(rb, b.trans), dbu)
        check(f'device_{number}_placed_all_layer_translation_xor_zero',
              all(value == 0 for value in placed_delta.values()), xor_um2=placed_delta)
        left_edge = rb[(11, 0)].bbox().left * dbu + b.trans.disp.x * dbu
        left_edges.append(left_edge)
        origins.append([b.trans.disp.x*dbu, b.trans.disp.y*dbu])
        check(f'device_{number}_left_bend_edge_at_minus_900', abs(left_edge + 900) < 1e-9,
              left_bend_edge_um=left_edge)
        for layer in (15, 16):
            count = len(list(rb[(layer, 0)].merged().each()))
            check(f'device_{number}_layer_{layer}_connected', count == 1, components=count)
    check('unchanged_y_positions_and_600um_pitch',
          np.allclose(np.asarray(origins)[:, 1], (-700, -100, 500), atol=1e-9, rtol=0),
          device_origins_um=origins, left_bend_edges_um=left_edges)
    delta = differences(global_only(before, old), global_only(after, new), dbu)
    check('chip_outline_and_global_markers_unchanged', all(value == 0 for value in delta.values()), xor_um2=delta)
    REPORT['status'] = 'PASS' if not REPORT['failures'] else 'FAIL'
    REPORT['baseline_gds'] = args.baseline.name
    (OUT / 'validation_report.json').write_text(json.dumps(REPORT, indent=2))
    if REPORT['failures']:
        raise AssertionError(', '.join(REPORT['failures']))
    preview(after, new)
    print('PLACEMENT_ONLY_ALIGNMENT_VALIDATION_PASS', flush=True)


if __name__ == '__main__':
    main()
