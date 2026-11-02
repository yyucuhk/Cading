"""Validate symmetric 2611 release windows and bounded Layer-14 hole fills.

Run `python validate_2611_windows.py baseline` before a mask change, and
`python validate_2611_windows.py` after it. The ordinary all-layer sweep
validator is run separately. Geometry comparisons include exported GDS.
"""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
import sys

import gdsfactory as gf
import klayout.db as kdb
import numpy as np
from validate_2611_sweep import snapshot, file_snapshot, box_region, compare, export

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'validation_2611_windows'
OUT.mkdir(exist_ok=True)
REPORT = {'checks': [], 'failures': []}
TOL = 0.005


def load_source():
    spec = importlib.util.spec_from_file_location('layout2611_windows', ROOT / '2611_tapered_yue_dual_resist.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    gf.gpdk.PDK.activate()
    return module


def check(name, passed, **details):
    REPORT['checks'].append({'name': name, 'passed': bool(passed), **details})
    if not passed:
        REPORT['failures'].append(name)
    (OUT / 'validation_report.json').write_text(json.dumps(REPORT, indent=2))
    print(json.dumps(REPORT['checks'][-1]), flush=True)


def read_snapshot(path):
    layout = kdb.Layout()
    layout.read(str(path))
    return file_snapshot(layout, layout.top_cell()), layout.dbu


def holes(region):
    small, large = kdb.Region(), kdb.Region()
    bounds = []
    for polygon in region.merged().each():
        for index in range(polygon.holes()):
            hole = kdb.Polygon(list(polygon.each_point_hole(index)))
            bounds.append(hole.bbox())
            if hole.area() < 1000000:  # <1 um^2 at the verified 1 nm DBU
                small.insert(hole)
            else:
                large.insert(hole)
    return small, large, bounds


def render_region(regions, title, filename, bounds, dbu):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import PathPatch
    fig, ax = plt.subplots(figsize=(13, 5))
    roi = box_region(bounds, dbu)
    region = (regions[(14, 0)] & roi).merged()
    for polygon in region.each():
        vertices, codes = [], []
        loops = [list(polygon.each_point_hull())]
        loops.extend(list(polygon.each_point_hole(i)) for i in range(polygon.holes()))
        for loop in loops:
            if not loop:
                continue
            points = [(p.x * dbu, p.y * dbu) for p in loop]
            vertices.extend(points + [points[0]])
            codes.extend([MplPath.MOVETO] + [MplPath.LINETO] * (len(points) - 1) + [MplPath.CLOSEPOLY])
        if vertices:
            ax.add_patch(PathPatch(MplPath(vertices, codes), alpha=0.3))
    for layer in (14, 10):
        points = []
        for polygon in (regions.get((layer, 0), kdb.Region()) & roi).merged().each():
            loops = [list(polygon.each_point_hull())]
            loops.extend(list(polygon.each_point_hole(i)) for i in range(polygon.holes()))
            for loop in loops:
                pp = [(p.x * dbu, p.y * dbu) for p in loop]
                if pp:
                    points.extend(pp + [pp[0], (np.nan, np.nan)])
        if points:
            pp = np.asarray(points)
            ax.plot(pp[:, 0], pp[:, 1], linewidth=0.8, label=f'Layer {layer}')
    x0, y0, x1, y1 = bounds
    ax.set(xlim=(x0, x1), ylim=(y0, y1), xlabel='Device-local x (um)', ylabel='Device-local y (um)', title=title)
    ax.set_aspect('equal', adjustable='box')
    ax.legend(loc='upper right')
    fig.tight_layout()
    fig.savefig(OUT / filename, dpi=180)
    plt.close(fig)


def main():
    module = load_source()
    lengths = tuple(module.MEANDER_STRAIGHT_LENGTHS)
    if len(sys.argv) > 1 and sys.argv[1] == 'baseline':
        for i, length in enumerate(lengths, 1):
            export(module.make_device(length), OUT / f'before_{i}.gds')
        print('WINDOW_BASELINES_CAPTURED')
        return

    baseline_available = all((OUT / f'before_{i}.gds').exists() for i in (1, 2, 3))
    snapshots, devices = [], []
    for i, length in enumerate(lengths, 1):
        device = module.make_device(length)
        devices.append(device)
        new = snapshot(device)
        snapshots.append(new)
        dbu = float(device.kcl.dbu)
        assert abs(dbu - 0.001) < 1e-12
        export(device, OUT / f'after_{i}.gds')
        reopened, _ = read_snapshot(OUT / f'after_{i}.gds')
        delta = compare(new, reopened, dbu)
        check(f'device_{i}_export_preserves_all_layers', all(a <= TOL for a in delta.values()), xor_um2=delta)

        tip_boxes = kdb.Region()
        block_zone = kdb.Region()
        for cx, cy in device.info['tip_centers']:
            tip_boxes += box_region((cx - 3.1, cy - 12, cx + 3.1, cy + 12), dbu)
            for sx in (-1, 1):
                for sy in (-1, 1):
                    wx, wy = cx + sx * 1.25, cy + sy * 7.3
                    block_zone += box_region((wx - 1, wy - 4, wx + 1, wy + 4), dbu)
            block_zone -= box_region((cx - 3, cy - 3.5, cx + 3, cy + 3.5), dbu)
        for layer in (2, 14):
            overlap = reopened[(layer, 0)] & block_zone
            check(f'device_{i}_layer_{layer}_no_companion_blocks_at_either_tip', overlap.area() == 0,
                  unwanted_block_area_um2=float(overlap.area() * dbu**2))

        small, pads, hole_bounds = holes(reopened[(14, 0)])
        check(f'device_{i}_layer14_no_small_center_holes', small.is_empty(),
              remaining_small_hole_area_um2=float(small.area() * dbu**2),
              remaining_hole_bboxes_um=[[b.left*dbu,b.bottom*dbu,b.right*dbu,b.top*dbu] for b in hole_bounds])
        check(f'device_{i}_fills_both_matching_taper_junctions',
              device.info['layer14_taper_holes_filled'] == 2,
              holes_filled=int(device.info['layer14_taper_holes_filled']),
              fill_area_um2=float(device.info['layer14_taper_fill_area_um2']))

        # Complete 180-degree end-to-end window symmetry, excluding label text.
        right_roi = box_region((150.001, -1000, 10000, -25), dbu)
        label = box_region(device.info['label_bbox_um'], dbu)
        for layer in (2, 14):
            reg = reopened[(layer, 0)]
            if layer == 2:
                reg = reg - label
            asymmetric = (reg ^ reg.transformed(kdb.Trans(kdb.Trans.R180))) & right_roi
            check(f'device_{i}_layer_{layer}_two_complete_ends_symmetric',
                  asymmetric.area() * dbu**2 <= TOL,
                  rotated_xor_um2=float(asymmetric.area() * dbu**2))

        if baseline_available:
            old, _ = read_snapshot(OUT / f'before_{i}.gds')
            unchanged = {layer: old[layer] ^ reopened[layer] for layer in old if layer not in ((2,0),(14,0))}
            check(f'device_{i}_all_non_window_layers_unchanged', all(r.area() == 0 for r in unchanged.values()),
                  xor_um2={f'{l[0]}/{l[1]}':float(r.area()*dbu**2) for l,r in unchanged.items()})
            old_small, old_pads, _ = holes(old[(14, 0)])
            check(f'device_{i}_previous_center_holes_now_fully_exposed',
                  (old_small - reopened[(14,0)]).is_empty(),
                  formerly_unexposed_area_um2=float(old_small.area()*dbu**2))
            check(f'device_{i}_support_pad_openings_preserved', (old_pads ^ pads).area() == 0,
                  pad_hole_difference_um2=float((old_pads ^ pads).area()*dbu**2))
            for layer in (2, 14):
                allowed = tip_boxes.dup()
                if layer == 14:
                    allowed += old_small
                unexpected = (old[(layer,0)] ^ reopened[(layer,0)]) - allowed
                check(f'device_{i}_layer_{layer}_changes_confined_to_requested_regions',
                      unexpected.area()*dbu**2 <= TOL,
                      unexpected_difference_um2=float(unexpected.area()*dbu**2))
        if i == 1:
            if baseline_available:
                render_region(old, 'Before: left output window and remaining blocks', 'before_left.png', (-265,182,-175,218), dbu)
                render_region(old, 'Before: right taper junction center hole', 'before_right.png', (175,-218,265,-182), dbu)
            render_region(reopened, 'After: left output window, Layer 10 tips retained', 'after_left.png', (-265,182,-175,218), dbu)
            render_region(reopened, 'After: right output window, Layer 10 tips retained', 'after_right.png', (175,-218,265,-182), dbu)

    # Actual output is regenerated by the ordinary sweep validator first.
    final = kdb.Layout()
    final.read(str(ROOT / 'gds/2611_tapered_yue_dual_resist_L10_L12.gds'))
    instances = [inst for inst in final.top_cell().each_inst()
                 if final.find_layer(51,0) is not None
                 and not kdb.Region(inst.cell.begin_shapes_rec(final.find_layer(51,0))).is_empty()]
    instances.sort(key=lambda inst: inst.trans.disp.y)
    check('export_contains_exactly_three_complete_devices', len(instances)==3)
    for i, instance in enumerate(instances, 1):
        delta = compare(snapshots[i-1], file_snapshot(final, instance.cell), final.dbu)
        check(f'export_device_{i}_matches_validated_window_geometry', all(a<=TOL for a in delta.values()), xor_um2=delta)
    positions = [[inst.trans.disp.x*final.dbu,inst.trans.disp.y*final.dbu] for inst in instances]
    left_bend_edges = [
        kdb.Region(inst.cell.begin_shapes_rec(final.layer(11, 0))).bbox().left
        * final.dbu + inst.trans.disp.x * final.dbu
        for inst in instances
    ]
    check('export_keeps_aligned_left_bends_and_600um_pitch',
          len(positions)==3 and
          np.allclose(left_bend_edges,module.LEFT_BEND_ALIGNMENT_X,atol=1e-9,rtol=0) and
          np.allclose(np.diff(np.asarray(positions)[:,1]),600,atol=1e-9,rtol=0),
          positions_um=positions, left_bend_edges_um=left_bend_edges)
    REPORT['baseline_comparison_performed'] = baseline_available
    REPORT['status'] = 'PASS' if not REPORT['failures'] else 'FAIL'
    (OUT / 'validation_report.json').write_text(json.dumps(REPORT,indent=2))
    if REPORT['failures']:
        raise AssertionError(', '.join(REPORT['failures']))
    print('SYMMETRIC_WINDOWS_VALIDATION_PASS', flush=True)


if __name__ == '__main__':
    main()
