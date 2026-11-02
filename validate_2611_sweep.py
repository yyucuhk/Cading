"""Regression checks for the complete 2611 device factory.

Run: python validate_2611_sweep.py
An optional validation_2611/baseline_bottom.gds is compared against device 1.
Only the device text label is excluded from between-length endpoint comparisons.
All geometry tests operate on KLayout regions, not sampled cross-sections.
"""
from __future__ import annotations

from functools import partial
import importlib.util
import json
from pathlib import Path
import traceback

import gdsfactory as gf
import klayout.db as kdb
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'validation_2611'
OUT.mkdir(exist_ok=True)
REPORT: dict = {'checks': [], 'failures': []}
# 0.005 um^2 permits only small GDS-fracture/grid differences over entire masks.
# In-memory repeatability and both-end comparisons require EXACT zero XOR.
SERIAL_TOL_UM2 = 0.005


def check(name, passed, **details):
    item = {'name': name, 'passed': bool(passed), **details}
    REPORT['checks'].append(item)
    if not passed:
        REPORT['failures'].append(name)
    (OUT / 'validation_report.json').write_text(json.dumps(REPORT, indent=2))
    print(json.dumps(item), flush=True)


def snapshot(component):
    return {
        tuple(map(int, layer)): kdb.Region(
            component.kdb_cell.begin_shapes_rec(gf.get_layer(tuple(layer)))).merged()
        for layer in component.layers
    }


def file_snapshot(layout, cell):
    return {
        (layout.get_info(index).layer, layout.get_info(index).datatype):
            kdb.Region(cell.begin_shapes_rec(index)).merged()
        for index in layout.layer_indices()
        if not kdb.Region(cell.begin_shapes_rec(index)).is_empty()
    }


def box_region(bounds, dbu):
    return kdb.Region(kdb.Box(*[int(round(float(value) / dbu)) for value in bounds]))


def compare(a, b, dbu, roi=None, exclude=None):
    differences = {}
    for layer in sorted(set(a) | set(b)):
        ra = a.get(layer, kdb.Region())
        rb = b.get(layer, kdb.Region())
        delta = ra ^ rb
        if roi is not None:
            delta &= roi
        if exclude is not None and layer == (2, 0):
            delta -= exclude
        differences[f'{layer[0]}/{layer[1]}'] = float(delta.area() * dbu**2)
    return differences


def export(component, path):
    opts = gf.kdb.SaveLayoutOptions()
    opts.gds2_max_vertex_count = 2000
    component.write_gds(path, save_options=opts, with_metadata=False)


def render(path, prefix, length):
    # Render actual exported polygon contours. One axes per image.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ly = kdb.Layout()
    ly.read(str(path))
    cell = ly.top_cell()
    dbu = ly.dbu
    views = {
        'full': (-max(450, length / 2 + 160), max(450, length / 2 + 160), -365, 365),
        'left_end': (-440, -135, 140, 265),
        'right_end': (135, 440, -265, -140),
    }
    for view, (x0, x1, y0, y1) in views.items():
        fig, ax = plt.subplots(figsize=(14, 5))
        clip = box_region((x0, y0, x1, y1), dbu)
        for layer in (2, 10, 11, 13, 14, 15, 16):
            region = kdb.Region(cell.begin_shapes_rec(ly.layer(layer, 0))) & clip
            points = []
            for polygon in region.each():
                hull = [(p.x * dbu, p.y * dbu) for p in polygon.each_point_hull()]
                if hull:
                    points.extend(hull + [hull[0], (np.nan, np.nan)])
                for index in range(polygon.holes()):
                    hole = [(p.x * dbu, p.y * dbu) for p in polygon.each_point_hole(index)]
                    if hole:
                        points.extend(hole + [hole[0], (np.nan, np.nan)])
            if points:
                p = np.asarray(points)
                ax.plot(p[:, 0], p[:, 1], linewidth=0.7, label=f'Layer {layer}')
        ax.set(xlim=(x0, x1), ylim=(y0, y1), xlabel='Local x (um)', ylabel='Local y (um)',
               title=f'Device {prefix}: meander straight {length:g} um / {view.replace("_", " ")}')
        ax.set_aspect('equal', adjustable='box')
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1))
        fig.tight_layout()
        fig.savefig(OUT / f'device_{prefix}_{view}.png', dpi=180)
        plt.close(fig)


def main():
    spec = importlib.util.spec_from_file_location('layout2611', ROOT / '2611_tapered_yue_dual_resist.py')
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    gf.gpdk.PDK.activate()
    lengths = tuple(m.MEANDER_STRAIGHT_LENGTHS)
    check('requested_sweep', lengths == (200.0, 800.0, 1400.0), lengths_um=list(lengths))
    check('vertical_pitch', m.DEVICE_Y_PITCH == 600.0, pitch_um=m.DEVICE_Y_PITCH)
    devices, snapshots = [], []
    for index, length in enumerate(lengths):
        device = m.make_device(length)
        dbu = device.kcl.dbu
        devices.append(device)
        snapshots.append(snapshot(device))
        export(device, OUT / f'device_{index + 1}.gds')
        check(f'device_{index + 1}_required_layers',
              {(2, 0), (10, 0), (11, 0), (13, 0), (14, 0), (15, 0), (16, 0)} <= set(snapshots[-1]),
              layers=[list(layer) for layer in sorted(snapshots[-1])])
        for previous in range(index):
            delta = compare(snapshots[previous], snapshot(devices[previous]), dbu)
            check(f'build_{index + 1}_does_not_mutate_device_{previous + 1}',
                  all(area == 0 for area in delta.values()), xor_um2=delta)

    # Full-area, every-layer comparisons of BOTH complete end assemblies.
    # The meander occupies the opposite half-planes; no small tip patches.
    left_roi = box_region((-10000, 0, -m.S_BEND_HALF_LENGTH, 1000), dbu)
    right_roi = box_region((m.S_BEND_HALF_LENGTH, -1000, 10000, 0), dbu)
    for index in (1, 2):
        label_exclusion = box_region(devices[0].info['label_bbox_um'], dbu)
        label_exclusion += box_region(devices[index].info['label_bbox_um'], dbu)
        for side, roi in (('left', left_roi), ('right', right_roi)):
            delta = compare(snapshots[0], snapshots[index], dbu, roi, label_exclusion)
            check(f'device_{index + 1}_{side}_complete_assembly_equals_bottom',
                  all(area == 0 for area in delta.values()), xor_um2=delta)

    # A repeated call after the others must match: this detects shared-cell
    # mutations even when earlier output components have been flattened.
    for index in (2, 1, 0):
        repeated = m.make_device(lengths[index])
        delta = compare(snapshots[index], snapshot(repeated), dbu)
        check(f'device_{index + 1}_build_order_independent',
              all(area == 0 for area in delta.values()), xor_um2=delta)

    baseline_path = OUT / 'baseline_bottom.gds'
    if baseline_path.exists():
        baseline = kdb.Layout()
        baseline.read(str(baseline_path))
        delta = compare(file_snapshot(baseline, baseline.top_cell()), snapshots[0], dbu)
        check('bottom_preserved_against_pre_refactor_reference',
              all(area <= SERIAL_TOL_UM2 for area in delta.values()), xor_um2=delta,
              serialization_tolerance_um2=SERIAL_TOL_UM2)

    ridge_areas = []
    for index, (device, regions, length) in enumerate(zip(devices, snapshots, lengths, strict=True), 1):
        centers = np.asarray(device.info['connection_centers'], dtype=float)
        orientations = device.info['connection_orientations']
        check(f'device_{index}_fixed_connection_centers',
              np.allclose(centers, np.asarray(devices[0].info['connection_centers']), rtol=0, atol=1e-12),
              centers_um=centers.tolist())
        angle = np.radians(float(orientations[0]))
        rotation = np.array([[np.cos(angle), np.sin(angle)], [-np.sin(angle), np.cos(angle)]])
        local_delta = (centers[1] - centers[0]) @ rotation
        half_lead = (local_delta[0] + length) / 2
        turn = 1 if local_delta[1] > 0 else -1
        path = gf.Path()
        path.append(gf.path.straight(length=half_lead))
        path.append(gf.path.euler(radius=m.S_BEND_COMPACT_180_RADIUS,
                    angle=turn * 180, p=m.S_BEND_EULER_P, use_eff=True))
        path.append(gf.path.straight(length=length))
        path.append(gf.path.euler(radius=m.S_BEND_COMPACT_180_RADIUS,
                    angle=-turn * 180, p=m.S_BEND_EULER_P, use_eff=True))
        path.append(gf.path.straight(length=half_lead))
        path.rotate(float(orientations[0]))
        path.move(centers[0])
        check(f'device_{index}_meander_reaches_both_ports',
              np.allclose(path.points[0], centers[0], atol=1e-8) and
              np.allclose(path.points[-1], centers[1], atol=1e-8),
              middle_straight_um=length, outer_lead_um=float(half_lead))

        tether = m.phononic_tether(width=m.REVIEW_TETHER_CONNECTION_WIDTH,
            segments=m.PHONONIC_SUPPORT_SEGMENTS, gaps=m.PHONONIC_SUPPORT_GAPS,
            overshoot=m.PHONONIC_SUPPORT_OVERSHOOT)
        tether_length = float(tether.ports['a2'].center[0] - tether.ports['a1'].center[0])
        for layer, width in ((15, m.RIDGE_EXTENSION_END_WIDTH),
                             (14, m.S_BEND_WIDTH + 2 * tether_length)):
            expected = gf.path.extrude(path, cross_section=partial(
                gf.cross_section.cross_section, width=width, layer=(layer, 0)))
            expected_region = snapshot(expected)[(layer, 0)]
            missing = expected_region - regions[(layer, 0)]
            area = missing.area() * dbu**2
            check(f'device_{index}_layer_{layer}_entire_meander_covered',
                  area <= SERIAL_TOL_UM2, missing_um2=float(area))
        for layer in (15, 16):
            components = list(regions[(layer, 0)].merged().each())
            check(f'device_{index}_layer_{layer}_connected', len(components) == 1,
                  connected_components=len(components))
        label_roi = box_region(device.info['label_bbox_um'], dbu)
        uncovered = (regions[(2, 0)] - regions[(14, 0)]) - label_roi
        check(f'device_{index}_layer14_covers_layer2_mask',
              uncovered.area() * dbu**2 <= SERIAL_TOL_UM2,
              missing_um2=float(uncovered.area() * dbu**2), excludes='device label only')
        ridge_areas.append(float(regions[(15, 0)].area() * dbu**2))
        render(OUT / f'device_{index}.gds', index, length)
    increments = np.diff(ridge_areas)
    check('ridge_area_matches_only_straight_length_sweep',
          np.allclose(increments, 2 * np.diff(lengths) * m.RIDGE_EXTENSION_END_WIDTH, atol=0.01, rtol=0),
          ridge_areas_um2=ridge_areas, increments_um2=increments.tolist())

    # Export the real top-level GDS, then inspect each complete device cell
    # and each instance placement after serialization.
    top = m.build_layout()
    output_path = ROOT / 'gds' / '2611_tapered_yue_dual_resist_L10_L12.gds'
    output_path.parent.mkdir(exist_ok=True)
    export(top, output_path)
    final = kdb.Layout()
    final.read(str(output_path))
    names = list(top.info['device_cells'])
    actual_positions = []
    for index, name in enumerate(names):
        cell = final.cell(name)
        check(f'export_device_{index + 1}_present', cell is not None, cell_name=name)
        if cell is not None:
            delta = compare(snapshots[index], file_snapshot(final, cell), dbu)
            check(f'export_device_{index + 1}_all_layers_preserved',
                  all(area <= SERIAL_TOL_UM2 for area in delta.values()), xor_um2=delta)
        instances = [i for i in final.top_cell().each_inst() if i.cell.name == name]
        check(f'export_device_{index + 1}_one_instance', len(instances) == 1, instance_count=len(instances))
        if len(instances) == 1:
            transform = instances[0].trans
            actual_positions.append([transform.disp.x * dbu, transform.disp.y * dbu])
    if len(actual_positions) == 3:
        check('export_actual_600um_y_pitch_and_same_x',
              np.allclose(np.diff(np.asarray(actual_positions)[:, 1]), 600.0, rtol=0, atol=1e-9)
              and np.ptp(np.asarray(actual_positions)[:, 0]) == 0,
              device_origins_um=actual_positions)
    REPORT['output_file'] = str(output_path.relative_to(ROOT))
    REPORT['status'] = 'PASS' if not REPORT['failures'] else 'FAIL'
    (OUT / 'validation_report.json').write_text(json.dumps(REPORT, indent=2))
    if REPORT['failures']:
        raise AssertionError('Failed checks: ' + ', '.join(REPORT['failures']))
    print('ALL_LAYER_SWEEP_VALIDATION_PASS', flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception:
        REPORT['status'] = 'FAIL'
        REPORT['exception'] = traceback.format_exc()
        (OUT / 'validation_report.json').write_text(json.dumps(REPORT, indent=2))
        raise
