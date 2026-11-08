"""Verify only the requested layer policy changed in the six-device GDS.

The audit parses GDSII records directly, using Shapely rather than the CAD
builder or KLayout. Supply the previous comparison GDS with --baseline.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from validate_2611_straight_output_raw import read, union

STEM = '2611_design1_L6000_P10_straight_output_comparison'


def validate(root: Path, baseline: Path) -> dict:
    current = root / 'gds' / (STEM + '.gds')
    old, ot = read(baseline)
    new, nt = read(current)
    checks = []
    def check(name, ok, **detail):
        checks.append(dict(name=name, passed=bool(ok), **detail))
        if not ok:
            raise AssertionError((name, detail))
    check('same_top_cell', ot == nt)
    check('same_cell_inventory', set(old) == set(new))
    for name, cell in new.items():
        check(name + '_no16_or51', not ({16, 51} & set(cell['polys'])))
        check(name + '_references_and_transforms_unchanged',
              sorted(cell['refs']) == sorted(old[name]['refs']))
        expected = set(old[name]['polys']) - {16, 51}
        if name == 'FRAME_AND_MARKERS':
            expected = (expected - {2}) | {15}
        check(name + '_exact_layer_inventory', set(cell['polys']) == expected)
        for layer in sorted(expected):
            source = 2 if name == 'FRAME_AND_MARKERS' and layer == 15 else layer
            old_polys = sorted(tuple(p) for p in old[name]['polys'].get(source, []))
            new_polys = sorted(tuple(p) for p in cell['polys'].get(layer, []))
            if old_polys == new_polys:
                area = 0.0  # identical raw boundary coordinates, stronger than area equality
            else:
                a, b = union(old[name], source), union(cell, layer)
                area = a.symmetric_difference(b).area
            check(name + f'_L{source}_to_L{layer}_exact', area < 1e-7,
                  xor_um2=area)
    refs = new[nt]['refs']
    devices = [(n, xy) for n, xy in refs if n.startswith('D1_L')]
    labels = [(n, xy) for n, xy in refs if n.startswith('LABEL_')]
    check('six_devices_six_labels_one_frame', len(devices) == 6 and len(labels) == 6 and len(refs) == 13)
    check('labels_still_on_layer2', all(set(new[n]['polys']) == {2} for n, _ in labels))
    check('600um_row_pitch', sorted(xy[1] for _, xy in devices) == [-700000,-100000,500000,1100000,1700000,2300000])
    for name, xy in devices:
        check(name + '_left_bend_aligned', abs(min(x for p in new[name]['polys'][11] for x, y in p)*.001 + xy[0]*.001 + 900) < 1e-7)
    frame = new['FRAME_AND_MARKERS']
    moved = union(frame, 15)
    previous = union(old['FRAME_AND_MARKERS'], 2)
    check('all_marker_polygons_moved_not_duplicated',
          len(frame['polys'][15]) == len(old['FRAME_AND_MARKERS']['polys'][2]) and not frame['polys'].get(2))
    check('marker_bounds_exactly_unchanged', moved.bounds == previous.bounds)
    counts = {str(layer): sum(len(c['polys'].get(layer, [])) for c in old.values()) for layer in (16,51)}
    result = dict(status='PASS', check_count=len(checks), checks=checks,
                  scope='six-device L6000 straight-output comparison only',
                  removed_layers=[16,51], baseline_removed_polygon_counts=counts,
                  marker_source_layer=2, marker_target_layer=15,
                  moved_marker_polygons=len(frame['polys'][15]),
                  marker_area_um2=moved.area,
                  layer2_release_windows_and_labels_preserved=True,
                  device_count=len(devices), label_count=len(labels),
                  baseline_sha256=hashlib.sha256(baseline.read_bytes()).hexdigest(),
                  gds_sha256=hashlib.sha256(current.read_bytes()).hexdigest(),
                  method='Independent raw GDSII parser and Shapely; no generator or KLayout import',
                  optical_or_mechanical_simulation_performed=False)
    (root / 'gds/2611_layer_cleanup_validation.json').write_text(json.dumps(result, indent=2))
    print('PASS', len(checks), 'independent layer-cleanup checks')
    return result


def render_marker(root: Path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cells, top = read(root / 'gds' / (STEM + '.gds'))
    frame = cells['FRAME_AND_MARKERS']
    dx, dy = dict(cells[top]['refs'])['FRAME_AND_MARKERS']
    fx, fy = dx*.001, dy*.001
    xmin, ymin, _, _ = union(frame, 0).bounds
    fig, ax = plt.subplots(figsize=(8, 7))
    for layer in (2, 10, 15):
        xs, ys = [], []
        for poly in frame['polys'].get(layer, []):
            if any(xmin-120 <= x*.001 <= xmin+120 and ymin-120 <= y*.001 <= ymin+120 for x, y in poly):
                xs.extend([x*.001 + fx for x, _ in poly] + [float('nan')])
                ys.extend([y*.001 + fy for _, y in poly] + [float('nan')])
        if xs:
            ax.plot(xs, ys, label=f'Layer {layer}', linewidth=1.2)
    ax.set(xlim=(xmin+fx-120,xmin+fx+120), ylim=(ymin+fy-120,ymin+fy+120),
           xlabel='X (um)', ylabel='Y (um)', title='Alignment markers: Layer 2 moved to Layer 15')
    ax.set_aspect('equal'); ax.legend(); ax.grid(alpha=.2); fig.tight_layout()
    out=root/'validation_2611_layer_cleanup';out.mkdir(exist_ok=True)
    fig.savefig(out/'markers_layer15.png', dpi=180);plt.close(fig)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--render',action='store_true')
    args=p.parse_args();validate(args.root,args.baseline)
    if args.render:render_marker(args.root)
