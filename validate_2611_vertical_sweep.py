"""Check actual exported vertical layout; optionally compare approved old devices.

Run after 2611_build_combined_sweep.py. Optional baseline inputs:
  validation_2611_vertical/before.gds
  validation_2611_vertical/before.json
  validation_2611_vertical/protected_hashes.json
Original text may overlap windows: compare old = new UNION old text, never cut out a text box.
These are geometric regression checks, not foundry DRC or optical validation.
"""
from pathlib import Path
import hashlib
import json
import math

import gdsfactory as gf
import klayout.db as kdb

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'validation_2611_vertical'
STEM = ROOT / 'gds/2611_design1_design2_L1400_18device_sweep'
REPORT = {'checks': [], 'failures': []}


def check(name, passed, **details):
    REPORT['checks'].append(dict(name=name, passed=bool(passed), **details))
    if not passed:
        REPORT['failures'].append(name)


def regions(layout, cell):
    result = {}
    for index in layout.layer_indexes():
        info = layout.get_info(index)
        reg = kdb.Region(cell.begin_shapes_rec(index)).merged()
        if not reg.is_empty():
            result[(info.layer, info.datatype)] = reg
    return result


def text_polygons(text, size, x=0.0, y=0.0, dbu=0.001):
    component = gf.components.text(text=text, size=size, layer=(2, 0))
    transform = kdb.Trans(round(x/dbu), round(y/dbu))
    polys = []
    for shape in component.kdb_cell.shapes(gf.get_layer((2, 0))).each():
        if shape.is_polygon() or shape.is_box():
            poly = kdb.Polygon(shape.box) if shape.is_box() else shape.polygon
            polys.append(poly.transformed(transform))
    return polys


def original_length_text_region(entry, dbu):
    """Reconstruct old text on its original grid for union-based comparison.

    Original text sometimes overlapped a physical window. The correct
    baseline relation is new_mask UNION old_text == approved_old_mask,
    not new_mask == approved_old_mask MINUS old_text.
    """
    total = entry['meander_length_um']
    text = (f"ML{total:.1f}_TL{entry['taper_total_um']:g}" if entry['design'] == 1
            else f"DESIGN 2_L{entry['straight_um']:g}_ML{total:.1f}")
    tip = min(entry['tip_centers'], key=lambda point: point[0])
    reg = kdb.Region()
    for poly in text_polygons(text, 3.0, tip[0]+5.0, tip[1]+10.0, dbu):
        reg.insert(poly)
    return reg.merged()


def run():
    gf.gpdk.PDK.activate()
    OUT.mkdir(exist_ok=True)
    manifest = json.loads(STEM.with_suffix('.json').read_text())
    ly = kdb.Layout()
    ly.read(str(STEM.with_suffix('.gds')))
    top = ly.top_cell()
    dbu = ly.dbu
    instances = list(top.each_inst())
    devices = {i.cell.name: i for i in instances if i.cell.name.startswith(('D1_', 'D2_'))}
    labels = {i.cell.name[6:]: i for i in instances if i.cell.name.startswith('LABEL_')}
    frames = [i for i in instances if i.cell.name == 'SWEEP_FRAME_AND_GLOBAL_MARKERS']
    expected = {d['case_id'] for d in manifest['devices']}
    check('one_top_cell', len(ly.top_cells()) == 1)
    check('no_temporary_text_layer_exported', all(ly.get_info(i).layer != 998
          or kdb.Region(top.begin_shapes_rec(i)).is_empty() for i in ly.layer_indexes()))
    check('exactly_18_device_cells', len(devices) == 18 and set(devices) == expected)
    check('exactly_one_label_for_every_device', len(labels) == 18 and set(labels) == expected)
    check('one_frame_and_no_obsolete_label_group', len(frames) == 1 and len(instances) == 37)
    check('single_column_manifest', manifest['placement_mode'] == 'single_vertical_column_design1_then_design2'
          and manifest['design_x_offset_um'] == 0)
    check('no_gds_text_elements_outside_polygon_labels',
          not any(shape.is_text() for cell in ly.each_cell()
                  for layer in ly.layer_indexes() for shape in cell.shapes(layer).each()))

    old_layout = None
    old_devices = {}
    if (OUT/'before.gds').exists():
        old_layout = kdb.Layout()
        old_layout.read(str(OUT/'before.gds'))
        old_devices = {i.cell.name: i.cell for i in old_layout.top_cell().each_inst()
                       if i.cell.name.startswith(('D1_', 'D2_'))}
        old_manifest = json.loads((OUT/'before.json').read_text())
        old_entries = {e['case_id']: e for e in old_manifest['devices']}

    physical = []
    global_labels = []
    for index, entry in enumerate(manifest['devices']):
        ident = entry['case_id']
        inst = devices[ident]
        local = regions(ly, inst.cell)
        origin = [inst.trans.disp.x*dbu, inst.trans.disp.y*dbu]
        check(ident+'_row_order', entry['global_row'] == index and
              entry['design'] == (1 if index < 9 else 2))
        check(ident+'_actual_y_position', abs(origin[1]-(-700+600*index)) < 1e-8,
              y_um=origin[1])
        left = local[(11,0)].bbox().left*dbu + origin[0]
        check(ident+'_actual_left_bend_x', abs(left+900) < 1e-8, x_um=left)
        check(ident+'_only_translation', inst.trans.angle == 0 and not inst.trans.is_mirror())
        check(ident+'_three_equal_1400um_straights',
              all(abs(v-1400) < 1e-8 for v in entry['measured_straights_um']))
        check(ident+'_old_length_glyph_removal_recorded',
              entry['removed_length_label']['removed_polygon_count'] > 0)

        li = labels[ident]
        actual_label = regions(ly, li.cell)
        label = f"D{entry['design']}_T{entry['taper_total_um']:g}_P{entry['bend_period_um']:g}"
        expected_text = kdb.Region()
        for poly in text_polygons(label, entry['label_size_um'], dbu=dbu):
            expected_text.insert(poly)
        difference = (actual_label[(2,0)] ^ expected_text).area()*dbu**2
        check(ident+'_label_only_design_taper_period', set(actual_label) == {(2,0)} and
              entry['label'] == label and difference == 0, text=label, xor_um2=difference)
        expected_origin = [-700.0, origin[1]+250.0]
        actual_origin = [li.trans.disp.x*dbu, li.trans.disp.y*dbu]
        check(ident+'_label_position', actual_origin == expected_origin)
        global_labels.append((ident, actual_label[(2,0)].transformed(li.trans)))
        physical.append((ident, {layer: reg.transformed(inst.trans)
                               for layer, reg in local.items() if layer != (51,0)}))

        if old_layout is not None:
            baseline = regions(old_layout, old_devices[ident])
            original_text = original_length_text_region(old_entries[ident], dbu)
            restored = dict(local)
            restored[(2,0)] = (local[(2,0)] + original_text).merged()
            outside = ((baseline[(2,0)] ^ local[(2,0)]) - original_text).area()*dbu**2
            check(ident+'_no_changes_outside_original_length_text', outside < 0.005,
                  outside_old_text_xor_um2=float(outside))
            delta = {f'{n}/{d}': float((baseline.get((n,d), kdb.Region()) ^
                                     restored.get((n,d), kdb.Region())).area()*dbu**2)
                     for n,d in sorted(set(baseline)|set(local))}
            check(ident+'_approved_geometry_recovered_by_adding_only_old_text',
                  all(v < 0.005 for v in delta.values()), xor_um2=delta)
            for key in ('stage_lengths_um','layer10_widths_um','layer15_widths_um',
                        'layer15_lengths_um','bend_pnc_measurements','tip_centers','straight_um'):
                check(ident+'_unchanged_'+key, entry[key] == old_entries[ident][key])

    overlap_count = 0
    for i, (name, a) in enumerate(physical):
        for other_name, b in physical[i+1:]:
            for layer in set(a)&set(b):
                if a[layer].bbox().overlaps(b[layer].bbox()):
                    area = (a[layer]&b[layer]).area()*dbu**2
                    if area > 0.001:
                        overlap_count += 1
    check('no_interdevice_fabrication_overlaps', overlap_count == 0, overlaps=overlap_count)
    for name, text in global_labels:
        check(name+'_label_clear_of_all_device_masks', all(
            (text&reg).is_empty() for _, masks in physical for reg in masks.values()))
    check('all_labels_mutually_separate', all((a&b).is_empty()
          for i, (_,a) in enumerate(global_labels) for _,b in global_labels[i+1:]))
    frame_regions = regions(ly, frames[0].cell)
    frame_box = frame_regions[(0,0)].transformed(frames[0].trans)
    frame_overlap = []
    for name, masks in physical:
        for layer, reg in masks.items():
            if not (reg-frame_box).is_empty():
                frame_overlap.append([name, str(layer)])
    check('frame_encloses_all_fabrication_geometry', not frame_overlap, outside=frame_overlap)
    # Global overlay marks must be outside every device and every single label.
    marks_clear = True
    for layer in ((2,0),(10,0)):
        marks = frame_regions.get(layer,kdb.Region()).transformed(frames[0].trans)
        for _, masks in physical:
            marks_clear &= (marks & masks.get(layer,kdb.Region())).is_empty()
        if layer == (2,0):
            for _, text in global_labels:
                marks_clear &= (marks & text).is_empty()
    check('global_alignment_marks_clear_of_devices_and_labels', marks_clear)
    if (OUT/'protected_hashes.json').exists():
        for name, expected_hash in json.loads((OUT/'protected_hashes.json').read_text()).items():
            check('unchanged_standalone_'+name, hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected_hash)
    REPORT.update(status='FAIL' if REPORT['failures'] else 'PASS',
                  check_count=len(REPORT['checks']), baseline_comparison_performed=old_layout is not None,
                  frame_bounds_um=manifest['frame_bounds_um'], label_format='D{design}_T{taper_um}_P{bend_period_um}',
                  y_origins_um=[e['origin_um'][1] for e in manifest['devices']],
                  gds_sha256=hashlib.sha256(STEM.with_suffix('.gds').read_bytes()).hexdigest())
    (OUT/'validation_report.json').write_text(json.dumps(REPORT,indent=2))
    (ROOT/'gds/2611_vertical_sweep_validation.json').write_text(json.dumps(REPORT,indent=2))
    print(REPORT['status'], REPORT['check_count'], 'checks', REPORT['failures'])
    if REPORT['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    run()
