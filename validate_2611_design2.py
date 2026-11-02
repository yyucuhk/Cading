"""Check Design 2 geometry, inherited masks, and preservation of Design 1.

python validate_2611_design2.py baseline
python validate_2611_design2.py
The 840 nm quantity in this draft is CENTERLINE offset, not edge gap.
"""
from __future__ import annotations
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import types

import gdsfactory as gf
import klayout.db as kdb
import numpy as np
from validate_2611_sweep import snapshot, file_snapshot, box_region, compare, export

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'validation_2611_design2'
OUT.mkdir(exist_ok=True)
REPORT = {'checks': [], 'failures': [], 'separation_definition': 'centerline_to_centerline',
          'assumptions': ['40 um inherited first ridge transition',
                          '100 um applies to final Layer-10 1.0 to 0.2 um taper',
                          '840 nm is final centerline offset, not edge-to-edge gap']}
TOL = 0.005


def check(name, passed, **details):
    record = {'name': name, 'passed': bool(passed), **details}
    REPORT['checks'].append(record)
    if not passed:
        REPORT['failures'].append(name)
    (OUT / 'validation_report.json').write_text(json.dumps(REPORT, indent=2))
    print(json.dumps(record), flush=True)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def read_gds(path):
    layout = kdb.Layout()
    layout.read(str(path))
    return layout, file_snapshot(layout, layout.top_cell())


def fingerprint(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def baseline():
    gf.gpdk.PDK.activate()
    source_path = ROOT / '2611_tapered_yue_dual_resist.py'
    source = source_path.read_text()
    assert source.count('SLAB_TAPER_LENGTH = 40.0') == 1
    reference = types.ModuleType('design2_reference100')
    reference.__file__ = str(source_path)
    sys.modules[reference.__name__] = reference
    exec(compile(source.replace('SLAB_TAPER_LENGTH = 40.0', 'SLAB_TAPER_LENGTH = 100.0'),
                 str(source_path), 'exec'), reference.__dict__)
    device = reference.make_device(800.0)
    export(device, OUT / 'design1_with_100um_outer_taper.gds')
    data = {'label_bbox_um': list(device.info['label_bbox_um']),
            'files': {str(p.relative_to(ROOT)): fingerprint(p) for p in (
                source_path, ROOT / 'gds/2611_tapered_yue_dual_resist_L10_L12.gds')}}
    (OUT / 'baseline.json').write_text(json.dumps(data, indent=2))
    print('DESIGN2_REFERENCE_BASELINE_CAPTURED')


def outlines(region, dbu):
    points = []
    for polygon in region.merged().each():
        loops = [list(polygon.each_point_hull())]
        loops += [list(polygon.each_point_hole(i)) for i in range(polygon.holes())]
        for loop in loops:
            pp = [(p.x * dbu, p.y * dbu) for p in loop]
            if pp:
                points += pp + [pp[0], (np.nan, np.nan)]
    return np.asarray(points)


def render(regions, dbu, bounds, layers, title, filename, equal=True):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(14, 5.5))
    roi = box_region(bounds, dbu)
    for layer in layers:
        region = regions.get((layer, 0), kdb.Region()) & roi
        pp = outlines(region, dbu)
        if len(pp):
            ax.plot(pp[:, 0], pp[:, 1], linewidth=0.9, label=f'Layer {layer}')
    x0, y0, x1, y1 = bounds
    ax.set(xlim=(x0, x1), ylim=(y0, y1), xlabel='Local x (um)', ylabel='Local y (um)', title=title)
    if equal:
        ax.set_aspect('equal', adjustable='box')
    else:
        ax.text(0.01, 0.02, 'Transverse scale expanded; coordinates are physical dimensions',
                transform=ax.transAxes, fontsize=9)
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1))
    fig.tight_layout()
    fig.savefig(OUT / filename, dpi=180)
    plt.close(fig)


def main():
    gf.gpdk.PDK.activate()
    b = load(ROOT / '2611_tapered_yue_dual_resist.py', 'unchanged_design1')
    design1 = b.make_device(800.0)
    old_d1 = snapshot(design1)
    m = load(ROOT / '2611_tapered_yue_design2.py', 'layout_design2')
    device = m.make_device(800.0)
    regions = snapshot(device)
    dbu = float(device.kcl.dbu)
    check('one_800um_design2', tuple(m.MEANDER_STRAIGHT_LENGTHS) == (800.0,)
          and device.info['design'] == 'Design 2')
    check('requested_converter_parameters', m.SLAB_TAPER_LENGTH == 100.0
          and m.RIDGE_TAPER_LENGTH == 40.0 and m.RIDGE_EXTENSION_START_WIDTH == 0.6
          and m.RIDGE_CONVERTER_LENGTH == 60.0 and m.RIDGE_CONVERTER_END_WIDTH == 0.2
          and m.RIDGE_CONVERTER_Y_SHIFT == 0.840,
          layer10_final_taper_um=m.SLAB_TAPER_LENGTH,
          ridge_widths_um=[1.0,0.6,0.2], stage_lengths_um=[40.0,60.0],
          centerline_shift_um=0.840)
    delta = compare(old_d1, snapshot(design1), dbu)
    check('design2_does_not_mutate_design1_cells', all(a == 0 for a in delta.values()), xor_um2=delta)
    delta = compare(old_d1, snapshot(b.make_device(800.0)), dbu)
    check('design1_rebuild_unchanged_after_design2', all(a == 0 for a in delta.values()), xor_um2=delta)
    check('required_layers', {(2,0),(10,0),(11,0),(13,0),(14,0),(15,0),(16,0),(51,0)} <= set(regions),
          layers=[list(layer) for layer in sorted(regions)])
    check('no_construction_layer_exported', (999,0) not in regions)

    base_info = json.loads((OUT / 'baseline.json').read_text())
    ref_layout, ref = read_gds(OUT / 'design1_with_100um_outer_taper.gds')
    label_roi = box_region(device.info['label_bbox_um'], dbu)
    label_roi += box_region(base_info['label_bbox_um'], dbu)
    delta = compare(ref, regions, dbu, exclude=label_roi)
    for layer in (2,10,11,13,14,51):
        area = delta[f'{layer}/0']
        check(f'layer_{layer}_inherited_from_100um_slab_reference', area <= TOL,
              difference_um2=area, excludes='label text only' if layer==2 else 'nothing')
    ticks = m.alignment_tip_region(device.info['tip_centers'], device.info['connection_orientations'], dbu)
    optical = regions[(15,0)] - ticks
    check('ridge_and_both_converters_connected', len(list(optical.merged().each())) == 1,
          connected_components=len(list(optical.merged().each())))
    check('all_four_alignment_tips_retained_on_layer15', (ticks-regions[(15,0)]).is_empty()
          and (ticks & regions[(10,0)]).is_empty(), count=len(list(ticks.each())))
    check('layer16_protection_connected', len(list(regions[(16,0)].merged().each())) == 1)
    check('layer16_covers_offset_ridge', (optical-regions[(16,0)]).area()*dbu**2 <= TOL,
          uncovered_um2=float((optical-regions[(16,0)]).area()*dbu**2))
    check('layer14_covers_complete_offset_ridge', (optical-regions[(14,0)]).area()*dbu**2 <= TOL,
          uncovered_um2=float((optical-regions[(14,0)]).area()*dbu**2))
    missing_l2 = (regions[(2,0)]-regions[(14,0)])-label_roi
    check('layer14_still_covers_layer2_windows', missing_l2.area()*dbu**2 <= TOL,
          missing_um2=float(missing_l2.area()*dbu**2))

    allowed = box_region((194.999,-203,295.001,-197), dbu)
    allowed += box_region((-295.001,197,-194.999,203), dbu)
    ridge_delta = (ref[(15,0)] ^ regions[(15,0)])-allowed
    check('meander_and_ridge_outside_converters_unchanged', ridge_delta.area()*dbu**2 <= TOL,
          unexpected_difference_um2=float(ridge_delta.area()*dbu**2))
    right_end = box_region((150.001,-1000,10000,-25), dbu)
    for layer in (10,11,13,14,15,16,51):
        reg = regions[(layer,0)]
        delta_r = (reg ^ reg.transformed(kdb.Trans(kdb.Trans.R180))) & right_end
        check(f'layer_{layer}_both_end_assemblies_rotationally_symmetric', delta_r.area()*dbu**2 <= TOL,
              rotated_xor_um2=float(delta_r.area()*dbu**2))

    coupler = m.design2_two_layer_tapered_fiber(add_support_tips=False,
                add_phononic_taper_support=True)
    cp = snapshot(coupler)
    section_data = []
    for x, expected_center, expected_width in ((45.001,0.0,1.0),(85.001,0.0,0.6),(144.999,0.840,0.2)):
        xi = int(round(x/dbu))
        section = cp[(15,0)] & kdb.Region(kdb.Box(xi-1,-3000,xi+1,3000))
        bb = section.bbox()
        center = (bb.bottom+bb.top)*dbu/2
        width = bb.height()*dbu
        ok = abs(center-expected_center) <= 0.002 and abs(width-expected_width) <= 0.002
        check(f'ridge_section_x_{x:g}', ok, x_um=x, center_y_um=center, vertical_width_um=width)
        section_data.append({'x_um':x,'center_y_um':center,'width_um':width})
    xend = m.PNC_TAPER_LENGTH+m.RIDGE_TAPER_LENGTH+m.RIDGE_CONVERTER_LENGTH
    xi = int(round(xend/dbu))
    section10 = cp[(10,0)] & kdb.Region(kdb.Box(xi-1,-3000,xi,3000))
    section15 = cp[(15,0)] & kdb.Region(kdb.Box(xi-1,-3000,xi,3000))
    gap = (section15.bbox().bottom-section10.bbox().top)*dbu
    check('terminal_layers_are_laterally_separated', gap > 0,
          terminal_centerline_offset_um=m.RIDGE_CONVERTER_Y_SHIFT,
          terminal_same_x_vertical_edge_gap_um=float(gap))
    REPORT['terminal_vertical_edge_gap_um'] = float(gap)
    REPORT['ridge_sections'] = section_data

    top = m.build_layout()
    path = ROOT / 'gds/2611_tapered_yue_design2.gds'
    export(top, path)
    final, _ = read_gds(path)
    names = list(top.info['device_cells'])
    check('export_contains_one_device', len(names)==1 and top.info['device_count']==1)
    cell = final.cell(names[0])
    reopened = file_snapshot(final, cell)
    delta = compare(regions, reopened, dbu)
    check('export_preserves_every_layer', all(a<=TOL for a in delta.values()), xor_um2=delta)
    insts = [inst for inst in final.top_cell().each_inst() if inst.cell.name==names[0]]
    check('exactly_one_device_instance', len(insts)==1)
    if len(insts)==1:
        tr=insts[0].trans
        edge=reopened[(11,0)].bbox().left*dbu+tr.disp.x*dbu
        check('left_bend_alignment_preserved', abs(edge+900.0)<1e-9, left_pnc_edge_um=edge)
        REPORT['device_origin_um']=[tr.disp.x*dbu,tr.disp.y*dbu]
    for rel, digest in base_info['files'].items():
        check(f'original_design1_file_unchanged_{Path(rel).name}', fingerprint(ROOT/rel)==digest)
    render(reopened, dbu, (-1600,-375,1600,375), (2,10,11,13,14,15,16),
           'Design 2: one 800 um meander, complete two-ended device', 'device_overview.png')
    render(cp, dbu, (35,-1.6,190,1.6), (10,15),
           'Design 2 coupler: Layer 10 final taper 100 um; Layer 15 1.0 -> 0.6 -> 0.2 um',
           'two_layer_coupler.png', equal=False)
    render(cp, dbu, (75,-0.75,155,1.25), (10,15),
           '60 um converter: 600 -> 200 nm, final centerline offset 840 nm',
           'converter_detail.png', equal=False)
    render(reopened, dbu, (-310,198.5,-180,201.5), (10,15,16),
           'Design 2 opposite end: same converter under 180 degree rotation',
           'opposite_converter.png', equal=False)
    REPORT['status']='PASS' if not REPORT['failures'] else 'FAIL'
    REPORT['output_file']=str(path.relative_to(ROOT))
    (OUT/'validation_report.json').write_text(json.dumps(REPORT,indent=2))
    (ROOT/'gds/2611_design2_validation.json').write_text(json.dumps(REPORT,indent=2))
    if REPORT['failures']:
        raise AssertionError(', '.join(REPORT['failures']))
    print('DESIGN2_GEOMETRY_VALIDATION_PASS', flush=True)


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='baseline':
        baseline()
    else:
        main()
