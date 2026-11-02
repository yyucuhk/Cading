"""Validate equal three-straight meanders in Design 1 and Design 2.

python validate_2611_equal_straights.py capture  # BEFORE a geometry migration
python validate_2611_equal_straights.py          # validate and export both GDSs

Optional baseline comparisons use entire end assemblies after port translation,
not a handful of cross-sections. The actual path sent to the PNC builder is
captured, compared to S(L)-Euler180-S(L)-Euler180-S(L), and measured directly.
"""
from __future__ import annotations
from pathlib import Path
from functools import partial
import importlib.util
import json
import sys
import traceback

import gdsfactory as gf
import klayout.db as kdb
import numpy as np
from validate_2611_sweep import snapshot, file_snapshot, compare, export, box_region

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'validation_2611_equal_straights'
BASE = OUT / 'baseline'
OUT.mkdir(exist_ok=True)
MODELS = {
    'design1': ('2611_tapered_yue_dual_resist.py', 'phononic_meander_shield'),
    'design2': ('2611_tapered_yue_design2.py', 'design2_phononic_meander_shield'),
}
REPORT = {'checks': [], 'failures': [], 'lengths': []}
TOL = 0.005  # um^2: GDS polygon fracture/grid tolerance, not a physical gap budget


def check(name, condition, **details):
    row = {'name': name, 'passed': bool(condition), **details}
    REPORT['checks'].append(row)
    if not condition:
        REPORT['failures'].append(name)
    (OUT / 'validation_report.json').write_text(json.dumps(REPORT, indent=2))
    print(json.dumps(row), flush=True)


def load(key, filename):
    spec = importlib.util.spec_from_file_location('equal_' + key, ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    gf.gpdk.PDK.activate()
    return module


def read(path):
    ly = kdb.Layout()
    ly.read(str(path))
    return ly, file_snapshot(ly, ly.top_cell())


def metadata(device):
    keys = ('connection_centers', 'connection_orientations', 'tip_centers',
            'label_bbox_um', 'total_meander_length_um', 'meander_straight_length_um')
    return {key: device.info[key] for key in keys}


def instances(ly):
    idx = ly.find_layer(51, 0)
    return sorted([inst for inst in ly.top_cell().each_inst()
                   if idx is not None and not kdb.Region(inst.cell.begin_shapes_rec(idx)).is_empty()],
                  key=lambda inst: inst.trans.disp.y)


def chip_regions(ly):
    device_indices = {i.cell_index for i in instances(ly)}
    result = {}
    for idx in ly.layer_indices():
        region = kdb.Region(ly.top_cell().shapes(idx))
        for inst in ly.top_cell().each_inst():
            if inst.cell_index not in device_indices:
                region += kdb.Region(inst.cell.begin_shapes_rec(idx)).transformed(inst.trans)
        info = ly.get_info(idx)
        result[(info.layer, info.datatype)] = region.merged()
    return result


def label_region(module, key, info, dbu):
    L = float(info['meander_straight_length_um'])
    total = float(info['total_meander_length_um'])
    text = (f'DESIGN 2_L{L:g}_ML{total:.1f}' if key == 'design2' else
            f'ML{total:.1f}_TL{module.POST_PNC_TAPER_LENGTH + module.SLAB_TAPER_LENGTH:g}')
    label = gf.components.text(text=text, size=module.DETAIL_LABEL_SIZE, layer=module.MARKER_LAYER)
    tip = min(info['tip_centers'], key=lambda p: p[0])
    shift = kdb.Trans(round((tip[0] + module.DETAIL_LABEL_TIP_X_OFFSET) / dbu),
                      round((tip[1] + module.DETAIL_LABEL_TIP_Y_OFFSET) / dbu))
    return snapshot(label)[(2, 0)].transformed(shift)


def expected_path(module, L, centers, orientations):
    turn = 1 if centers[1, 1] > centers[0, 1] else -1
    p = gf.Path()
    p.append(gf.path.straight(length=L))
    p.append(gf.path.euler(radius=module.S_BEND_COMPACT_180_RADIUS,
                          angle=turn*180, p=module.S_BEND_EULER_P, use_eff=True))
    p.append(gf.path.straight(length=L))
    p.append(gf.path.euler(radius=module.S_BEND_COMPACT_180_RADIUS,
                          angle=-turn*180, p=module.S_BEND_EULER_P, use_eff=True))
    p.append(gf.path.straight(length=L))
    p.rotate(float(orientations[0]))
    p.move(centers[0])
    return p


def preview(path, module, key, L):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ly, r = read(path)
    fig, ax = plt.subplots(figsize=(13, 6))
    bounds = (-L/2-500, -355, L/2+500, 355)
    roi = box_region(bounds, ly.dbu)
    for layer in (2, 10, 11, 14, 15, 16):
        points = []
        for poly in (r.get((layer, 0), kdb.Region()) & roi).merged().each():
            loops = [list(poly.each_point_hull())]
            loops += [list(poly.each_point_hole(i)) for i in range(poly.holes())]
            for loop in loops:
                pp = [(q.x*ly.dbu, q.y*ly.dbu) for q in loop]
                if pp:
                    points += pp + [pp[0], (np.nan, np.nan)]
        if points:
            pp = np.asarray(points)
            ax.plot(pp[:,0], pp[:,1], linewidth=0.6, label=f'Layer {layer}')
    for y in (-200, 0, 200):
        ax.text(0, y+18, f'{L:g} um straight', ha='center', fontsize=10)
    ax.set(xlim=(bounds[0], bounds[2]), ylim=(bounds[1], bounds[3]),
           xlabel='Device-local X (um)', ylabel='Device-local Y (um)',
           title=f'{key.title()}: three equal {L:g} um straights; intact all-layer ends')
    ax.set_aspect('equal', adjustable='box')
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1))
    fig.tight_layout()
    fig.savefig(OUT / f'{key}_equal_straights.png', dpi=180)
    plt.close(fig)


def run(capture=False):
    if capture:
        BASE.mkdir(exist_ok=True)
    baseline_used = True
    for key, (filename, shield_name) in MODELS.items():
        module = load(key, filename)
        lengths = tuple(dict.fromkeys([*map(float, module.MEANDER_STRAIGHT_LENGTHS), 200., 800., 1400.]))
        devices, saved_regions = [], []
        for L in lengths:
            tag = f'{key}_L{L:g}'
            if capture:
                d = module.make_device(L)
                export(d, BASE / f'{tag}.gds')
                (BASE / f'{tag}.json').write_text(json.dumps(metadata(d), indent=2))
                continue
            # Capture the actual uncut path delivered by tapered_fiber_inout.
            captured = []
            original_shield = getattr(module, shield_name)
            def recording_shield(path, **kwargs):
                captured.append(np.asarray(path.points).copy())
                return original_shield(path=path, **kwargs)
            recording_shield.__name__ = f'{key}_capture_L{L:g}'
            setattr(module, shield_name, recording_shield)
            try:
                d = module.make_device(L)
            finally:
                setattr(module, shield_name, original_shield)
            r = snapshot(d)
            devices.append(d)
            saved_regions.append(r)
            dbu = float(d.kcl.dbu)
            info = metadata(d)
            centers = np.asarray(info['connection_centers'], dtype=float)
            orientations = info['connection_orientations']
            check(tag + '_actual_path_captured', len(captured) == 1, paths=len(captured))
            if len(captured) != 1:
                raise AssertionError('Could not independently measure the native generated path')
            actual = captured[0]
            segments = np.diff(actual, axis=0)
            long_horizontal = np.abs(segments[np.abs(segments[:,1]) < 1e-9,0])
            measured = long_horizontal[long_horizontal > L/4].tolist()
            check(tag + '_three_actual_straights_equal', len(measured) == 3 and
                  np.allclose(measured, [L]*3, rtol=0, atol=1e-8), measured_straights_um=measured)
            check(tag + '_ports_follow_equal_length',
                  np.allclose(centers, [[-L/2, module.S_BEND_VERTICAL_HALF_SPAN],
                                        [L/2, -module.S_BEND_VERTICAL_HALF_SPAN]], atol=1e-8, rtol=0),
                  connection_centers_um=centers.tolist())
            p = expected_path(module, L, centers, orientations)
            check(tag + '_entire_path_equals_three_straights_two_bends',
                  actual.shape == p.points.shape and np.allclose(actual, p.points, atol=1e-8, rtol=0))
            check(tag + '_path_ends_at_both_ports',
                  np.allclose(actual[[0,-1]], centers, rtol=0, atol=1e-8))
            bend = float(gf.path.euler(radius=module.S_BEND_COMPACT_180_RADIUS, angle=180,
                                      p=module.S_BEND_EULER_P, use_eff=True).length())
            total = 3*L + 2*bend
            measured_total = float(np.linalg.norm(segments, axis=1).sum())
            check(tag + '_total_length_formula_and_metadata',
                  abs(float(info['total_meander_length_um'])-total)<1e-8 and
                  abs(measured_total-total)<0.003 and
                  list(d.info['meander_straight_segment_lengths_um']) == [L]*3,
                  bend180_um=bend, formula_total_um=total, actual_path_total_um=measured_total)
            REPORT['lengths'].append(dict(design=key, common_straight_um=L,
                straights_um=measured, each_euler_bend_um=bend, total_meander_um=total))
            exact_label = label_region(module, key, info, dbu)
            check(tag + '_correct_length_label_present', (exact_label-r[(2,0)]).is_empty())
            check(tag + '_all_layers_retained', set(r) == {(2,0),(10,0),(11,0),(13,0),(14,0),(15,0),(16,0),(51,0)})
            tips = module.alignment_tip_region(info['tip_centers'], orientations, dbu)
            optical = (r[(15,0)]-tips).merged()
            check(tag + '_four_tips_on_layer15_only', len(list(tips.each()))==4 and
                  (tips-r[(15,0)]).is_empty() and (tips&r[(10,0)]).is_empty())
            check(tag + '_ridge_and_both_ends_connected', len(list(optical.each()))==1)
            check(tag + '_protection_connected', len(list(r[(16,0)].merged().each()))==1)
            for layer in (14,16):
                missing = float((optical-r[(layer,0)]).area()*dbu**2)
                check(tag+f'_layer{layer}_covers_complete_optical_ridge', missing<=TOL, missing_um2=missing)
            check(tag + '_layer14_covers_layer2',
                  ((r[(2,0)]-r[(14,0)])-exact_label).area()*dbu**2<=TOL)
            support = module.phononic_tether(width=module.REVIEW_TETHER_CONNECTION_WIDTH,
                segments=module.PHONONIC_SUPPORT_SEGMENTS, gaps=module.PHONONIC_SUPPORT_GAPS,
                overshoot=module.PHONONIC_SUPPORT_OVERSHOOT)
            span = float(support.ports['a2'].center[0]-support.ports['a1'].center[0])
            for layer, width in ((15,module.RIDGE_MEANDER_WIDTH),(14,module.S_BEND_WIDTH+2*span)):
                expected = gf.path.extrude(p, cross_section=partial(gf.cross_section.cross_section,
                                                               width=width, layer=(layer,0)))
                missing = float((snapshot(expected)[(layer,0)]-r[(layer,0)]).area()*dbu**2)
                check(tag+f'_layer{layer}_whole_meander_and_all_junctions_covered', missing<=TOL, missing_um2=missing)
            before = BASE / f'{tag}.gds'
            if before.exists():
                old_ly, old = read(before)
                old_info = json.loads((BASE / f'{tag}.json').read_text())
                old_text = label_region(module,key,old_info,dbu)
                for side in range(2):
                    delta = centers[side]-np.asarray(old_info['connection_centers'][side])
                    trans = kdb.Trans(round(delta[0]/dbu),round(delta[1]/dbu))
                    roi = box_region((-10000,25,-L/2,1000) if side==0 else
                                     (L/2,-1000,10000,-25),dbu)
                    errors = {}
                    for layer in sorted(set(old)|set(r)):
                        diff = (old.get(layer,kdb.Region()).transformed(trans)^r.get(layer,kdb.Region())) & roi
                        if layer == (2,0):
                            diff -= old_text.transformed(trans)+exact_label
                        errors[f'{layer[0]}/{layer[1]}'] = float(diff.area()*dbu**2)
                    check(tag+f'_complete_end_{side+1}_every_layer_unchanged_after_translation',
                          all(v<=TOL for v in errors.values()), translation_um=delta.tolist(), xor_um2=errors)
                check(tag+'_left_bend_reference_unchanged',
                      old[(11,0)].bbox().left==r[(11,0)].bbox().left)
            else:
                baseline_used = False
            export(d, OUT / f'{tag}.gds')
            ly, reopened = read(OUT / f'{tag}.gds')
            errors = compare(r,reopened,dbu)
            check(tag+'_gds_preserves_every_layer',all(v<=TOL for v in errors.values()),xor_um2=errors)
            for i, previous in enumerate(devices):
                errors = compare(saved_regions[i],snapshot(previous),dbu)
                check(tag+f'_does_not_mutate_device_{i+1}',all(v==0 for v in errors.values()))
            if L == float(module.MEANDER_STRAIGHT_LENGTHS[0]):
                preview(OUT/f'{tag}.gds',module,key,L)
        if capture:
            export(module.build_layout(),BASE/f'{key}_layout.gds')
            continue
        for L, saved in zip(lengths,saved_regions,strict=True):
            repeated = module.make_device(L)
            check(f'{key}_L{L:g}_build_order_independent',
                  all(v==0 for v in compare(saved,snapshot(repeated),dbu).values()))
        for bad in (0.,-1.,float('nan'),float('inf')):
            try:
                module.make_device(bad)
                ok=False
            except ValueError:
                ok=True
            check(key+'_rejects_'+str(bad),ok)
        top = module.build_layout()
        outpath = ROOT/'gds'/module.OUTPUT_GDS_FILENAME
        export(top,outpath)
        ly,_ = read(outpath)
        placed = instances(ly)
        check(key+'_exported_device_count',len(placed)==len(module.MEANDER_STRAIGHT_LENGTHS))
        positions=[]
        for i, inst in enumerate(placed):
            rr=file_snapshot(ly,inst.cell)
            position=[inst.trans.disp.x*ly.dbu,inst.trans.disp.y*ly.dbu]
            positions.append(position)
            check(key+f'_placed_device_{i+1}_left_bend_alignment',
                  abs(position[0]+rr[(11,0)].bbox().left*ly.dbu-module.LEFT_BEND_ALIGNMENT_X)<1e-8)
            check(key+f'_placed_device_{i+1}_y_unchanged',
                  abs(position[1]-(module.BOTTOM_DEVICE_Y+i*module.DEVICE_Y_PITCH))<1e-8)
        before=BASE/f'{key}_layout.gds'
        if before.exists():
            oldly,_=read(before)
            errors=compare(chip_regions(oldly),chip_regions(ly),ly.dbu)
            check(key+'_chip_and_global_markers_unchanged',all(v==0 for v in errors.values()),xor_um2=errors)
        saved_name,saved_lengths=module.TOP_CELL_NAME,module.MEANDER_STRAIGHT_LENGTHS
        for sweep in ((800.,),(200.,800.,1400.),(200.,500.,800.,1100.,1400.)):
            module.TOP_CELL_NAME=saved_name+f'_test_{len(sweep)}'
            module.MEANDER_STRAIGHT_LENGTHS=sweep
            trial=module.build_layout()
            check(key+f'_supports_{len(sweep)}_devices_with_600um_pitch',
                  trial.info['device_count']==len(sweep) and
                  np.allclose(trial.info['device_y_positions'],
                    module.BOTTOM_DEVICE_Y+np.arange(len(sweep))*module.DEVICE_Y_PITCH,rtol=0,atol=1e-8))
        module.TOP_CELL_NAME,module.MEANDER_STRAIGHT_LENGTHS=saved_name,saved_lengths
    if capture:
        print('EQUAL_STRAIGHT_BASELINE_CAPTURED',flush=True)
        return
    REPORT['baseline_comparison_performed']=baseline_used
    REPORT['status']='PASS' if not REPORT['failures'] else 'FAIL'
    (OUT/'validation_report.json').write_text(json.dumps(REPORT,indent=2))
    if REPORT['failures']:
        raise AssertionError(', '.join(REPORT['failures']))
    (ROOT/'gds/2611_equal_straights_validation.json').write_text(json.dumps(REPORT,indent=2))
    print('EQUAL_STRAIGHTS_VALIDATION_PASS',flush=True)


if __name__=='__main__':
    try:
        run(capture=len(sys.argv)>1 and sys.argv[1]=='capture')
    except Exception:
        REPORT['status']='FAIL'
        REPORT['exception']=traceback.format_exc()
        (OUT/'validation_report.json').write_text(json.dumps(REPORT,indent=2))
        raise
