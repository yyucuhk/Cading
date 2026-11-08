"""One-use guarded migration; removed after successful build and validation."""
from pathlib import Path
import hashlib

p = Path('2611_design1_L6000_straight_output.py')
s = p.read_text()
b = p.read_bytes()
assert hashlib.sha1(f'blob {len(b)}\0'.encode()+b).hexdigest() == '07bef57240245ca8d30ef5fa6b2ad1f57a871224'

def replace(old, new):
    global s
    assert s.count(old) == 1, old[:100]
    s = s.replace(old, new)

replace('INCLUDE_REFERENCE_DEVICES = True     # originals first, new variants above',
        'INCLUDE_REFERENCE_DEVICES = True     # two-coupler set first, straight outputs above\nTWO_COUPLER_END_WINDOWS_ONLY = True   # remove L11; keep L14 outward of each L10 entrance')
helper = '''def finalize_two_coupler_windows(layout, entries):
    """Restrict ONLY the two-coupler cells, after all placements are fixed.

    L10's actual polygon edge locates each entrance. Transverse half-planes
    isolate the complete terminal assemblies from the nearby Euler bends.
    L14 is intersected, never expanded or redrawn; L11 is removed in every
    datatype. Other layers and all straight-output cells are read-only.
    """
    import klayout.db as kdb
    policies = {}
    if not TWO_COUPLER_END_WINDOWS_ONLY:
        return policies
    for e in entries:
        if e['mode'] != 'original_two_tapers':
            continue
        cell = layout.cell(e['case_id'])
        require(not list(cell.each_inst()), 'Expected independently owned flat device cell')
        r10 = kdb.Region(cell.shapes(layout.layer(10, 0))).merged()
        r14 = kdb.Region(cell.shapes(layout.layer(14, 0))).merged()
        require(not r10.is_empty() and not r14.is_empty(), 'Missing taper or release geometry')
        b = r14.bbox() + r10.bbox()
        xmin, ymin, xmax, ymax = (b.left*DBU-1, b.bottom*DBU-1,
                                 b.right*DBU+1, b.top*DBU+1)
        keep = kdb.Region()
        ends = []
        for side, port in (('input', e['input_port_local_um']),
                           ('output', e['output_port_local_um'])):
            px, py = map(float, port)
            require((px < 0 < py) or (py < 0 < px), 'Unexpected terminal orientation')
            split_y = py / 8.0
            search_bounds = ((xmin, split_y, px, ymax) if px < 0 else
                             (px, ymin, xmax, split_y))
            taper = (r10 & box(search_bounds)).merged()
            require(not taper.is_empty(), 'No L10 taper found at '+side)
            start_x = (taper.bbox().right if px < 0 else taper.bbox().left)*DBU
            require((start_x-px)*(-1 if px < 0 else 1) > 0, 'Entrance is not outward of port')
            limits = ((xmin, split_y, start_x, ymax) if px < 0 else
                      (start_x, ymin, xmax, split_y))
            terminal = box(limits)
            require(not (r14 & terminal).is_empty(), 'Empty retained '+side+' window')
            keep |= terminal
            ends.append(dict(side=side, l10_start_local_x_um=start_x,
                             center_y_um=py, keep_box_local_um=list(limits)))
        before_area = r14.area()*DBU**2
        removed_l11 = 0
        for li in tuple(layout.layer_indexes()):
            info = layout.get_info(li)
            if info.layer == 11:
                removed_l11 += cell.shapes(li).size()
                cell.shapes(li).clear()
            elif info.layer == 14:
                original = kdb.Region(cell.shapes(li)).merged()
                retained = (original & keep).merged()
                cell.shapes(li).clear()
                cell.shapes(li).insert(retained)
        after = kdb.Region(cell.shapes(layout.layer(14, 0))).merged()
        policies[e['case_id']] = dict(removed_layer=11, removed_l11_shapes=removed_l11,
            layer14_before_area_um2=before_area, layer14_after_area_um2=after.area()*DBU**2,
            layer14_region='Both existing end windows from the L10 entrance outward', ends=ends)
        e['final_pnc_present'] = False
        e['final_layer14_scope'] = 'two_terminal_windows_only'
    return policies


'''
replace('# %% [02] Complete all-layer straight-output factory', helper+'# %% [02] Complete all-layer straight-output factory')
replace('    layer_policy = finalize_output_layers(ly, fc)',
        '    layer_policy = finalize_output_layers(ly, fc)\n    two_coupler_policy = finalize_two_coupler_windows(ly, entries)')
replace("        check(checks, e['case_id'] + '_exported_layer_set', set(rr) == set(expected))", '''        if e['case_id'] in two_coupler_policy:
            expected = {layer: r for layer, r in expected.items() if layer[0] != 11}
            keep = kdb.Region()
            for end in two_coupler_policy[e['case_id']]['ends']:
                keep |= box(end['keep_box_local_um'])
            expected = {layer: ((r & keep).merged() if layer[0] == 14 else r)
                        for layer, r in expected.items()}
            expected = {layer: r for layer, r in expected.items() if not r.is_empty()}
            check(checks, e['case_id']+'_no_layer11_on_two_coupler_device',
                  all(layer[0] != 11 for layer in rr))
            check(checks, e['case_id']+'_layer14_only_outward_end_windows',
                  (rr[(14,0)] - keep).is_empty())
        check(checks, e['case_id'] + '_exported_layer_set', set(rr) == set(expected))''')
replace("        check(checks,e['case_id']+'_actual_x_alignment',\n              abs((rr[(11,0)].bbox().left+inst.trans.disp.x)*DBU-LEFT_BEND_X)<1e-8)",
        "        check(checks,e['case_id']+'_original_left_bend_reference_placement',\n              abs((snapshots[e['case_id']][(11,0)].bbox().left+inst.trans.disp.x)*DBU-LEFT_BEND_X)<1e-8)")
replace('                final_layer_policy=layer_policy,',
        '                final_layer_policy=layer_policy,\n                two_coupler_layer_policy=two_coupler_policy,')
replace("worker_checks_scope='Unfiltered construction/reference masks; final export checks exclude Layers 16 and 51'",
        "worker_checks_scope='Unfiltered source masks; final checks enforce global removal of L16/L51 and selective two-coupler L11 removal/L14 cropping'")
replace("'Design 1: original devices below; three straight-output variants above'",
        "'Design 1: two-coupler devices without PNC below; unchanged straight-output variants above'")
replace('Other saved designs are not modified. Final output omits Layers 16 and 51;',
        'Other saved designs are not modified. The two-coupler set also omits Layer 11\nand keeps Layer 14 only outward from both Layer 10 entrances. Straight-output\nvariants are unchanged. Final output omits Layers 16 and 51;')
p.write_text(s)

p = Path('validate_2611_straight_output_raw.py')
s = p.read_text()
old = "        check(n+'_six_output_layers',set(new[n]['polys'])=={2,10,11,13,14,15})"
assert old in s
s = s.replace(old, "        check(n+'_five_output_layers',set(new[n]['polys'])=={2,10,13,14,15})")
s = s.replace("            ar=a.symmetric_difference(b).area\n            check(n+f'_reference_L{l}_unchanged',ar<1e-7,xor_um2=ar)\n            ar=b.symmetric_difference(c).difference(roi).area", "            expected = (a.intersection(box(-10000,25,-3045,10000).union(box(3045,-10000,10000,-25))) if l == 14 else a)\n            if l == 11: expected = Polygon()\n            ar=expected.symmetric_difference(b).area\n            check(n+f'_two_coupler_L{l}_matches_policy',ar<1e-7,xor_um2=ar)\n            ar=a.symmetric_difference(c).difference(roi).area")
s = s.replace("x=union(new[c],11).bounds[0]+devrefs[c][0]*.001", "x=union(old[n] if c == n else new[c],11).bounds[0]+devrefs[c][0]*.001")
s = s.replace("check(c+'_left_bend_minus900'", "check(c+'_original_left_bend_reference_minus900'")
p.write_text(s)

p=Path('DESIGN1_STRAIGHT_OUTPUT_README.txt')
s=p.read_text()
s += '''\n\nTWO-COUPLER LAYER POLICY (current)\nOnly the three two-coupler devices in this comparison omit Layer 11.\nLayer 14 is retained only at the two existing terminal assemblies, starting\nat the actual Layer-10 taper entrance (local X=-3045 and +3045 um).\nBoth outward channels and their existing pads/openings are preserved.\nThe meander and preceding 15-to-5 um PNC-transition openings are removed\nfrom Layer 14; no other device layer is edited. All three straight-output\nvariants remain identical to the preceding comparison.\nThe global removal of Layers 16/51 and L2-to-L15 marker move remain active.\nDevice labels are intentionally unchanged; P10 in the two-coupler labels\nis a legacy source identifier, NOT a remaining Layer-11 pattern.\nThe original three-device source GDS remains unchanged because it is the\nshared source for constructing the straight-output variants.\nControl: TWO_COUPLER_END_WINDOWS_ONLY in Section 01 of the builder.\nCurrent exported-mask audit: gds/2611_two_coupler_windows_validation.json.\nThe old layer-cleanup report is superseded by this current report.\n'''
p.write_text(s)

import json
Path('gds/2611_layer_cleanup_validation.json').write_text(json.dumps(dict(
    status='SUPERSEDED', superseded_by='gds/2611_two_coupler_windows_validation.json',
    note='The L16/L51 removal and marker move remain active; current two-coupler cells additionally omit L11 and crop L14.'),indent=2))
print('Selective two-coupler source changes applied')
