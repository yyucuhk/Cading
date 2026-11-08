"""Guarded migration for the six-device Design-1 comparison only."""
from pathlib import Path
import ast
import hashlib

root = Path.cwd()
p = root / '2611_design1_L6000_straight_output.py'
s = p.read_text()
assert hashlib.sha1(b'blob ' + str(len(p.read_bytes())).encode() + b'\0' + p.read_bytes()).hexdigest() == 'ddd5d5b9a6febc3483959af9fdde7d6cdf8f47d0', 'Builder changed; reconcile before patching'

def replace(old, new):
    global s
    assert s.count(old) == 1, repr(old[:110])
    s = s.replace(old, new)

replace('Run from the Cading root. Existing source/GDS files are never overwritten.',
        'Run from the Cading root. Only this comparison output is overwritten.\nOther saved designs are not modified. Final output omits Layers 16 and 51;\nonly global alignment markers move from Layer 2 to Layer 15. Device Layer-2\nrelease windows and the six text labels stay on Layer 2.')
replace('FAB_LAYERS = (2, 10, 11, 13, 14, 15, 16)',
'''FAB_LAYERS = (2, 10, 11, 13, 14, 15, 16)  # pre-export construction checks
REMOVED_OUTPUT_LAYERS = (16, 51)  # remove all datatypes from every output cell
ALIGNMENT_MARKER_SOURCE_LAYER = (2, 0)
ALIGNMENT_MARKER_TARGET_LAYER = (15, 0)  # global alignment markers ONLY''')
anchor = '# %% [02] Complete all-layer straight-output factory'
helper = '''def finalize_output_layers(layout, frame_cell):
    """Apply only the requested layer changes to the independently owned layout.

    The frame cell contains chip outline and alignment marks, never device
    windows or labels. Delete excluded layers across every cell, including
    unreferenced cells. Never remap all of Layer 2: it also carries release
    windows and device labels. Original source/baseline layouts are read-only.
    """
    import klayout.db as kdb
    require(not list(frame_cell.each_inst()), 'Expected an owned flat frame cell')
    source = layout.layer(*ALIGNMENT_MARKER_SOURCE_LAYER)
    target = layout.layer(*ALIGNMENT_MARKER_TARGET_LAYER)
    moved = kdb.Region(frame_cell.shapes(source)).merged()
    require(not moved.is_empty(), 'No source alignment markers were found')
    existing = kdb.Region(frame_cell.shapes(target)).merged()
    frame_cell.shapes(target).clear()
    frame_cell.shapes(target).insert((existing | moved).merged())
    frame_cell.shapes(source).clear()
    removed_shapes = {}
    for index in tuple(layout.layer_indexes()):
        info = layout.get_info(index)
        if info.layer in REMOVED_OUTPUT_LAYERS:
            key = str(info.layer) + '/' + str(info.datatype)
            removed_shapes[key] = sum(cell.shapes(index).size() for cell in layout.each_cell())
            layout.delete_layer(index)
    return dict(removed_layers=list(REMOVED_OUTPUT_LAYERS),
                removed_shape_counts=removed_shapes,
                marker_source_layer=list(ALIGNMENT_MARKER_SOURCE_LAYER),
                marker_target_layer=list(ALIGNMENT_MARKER_TARGET_LAYER),
                marker_cell=frame_cell.name,
                moved_marker_polygon_count=moved.count(),
                moved_marker_area_um2=moved.area() * layout.dbu ** 2,
                layer2_windows_and_labels_preserved=True)


'''
replace(anchor, helper + anchor)
replace('    write_layout(ly,OUTPUT)\n    final,top=read_gds(OUTPUT)',
'''    # Apply the layer policy only after placement and frame size are fixed.
    # Internal construction checks above still inspect the approved reference
    # masks; all checks below inspect the final deliverable without L16/L51.
    layer_policy = finalize_output_layers(ly, fc)
    final_frame = regions(ly, fc)
    check(checks, 'global_markers_moved_2_to_15_exactly',
          (final_frame[ALIGNMENT_MARKER_TARGET_LAYER] ^ fr[ALIGNMENT_MARKER_SOURCE_LAYER]).is_empty())
    check(checks, 'no_marker_copies_left_in_frame_layer2', (2, 0) not in final_frame)
    target_marks = final_frame[ALIGNMENT_MARKER_TARGET_LAYER].transformed(ft)
    check(checks, 'moved_markers_clear_of_devices_and_labels',
          all((target_marks & p[15]).is_empty() for p in placed)
          and all((target_marks & t).is_empty() for t in texts))
    write_layout(ly,OUTPUT)
    final,top=read_gds(OUTPUT)
    exported_layers = {final.get_info(i).layer for i in final.layer_indexes()}
    check(checks, 'layers16_and51_absent_from_entire_export',
          not (exported_layers & set(REMOVED_OUTPUT_LAYERS)))''')
replace("        for layer,r in snapshots[e['case_id']].items():\n            check(checks,e['case_id']+f'_saved_L{layer[0]}_exact',(r^rr[layer]).is_empty())",
'''        expected = {layer: r for layer, r in snapshots[e['case_id']].items()
                    if layer[0] not in REMOVED_OUTPUT_LAYERS}
        check(checks, e['case_id'] + '_exported_layer_set', set(rr) == set(expected))
        for layer,r in expected.items():
            check(checks,e['case_id']+f'_saved_L{layer[0]}_exact',(r^rr[layer]).is_empty())''')
replace('extension_layers=[2,11,14,15,16],frame_bounds_um=[xmin,ymin,xmax,ymax],',
'''extension_layers=[2,11,14,15],frame_bounds_um=[xmin,ymin,xmax,ymax],
                final_layer_policy=layer_policy,
                worker_checks_scope='Unfiltered construction/reference masks; final export checks exclude Layers 16 and 51', ''')
s=s.replace('layers=(2,10,11,14,15,16)', 'layers=(2,10,11,14,15)')
replace('\n\n# %% [06] Entry point: one isolated process per width',
'''    xmin, ymin, _, _ = result['frame_bounds_um']
    view(top, 'alignment_markers_layer15',
         (xmin-120, ymin-120, xmin+120, ymin+120),
         'Alignment marks moved from Layer 2 to Layer 15; Layer 10 unchanged',
         layers=(2,10,15))


# %% [06] Entry point: one isolated process per width''')
ast.parse(s)
p.write_text(s)

# The independent validator now checks the six retained device layers.
p = root/'validate_2611_straight_output_raw.py'
s = p.read_text()
s=s.replace("'_eight_layers'", "'_six_output_layers'")
s=s.replace('{2,10,11,13,14,15,16,51}', '{2,10,11,13,14,15}')
s=s.replace('for l in (2,10,11,13,14,15,16,51):', 'for l in (2,10,11,13,14,15):')
old='''        protection=union(new[v],16)
        check(v+'_protection_connected',protection.geom_type=='Polygon')
        body=max(geoms,key=lambda p:p.area)
        check(v+'_ridge_covered_by_protection',body.difference(protection).area<1e-7)'''
new='''        check(v+'_protection_and_fiber_layers_removed',
              not new[v]['polys'].get(16) and not new[v]['polys'].get(51))
        body=max(geoms,key=lambda p:p.area)'''
assert old in s
s=s.replace(old,new)
old="    refs=new[top]['refs']"
new="""    check('layers16_and51_absent_in_all_cells',
          all(not c['polys'].get(16) and not c['polys'].get(51) for c in new.values()))
    check('global_marker_cell_on15_not2',
          bool(new['FRAME_AND_MARKERS']['polys'].get(15))
          and not new['FRAME_AND_MARKERS']['polys'].get(2))
    refs=new[top]['refs']"""
assert s.count(old)==1
s=s.replace(old,new)
ast.parse(s)
p.write_text(s)

# Current editing notes, not a contradictory historical instruction.
p=root/'DESIGN1_STRAIGHT_OUTPUT_README.txt'
s=p.read_text()
s=s.replace('New output (existing layouts are not overwritten):',
            'Output (this file is overwritten; other saved layouts are unchanged):')
s=s.replace('release windows and fiber guide.', 'release windows. Its Layer-51 fiber guide is omitted from this comparison.')
s=s.replace('with Layer-16 protection. To retain the existing suspended-waveguide support',
            'without Layer-16 protection in the exported file. To retain the existing suspended-waveguide support')
s=s.replace('connected ridge/protection, mask', 'connected ridge and release-mask')
s += '''\nOUTPUT LAYER CLEANUP\n====================\nLayers 16 and 51 are omitted everywhere in the six-device comparison GDS.\nGlobal alignment markers formerly on Layer 2 are moved (not copied) to\nLayer 15 at exactly the same positions, with the same shapes. Layer-10\nmarkers are retained. Device Layer-2 release windows and all six device\nlabels remain unchanged on Layer 2. Alignment tips already on Layer 15\nremain unchanged. No waveguide, PNC, taper, placement, or frame changes.\n\nThe policy is explicit in Section 01 and finalize_output_layers(). It runs\nafter all construction/placement, on an independent KLayout layout, before\nwriting the deliverable. Intermediate worker masks may include Layer 16\nfor baseline checks; they are not the final export.\n\nThis change applies ONLY to the comparison builder/output. The source\nthree-device width sweep, standalone designs and 18-device sweep are not\nmodified. The original baseline GDS must retain its manifest hash.\n'''
p.write_text(s)
print('Patched comparison builder, independent validator and editing notes')
