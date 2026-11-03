"""Keep original length text off all mask layers before it can merge with them."""
from pathlib import Path
import ast
ROOT=Path(__file__).resolve().parent
p=ROOT/'2611_build_combined_sweep.py'
s=p.read_text()
anchor="    ast.fix_missing_locations(tree)\n    module = types.ModuleType"
replacement='''    # Isolate ONLY the final standalone length label. It is drawn after all
    # exposure/protection masks. Original text can overlap PNC windows and
    # flatten() merges overlapping polygons, so post-hoc glyph subtraction
    # cannot safely distinguish text from the underlying physical opening.
    # Give that one text call a temporary layer before construction, then
    # remove only the temporary layer from the independent finished cell.
    label_calls = 0
    for function in tree.body:
        if isinstance(function, ast.FunctionDef) and function.name == 'make_device':
            for node in ast.walk(function):
                if isinstance(node, ast.Call) and ast.unparse(node.func) == 'gf.components.text':
                    layer_kw = next((k for k in node.keywords if k.arg == 'layer'), None)
                    require(layer_kw is not None and ast.unparse(layer_kw.value) == 'MARKER_LAYER',
                            'Unexpected standalone label construction')
                    layer_kw.value = ast.parse('(998, 0)', mode='eval').body
                    label_calls += 1
    require(label_calls == 1, 'Expected exactly one standalone length-label call')
    ast.fix_missing_locations(tree)
    module = types.ModuleType'''
assert s.count(anchor)==1
s=s.replace(anchor,replacement)
a=s.index('def remove_original_length_label(')
b=s.index('# %% [03]',a)
s=s[:a]+'''def remove_original_length_label(component, module):
    """Remove only temporary length text, never any fabrication-mask polygon.

    configured_module places the original final length label on 998/0.
    This layer is isolated before flattening can merge text into a release
    window. Standalone source files still draw their usual labels unchanged.
    """
    import gdsfactory as gf
    import klayout.db as kdb
    component.flatten()
    layer = gf.get_layer((998, 0))
    text = kdb.Region(component.kdb_cell.begin_shapes_rec(layer)).merged()
    require(not text.is_empty(), 'Isolated original length label is missing')
    count = int(text.count())
    area = float(text.area()) * float(component.kcl.dbu)**2
    component.remove_layers([(998, 0)], recursive=False)
    require(kdb.Region(component.kdb_cell.begin_shapes_rec(layer)).is_empty(),
            'Temporary text layer was not cleared')
    component.info['length_label_removed_for_combined_sweep'] = True
    return dict(removed_polygon_count=count, removed_area_um2=area,
                method='isolate_final_length_text_on_998_before_flatten_then_remove_998')


'''+s[b:]
s=s.replace('All geometry functions and expressions are left untouched.',
            'All physical geometry functions and expressions are left untouched.\n    Only the final standalone text call is isolated on temporary layer 998.')
s=s.replace('length text by exact glyph matching; the old extra sweep ID is gone.',
            'length text via its isolated temporary layer; the old sweep ID is gone.')
s=s.replace('Original length-label glyphs are erased individually, never subtracted from mask geometry.',
            'Original length text is isolated on temporary 998/0 before flattening, then removed. Fabrication masks are never subtracted.')
ast.parse(s);p.write_text(s)

p=ROOT/'validate_2611_vertical_sweep.py';s=p.read_text()
a=s.index('def baseline_without_separate_length_glyphs(')
b=s.index('\n\ndef run():',a)
s=s[:a]+'''def original_length_text_region(entry, dbu):
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
'''+s[b:]
s=s.replace("            baseline = baseline_without_separate_length_glyphs(old_layout, old_devices[ident], old_entries[ident])",
'''            baseline = regions(old_layout, old_devices[ident])
            original_text = original_length_text_region(old_entries[ident], dbu)
            restored = dict(local)
            restored[(2,0)] = (local[(2,0)] + original_text).merged()
            outside = ((baseline[(2,0)] ^ local[(2,0)]) - original_text).area()*dbu**2
            check(ident+'_no_changes_outside_original_length_text', outside < 0.005,
                  outside_old_text_xor_um2=float(outside))''')
s=s.replace("local.get((n,d), kdb.Region())).area()*dbu**2)",
            "restored.get((n,d), kdb.Region())).area()*dbu**2)")
s=s.replace("'_every_nonlabel_polygon_preserved'", "'_approved_geometry_recovered_by_adding_only_old_text'")
s=s.replace('Only the old standalone length-label polygons may differ inside device cells.',
            'Original text may overlap windows: compare old = new UNION old text, never cut out a text box.')
s=s.replace("    check('one_top_cell', len(ly.top_cells()) == 1)",
'''    check('one_top_cell', len(ly.top_cells()) == 1)
    check('no_temporary_text_layer_exported', all(ly.get_info(i).layer != 998
          or kdb.Region(top.begin_shapes_rec(i)).is_empty() for i in ly.layer_indexes()))''')
ast.parse(s);p.write_text(s)

p=ROOT/'COMBINED_SWEEP_README.md';s=p.read_text()
s=s.replace('Original length glyphs are removed individually\nby exact polygon matching, never by cutting a rectangle into a mask.',
'''The final standalone length-label call is isolated on temporary layer 998
in the in-memory source before geometry is built. It is then removed from
the independent device. This prevents old text from merging into release
windows; neither glyph areas nor text bounding boxes are cut out of masks.''')
s=s.replace('it also compares every local device layer against that baseline with only\nthe individually identified original text glyphs removed. This proves that\nstacking and label cleanup did not alter the optical or fabrication shapes.',
'''it also compares every local device layer against that baseline. For Layer 2,
new_mask UNION original_text must recover the approved old mask; this is
important where the old text overlapped real release windows. Every other
layer must match directly. Physical source-function bodies are unchanged;
only the single final length-label layer is redirected in the memory copy.''')
p.write_text(s)
print('Safe pre-flatten label isolation and overlap-aware baseline comparison applied.')
