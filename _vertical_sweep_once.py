"""Guarded one-time update; the device geometry generators are not edited."""
from pathlib import Path
import ast

ROOT = Path(__file__).resolve().parent
path = ROOT / '2611_build_combined_sweep.py'
s = path.read_text()

def replace(old, new):
    global s
    if s.count(old) != 1:
        raise RuntimeError(f'Expected one migration anchor: {old[:100]!r}; found {s.count(old)}')
    s = s.replace(old, new)

replace('Build two 3x3 CAD sweeps in one GDS, without sharing cached device cells.',
        'Build 18 vertically stacked CAD devices with one parameter label each.')
replace('DESIGN_X_OFFSET = 5000.0  # Design 2 relative to Design 1',
        'DESIGN_ORDER = (1, 2)  # bottom to top: nine Design-1, then nine Design-2 devices\nDESIGN_X_OFFSET = 0.0  # retained in manifest; both designs share the same bend X')
replace('LEFT_BEND_X = -900.0  # leftmost L11 PNC edge in the Design-1 column',
        'LEFT_BEND_X = -900.0  # common leftmost Layer-11 PNC bend edge for ALL devices')
replace('LABEL_SIZE = 5.0',
        'LABEL_SIZE = 5.0\nLABEL_X_OFFSET = 200.0  # from the aligned left PNC bend edge\nLABEL_Y_OFFSET = 250.0  # from the device origin, clear of the optical route\n# Single GDS label format: D1_T80_P2 (all lengths/periods in um).')
replace("    require(STRAIGHT_LENGTH > 0 and ROW_PITCH > 0 and DESIGN_X_OFFSET > 0,\n            'Lengths and spacing must be positive')",
        "    require(STRAIGHT_LENGTH > 0 and ROW_PITCH > 0, 'Lengths and spacing must be positive')\n    require(DESIGN_X_OFFSET == 0 and DESIGN_ORDER == (1, 2),\n            'This layout uses one vertical column: Design 1 below Design 2')")
replace('    for design in (1, 2):', '    for design in DESIGN_ORDER:')
replace("                               straight_um=STRAIGHT_LENGTH,\n                               case_id=", "                               straight_um=STRAIGHT_LENGTH,\n                               global_row=len(result),\n                               case_id=")

helper = '''def sweep_label(case):
    """Exactly the requested three quantities; widths/lengths are in um."""
    return f"D{case['design']}_T{case['taper_total_um']:g}_P{case['bend_period_um']:g}"


def remove_original_length_label(component, module):
    """Erase ONLY individually matched original text polygons on an owned cell.

    Never subtract a text bounding box or a merged text region from a mask:
    a label may overlap an exposure polygon. Deleting the separate original
    glyph shapes preserves any underlying mask area exactly.
    """
    import gdsfactory as gf
    import klayout.db as kdb
    component.flatten()
    total = float(component.info['total_meander_length_um'])
    length = float(component.info['meander_straight_length_um'])
    if getattr(module, 'DESIGN_ID', None) == 'Design 2':
        text = f'DESIGN 2_L{length:g}_ML{total:.1f}'
    else:
        text = f'ML{total:.1f}_TL{module.POST_PNC_TAPER_LENGTH + module.SLAB_TAPER_LENGTH:g}'
    tip = min(component.info['tip_centers'], key=lambda point: float(point[0]))
    expected_component = gf.Component()
    expected_component.add_ref(gf.components.text(
        text=text, size=module.DETAIL_LABEL_SIZE, layer=module.MARKER_LAYER
    )).move((float(tip[0]) + module.DETAIL_LABEL_TIP_X_OFFSET,
             float(tip[1]) + module.DETAIL_LABEL_TIP_Y_OFFSET))
    expected_component.flatten()
    layer_index = gf.get_layer(module.MARKER_LAYER)
    expected = []
    for shape in expected_component.kdb_cell.shapes(layer_index).each():
        if shape.is_polygon() or shape.is_box():
            expected.append(kdb.Polygon(shape.box) if shape.is_box() else shape.polygon)
    require(bool(expected), 'Original length label was unexpectedly empty')
    shapes = component.kdb_cell.shapes(layer_index)
    matched = []
    remaining = list(expected)
    for shape in shapes.each():
        if not (shape.is_polygon() or shape.is_box()):
            continue
        polygon = kdb.Polygon(shape.box) if shape.is_box() else shape.polygon
        for i, target in enumerate(remaining):
            if polygon.bbox() == target.bbox() and (kdb.Region(polygon) ^ kdb.Region(target)).is_empty():
                matched.append(shape)
                remaining.pop(i)
                break
    require(not remaining and len(matched) == len(expected),
            f'Could not identify all separate original text glyphs: {text}')
    for shape in matched:
        shapes.erase(shape)
    component.info['length_label_removed_for_combined_sweep'] = True
    return dict(text=text, removed_polygon_count=len(matched),
                method='erase_exact_separate_glyph_shapes_not_mask_subtraction')


'''
replace('# %% [03] Build and measure one COMPLETE all-layer device (child process)',
        helper + '# %% [03] Build and measure one COMPLETE all-layer device (child process)')
replace("    d = m.make_device(case['straight_um'])\n    d.flatten()", "    d = m.make_device(case['straight_um'])\n    removed_label = remove_original_length_label(d, m)\n    d.flatten()")
replace("    data = dict(case, source=SOURCES[case['design']],",
        "    data = dict(case, source=SOURCES[case['design']],\n                removed_length_label=removed_label,")
replace("    labels = gf.Component('SWEEP_CASE_LABELS')", "    label_instances = []")
replace("        x = LEFT_BEND_X + (case['design']-1)*DESIGN_X_OFFSET - meta['left_pnc_x_um']\n        y = FIRST_ROW_Y + case['row']*ROW_PITCH",
        "        x = LEFT_BEND_X - meta['left_pnc_x_um']\n        y = FIRST_ROW_Y + case['global_row']*ROW_PITCH")
replace("        # Keep the existing length label; add explicit sweep ID above each row.\n        labels.add_ref(gf.components.text(\n            text=f'{ident}_L{STRAIGHT_LENGTH:g}', size=LABEL_SIZE, layer=(2,0)\n        )).move((LEFT_BEND_X+(case['design']-1)*DESIGN_X_OFFSET+200, y+250))",
'''        # ONE independent label per device. The worker has removed its old
        # length text by exact glyph matching; the old extra sweep ID is gone.
        label_text = sweep_label(case)
        label = gf.components.text(text=label_text, size=LABEL_SIZE, layer=(2, 0))
        label_region = kdb.Region(label.kdb_cell.begin_shapes_rec(gf.get_layer((2, 0))))
        label_cell = ly.create_cell('LABEL_' + ident)
        label_cell.shapes(ly.layer(2, 0)).insert(label_region)
        label_trans = kdb.Trans(round((LEFT_BEND_X + LABEL_X_OFFSET)/DBU),
                               round((y + LABEL_Y_OFFSET)/DBU))
        top.insert(kdb.CellInstArray(label_cell.cell_index(), label_trans))
        meta['label'] = label_text
        meta['label_layer'] = [2, 0]
        meta['label_origin_um'] = [label_trans.disp.x*DBU, label_trans.disp.y*DBU]
        meta['label_size_um'] = LABEL_SIZE
        label_instances.append((ident, label_region.transformed(label_trans)))''')
old = '''    for i in range(len(entries)//2):
        a=placed[i][0]; b=placed[i+len(entries)//2][0]
        require(abs(b['origin_um'][0]-a['origin_um'][0]-DESIGN_X_OFFSET)<.0011,
                'Design column displacement is not 5000 um')
        require(b['origin_um'][1]==a['origin_um'][1], 'Design columns have unequal row coordinates')'''
new = '''    for (_, a), (_, b) in zip(placed, placed[1:]):
        pass  # Individual fabrication-layer intersections were tested above.
    origins = [meta['origin_um'][1] for meta, _ in placed]
    require(all(abs(b-a-ROW_PITCH) < 1e-8 for a,b in zip(origins, origins[1:])),
            'The full vertical column must have uniform Y pitch')
    require(all(abs(meta['left_bend_reference_x_um']-LEFT_BEND_X) < 1e-8
                for meta, _ in placed), 'Left bends are not aligned')
    for ident, text_region in label_instances:
        for meta, masks in placed:
            require(all((text_region & reg).is_empty() for reg in masks.values()),
                    f'Label {ident} touches fabrication geometry of {meta["case_id"]}')
        fabrication_bbox += text_region.bbox()
    for i, (name, reg) in enumerate(label_instances):
        for other_name, other in label_instances[i+1:]:
            require((reg & other).is_empty(), f'Labels overlap: {name}, {other_name}')'''
new = new.replace("    for (_, a), (_, b) in zip(placed, placed[1:]):\n        pass  # Individual fabrication-layer intersections were tested above.\n", '')
replace(old, new)
replace("    ymin=min(-MIN_FRAME_HEIGHT/2, math.floor((fabrication_bbox.bottom*DBU-FRAME_MARGIN)/100)*100)\n    ymax=max(MIN_FRAME_HEIGHT/2, math.ceil((fabrication_bbox.top*DBU+FRAME_MARGIN)/100)*100)",
'''    ymin=math.floor((fabrication_bbox.bottom*DBU-FRAME_MARGIN)/100)*100
    ymax=math.ceil((fabrication_bbox.top*DBU+FRAME_MARGIN)/100)*100
    if ymax-ymin < MIN_FRAME_HEIGHT:
        middle=(ymin+ymax)/2
        ymin=math.floor((middle-MIN_FRAME_HEIGHT/2)/100)*100
        ymax=math.ceil((middle+MIN_FRAME_HEIGHT/2)/100)*100''')
replace("    m.add_chip_outline_and_markers(frame, len(entries)//2)\n    frame.add_ref(labels).move((-(xmax+xmin)/2, -(ymax+ymin)/2))",
        "    m.add_chip_outline_and_markers(frame, len(entries))")
replace("'SWEEP_FRAME_AND_LABELS'", "'SWEEP_FRAME_AND_GLOBAL_MARKERS'")
replace("    require(len(actual)==len(entries)==18,'Incorrect exported device count')",
'''    require(len(actual)==len(entries)==18,'Incorrect exported device count')
    actual_labels=[i for i in ex_top.each_inst() if i.cell.name.startswith('LABEL_')]
    require(len(actual_labels)==len(entries), 'There must be exactly one label per device')''')
replace("                  row_pitch_um=ROW_PITCH, design_x_offset_um=DESIGN_X_OFFSET,",
'''                  row_pitch_um=ROW_PITCH, design_x_offset_um=0.0,
                  placement_mode='single_vertical_column_design1_then_design2',
                  first_row_y_um=FIRST_ROW_Y, left_bend_x_um=LEFT_BEND_X,
                  labels_per_device=1, label_format='D{design}_T{taper_um}_P{bend_period_um}',
                  design_y_offset_um=len(entries)//2*ROW_PITCH,''')
replace("                         'Layer 51 is a nonfabrication fiber guide and may cross columns. No fabrication layer is clipped.'])",
        "                         'Layer 51 is a retained nonfabrication fiber guide. No fabrication layer is clipped.',\n                         'Original length-label glyphs are erased individually, never subtracted from mask geometry.'])")
replace("    fields=['case_id','design','row','bend_period_um','taper_total_um','stage_lengths_um',",
        "    fields=['case_id','design','row','global_row','label','bend_period_um','taper_total_um','stage_lengths_um',")
replace("    fig,ax=plt.subplots(figsize=(13,10))", "    fig,ax=plt.subplots(figsize=(10,19))")
replace("        ax.text(x+200,y+255,e['case_id'],fontsize=7)\n    ax.set(xlim=(-3200,8200),ylim=(-1200,4550),xlabel='X (um)',ylabel='Y (um)',\n           title='Design 1 (left) / Design 2 (right): 9 cases each, 5000 um column offset')",
'''        ax.text(x+LABEL_X_OFFSET,y+LABEL_Y_OFFSET+15,e['label'],fontsize=7)
    frame=manifest['frame_bounds_um']
    ax.set(xlim=(frame[0],frame[2]),ylim=(FIRST_ROW_Y-450,
           FIRST_ROW_Y+(len(manifest['devices'])-1)*ROW_PITCH+450),
           xlabel='X (um)',ylabel='Y (um)',
           title='One vertical column: Design 1 (bottom 9), Design 2 (top 9)')''')
ast.parse(s)
path.write_text(s)

p = ROOT/'validate_2611_combined_sweep.py'
v = p.read_text()
v = v.replace("abs(left-(-900+(e['design']-1)*5000))", "abs(left-manifest['left_bend_x_um'])")
a=v.index('        # Full end-assembly symmetry:')
b=v.index('        roi=box(',a)
v=v[:a] + '        # No device-level length text remains; compare the full end masks.\n' + v[b:]
v=v.replace('            if n==2: delta-=text\n','')
v=v.replace("    check('matching_rows_5000um_apart',all(abs(b[1]-a[1]-5000)<1e-8 and b[2]==a[2]\n           for a,b in zip(rows[1],rows[2])))",
"    check('matching_design_rows_same_x_5400um_y_offset',all(abs(b[1]-a[1])<1e-8 and abs(b[2]-a[2]-5400)<1e-8\n           for a,b in zip(rows[1],rows[2])))\n    ys=sorted(i.trans.disp.y*dbu for i in instances.values())\n    check('all_18_rows_600um_pitch',all(abs(b-a-600)<1e-8 for a,b in zip(ys,ys[1:])))")
ast.parse(v)
p.write_text(v)

readme = ROOT/'COMBINED_SWEEP_README.md'
d=readme.read_text()
# Replace obsolete placement/label documentation; preserve taper and PNC definitions.
d=d.replace('DESIGN_X_OFFSET = 5000.0','DESIGN_X_OFFSET = 0.0')
a=d.index('Design 1 is the left column;')
b=d.index('The unchanged PNC sampler',a)
d=d[:a]+'''## Single vertical column and one label

Design 1 occupies rows 1-9 at origin Y=-700,-100,...,4100 um.
Design 2 occupies rows 10-18 at origin Y=4700,5300,...,9500 um.
All consecutive rows, including the design boundary, have 600 um pitch.
ALL leftmost PNC bend edges align at X=-900 um. There is NO 5000 um X shift.
The table above gives rows within Design 1; add 5400 um to its Y coordinates
for corresponding Design-2 cases. Device geometry is translated as a whole.

Exactly ONE label per device is drawn on Layer 2, for example `D2_T100_P3`:
Design 2, total post-PNC Layer-10 taper 100 um, nominal bend PNC period 3 um.
No row number, meander length, ML/TL label, or second label is drawn.
Each label is a separate `LABEL_<case_id>` cell; the physical device cell
contains no length text. Original length glyphs are removed individually
by exact polygon matching, never by cutting a rectangle into a mask.
Editable label size and offsets are in Section 01 of the combined builder.
The individual standalone Design-1/Design-2 scripts and GDSs are unchanged.

'''+d[b:]
a=d.index('The combined frame encompasses fabrication geometry:')
b=d.index('These checks are geometric',a)
d=d[:a]+'''The reference frame and global overlay markers are rebuilt to enclose the
vertical column, with 500 um margin and a minimum height of 10000 um.
Actual bounds are recorded in the JSON manifest; the 18-row layout is taller
than the old two-column layout. Layer-51 fiber guides are retained as
nonfabrication geometry and are not used to set exposure-frame extents.
All fabrication-mask pairwise overlaps and label-to-mask overlaps are checked.

'''+d[b:]
d=d.replace('compares full end-assembly symmetry, excluding only exact label polygons.',
            'compares full end-assembly symmetry; device cells now have no text.')
d+='''
## Independent placement/label regression

Run `python validate_2611_vertical_sweep.py` after building. When an approved
pre-update GDS and manifest are present in `validation_2611_vertical/before.*`,
it also compares every local device layer against that baseline with only
the individually identified original text glyphs removed. This proves that
stacking and label cleanup did not alter the optical or fabrication shapes.
'''
readme.write_text(d)
print('Updated builder, exported-GDS checks, and documentation; standalone generators untouched.')
