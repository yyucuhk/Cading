"""One-use guarded migration; run in the repository root."""
from pathlib import Path
import ast
import hashlib

p = Path('2611_design1_L6000_straight_output.py')
s = p.read_text()
blob = hashlib.sha1(b'blob ' + str(len(p.read_bytes())).encode() + b'\0' + p.read_bytes()).hexdigest()
if blob != 'f5450677ea7bcbc5798988d073edf246f50ef269':
    raise RuntimeError('Builder changed since the approved baseline: ' + blob)
original = s

def replace(old, new, count=1):
    global s
    found = s.count(old)
    if found != count:
        raise RuntimeError(f'Expected {count}, found {found}: {old[:90]}')
    s = s.replace(old, new)

replace('Layer 14 alone uses a 1 um taper-window expansion. Its mask is rebuilt\nwith the original mask construction and an isolated margin parameter;\nLayer 2 retains its original shape inside the terminal windows. Long\nchannels and pads keep their original definitions. All other geometry',
'''Both release layers use a 1 um taper-window expansion, with independent
controls. Their retained masks begin 10 um outward from each Layer-10
entrance, not at a 10 um waveguide width or at a global X coordinate.
Long channels and pads keep their original definitions. All other geometry''')
replace('actual Layer-10 entrances: two windows on two-coupler devices and one',
        'a configurable offset after Layer-10 entrances: two windows on two-coupler devices and one')
replace('LAYER14_WINDOW_MARGIN = 1.0           # per-side expansion; independent of original L2 margin',
'''LAYER2_WINDOW_MARGIN = 1.0            # per-side taper-window expansion on Layer 2
LAYER14_WINDOW_MARGIN = 1.0           # independent Layer 14 expansion
RELEASE_WINDOW_START_OFFSET = 10.0    # um OUTWARD after L10 starts; applies to BOTH masks''')
s = s.replace('layer14_window_file', 'release_window_file')
s = s.replace('build_layer14_window_candidate', 'build_release_window_candidate')
replace('"""Crop L2/L14 in both sets and independently narrow only L14.',
        '"""Narrow L2/L14 and start both masks after the L10 entrance.')
replace("        candidate = regions(ml,mc)[(14,0)]", "        candidates = regions(ml,mc)")
replace("            limits = (xmin,split_y,sx,ymax) if px<0 else (sx,ymin,xmax,split_y)",
'''            outward = -1.0 if px < 0 else 1.0
            start_x = sx + outward * RELEASE_WINDOW_START_OFFSET
            require(RELEASE_WINDOW_START_OFFSET < taper.bbox().width()*DBU,
                    'Release start lies beyond the Layer10 terminal assembly')
            limits = (xmin,split_y,start_x,ymax) if px<0 else (start_x,ymin,xmax,split_y)''')
replace("            ends.append(dict(side=side, l10_start_local_x_um=sx,\n                             center_y_um=py, keep_box_local_um=list(limits)))",
'''            ends.append(dict(side=side, l10_start_local_x_um=sx,
                             release_start_local_x_um=start_x,
                             start_offset_um=RELEASE_WINDOW_START_OFFSET,
                             center_y_um=py, keep_box_local_um=list(limits)))''')
replace('''        # Intersection with the approved mask forbids widening anywhere.
        # Only L14 uses the fresh 1 um-margin candidate; L2 is a pure crop.
        retained2 = (r2 & keep).merged()
        retained14 = (r14 & candidate & keep).merged()''',
'''        # Rebuild only release-mask candidates; intersecting with the original
        # masks prevents widening. All optical polygons remain untouched.
        retained2 = (r2 & candidates[(2,0)] & keep).merged()
        retained14 = (r14 & candidates[(14,0)] & keep).merged()''')
replace('''            layer14_margin_um=LAYER14_WINDOW_MARGIN,
            release_window_file=e['release_window_file'])''',
'''            layer2_margin_um=LAYER2_WINDOW_MARGIN,
            layer14_margin_um=LAYER14_WINDOW_MARGIN,
            start_offset_um=RELEASE_WINDOW_START_OFFSET,
            release_window_file=e['release_window_file'])''')
replace("        e['final_layer2_scope']=e['final_layer14_scope']='L10_entrances_outward_only'\n        e['layer14_margin_um']=LAYER14_WINDOW_MARGIN",
'''        e['final_layer2_scope']=e['final_layer14_scope']='L10_entrances_plus_offset_outward_only'
        e['release_start_offset_um']=RELEASE_WINDOW_START_OFFSET
        e['layer2_margin_um']=LAYER2_WINDOW_MARGIN
        e['layer14_margin_um']=LAYER14_WINDOW_MARGIN''')
start = s.index('def build_release_window_candidate(')
end = s.index('\n\n# %% [02]', start)
s = s[:start] + '''def build_release_window_candidate(base, width):
    """Create isolated L2/L14 candidates using the existing mask construction.

    Other candidate layers, including labels and optical shapes, are discarded.
    Final intersection preserves the old channels, pads and terminal geometry.
    Each worker owns its configuration and uses a separate cached-cell namespace.
    """
    import gdsfactory as gf
    import klayout.db as kdb
    require(math.isfinite(RELEASE_WINDOW_START_OFFSET) and RELEASE_WINDOW_START_OFFSET >= 0,
            'Release start offset must be finite and non-negative')
    candidates = {}
    margins = {2: LAYER2_WINDOW_MARGIN, 14: LAYER14_WINDOW_MARGIN}
    for margin in sorted(set(margins.values())):
        require(math.isfinite(margin) and 0 < margin <= 3.0,
                'Release margins must be finite, positive and at most the original 3 um')
        mask_model = base.configured_design(width, {'POSITIVE_WINDOW_MARGIN': margin})
        mask_model.validate_settings()
        component = mask_model.make_device(STRAIGHT_LENGTH)
        for layer, requested in margins.items():
            if requested != margin:
                continue
            region = kdb.Region(component.kdb_cell.begin_shapes_rec(gf.get_layer((layer,0)))).merged()
            require(not region.is_empty(), f'No Layer{layer} mask candidate generated')
            candidates[(layer,0)] = region
    name = ('RELEASE_M2_' + str(LAYER2_WINDOW_MARGIN).replace('.', 'p')
            + '_M14_' + str(LAYER14_WINDOW_MARGIN).replace('.', 'p')
            + '_' + base.case_id(width))
    ly = kdb.Layout(); ly.dbu = DBU
    insert_regions(ly, name, candidates)
    write_layout(ly, WORK/(name+'.gds'))
    return name+'.gds'
''' + s[end:]
replace('''                  release_window_file=candidate_file, layer14_margin_um=LAYER14_WINDOW_MARGIN,''',
'''                  release_window_file=candidate_file,
                  layer2_margin_um=LAYER2_WINDOW_MARGIN,
                  layer14_margin_um=LAYER14_WINDOW_MARGIN,
                  release_start_offset_um=RELEASE_WINDOW_START_OFFSET,''')
replace('''            candidate=regions(ml,mc)[(14,0)]
            expected={layer: ((r & keep).merged() if layer[0]==2 else
                             (r & keep & candidate).merged() if layer[0]==14 else r)
                      for layer,r in expected.items()}''',
'''            candidates=regions(ml,mc)
            expected={layer: ((r & keep & candidates[layer]).merged()
                             if layer[0] in (2,14) else r)
                      for layer,r in expected.items()}''')
replace('''                layer14_window_margin_um=LAYER14_WINDOW_MARGIN,
                worker_checks_scope='Pre-export source geometry checks; final checks remove L11/L16/L51 and crop L2/L14 on BOTH sets with isolated 1 um L14 margin', ''',
'''                layer2_window_margin_um=LAYER2_WINDOW_MARGIN,
                layer14_window_margin_um=LAYER14_WINDOW_MARGIN,
                release_window_start_offset_um=RELEASE_WINDOW_START_OFFSET,
                worker_checks_scope='Pre-export source geometry; final checks remove L11/L16/L51, narrow both masks and apply a 10 um outward start offset', ''')
replace("'Input: original optical taper, terminal-only release masks, 1 um L14 expansion'",
        "'Input: unchanged optical taper; both masks 1 um expansion, start offset 10 um'")
ast.parse(s)
a = {n.name:ast.dump(n,include_attributes=False) for n in ast.parse(original).body if isinstance(n,ast.FunctionDef)}
b = {n.name:ast.dump(n,include_attributes=False) for n in ast.parse(s).body if isinstance(n,ast.FunctionDef)}
for name in ['make_straight_output','load_base','finalize_output_layers','read_gds','regions','insert_regions','write_layout']:
    assert a[name] == b[name], name
p.write_text(s)
print('Updated only release masks, metadata and related export expectations')
