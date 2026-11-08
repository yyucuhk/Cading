"""Guarded one-use migration for the six-device comparison only."""
from pathlib import Path
import ast
import hashlib
ROOT=Path.cwd()
p=ROOT/'2611_design1_L6000_straight_output.py'
s=p.read_text()
blob=hashlib.sha1(b'blob '+str(len(s.encode())).encode()+b'\0'+s.encode()).hexdigest()
assert blob=='500f5f9ad10c0e610e4d0ddced57c366ac9791c6',blob
original=s

def replace(old,new):
    global s
    assert s.count(old)==1,(old[:150],s.count(old))
    s=s.replace(old,new)

end=s.index('"""',3)+3
s='''"""Design 1 L6000: two-coupler and 2000 um straight-output comparison.

The final six-device output omits Layers 11, 16 and 51 everywhere.
Both Layer 2 and Layer 14 release masks are retained only outward from
actual Layer-10 entrances: two windows on two-coupler devices and one
input window on straight-output variants. No release mask remains along
the meander, preceding transition, or the taper-free output extension.

Layer 14 alone uses a 1 um taper-window expansion. Its mask is rebuilt
with the original mask construction and an isolated margin parameter;
Layer 2 retains its original shape inside the terminal windows. Long
channels and pads keep their original definitions. All other geometry
comes from the approved saved polygons, not a new reconstruction.
Labels remain on Layer 2; global alignment marks remain on Layer 15.
These are geometry checks, not optical, mechanical or fabrication qualification.
"""'''+s[end:]
replace('TWO_COUPLER_END_WINDOWS_ONLY = True   # remove L11; keep L14 outward of each L10 entrance',
'''COUPLER_END_WINDOWS_ONLY = True       # applies to BOTH device sets; crops L2 and L14
LAYER14_WINDOW_MARGIN = 1.0           # per-side expansion; independent of original L2 margin''')
replace('REMOVED_OUTPUT_LAYERS = (16, 51)  # remove all datatypes from every output cell',
        'REMOVED_OUTPUT_LAYERS = (11, 16, 51)  # all datatypes, every final output cell')
replace('layer2_windows_and_labels_preserved=True)',
        'layer2_labels_preserved=True, release_masks_handled_separately=True)')
a=s.index('def finalize_two_coupler_windows(');b=s.index('# %% [02]',a)
s=s[:a]+'''def finalize_coupler_windows(layout, entries):
    """Crop L2/L14 in both sets and independently narrow only L14.

    Optical/device cells are flat and privately owned. Labels and frame are
    separate cells and are not edited. The lower-right straight output has
    no L10 taper, so no release window is retained at that end.
    """
    import klayout.db as kdb
    policies = {}
    if not COUPLER_END_WINDOWS_ONLY:
        return policies
    for e in entries:
        cell = layout.cell(e['case_id'])
        require(not list(cell.each_inst()), 'Expected an owned flat device cell')
        rr = regions(layout, cell)
        require(all((n,0) in rr for n in (2,10,14)), 'Missing device masks/taper')
        require(all(d == 0 for n,d in rr if n in (2,14)), 'Unexpected mask datatype')
        r10, r2, r14 = (rr[(n,0)] for n in (10,2,14))
        ml, mc = read_gds(WORK/e['layer14_window_file'])
        candidate = regions(ml,mc)[(14,0)]
        bounds = r14.bbox() + r2.bbox() + r10.bbox()
        xmin,ymin,xmax,ymax = (bounds.left*DBU-1, bounds.bottom*DBU-1,
                              bounds.right*DBU+1, bounds.top*DBU+1)
        keep, ends = kdb.Region(), []
        for side, port in (('input',e['input_port_local_um']),
                           ('output',e['output_port_local_um'])):
            px,py = map(float,port)
            require((px < 0 < py) or (py < 0 < px), 'Unexpected end orientation')
            split_y = py/8.0
            search = box((xmin,split_y,px,ymax) if px < 0 else (px,ymin,xmax,split_y))
            taper = (r10 & search).merged()
            if taper.is_empty():
                require(side == 'output' and e['mode']=='straight_output',
                        'Unexpected missing L10 taper at '+side)
                continue
            sx = (taper.bbox().right if px < 0 else taper.bbox().left)*DBU
            require((sx-px)*(-1 if px<0 else 1)>0, 'L10 entrance must be outward')
            limits = (xmin,split_y,sx,ymax) if px<0 else (sx,ymin,xmax,split_y)
            keep |= box(limits)
            ends.append(dict(side=side, l10_start_local_x_um=sx,
                             center_y_um=py, keep_box_local_um=list(limits)))
        require(len(ends)==(1 if e['mode']=='straight_output' else 2), 'Wrong end count')
        # Intersection with the approved mask forbids widening anywhere.
        # Only L14 uses the fresh 1 um-margin candidate; L2 is a pure crop.
        retained2 = (r2 & keep).merged()
        retained14 = (r14 & candidate & keep).merged()
        for layer,r in ((2,retained2),(14,retained14)):
            require(not r.is_empty(), 'Empty retained mask')
            cell.shapes(layout.layer(layer,0)).clear()
            cell.shapes(layout.layer(layer,0)).insert(r)
        policies[e['case_id']] = dict(ends=ends,
            removed_layer=11, layer2_before_area_um2=r2.area()*DBU**2,
            layer2_after_area_um2=retained2.area()*DBU**2,
            layer14_before_area_um2=r14.area()*DBU**2,
            layer14_after_area_um2=retained14.area()*DBU**2,
            layer14_margin_um=LAYER14_WINDOW_MARGIN,
            layer14_window_file=e['layer14_window_file'])
        e['final_pnc_present']=False
        e['final_release_window_count']=len(ends)
        e['final_layer2_scope']=e['final_layer14_scope']='L10_entrances_outward_only'
        e['layer14_margin_um']=LAYER14_WINDOW_MARGIN
    return policies


def build_layer14_window_candidate(base, width):
    """Use the original mask code with an isolated L14-only margin override.

    This creates a separate candidate mask; all its other layers are discarded.
    The original L2 polygons and all optical polygons are never replaced.
    Input/output positions and channel/pad definitions stay unchanged.
    """
    import gdsfactory as gf
    import klayout.db as kdb
    require(math.isfinite(LAYER14_WINDOW_MARGIN) and 0 < LAYER14_WINDOW_MARGIN <= 3.0,
            'L14 margin must be finite, positive and at most the original 3 um')
    mask_model=base.configured_design(width, {'POSITIVE_WINDOW_MARGIN': LAYER14_WINDOW_MARGIN})
    mask_model.validate_settings()
    component=mask_model.make_device(STRAIGHT_LENGTH)
    r14=kdb.Region(component.kdb_cell.begin_shapes_rec(gf.get_layer((14,0)))).merged()
    require(not r14.is_empty(), 'No Layer14 candidate was generated')
    name='L14_M'+str(LAYER14_WINDOW_MARGIN).replace('.','p')+'_'+base.case_id(width)
    ly=kdb.Layout();ly.dbu=DBU
    insert_regions(ly,name,{(14,0):r14})
    write_layout(ly,WORK/(name+'.gds'))
    return name+'.gds'


'''+s[b:]
replace("    report = dict(width_um=width, reference_id=ref_id, variant_id=variant_id,",
        "    candidate_file=build_layer14_window_candidate(base,width)\n    report = dict(width_um=width, reference_id=ref_id, variant_id=variant_id,\n                  layer14_window_file=candidate_file, layer14_margin_um=LAYER14_WINDOW_MARGIN,")
replace("                       input_port_local_um=meta['input_port_um'],output_port_local_um=meta['output_port_um'])",
        "                       layer14_window_file=meta['layer14_window_file'],\n                       input_port_local_um=meta['input_port_um'],output_port_local_um=meta['output_port_um'])")
replace('    two_coupler_policy = finalize_two_coupler_windows(ly, entries)',
        '    coupler_policy = finalize_coupler_windows(ly, entries)')
a=s.index("        if e['case_id'] in two_coupler_policy:");b=s.index("        check(checks, e['case_id'] + '_exported_layer_set'",a)
s=s[:a]+'''        if e['case_id'] in coupler_policy:
            keep=kdb.Region()
            for end in coupler_policy[e['case_id']]['ends']:
                keep |= box(end['keep_box_local_um'])
            ml,mc=read_gds(WORK/e['layer14_window_file'])
            candidate=regions(ml,mc)[(14,0)]
            expected={layer: ((r & keep).merged() if layer[0]==2 else
                             (r & keep & candidate).merged() if layer[0]==14 else r)
                      for layer,r in expected.items()}
            expected={layer:r for layer,r in expected.items() if not r.is_empty()}
            check(checks,e['case_id']+'_no_layer11',all(layer[0]!=11 for layer in rr))
            for n in (2,14):
                check(checks,e['case_id']+f'_L{n}_terminal_windows_only',
                      (rr[(n,0)]-keep).is_empty())
'''+s[b:]
replace("                extension_layers=[2,11,14,15],frame_bounds_um=[xmin,ymin,xmax,ymax],",
        "                extension_layers=[15],frame_bounds_um=[xmin,ymin,xmax,ymax],")
replace('                two_coupler_layer_policy=two_coupler_policy,',
        '                coupler_layer_policy=coupler_policy,\n                layer14_window_margin_um=LAYER14_WINDOW_MARGIN,')
replace("worker_checks_scope='Unfiltered source masks; final checks enforce global removal of L16/L51 and selective two-coupler L11 removal/L14 cropping'",
        "worker_checks_scope='Pre-export source geometry checks; final checks remove L11/L16/L51 and crop L2/L14 on BOTH sets with isolated 1 um L14 margin'")
replace("'layers16_and51_absent_from_entire_export'","'layers11_16_51_absent_from_entire_export'")
replace("'Design 1: two-coupler devices without PNC below; unchanged straight-output variants above'",
        "'Design 1: no PNCs in either set; release windows only at existing L10 couplers'")
replace("'No output taper: constant-width ridge and 10 um PNC supports (axes expanded)'",
        "'Straight output: unchanged ridge; PNC and release masks removed'")
replace("'Input preserved: 15-to-5 um PNC transition, dense end group, original optical taper'",
        "'Input: original optical taper, terminal-only release masks, 1 um L14 expansion'")
s=s.replace('final deliverable without L16/L51','final deliverable without L11/L16/L51')
ast.parse(s)
def fn(text,n):
    return ast.dump(next(x for x in ast.parse(text).body if isinstance(x,ast.FunctionDef) and x.name==n),include_attributes=False)
for n in ('make_straight_output','load_base','regions','read_gds','write_layout'):
    assert fn(original,n)==fn(s,n),n
p.write_text(s)
print('Updated comparison builder only; optical construction functions unchanged')
