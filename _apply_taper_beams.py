from pathlib import Path
import ast,hashlib
root=Path.cwd()
p=root/'2611_design1_L6000_straight_output.py'
s=p.read_text()
assert hashlib.sha256(p.read_bytes()).hexdigest()=='1c9afcfe2be5c45f308a77321cfa2fce952278558bb3e6b658e6f2b7371171ef'
def replace(a,b):
    global s
    assert s.count(a)==1,(a,s.count(a))
    s=s.replace(a,b,1)
replace('L11/L16/L51 are absent; global marks are on L15 and labels on L2.',
        'L16/L51 and the old PNC patterns are absent. L11 contains only local plain\nbeam supports at taper entrances; global marks are on L15 and labels on L2.')
replace('REMOVED_OUTPUT_LAYERS = (11, 16, 51)  # all datatypes, every final output cell',
'''REMOVED_OUTPUT_LAYERS = (16, 51)       # absent from every final output cell
CLEARED_SOURCE_PNC_LAYERS = (11,)     # remove all old PNCs before adding local beams
LOCAL_TAPER_BEAM_LAYER = (11, 0)      # same solid-support mask as original PNCs
LOCAL_TAPER_BEAM_WIDTH = 0.5          # longitudinal width of each transverse beam
LOCAL_TAPER_BEAM_PITCH = 1.0          # center-to-center, measured along the taper
LOCAL_TAPER_BEAM_COUNT = 3            # three stations; TWO beams per station
LOCAL_TAPER_BEAM_START = 0.0          # first center aligned with the L10 entrance
LOCAL_TAPER_ANCHOR_SPAN = 5.096       # original PNC anchor-to-anchor span
LOCAL_TAPER_ANCHOR_OVERLAP = 0.5      # overlap into slab and outer anchor''')
replace('if info.layer in REMOVED_OUTPUT_LAYERS:',
        'if info.layer in (*REMOVED_OUTPUT_LAYERS, *CLEARED_SOURCE_PNC_LAYERS):')
replace('return dict(removed_layers=list(REMOVED_OUTPUT_LAYERS),',
        'return dict(removed_layers=list(REMOVED_OUTPUT_LAYERS),\n                cleared_source_pnc_layers=list(CLEARED_SOURCE_PNC_LAYERS),')
replace('if key[0] in REMOVED_OUTPUT_LAYERS:',
        'if key[0] in (*REMOVED_OUTPUT_LAYERS, *CLEARED_SOURCE_PNC_LAYERS):')
replace('    layer_policy=finalize_output_layers(ly,fc)\n',
        '    layer_policy=finalize_output_layers(ly,fc)\n    beam_policy=add_local_taper_beams(ly,entries,policies)\n')
replace('for n in (2,10,13,14,15))','for n in (2,10,11,13,14,15))')
replace('final_layer_policy=layer_policy,coupler_layer_policy=policies,checks=checks,workers=workers,',
        'final_layer_policy=layer_policy,coupler_layer_policy=policies,local_taper_beams=beam_policy,checks=checks,workers=workers,')
replace('# %% [02] Complete all-layer straight-output factory',
        (root/'_local_taper_beams_function.txt').read_text()+'# %% [02] Complete all-layer straight-output factory')
ast.parse(s);p.write_text(s)
p=root/'validate_2611_four_sets.py';s=p.read_text()
assert hashlib.sha256(p.read_bytes()).hexdigest()=='1e9a922e1a02a38687dca323adeaf8f535ca47ff7784067c409b8e487a7cc6a5'
s=s.replace('not ({11,16,51}&set(c[\'polys\']))','not ({16,51}&set(c[\'polys\']))')
s=s.replace('                if name in devices and layer == 14:continue\n','')
a=s.index('    if old:\n        check(\'all_cell_names_preserved\'')
b=s.index('    for name,pos in devices.items():',a)
s=s[:a]+'''    if old:
        check('all_cell_names_preserved',set(cells)==set(old))
        check('all_cell_references_preserved',all(cells[n]['refs']==old[n]['refs'] for n in cells))
'''+s[b:]
s=s.replace("set(cells[name]['polys'])=={2,10,13,14,15}","set(cells[name]['polys'])=={2,10,11,13,14,15}")
s=s.replace('for layer in (2,10,13,14,15):','for layer in (2,10,11,13,14,15):')
s=s.replace('for l in (2,10,13,14,15)])','for l in (2,10,11,13,14,15)])')
s=s.replace('removed_layers=[11,16,51]','removed_layers=[16,51]')
needle="    report=dict(status='PASS',check_count=len(checks),checks=checks,measured_entrances=measured,"
assert s.count(needle)==1
s=s.replace(needle,"    from validate_2611_taper_beams import validate as validate_beams\n    beam_report=validate_beams(root,baseline)\n    checks.extend(beam_report['checks'])\n"+needle)
ast.parse(s);p.write_text(s)
p=root/'DESIGN1_STRAIGHT_OUTPUT_README.txt';s=p.read_text()
s=s.replace('Layers 11,16,51 remain absent', 'Layers 16,51 and the old PNC pattern remain absent')
s+='''\nLOCAL PLAIN TAPER SUPPORTS\nAt each retained taper entrance: 3 paired beam stations at 0, 1 and 2 um\noutward from the Layer-10 start. Each transverse beam is 0.5 um wide.\nLayer 11 contains ONLY these 108 plain rectangular beams (18 ends,\n3 stations per end, 2 sides); no square masses or meander PNCs are restored.\nThe anchor-to-anchor span is the original PNC span, 5.096 um, with 0.5 um\noverlap into each anchor. Outer anchors extend beyond both release masks.\nThe bare 2000 um outputs receive no tethers. Layer 2 retains its 15 um\nentrance, 2.5 um downstream margins and 0.5 um interface offset. Layer 14\nretains its 8 um entrance, shared L10 start and 1 um downstream margins.\nThese geometry checks do not establish mechanical strength or optical loss.\n'''
p.write_text(s)
print('Applied local beam support helper and independent validator integration.')
