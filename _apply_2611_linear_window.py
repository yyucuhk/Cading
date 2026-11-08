from pathlib import Path
import ast, hashlib
p=Path('2611_design1_L6000_straight_output.py');s=p.read_text()
blob=hashlib.sha1(b'blob '+str(len(p.read_bytes())).encode()+b'\0'+p.read_bytes()).hexdigest()
assert blob=='89ad3332a9c111148d981a7b75ea825b818b9b52',blob
s=s.replace('LAYER14_WINDOW_MARGIN = 1.0           # independent Layer 14 expansion',
'''LAYER14_WINDOW_MARGIN = 1.0           # retained upstream taper-window clearance
LAYER14_LINEAR_START_FROM_L10 = 58.847 # gives global output X=5346 um in this layout
LAYER14_GLUE_DIAMETER = 300.0         # preserve both existing circular windows''',1)
old='    beam_policy=add_local_taper_beams(ly,entries,policies)'
new=old+'''
    from linear_layer14_window import linearize_layer14_windows
    linear_policy=linearize_layer14_windows(ly,entries,policies,
        start_distance_um=LAYER14_LINEAR_START_FROM_L10,
        circle_diameter_um=LAYER14_GLUE_DIAMETER)'''
assert s.count(old)==1;s=s.replace(old,new,1)
old='final_layer_policy=layer_policy,coupler_layer_policy=policies,local_taper_beams=beam_policy,checks=checks,workers=workers,'
new='final_layer_policy=layer_policy,coupler_layer_policy=policies,local_taper_beams=beam_policy,linear_layer14_windows=linear_policy,checks=checks,workers=workers,'
assert s.count(old)==1;s=s.replace(old,new,1)
s=s.replace('Requires the existing source modules and source width-sweep GDS/manifest.',
'''After 58.847 um from each Layer-10 entrance, only L14 continues as a
straight-sided channel through the two existing 300 um glue circles.
The intermediate L14 frame is omitted; Layer 2 and optics are unchanged.
Requires the existing source modules and source width-sweep GDS/manifest.''',1)
ast.parse(s);p.write_text(s)
p=Path('validate_2611_four_sets.py');s=p.read_text()
old='                for dist in (50,60,70,75):\n'
new='''                for dist in (50,60,70,75):
                    # The requested linear L14 continuation replaces constant
                    # clearance here; its edges receive an independent full audit.
                    if layer==14 and meta.get('linear_layer14_windows',{}).get('enabled') and dist>=58.847:
                        continue
'''
assert s.count(old)==1;s=s.replace(old,new,1);ast.parse(s);p.write_text(s)
p=Path('DESIGN1_STRAIGHT_OUTPUT_README.txt')
p.write_text(p.read_text()+'''
LINEAR LAYER 14 CONTINUATION (current revision)
Only L14 changes. Keep its complete upstream geometry through 58.847 um
from the L10 entrance (output global X=5346 um; mirrored input X=-861.694 um).
Measure and retain the existing width there (2.644 um on the 1 nm grid).
Replace all subsequent channel/rectangular-frame fragments with ONE linear
trapezoid to the unchanged far-end cross-section, plus the two 300 um glue
circles at their existing centers. No extra intermediate L14 frame remains.
The 1 um margin applies to the preserved upstream taper, not the widening
linear continuation. Circle outlines use the original 2.5 degree sampling.
L2, L10, L11 beams, L13, L15, labels, markers and placements are unchanged.
Controls: LAYER14_LINEAR_START_FROM_L10 and LAYER14_GLUE_DIAMETER in Section 01.
Audit against the preceding GDS: python validate_2611_linear_window.py --baseline <previous.gds> --render
The linear channel is width-continuous at the join; no optical adiabaticity
or tangent-continuity claim is made for this exposure-mask modification.
''')
print('Patched isolated L14 continuation and corresponding audit expectations')
