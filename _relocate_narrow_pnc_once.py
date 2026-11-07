"""One-time guarded correction: dense PNCs belong at the 3 um-wide end."""
from pathlib import Path
import ast


def once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('Expected one match: ' + old[:100])
    return text.replace(old, new, 1)


for name, function in [
    ('2611_tapered_yue_dual_resist.py', 'phononic_taper_shield'),
    ('2611_tapered_yue_design2.py', 'design2_phononic_taper_shield'),
]:
    path = Path(name)
    original = path.read_text()
    code = once(original,
        '# local longitudinal center-to-center pitch at BOTH taper starts',
        '# longitudinal center-to-center pitch at BOTH narrow PNC ends')
    code = once(code,
        '# five intervals / six stations across the first 5 um; 0 disables',
        '# five intervals / six stations in the LAST 5 um, ending at 3 um width; 0 disables')
    begin = code.index('    dense_end = connection_pitch * connection_periods\n')
    end = code.index('    if len(support_positions) > 1', begin)
    replacement = '''    dense_span = connection_pitch * connection_periods
    if dense_span > cutoff_x + 1e-9:
        raise ValueError("Dense connection segment must fit inside the existing PNC taper")
    # Local +x points from the wide meander toward the narrow optical taper.
    # Anchor reinforcement to cutoff_x (the 3 um-wide PNC end), NOT x=0.
    # Five 1 um periods occupy x=40,41,42,43,44,45 um for a 45 um taper.
    # Restore the regular grid at the wide end and retain its original phase
    # everywhere before the final dense segment. Never shift optical geometry.
    dense_start = max(0.0, cutoff_x - dense_span)
    if connection_periods:
        dense_positions = cutoff_x - np.arange(connection_periods, -1, -1, dtype=float) * connection_pitch
        support_positions = np.concatenate((
            regular_positions[regular_positions < dense_start - 1e-9], dense_positions
        ))
    else:
        support_positions = regular_positions
'''
    code = code[:begin] + replacement + code[end:]
    code = once(code,
        '    component.info["connection_pnc_span_um"] = float(dense_end)',
        '    component.info["connection_pnc_span_um"] = float(dense_span)\n'
        '    component.info["connection_pnc_region"] = "narrow_pnc_end"\n'
        '    component.info["connection_pnc_start_um"] = float(dense_start)\n'
        '    component.info["connection_pnc_end_um"] = float(cutoff_x)')
    def functions(text):
        return {n.name: ast.dump(n, include_attributes=False) for n in ast.parse(text).body
                if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))}
    before, after = functions(original), functions(code)
    changed = [key for key in before if before[key] != after[key]]
    if changed != [function]:
        raise ValueError(f'Unexpected changed geometry functions in {name}: {changed}')
    compile(code, str(path), 'exec')
    path.write_text(code)
    print(name, 'only drawing function changed:', changed)

for name in ('2611_build_combined_sweep.py', '2611_design1_L6000_P10_width_sweep.py'):
    path = Path(name)
    code = path.read_text()
    anchor = '                connection_pnc_span_um=float(m.PNC_CONNECTION_PITCH*m.PNC_CONNECTION_PERIODS),'
    code = once(code, anchor, anchor + '\n'
        '                connection_pnc_region="narrow_pnc_end",\n'
        '                connection_pnc_start_um=float(m.PNC_TAPER_LENGTH-m.PNC_CONNECTION_PITCH*m.PNC_CONNECTION_PERIODS),\n'
        '                connection_pnc_end_um=float(m.PNC_TAPER_LENGTH),')
    path.write_text(code)

path = Path('validate_2611_taper_pitch.py')
code = path.read_text()
code = once(code,
    'Only Layer 11 inside the first 5 um of each existing PNC taper may change.\nFive 1 um intervals are checked, along with the unchanged surrounding grid.',
    'Correct the mistaken wide-end reinforcement and place five 1 um intervals\nin the LAST 5 um, at the 3 um-wide end of each PNC taper. Only Layer 11 in\nthe old and new reinforcement neighborhoods may change. The rest is checked.')
code = once(code, 'DENSE_SPAN = DENSE_PITCH * DENSE_PERIODS',
    'DENSE_SPAN = DENSE_PITCH * DENSE_PERIODS\nDENSE_END = 45.0\nDENSE_START = DENSE_END - DENSE_SPAN')
code = once(code, '''        dense = DENSE_PITCH * np.arange(DENSE_PERIODS + 1)
        expected = np.r_[np.arange(-3*pitch, -.0001, pitch), dense,
                         regular[regular > DENSE_SPAN + 1e-9]]
        allowed_global = kdb.Region()
        local_allowed = box((L/2-.4, -230, L/2+DENSE_SPAN+.4, -170), dbu)
        local_allowed += local_allowed.transformed(kdb.Trans(kdb.Trans.R180))''',
'''        dense = DENSE_START + DENSE_PITCH * np.arange(DENSE_PERIODS + 1)
        expected = np.r_[np.arange(-3*pitch, -.0001, pitch),
                         regular[regular < DENSE_START - 1e-9], dense]
        old_dense_global = kdb.Region()
        new_dense_global = kdb.Region()
        old_dense_local = box((L/2-.4, -230, L/2+DENSE_SPAN+.4, -170), dbu)
        old_dense_local += old_dense_local.transformed(kdb.Trans(kdb.Trans.R180))
        new_dense_local = box((L/2+DENSE_START-.4, -230, L/2+DENSE_END+.4, -170), dbu)
        new_dense_local += new_dense_local.transformed(kdb.Trans(kdb.Trans.R180))''')
code = once(code, '            allowed_global += local_allowed.transformed(inst.trans)',
    '            old_dense_global += old_dense_local.transformed(inst.trans)\n'
    '            new_dense_global += new_dense_local.transformed(inst.trans)')
code = once(code,
    '                    dense_centers = centers[(centers >= -.0011) & (centers <= DENSE_SPAN+.0011)]',
    '                    dense_centers = centers[(centers >= DENSE_START-.0011) & (centers <= DENSE_END+.0011)]')
code = once(code, "_exactly_five_1um_periods',", "_narrow_end_five_1um_periods',")
anchor = "                          sum(abs(centers) < .0011) == 1)"
code = once(code, anchor, anchor + '''
                    check(f'{ident}_{end}_{side}_wide_end_regular_grid_restored',
                          np.array_equal(np.round(centers[(centers >= -.0011) & (centers <= DENSE_SPAN+.0011)],3),
                                         regular[regular <= DENSE_SPAN+1e-9]))
                    check(f'{ident}_{end}_{side}_one_PNC_at_3um_endpoint',
                          sum(abs(centers-DENSE_END) < .0011) == 1)''')
code = once(code,
    '                    dense_pitch_um=DENSE_PITCH,dense_periods=DENSE_PERIODS,dense_span_um=DENSE_SPAN,',
    '                    dense_pitch_um=DENSE_PITCH,dense_periods=DENSE_PERIODS,dense_span_um=DENSE_SPAN,\n'
    '                    dense_start_um=DENSE_START,dense_end_um=DENSE_END,region="narrow_pnc_end",')
code = once(code, '            outside = (delta-allowed_global).area()*dbu**2',
    '            outside = (delta-(old_dense_global | new_dense_global)).area()*dbu**2')
code = once(code, "'_PNC_changes_only_in_first_5um'", "'_PNC_changes_only_in_old_and_new_dense_regions'")
code = once(code,
    "            check(tag+'_no_existing_PNC_removed', (old[(11,0)]-new[(11,0)]).is_empty())",
    "            check(tag+'_only_mistaken_wide_end_additions_removed', ((old[(11,0)]-new[(11,0)])-old_dense_global).is_empty())\n"
    "            check(tag+'_new_supports_only_at_narrow_end', ((new[(11,0)]-old[(11,0)])-new_dense_global).is_empty())")
code = once(code,
    '\n                  dense_pitch_um=DENSE_PITCH,dense_periods=DENSE_PERIODS,dense_span_um=DENSE_SPAN,',
    '\n                  dense_pitch_um=DENSE_PITCH,dense_periods=DENSE_PERIODS,dense_span_um=DENSE_SPAN,\n'
    '                  dense_start_um=DENSE_START,dense_end_um=DENSE_END,region="narrow_pnc_end",')
code = once(code, "            ax.text(x0+.8,15,'Straight / taper junction',fontsize=9)",
    "            ax.text(x0+.8,15,'Wide PNC start: regular pitch restored',fontsize=9)\n"
    "            ax.axvline(x0+DENSE_END,linestyle=':',linewidth=1)\n"
    "            ax.annotate('3 um-wide end: five 1 um periods',\n"
    "                        xy=(x0+DENSE_END-2.5,6.5), xytext=(x0+15,12),\n"
    "                        arrowprops=dict(arrowstyle='->'),fontsize=9)")
code = once(code,
    "('5 um regular pitch' if state=='before' else 'five 1 um periods at connection; rest unchanged')",
    "('previous reinforcement at wrong (wide) end' if state=='before' else 'reinforcement moved to narrow PNC end')")
compile(code, str(path), 'exec')
path.write_text(code)

new_section = '''## Local input/output PNC reinforcement (current)

The reinforcement is at the NARROW, 3 um-wide END of the PNC taper, directly
before the optical taper, not at the wide meander-to-PNC transition.
For the unchanged 45 um PNC taper, five 1 um periods occupy the final 5 um:
six paired stations at local x = 40, 41, 42, 43, 44, 45 um.
Coordinates are measured outward from each meander connection port.
The normal whole-end rotation applies the same geometry at both ends.

The erroneously added wide-end stations at x=1,2,3,4 (and x=5 in the long
layout) are removed. The wide end and all regular stations before x=40 keep
their original pitch and phase. Short-layout stations are 0,5,10,...,35,
then 40,41,42,43,44,45 um. The long-layout stations are 0,10,20,30,
then 40,41,42,43,44,45 um. Straight and bend PNCs are unchanged.

Edit Section 01 / 1D in either underlying design generator:

```python
PNC_CONNECTION_PITCH = 1.0
PNC_CONNECTION_PERIODS = 5  # five intervals / six stations; 0 disables
```

The dense group is anchored to the narrow endpoint (cutoff_x); its start is
cutoff_x minus pitch times periods. Both derived builders inherit these
controls and record the region and start/end distances in their manifests.
The 45 um PNC envelope, optical taper dimensions, other layers, labels and
placements remain unchanged. Run validate_2611_taper_pitch.py after
regenerating all four GDSs. Capture a baseline before editing to compare all
non-PNC layers and confine changes to the old/new reinforcement regions.
Mechanical robustness has not been simulated; checks are geometric only.
'''
for name in ('CAD_EDITING_GUIDE.md', 'COMBINED_SWEEP_README.md', 'DESIGN1_L6000_WIDTH_SWEEP_README.md'):
    path = Path(name)
    text = path.read_text()
    anchor = '## Local input/output PNC reinforcement (current)'
    if text.count(anchor) != 1:
        raise ValueError('Missing documentation section: '+name)
    path.write_text(text.split(anchor)[0] + new_section)
print('Narrow-end correction applied; optical/mask/placement drawing functions preserved.')
