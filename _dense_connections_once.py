"""Guarded migration: add five 1 um periods only at both PNC-taper starts."""
from pathlib import Path
import ast

ROOT = Path(__file__).resolve().parent


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise RuntimeError(f'Expected one migration anchor: {old[:100]!r}')
    return text.replace(old, new, 1)


for filename, function in (
    ('2611_tapered_yue_dual_resist.py', 'phononic_taper_shield'),
    ('2611_tapered_yue_design2.py', 'design2_phononic_taper_shield'),
):
    path = ROOT / filename
    before = path.read_text()
    source = replace_once(before,
        'PNC_TAPER_TRANSITION_LENGTH = 45.0  # fixed 15 -> 3 um envelope; independent of PNC pitch',
        'PNC_TAPER_TRANSITION_LENGTH = 45.0  # fixed 15 -> 3 um envelope; independent of PNC pitch\n'
        'PNC_CONNECTION_PITCH = 1.0  # local longitudinal center-to-center pitch at BOTH taper starts\n'
        'PNC_CONNECTION_PERIODS = 5  # five intervals / six stations across the first 5 um; 0 disables')
    tree = ast.parse(source)
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == function)
    lines = source.splitlines(keepends=True)
    block = ''.join(lines[node.lineno - 1:node.end_lineno])
    block = replace_once(block,
        '    npoints: int = ADIABATIC_TAPER_POINTS,\n',
        '    npoints: int = ADIABATIC_TAPER_POINTS,\n'
        '    connection_pitch: float = PNC_CONNECTION_PITCH,\n'
        '    connection_periods: int = PNC_CONNECTION_PERIODS,\n')
    block = replace_once(block,
        '    support_positions = np.arange(0.0, cutoff_x + 1e-9, pitch)\n',
        '''    regular_positions = np.arange(0.0, cutoff_x + 1e-9, pitch)
    if not np.isfinite(connection_pitch) or connection_pitch <= 0:
        raise ValueError("connection_pitch must be finite and positive")
    if isinstance(connection_periods, bool) or not isinstance(connection_periods, int) or connection_periods < 0:
        raise ValueError("connection_periods must be a non-negative integer")
    dense_end = connection_pitch * connection_periods
    if dense_end > cutoff_x + 1e-9:
        raise ValueError("Dense connection segment must fit inside the existing PNC taper")
    # Local +x is outward from the meander at BOTH input/output ends.
    # Five 1 um periods mean stations 0,1,2,3,4,5 um, not a global pitch change.
    # Retain the original grid/phase after this segment: do not shift or stretch
    # any remaining PNC. With the approved 5/10 um regular grids this only adds
    # supports, and does not remove or move any original station.
    if connection_periods:
        dense_positions = np.arange(connection_periods + 1, dtype=float) * connection_pitch
        support_positions = np.concatenate((
            dense_positions, regular_positions[regular_positions > dense_end + 1e-9]
        ))
    else:
        support_positions = regular_positions
    if len(support_positions) > 1 and np.min(np.diff(support_positions)) <= cell_extent:
        raise ValueError("Selected PNC pitches would make neighboring support cells touch")
''')
    block = replace_once(block,
        '    component.info["support_positions"] = [float(value) for value in support_positions]\n',
        '    component.info["support_positions"] = [float(value) for value in support_positions]\n'
        '    component.info["connection_pnc_pitch_um"] = float(connection_pitch)\n'
        '    component.info["connection_pnc_periods"] = int(connection_periods)\n'
        '    component.info["connection_pnc_span_um"] = float(dense_end)\n')
    source = ''.join(lines[:node.lineno - 1]) + block + ''.join(lines[node.end_lineno:])
    old_functions = {n.name: ast.dump(n, include_attributes=False) for n in ast.parse(before).body if isinstance(n, ast.FunctionDef)}
    new_functions = {n.name: ast.dump(n, include_attributes=False) for n in ast.parse(source).body if isinstance(n, ast.FunctionDef)}
    changed = [name for name in old_functions if old_functions[name] != new_functions[name]]
    assert changed == [function], changed
    compile(source, filename, 'exec')
    path.write_text(source)
    print(filename, 'only changed geometry function:', changed)

for filename in ('2611_build_combined_sweep.py', '2611_design1_L6000_P10_width_sweep.py'):
    path = ROOT / filename
    source = path.read_text()
    source = replace_once(source,
        '                pnc_taper_length_um=float(m.PNC_TAPER_LENGTH),\n',
        '                pnc_taper_length_um=float(m.PNC_TAPER_LENGTH),\n'
        '                connection_pnc_pitch_um=float(m.PNC_CONNECTION_PITCH),\n'
        '                connection_pnc_periods=int(m.PNC_CONNECTION_PERIODS),\n'
        '                connection_pnc_span_um=float(m.PNC_CONNECTION_PITCH*m.PNC_CONNECTION_PERIODS),\n')
    compile(source, filename, 'exec')
    path.write_text(source)

path = ROOT / 'validate_2611_taper_pitch.py'
source = path.read_text()
source = replace_once(source, 'Only Layer 11 inside each existing 45 um PNC taper may change.',
    'Only Layer 11 inside the first 5 um of each existing PNC taper may change.\nFive 1 um intervals are checked, along with the unchanged surrounding grid.')
source = replace_once(source, "WORK = ROOT / 'validation_2611_taper_pitch'", "WORK = ROOT / 'validation_2611_taper_pitch'\nDENSE_PITCH = 1.0\nDENSE_PERIODS = 5\nDENSE_SPAN = DENSE_PITCH * DENSE_PERIODS")
source = replace_once(source,
    '        expected = np.r_[np.arange(-3*pitch, -.0001, pitch), np.arange(0, 45+1e-9, pitch)]',
    '''        regular = np.arange(0, 45+1e-9, pitch)
        dense = DENSE_PITCH * np.arange(DENSE_PERIODS + 1)
        expected = np.r_[np.arange(-3*pitch, -.0001, pitch), dense,
                         regular[regular > DENSE_SPAN + 1e-9]]''')
source = replace_once(source,
    '        local_allowed = box((L/2-.4, -230, L/2+45.4, -170), dbu)',
    '        local_allowed = box((L/2-.4, -230, L/2+DENSE_SPAN+.4, -170), dbu)')
source = replace_once(source,
    "check(f'{ident}_{end}_{side}_no_compressed_or_missing_interval',\n                          len(spacing) > 0 and np.allclose(spacing, pitch, atol=.0011, rtol=0),",
    "check(f'{ident}_{end}_{side}_correct_piecewise_intervals',\n                          len(spacing) == len(expected)-1 and np.allclose(spacing, np.diff(expected), atol=.0011, rtol=0),")
source = replace_once(source,
    "                    check(f'{ident}_{end}_{side}_one_PNC_at_interface',",
    "                    dense_centers = centers[(centers >= -.0011) & (centers <= DENSE_SPAN+.0011)]\n"
    "                    check(f'{ident}_{end}_{side}_exactly_five_1um_periods',\n"
    "                          len(dense_centers) == 6 and np.allclose(dense_centers, dense, atol=.0011, rtol=0),\n"
    "                          stations_um=dense_centers.tolist())\n"
    "                    check(f'{ident}_{end}_{side}_one_PNC_at_interface',")
source = replace_once(source,
    '                    pitch_um=pitch,taper_length_um=45.0,',
    '                    pitch_um=pitch,taper_length_um=45.0,\n'
    '                    dense_pitch_um=DENSE_PITCH,dense_periods=DENSE_PERIODS,dense_span_um=DENSE_SPAN,')
source = replace_once(source,
    "check(tag+'_PNC_changes_only_in_tapers',outside==0,\n                  outside_tapers_xor_um2=float(outside),",
    "check(tag+'_PNC_changes_only_in_first_5um',outside==0,\n                  outside_dense_region_xor_um2=float(outside),")
source = replace_once(source,
    "            check(tag+'_physical_PNC_update_present',not delta.is_empty())",
    "            check(tag+'_no_existing_PNC_removed', (old[(11,0)]-new[(11,0)]).is_empty())\n"
    "            check(tag+'_physical_PNC_update_present',not delta.is_empty())")
source = replace_once(source,
    "                  taper_ends_checked=2*count,",
    "                  taper_ends_checked=2*count,\n"
    "                  dense_pitch_um=DENSE_PITCH,dense_periods=DENSE_PERIODS,dense_span_um=DENSE_SPAN,")
source = replace_once(source,
    "title=f'Design {design}: {state} - taper PNC pitch '+('15 um' if state=='before' else '5 um, matching straight'))",
    "title=f'Design {design}: {state} - '+('5 um regular pitch' if state=='before' else 'five 1 um periods at connection; rest unchanged'))")
compile(source, 'validate_2611_taper_pitch.py', 'exec')
path.write_text(source)

note = '''

## Local input/output PNC reinforcement (current)

At both input/output straight-to-PNC-taper connections, only the first 5 um
on the taper side now uses five 1 um longitudinal center-to-center periods.
This means six paired PNC stations at local distances 0, 1, 2, 3, 4, 5 um.
The remaining PNCs retain their original positions, dimensions and grid phase.
Short-layout stations continue at 10, 15, ..., 45 um (regular pitch 5 um).
The 6000 um width-sweep stations continue at 10, 20, 30, 40 um (regular pitch
10 um); its 5-to-10 um transition gap is deliberately retained to avoid moving
any original station. Straight and bend PNC patterns are unchanged.

Edit Section 01 / 1D in either underlying design generator:

```python
PNC_CONNECTION_PITCH = 1.0
PNC_CONNECTION_PERIODS = 5  # intervals; six stations; 0 restores the regular grid
```

Both combined/long builders inherit these settings from the applicable source.
The dense segment is mirrored by the normal whole-end construction; no local
patch is applied to just one end. The 45 um envelope, optical tapers, windows,
other layers, labels and placements stay unchanged. The JSON manifests record
the separate connection pitch/count/span. Run `validate_2611_taper_pitch.py`
after regenerating all four outputs; capture an approved baseline before the
edit to enable exact polygon-difference checks. Increased support density is
verified geometrically; mechanical robustness has not been simulated.
'''
for filename in ('CAD_EDITING_GUIDE.md', 'COMBINED_SWEEP_README.md', 'DESIGN1_L6000_WIDTH_SWEEP_README.md'):
    p = ROOT / filename
    p.write_text(p.read_text() + note)
print('Dense PNC controls, case metadata, localized exported-GDS tests and documentation updated.')
