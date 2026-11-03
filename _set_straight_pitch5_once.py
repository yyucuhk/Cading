"""Guarded migration: exact 5 um straight pitch, independent of PNC taper pitch."""
from pathlib import Path
import ast

ROOT = Path(__file__).resolve().parent

def replace(text, old, new, count=1):
    assert text.count(old) == count, (old[:90], text.count(old), count)
    return text.replace(old, new)

for name in ('2611_tapered_yue_dual_resist.py', '2611_tapered_yue_design2.py'):
    path = ROOT / name
    text = path.read_text()
    text = replace(text,
        'REVIEW_SUPPORT_ANCHOR_PERIOD = 15.0  # requested PNC pitch on straight regions',
        'MEANDER_STRAIGHT_PNC_PERIOD = 5.0  # exact center-to-center pitch on all three straights\n'
        'REVIEW_SUPPORT_ANCHOR_PERIOD = 15.0  # PNC TAPER pitch only; preserves its 45 um length')
    shield = 'design2_phononic_meander_shield' if 'design2' in name else 'phononic_meander_shield'
    start = text.index('def ' + shield + '(')
    end = text.index('\n\n@gf.cell', start)
    block = text[start:end]
    block = replace(block, 'pitch: float = REVIEW_SUPPORT_ANCHOR_PERIOD,',
                    'pitch: float = MEANDER_STRAIGHT_PNC_PERIOD,')
    block = replace(block,
        '        usable_length = usable_end - usable_start\n',
        '''        # Straights use a fixed pitch, NOT a linspace fitted to the span.
        # Anchor the first straight backwards from its bend. Anchor the other
        # straights forwards from their starting bend. Shared boundary pairs
        # are added below exactly once. Any residual stays at an end clearance.
        if not is_bend:
            if start == 0:
                count = max(0, int(np.floor(
                    (section_length - cell_half_length) / requested_pitch + 1e-9
                )))
                path_positions = section_end_length - requested_pitch * np.arange(1, count + 1)
            else:
                clearance = requested_pitch if stop < len(lengths) else cell_half_length
                count = max(0, int(np.floor(
                    (section_length - clearance) / requested_pitch + 1e-9
                )))
                path_positions = section_start_length + requested_pitch * np.arange(1, count + 1)
            for path_position in np.sort(path_positions):
                add_pnc_pair(float(path_position))
            continue

        # Preserve the existing fitted-to-arc BEND sampler without changes.
        usable_length = usable_end - usable_start
''')
    block = replace(block, '    component.info["support_positions"] = support_positions\n',
                    '    component.info["straight_pnc_pitch_um"] = float(pitch)\n'
                    '    component.info["straight_pnc_sampling"] = "exact_pitch_bend_anchored"\n'
                    '    component.info["support_positions"] = support_positions\n')
    text = text[:start] + block + text[end:]
    text = replace(text, '    support_period = REVIEW_SUPPORT_ANCHOR_PERIOD\n',
                   '    support_period = MEANDER_STRAIGHT_PNC_PERIOD\n')
    text = replace(text, '    dimensions = {\n',
                   '    dimensions = {\n        "MEANDER_STRAIGHT_PNC_PERIOD": MEANDER_STRAIGHT_PNC_PERIOD,\n')
    ast.parse(text)
    path.write_text(text)

p = ROOT / '2611_build_combined_sweep.py'
s = p.read_text()
s = replace(s, 'BEND_PNC_PERIODS = (2.0, 3.0, 5.0)',
            'STRAIGHT_PNC_PERIOD = 5.0  # exact center-to-center pitch; independent of bend/taper PNC\n'
            'BEND_PNC_PERIODS = (2.0, 3.0, 5.0)')
s = replace(s, '                  REVIEW_SUPPORT_BEND_PERIOD=case[\'bend_period_um\'])',
            '                  REVIEW_SUPPORT_BEND_PERIOD=case[\'bend_period_um\'],\n'
            '                  MEANDER_STRAIGHT_PNC_PERIOD=STRAIGHT_PNC_PERIOD)')
s = replace(s, "        capture['bend_pitch'] = float(kwargs['bend_pitch'])",
            "        capture['bend_pitch'] = float(kwargs['bend_pitch'])\n"
            "        capture['straight_pitch'] = float(kwargs['pitch'])")
s = replace(s, '    bend_measurements = []\n',
'''    # Measure each actual straight population, including bend-boundary pitch.
    require(abs(capture['straight_pitch'] - STRAIGHT_PNC_PERIOD) < 1e-12,
            'Wrong straight PNC pitch sent to the shield')
    require(m.PNC_TAPER_LENGTH == 45.0 and m.REVIEW_SUPPORT_ANCHOR_PERIOD == 15.0,
            'Straight pitch must not change the PNC taper region')
    straight_measurements = []
    for a, b, is_bend in sections:
        if is_bend:
            continue
        positions = sorted(p for p in capture['positions']
                           if cumulative[a] + 1e-6 < p < cumulative[b] - 1e-6)
        diffs = np.diff(positions)
        require(len(diffs) > 1 and np.allclose(diffs, STRAIGHT_PNC_PERIOD, atol=1e-8, rtol=0),
                'Straight PNC spacing is not the exact requested pitch')
        if a > 0:
            require(abs(positions[0] - cumulative[a] - STRAIGHT_PNC_PERIOD) < 1e-8,
                    'Incorrect straight pitch next to its starting bend')
        if b < len(seg):
            require(abs(cumulative[b] - positions[-1] - STRAIGHT_PNC_PERIOD) < 1e-8,
                    'Incorrect straight pitch next to its ending bend')
        straight_measurements.append(dict(
            interior_pairs=len(positions), min_pitch_um=float(min(diffs)),
            max_pitch_um=float(max(diffs)),
            start_um=capture['points'][a].tolist(), end_um=capture['points'][b].tolist()))
    bend_measurements = []
''')
s = replace(s, 'straight_pnc_period_um=float(m.REVIEW_SUPPORT_ANCHOR_PERIOD),',
            'straight_pnc_period_um=float(m.MEANDER_STRAIGHT_PNC_PERIOD),\n'
            '                straight_pnc_measurements=straight_measurements,\n'
            '                taper_pnc_period_um=float(m.REVIEW_SUPPORT_ANCHOR_PERIOD),\n'
            '                pnc_taper_length_um=float(m.PNC_TAPER_LENGTH),')
s = replace(s, "    fields=['case_id','design','row','global_row','label','bend_period_um'",
            "    fields=['case_id','design','row','global_row','label','straight_pnc_period_um','bend_period_um'")
s = replace(s, "                  notes=['Bend period is the existing nominal parameter.",
            "                  notes=['Straight PNC center-to-center pitch is exactly 5 um, independent of the unchanged 15 um taper-PNC pitch and 45 um PNC taper.',\n"
            "                         'Bend period is the existing nominal parameter.")
ast.parse(s)
p.write_text(s)

for filename in ('CAD_EDITING_GUIDE.md', 'COMBINED_SWEEP_README.md'):
    p = ROOT / filename
    s = p.read_text()
    s = s.replace('The straight-region\nPNC nominal period remains 15 um.',
                  'The straight-meander PNC period is exactly 5 um. The separate\n45 um PNC taper still uses its existing 15 um period.')
    s = s.replace('REVIEW_SUPPORT_ANCHOR_PERIOD`: requested PNC pitch on straights',
                  'MEANDER_STRAIGHT_PNC_PERIOD`: exact PNC pitch on meander straights')
    s += '''\n## Independent straight-meander PNC pitch (current)\n\n`STRAIGHT_PNC_PERIOD = 5.0` in Section 01 of `2611_build_combined_sweep.py`\ncontrols the straight PNC spacing for all 18 devices. The two standalone\ngenerators use `MEANDER_STRAIGHT_PNC_PERIOD = 5.0` in Section 01 / 1D.\nThis is CENTER-TO-CENTER pitch, not edge-to-edge clearance. All three 1400 um\nstraights use an exact 5 um grid, with 279 interior PNC pairs per straight.\nOne shared PNC pair at each straight/bend junction is retained.\n\n`REVIEW_SUPPORT_ANCHOR_PERIOD = 15.0` is now explicitly the independent\nPNC-taper pitch. Its three-period 45 um transition and all optical taper\nlengths stay unchanged. Do not change it to edit straight-meander pitch.\nThe bend sweep remains nominal 2/3/5 um with the existing fitted-to-arc\nsampler. Labels still report design, post-PNC taper length and bend period\nonly; the common straight period is recorded separately in CSV/JSON.\n\nRun `python validate_2611_straight_pitch.py` after regenerating the combined\nsweep to measure every straight PNC population directly in the exported GDS.\n'''
    p.write_text(s)
print('Patched exact straight PNC pitch in both designs and combined builder.')
