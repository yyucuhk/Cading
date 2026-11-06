"""One-time guarded migration: match taper PNC pitch without moving masks."""
from pathlib import Path
import ast

ROOT = Path(__file__).resolve().parent

def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f'Expected one source match: {old[:100]!r}, found {text.count(old)}')
    return text.replace(old, new, 1)

for name in ('2611_tapered_yue_dual_resist.py', '2611_tapered_yue_design2.py'):
    path = ROOT / name
    before = path.read_text()
    after = replace_once(before,
        'REVIEW_SUPPORT_ANCHOR_PERIOD = 15.0  # PNC TAPER pitch only; preserves its 45 um length',
        'REVIEW_SUPPORT_ANCHOR_PERIOD = MEANDER_STRAIGHT_PNC_PERIOD  # taper follows straight pitch\n'
        'PNC_TAPER_TRANSITION_LENGTH = 45.0  # fixed 15 -> 3 um envelope; independent of PNC pitch')
    after = replace_once(after,
        'PNC_TAPER_LENGTH = 3 * REVIEW_SUPPORT_ANCHOR_PERIOD',
        'PNC_TAPER_LENGTH = PNC_TAPER_TRANSITION_LENGTH')
    after = replace_once(after,
        '# The PNC taper length is three straight periods; see Section 02.',
        '# PNC taper pitch follows the straight sections; its physical length stays fixed.')
    after = replace_once(after,
        '    support_positions = np.arange(pitch, cutoff_x + 1e-9, pitch)',
        '    # Continue the straight-section grid across the meander/taper interface.\n'
        '    # The meander stops one pitch before this port, so x=0 is not a duplicate.\n'
        '    # Keep the last full-pitch station; do not squeeze a final PNC onto the\n'
        '    # endpoint when cutoff_x is not an integer multiple of pitch.\n'
        '    support_positions = np.arange(0.0, cutoff_x + 1e-9, pitch)')
    after = replace_once(after,
        '    # Extend the PNC mask by half a cell at both longitudinal ends. The last\n'
        '    # cell remains centered exactly at the 3 um taper endpoint but is no\n'
        '    # longer cut in half by the end of the surrounding positive-resist band.',
        '    # Keep the approved continuous mask envelope, extended by half a cell\n'
        '    # at both longitudinal ends. PNC sampling does not change this mask;\n'
        '    # a cell on either endpoint is fully covered by the surrounding band.')
    # Only the taper sampler may change among the drawing functions.
    old_functions = {n.name: ast.dump(n, include_attributes=False) for n in ast.parse(before).body
                     if isinstance(n, ast.FunctionDef)}
    new_functions = {n.name: ast.dump(n, include_attributes=False) for n in ast.parse(after).body
                     if isinstance(n, ast.FunctionDef)}
    changed = [n for n in old_functions if old_functions[n] != new_functions[n]]
    expected = 'phononic_taper_shield' if 'dual_resist' in name else 'design2_phononic_taper_shield'
    if changed != [expected]:
        raise ValueError((name, changed, expected))
    path.write_text(after)

path = ROOT / '2611_build_combined_sweep.py'
s = path.read_text()
s = replace_once(s,
    '    require(m.PNC_TAPER_LENGTH == 45.0 and m.REVIEW_SUPPORT_ANCHOR_PERIOD == 15.0,\n'
    "            'Straight pitch must not change the PNC taper region')",
    '    require(m.PNC_TAPER_LENGTH == 45.0 and\n'
    '            m.REVIEW_SUPPORT_ANCHOR_PERIOD == STRAIGHT_PNC_PERIOD,\n'
    "            'Taper PNC pitch must match the straight pitch without changing its 45 um length')")
s = replace_once(s,
    'Straight PNC center-to-center pitch is exactly 5 um, independent of the unchanged 15 um taper-PNC pitch and 45 um PNC taper.',
    'Straight and taper PNC center-to-center pitches match. The taper grid starts at the connection port and its 45 um envelope is unchanged.')
path.write_text(s)

path = ROOT / '2611_design1_L6000_P10_width_sweep.py'
s = path.read_text()
s = replace_once(s,
    'The separate PNC-taper pitch stays 15 um. No existing source/GDS is rewritten.',
    'The PNC-taper pitch follows the 10 um straight pitch; its length stays 45 um.\n'
    'This builder does not rewrite the source or the other GDS outputs.')
s = replace_once(s,
    "    check(checks, 'separate_PNC_taper_unchanged', m.PNC_TAPER_LENGTH==45 and m.REVIEW_SUPPORT_ANCHOR_PERIOD==15)",
    "    check(checks, 'taper_PNC_matches_straight_without_changing_length',\n"
    '          m.PNC_TAPER_LENGTH==45 and m.REVIEW_SUPPORT_ANCHOR_PERIOD==STRAIGHT_PNC_PERIOD)')
s = replace_once(s,
    "                post_pnc_taper_total_um=80., label=f'D1_W{width:.2f}_T80_P{BEND_PNC_PERIOD:g}',",
    '                straight_pnc_period_um=float(m.MEANDER_STRAIGHT_PNC_PERIOD),\n'
    '                taper_pnc_period_um=float(m.REVIEW_SUPPORT_ANCHOR_PERIOD),\n'
    '                pnc_taper_length_um=float(m.PNC_TAPER_LENGTH),\n'
    "                post_pnc_taper_total_um=80., label=f'D1_W{width:.2f}_T80_P{BEND_PNC_PERIOD:g}',")
path.write_text(s)

path = ROOT / 'validate_2611_straight_pitch.py'
s = replace_once(path.read_text(),
    "entry['straight_pnc_period_um'] == 5.0 and entry['taper_pnc_period_um'] == 15.0",
    "entry['straight_pnc_period_um'] == 5.0 and entry['taper_pnc_period_um'] == 5.0")
path.write_text(s)

for name in ('CAD_EDITING_GUIDE.md', 'COMBINED_SWEEP_README.md', 'DESIGN1_L6000_WIDTH_SWEEP_README.md'):
    path = ROOT / name
    s = path.read_text()
    s = s.replace('has its own derived length (currently 45 um).',
                  'has its own fixed length (45 um), independent of its PNC pitch.')
    s = s.replace('45 um PNC taper still uses its existing 15 um period. Only bend PNC periods are swept.',
                  '45 um PNC taper now uses the same 5 um period as the straights. Only bend PNC periods are swept.')
    s = s.replace('`REVIEW_SUPPORT_ANCHOR_PERIOD` controls the separate PNC-taper pitch\n'
                  '(currently 15 um); the PNC taper length is three of those periods, or 45 um.',
                  '`REVIEW_SUPPORT_ANCHOR_PERIOD` is an alias of `MEANDER_STRAIGHT_PNC_PERIOD`.\n'
                  'The PNC taper uses the same pitch; `PNC_TAPER_TRANSITION_LENGTH` fixes its length at 45 um.')
    s = s.replace('`REVIEW_SUPPORT_ANCHOR_PERIOD = 15.0` is now explicitly the independent\n'
                  'PNC-taper pitch. Its three-period 45 um transition and all optical taper\n'
                  'lengths stay unchanged. Do not change it to edit straight-meander pitch.',
                  '`REVIEW_SUPPORT_ANCHOR_PERIOD = MEANDER_STRAIGHT_PNC_PERIOD` makes the\n'
                  'PNC taper follow the straight-section pitch. Its physical length is fixed\n'
                  'separately by `PNC_TAPER_TRANSITION_LENGTH = 45.0`; optical tapers do not move.')
    s = s.replace('The separate 45 um PNC TAPER transition retains its prior 15 um pitch.',
                  'The 45 um PNC TAPER transition now uses the same 10 um pitch as the straights.\n'
                  'Stations are at 0, 10, 20, 30, 40 um from its port. The final 5 um\n'
                  'is an unpopulated remainder, not a compressed last interval.')
    s += '\n## Taper-junction pitch correction\n\n'
    s += ('PNCs now continue on the straight-section grid through each 45 um PNC taper,\n'
          'including a pair at the connection port. The 1400 um designs use 5 um\n'
          'pitch on both straights and tapers; the separate 6000 um width sweep uses\n'
          '10 um on both. Bend samplers and all optical/release/protection geometry\n'
          'are unchanged. No extra cell is forced onto a non-grid taper endpoint.\n'
          'Run `python validate_2611_taper_pitch.py` after regenerating the four GDSs\n'
          'to check both taper ends and the junction-to-straight spacings. Historical\n'
          'baseline reports describe the previous geometry; the new report is\n'
          '`gds/2611_taper_pitch_validation.json`.\n')
    path.write_text(s)
print('Only taper-PNC sampling and its pitch/length controls changed in the drawing code.')
