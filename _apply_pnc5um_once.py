"""One-time, guarded migration of the approved source to an earlier 5 um cutoff."""
from pathlib import Path
import ast
ROOT = Path(__file__).resolve().parent

def replace(s, old, new, count=1):
    if s.count(old) != count:
        raise ValueError(f'Unexpected match count {s.count(old)} != {count}: {old[:100]}')
    return s.replace(old, new)

helper = '''def pnc_transition_cutoff_length(length, width1, width2, cutoff_width):
    """Find the requested width on the original quintic transition.

    length is the 45 um REFERENCE profile, not the active PNC length.
    All downstream geometry shares this root before database-grid rounding.
    """
    values = (length, width1, width2, cutoff_width)
    if not all(np.isfinite(value) for value in values):
        raise ValueError("PNC transition dimensions must be finite")
    if length <= 0 or not 0 < width2 <= cutoff_width < width1:
        raise ValueError("Require length > 0 and 0 < width2 <= cutoff < width1")
    if abs(cutoff_width - width2) < 1e-12:
        return float(length)
    lo, hi = 0.0, 1.0
    for _ in range(64):
        mid = (lo + hi) / 2
        weight = 6 * mid**5 - 15 * mid**4 + 10 * mid**3
        width = width1 + (width2 - width1) * weight
        if width > cutoff_width:
            lo = mid
        else:
            hi = mid
    return float(length * (lo + hi) / 2)


PNC_TAPER_LENGTH = pnc_transition_cutoff_length(
    PNC_TAPER_TRANSITION_LENGTH,
    SLAB_EXTENSION_END_WIDTH,
    PNC_TAPER_PROFILE_END_WIDTH,
    PNC_TAPER_SUPPORT_CUTOFF_WIDTH,
)'''
for name in ('2611_tapered_yue_dual_resist.py', '2611_tapered_yue_design2.py'):
    p = ROOT / name
    s = p.read_text()
    s = replace(s, 'PNC_TAPER_SUPPORT_CUTOFF_WIDTH = 3.0', 'PNC_TAPER_SUPPORT_CUTOFF_WIDTH = 5.0')
    s = replace(s, 'PNC_TAPER_TRANSITION_LENGTH = 45.0  # fixed 15 -> 3 um envelope; independent of PNC pitch',
                'PNC_TAPER_TRANSITION_LENGTH = 45.0  # REFERENCE 15 -> 3 um profile length\nPNC_TAPER_PROFILE_END_WIDTH = 3.0  # reference curve only; active PNC stops at 5 um')
    s = replace(s, 'PNC_TAPER_LENGTH = PNC_TAPER_TRANSITION_LENGTH', helper)
    prefix = 'design2_' if 'design2' in name else ''
    s = replace(s, f'def {prefix}phononic_taper_shield(\n    length: float = SLAB_EXTENSION_LENGTH,\n    width1: float = SLAB_EXTENSION_END_WIDTH,\n    width2: float = SLAB_EXTENSION_START_WIDTH,',
                f'def {prefix}phononic_taper_shield(\n    length: float = PNC_TAPER_TRANSITION_LENGTH,\n    width1: float = SLAB_EXTENSION_END_WIDTH,\n    width2: float = PNC_TAPER_PROFILE_END_WIDTH,')
    start = s.index('    # Solve the same quintic width profile used by ' + prefix + 'tapered_fiber_adiabatic.')
    stop = s.index('    x = np.linspace(0.0, cutoff_x, npoints)', start)
    s = s[:start] + '    # Share the exact root with the optical tapers and exposure masks.\n    cutoff_x = pnc_transition_cutoff_length(length, width1, width2, cutoff_width)\n\n' + s[stop:]
    s = replace(s, f'{prefix}phononic_taper_shield(\n                length=PNC_TAPER_LENGTH,\n                width1=SLAB_EXTENSION_END_WIDTH,\n                width2=PNC_TAPER_END_WIDTH,',
                f'{prefix}phononic_taper_shield(\n                length=PNC_TAPER_TRANSITION_LENGTH,\n                width1=SLAB_EXTENSION_END_WIDTH,\n                width2=PNC_TAPER_PROFILE_END_WIDTH,')
    s = replace(s, 'support_u = SUPPORT_TIP_POSITION / PNC_TAPER_LENGTH', 'support_u = SUPPORT_TIP_POSITION / PNC_TAPER_TRANSITION_LENGTH')
    s = replace(s, 'SLAB_EXTENSION_END_WIDTH + (PNC_TAPER_END_WIDTH - SLAB_EXTENSION_END_WIDTH) * support_weight',
                'SLAB_EXTENSION_END_WIDTH + (PNC_TAPER_PROFILE_END_WIDTH - SLAB_EXTENSION_END_WIDTH) * support_weight')
    s = replace(s, '    if add_support_tips:\n        for side in (-1, 1):',
                '    if add_support_tips:\n        if not 0 <= SUPPORT_TIP_POSITION <= PNC_TAPER_LENGTH:\n            raise ValueError("Optional support tips must lie inside the active PNC region")\n        for side in (-1, 1):')
    s = replace(s, '    component.info["pnc_taper_length"] = PNC_TAPER_LENGTH',
                '    component.info["pnc_taper_length"] = PNC_TAPER_LENGTH\n    component.info["pnc_reference_profile_length_um"] = PNC_TAPER_TRANSITION_LENGTH\n    component.info["pnc_cutoff_width_um"] = PNC_TAPER_SUPPORT_CUTOFF_WIDTH')
    s = s.replace('outward from the 3 um PNC endpoint', 'outward from the 5 um PNC endpoint')
    s = s.replace('ending at 3 um width', 'ending at 5 um width')
    s = s.replace('cutoff_x (the 3 um-wide PNC end)', 'cutoff_x (the 5 um-wide PNC end)')
    s = s.replace('Five 1 um periods occupy x=40,41,42,43,44,45 um for a 45 um taper.',
                  'Five 1 um periods occupy x=26.378626..31.378626 um at the shortened endpoint.')
    s = s.replace('Layer 10: 3 -> 1 over', 'Layer 10: 5 -> 1 over')
    s = s.replace('Layer 10:       3 -> 1.65', 'Layer 10:       5 -> 1.65')
    s = '# Current revision: active PNC stops at 5 um width; see PNC5UM_README.txt.\n' + s
    ast.parse(s); p.write_text(s)

p=ROOT/'2611_build_combined_sweep.py';s=p.read_text()
s=replace(s,'require(m.PNC_TAPER_LENGTH == 45.0 and','require(0 < m.PNC_TAPER_LENGTH < m.PNC_TAPER_TRANSITION_LENGTH and\n            m.PNC_TAPER_SUPPORT_CUTOFF_WIDTH == 5.0 and')
s=replace(s,"'Taper PNC pitch must match the straight pitch without changing its 45 um length'","'PNC pitch must match the straight pitch and stop at the 5 um cutoff'")
s=replace(s,"widths10 = [3., 1., .2] if case['design'] == 1 else [3., 1.65, .98, .2]","widths10 = [float(m.PNC_TAPER_END_WIDTH), 1., .2] if case['design'] == 1 else [float(m.PNC_TAPER_END_WIDTH), 1.65, .98, .2]")
s=replace(s,'                pnc_taper_length_um=float(m.PNC_TAPER_LENGTH),','                pnc_taper_length_um=float(m.PNC_TAPER_LENGTH),\n                pnc_cutoff_width_um=float(m.PNC_TAPER_SUPPORT_CUTOFF_WIDTH),\n                pnc_reference_profile_length_um=float(m.PNC_TAPER_TRANSITION_LENGTH),')
ast.parse(s);p.write_text(s)

p=ROOT/'2611_design1_L6000_P10_width_sweep.py';s=p.read_text()
a=s.index('        for layer,x,w in ((10,45.002,3.)');b=s.index("    check(checks, 'Design1_taper_lengths_preserved'",a)
s=s[:a]+'''        pnc_end = float(m.PNC_TAPER_LENGTH)
        first_end = pnc_end + float(m.RIDGE_TAPER_LENGTH)
        last_end = first_end + float(m.SLAB_TAPER_LENGTH)
        for layer,x,w in ((10,pnc_end+.002,float(m.PNC_TAPER_END_WIDTH)),
                          (10,first_end-.002,1.),(10,last_end-.002,.2),
                          (15,pnc_end+.002,width),(15,first_end-.002,.1)):
            x += STRAIGHT_LENGTH/2
            half_clip = max(2.0, w/2 + .02)
            section = (rr[(layer,0)] & box((x,-200-half_clip,x+DBU,-200+half_clip))).merged()
            measured = section.bbox().height()*DBU
            check(checks, f'end{side}_L{layer}_x{x:g}_taper_width', section.count()==1 and abs(measured-w)<.0041,
                  measured_width_um=measured, expected_width_um=w)
    check(checks, 'taper_PNC_stops_at_5um_and_matches_straight',
          0 < m.PNC_TAPER_LENGTH < m.PNC_TAPER_TRANSITION_LENGTH and
          m.PNC_TAPER_SUPPORT_CUTOFF_WIDTH == 5.0 and
          m.REVIEW_SUPPORT_ANCHOR_PERIOD == STRAIGHT_PNC_PERIOD)
'''+s[b:]
s=replace(s,'                pnc_taper_length_um=float(m.PNC_TAPER_LENGTH),','                pnc_taper_length_um=float(m.PNC_TAPER_LENGTH),\n                pnc_cutoff_width_um=float(m.PNC_TAPER_SUPPORT_CUTOFF_WIDTH),\n                pnc_reference_profile_length_um=float(m.PNC_TAPER_TRANSITION_LENGTH),')
s=s.replace('L10 3 -> 1 -> 0.2 over 40/40; L15 W -> 0.1 over 40 after the 45 um PNC.','L10 5 -> 1 -> 0.2 over 40/40; L15 W -> 0.1 over 40 after the shortened PNC.')
s=s.replace('The PNC-taper pitch follows the 10 um straight pitch; its length stays 45 um.','The PNC ends at the 5 um width point of the original 45 um reference profile.')
ast.parse(s);p.write_text(s)

p=ROOT/'validate_2611_combined_sweep.py';s=p.read_text()
s=replace(s,"x=e['straight_um']/2+45.","x=e['straight_um']/2+e['pnc_taper_length_um']")
s=replace(s,'                    cr=(r & kdb.Region(kdb.Box(ix,round(-202/dbu),ix+1,round(-198/dbu)))).merged()',
'''                    half_clip=max(2.0,float(width)/2+.02)
                    cr=(r & kdb.Region(kdb.Box(ix,round((-200-half_clip)/dbu),ix+1,round((-200+half_clip)/dbu)))).merged()''')
ast.parse(s);p.write_text(s)

note = '''PNC CUTOFF AT 5 UM - CURRENT GEOMETRY
====================================
The original 45 um quintic 15-to-3 um reference profile is truncated at
its 5 um width point: active PNC length = 31.378626347 um. The 45 um value
is a reference-curve scale, NOT the active PNC length.

Five 1 um reinforcement periods end at this new cutoff. Six paired centers
are at 26.378626347 through 31.378626347 um. Normal straight, bend, and earlier
taper PNC positions are preserved. Both input/output ends use this rule.

Layer 10 starts at width 5 um. Design 1 retains 5 -> 1 -> 0.2 um stages;
Design 2 retains 5 -> 1.65 -> 0.980 -> 0.200 um. Downstream stage lengths,
ridge widths, the 840 nm Design-2 shift and all sweep values are unchanged.
Layer 15's constant PNC section shortens. Layers 2/14 windows and Layer 16
protection follow the new endpoint. Tapers, tips, tethers, pads and channels
move 13.621373653 um toward the meander. No device is clipped or spliced.

Run the complete scripts in a fresh Python process from Cading:
  python 2611_tapered_yue_dual_resist.py
  python 2611_tapered_yue_design2.py
  python 2611_build_combined_sweep.py
  python 2611_design1_L6000_P10_width_sweep.py
  python validate_2611_combined_sweep.py
  python validate_2611_pnc5um.py

Use GDSFactory 9.45.0 / kfactory 2.6.1 / NumPy 2.4.6; KLayout 0.30.12
for the first three outputs and 0.30.13 for the separate 6000 um output.
Each generator writes its previous gds/ filename. Back up before rebuilding.
Old validate_2611_taper_pitch.py targets the retired 3 um cutoff and is not
the validator for this revision. Use validate_2611_pnc5um.py, which measures
the saved polygons. Report: gds/2611_pnc5um_validation.json.
Validation concerns geometry only; it does not demonstrate optical loss
or mechanical performance and does not replace fabrication DRC.
'''
(ROOT/'PNC5UM_README.txt').write_text(note)
for name in ('CAD_EDITING_GUIDE.md','COMBINED_SWEEP_README.md','DESIGN1_L6000_WIDTH_SWEEP_README.md'):
    p=ROOT/name
    p.write_text('> CURRENT REVISION: PNC terminates at 5 um width, at 31.378626347 um on the reference profile. See `PNC5UM_README.txt` and `validate_2611_pnc5um.py`. Earlier 45 um active-length / 3 um cutoff descriptions below are historical.\n\n'+p.read_text())
print('Guarded 5 um cutoff migration complete.')
