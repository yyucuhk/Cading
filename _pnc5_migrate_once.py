"""One-time guarded source migration; remove after validated output is committed."""
from pathlib import Path
import ast

ROOT = Path.cwd()
def once(s, old, new):
    assert s.count(old) == 1, (old[:120], s.count(old))
    return s.replace(old, new)

helper = '''def pnc_cutoff_distance(length, width1, profile_end_width, cutoff_width):
    """Truncate the original quintic profile, snapping the endpoint to 1 nm.

    A 45 um reference profile from 15 to 3 um reaches 5 um at 31.378626 um.
    Its retained CAD length is 31.379 um. This is not a refitted 15->5 taper
    stretched across the old full length. All following layers use this
    same endpoint; the endpoint transverse width is exactly cutoff_width.
    """
    values = np.asarray([length, width1, profile_end_width, cutoff_width], dtype=float)
    if not np.all(np.isfinite(values)) or length <= 0:
        raise ValueError("PNC profile dimensions must be finite and length positive")
    if not 0 < profile_end_width <= cutoff_width < width1:
        raise ValueError("Require 0 < profile_end_width <= cutoff_width < width1")
    if abs(cutoff_width - profile_end_width) < 1e-12:
        return round(float(length), 3)
    lo, hi = 0.0, 1.0
    for _ in range(64):
        u = (lo + hi) / 2
        w = width1 + (profile_end_width - width1) * (6*u**5 - 15*u**4 + 10*u**3)
        if w > cutoff_width:
            lo = u
        else:
            hi = u
    return round(float(length * (lo + hi) / 2), 3)


PNC_TAPER_LENGTH = pnc_cutoff_distance(
    PNC_TAPER_TRANSITION_LENGTH, SLAB_EXTENSION_END_WIDTH,
    PNC_TAPER_PROFILE_END_WIDTH, PNC_TAPER_SUPPORT_CUTOFF_WIDTH,
)'''
solver='''    # Solve the same quintic width profile used by tapered_fiber_adiabatic.
    # When the taper ends at the cutoff width, the PNC mask reaches the exact
    # taper endpoint and the Layer-15 taper can start at that same coordinate.
    if np.isclose(width2, cutoff_width):
        cutoff_x = length
    else:
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = (lo + hi) / 2
            weight = 6 * mid**5 - 15 * mid**4 + 10 * mid**3
            local_width = width1 + (width2 - width1) * weight
            if local_width > cutoff_width:
                lo = mid
            else:
                hi = mid
        cutoff_x = length * (lo + hi) / 2

    x = np.linspace(0.0, cutoff_x, npoints)'''
for name in ['2611_tapered_yue_dual_resist.py', '2611_tapered_yue_design2.py']:
    p=ROOT/name;s=p.read_text()
    s=once(s,'PNC_TAPER_SUPPORT_CUTOFF_WIDTH = 3.0','PNC_TAPER_SUPPORT_CUTOFF_WIDTH = 5.0')
    s=once(s,'PNC_TAPER_TRANSITION_LENGTH = 45.0  # fixed 15 -> 3 um envelope; independent of PNC pitch',
        'PNC_TAPER_TRANSITION_LENGTH = 45.0  # ORIGINAL reference profile length; actual retained length is derived\nPNC_TAPER_PROFILE_END_WIDTH = 3.0  # reference endpoint; PNC stops earlier at the 5 um cutoff')
    s=once(s,'PNC_TAPER_LENGTH = PNC_TAPER_TRANSITION_LENGTH',helper)
    s=once(s,'    length: float = SLAB_EXTENSION_LENGTH,\n    width1: float = SLAB_EXTENSION_END_WIDTH,\n    width2: float = SLAB_EXTENSION_START_WIDTH,\n    cutoff_width:',
        '    length: float = PNC_TAPER_TRANSITION_LENGTH,\n    width1: float = SLAB_EXTENSION_END_WIDTH,\n    width2: float = PNC_TAPER_PROFILE_END_WIDTH,\n    cutoff_width:')
    s=once(s,solver.replace('used by tapered_fiber_adiabatic.', 'used by design2_tapered_fiber_adiabatic.') if 'design2' in name else solver,'''    # Use the original profile up to the selected cutoff; all following
    # taper stages and masks share this exact grid-snapped endpoint.
    cutoff_x = pnc_cutoff_distance(length, width1, width2, cutoff_width)
    profile_samples = np.linspace(0.0, length, npoints)
    x = np.r_[profile_samples[profile_samples < cutoff_x - 1e-9], cutoff_x]''')
    s=once(s,'    widths = width1 + (width2 - width1) * weight\n\n    tether = phononic_tether(',
        '    widths = width1 + (width2 - width1) * weight\n    widths[-1] = cutoff_width  # exact transverse width at the snapped junction\n\n    tether = phononic_tether(')
    s=once(s,'                length=PNC_TAPER_LENGTH,\n                width1=SLAB_EXTENSION_END_WIDTH,\n                width2=PNC_TAPER_END_WIDTH,',
        '                length=PNC_TAPER_TRANSITION_LENGTH,\n                width1=SLAB_EXTENSION_END_WIDTH,\n                width2=PNC_TAPER_PROFILE_END_WIDTH,')
    s=once(s,'        support_width = width1 + (width2 - width1) * support_weight\n        cuts.add_ref',
        '        support_width = width1 + (width2 - width1) * support_weight\n        if abs(support_x - cutoff_x) < 1e-9:\n            support_width = cutoff_width\n        cuts.add_ref')
    s=once(s,'    support_u = SUPPORT_TIP_POSITION / PNC_TAPER_LENGTH',
        '    support_u = SUPPORT_TIP_POSITION / PNC_TAPER_TRANSITION_LENGTH')
    s=once(s,'SLAB_EXTENSION_END_WIDTH + (PNC_TAPER_END_WIDTH - SLAB_EXTENSION_END_WIDTH) * support_weight',
        'SLAB_EXTENSION_END_WIDTH + (PNC_TAPER_PROFILE_END_WIDTH - SLAB_EXTENSION_END_WIDTH) * support_weight')
    s=s.replace('outward from the 3 um PNC endpoint','outward from the 5 um PNC endpoint')
    s=s.replace('Layer 10: 3 -> 1 over','Layer 10: 5 -> 1 over').replace('Layer 10:       3 -> 1.65','Layer 10:       5 -> 1.65')
    s=s.replace('ending at 3 um width','ending at 5 um width')
    s=s.replace('the 3 um-wide PNC end','the selected narrow PNC end')
    s=s.replace('Five 1 um periods occupy x=40,41,42,43,44,45 um for a 45 um taper.',
                'Five 1 um periods occupy 26.379..31.379 um for the 5 um cutoff.')
    s=s.replace('everywhere before the final dense segment. Never shift optical geometry.',
                'everywhere before the final dense segment. Following masks use this endpoint.')
    s=s.replace('The preceding 45 um PNC transition and its 1 um ridge stay unchanged.',
                'The preceding PNC profile ends at the derived 5 um-wide cutoff.')
    s=s.replace('starts again at the 3 um endpoint, then tapers 3->1 um',
                'starts at the 5 um endpoint, then tapers 5->1 um')
    s=s.replace('Along +x, the 15->3 um PNC transition is excluded from Layer 10.',
                'Along +x, the retained 15->5 um PNC transition is excluded from Layer 10.')
    s=s.replace('    # The first transition is intentionally not drawn on Layer 10 because its\n    # width is greater than 3 um. Keep its longitudinal distance so all PNC,\n    # ridge, and downstream taper coordinates remain unchanged.',
                '    # Layer 10 starts at the shortened PNC endpoint (5 um wide).\n    # All ridge/taper/tether/channel landmarks follow the derived endpoint;\n    # post-PNC stage lengths and widths after the first stage are unchanged.')
    s=s.replace('3 um Layer-10 starting point','5 um Layer-10 starting point')
    s=s.replace('(the 3->1 um quintic taper)','(the 5->1 um quintic taper)')
    s=s.replace('FULL 15->3 um taper', 'retained 15->5 um taper')
    s=once(s,'    component.info["pnc_taper_length"] = PNC_TAPER_LENGTH',
        '    component.info["pnc_taper_length"] = PNC_TAPER_LENGTH\n    component.info["pnc_cutoff_width_um"] = PNC_TAPER_SUPPORT_CUTOFF_WIDTH\n    component.info["pnc_reference_profile_length_um"] = PNC_TAPER_TRANSITION_LENGTH')
    ast.parse(s);p.write_text(s)

p=ROOT/'2611_build_combined_sweep.py';s=p.read_text()
s=once(s,'m.PNC_TAPER_LENGTH == 45.0 and','abs(m.PNC_TAPER_LENGTH - 31.379) < 1e-9 and')
s=s.replace('Taper PNC pitch must match the straight pitch without changing its 45 um length',
            'Taper pitch must match straight pitch and end at the derived 5 um cutoff')
s=once(s,"widths10 = [3., 1., .2] if case['design'] == 1 else [3., 1.65, .98, .2]",
        "widths10 = [5., 1., .2] if case['design'] == 1 else [5., 1.65, .98, .2]")
s=once(s,'                pnc_taper_length_um=float(m.PNC_TAPER_LENGTH),',
        '                pnc_taper_length_um=float(m.PNC_TAPER_LENGTH),\n                pnc_cutoff_width_um=float(m.PNC_TAPER_SUPPORT_CUTOFF_WIDTH),')
ast.parse(s);p.write_text(s)

p=ROOT/'2611_design1_L6000_P10_width_sweep.py';s=p.read_text()
s=once(s,'for layer,x,w in ((10,45.002,3.),(10,84.998,1.),(10,124.998,.2),\n                          (15,45.002,width),(15,84.998,.1)):',
        'for layer,x,w in ((10,m.PNC_TAPER_LENGTH+.002,5.),\n                          (10,m.PNC_TAPER_LENGTH+39.998,1.),\n                          (10,m.PNC_TAPER_LENGTH+79.998,.2),\n                          (15,m.PNC_TAPER_LENGTH+.002,width),\n                          (15,m.PNC_TAPER_LENGTH+39.998,.1)):')
s=once(s,'box((x,-202,x+DBU,-198))','box((x,-204,x+DBU,-196))')
s=s.replace('taper_PNC_matches_straight_without_changing_length','taper_PNC_matches_straight_with_5um_cutoff')
s=once(s,'m.PNC_TAPER_LENGTH==45 and','abs(m.PNC_TAPER_LENGTH-31.379)<1e-9 and')
s=once(s,'                pnc_taper_length_um=float(m.PNC_TAPER_LENGTH),',
        '                pnc_taper_length_um=float(m.PNC_TAPER_LENGTH),\n                pnc_cutoff_width_um=float(m.PNC_TAPER_SUPPORT_CUTOFF_WIDTH),')
s=s.replace('L10 3 -> 1 -> 0.2 over 40/40; L15 W -> 0.1 over 40 after the 45 um PNC.',
            'L10 5 -> 1 -> 0.2 over 40/40; L15 W -> 0.1 over 40 after the 31.379 um PNC.')
s=s.replace('its length stays 45 um.', 'it now ends at 5 um width (31.379 um from its start).')
ast.parse(s);p.write_text(s)

p=ROOT/'validate_2611_combined_sweep.py';s=p.read_text()
s=once(s,"x=e['straight_um']/2+45.","x=e['straight_um']/2+e['pnc_taper_length_um']")
s=once(s,'round(-202/dbu),ix+1,round(-198/dbu)','round(-204/dbu),ix+1,round(-196/dbu)')
s=s.replace('transverse +/-2 um clip','transverse +/-4 um clip')
ast.parse(s);p.write_text(s)
print('Patched the two generators, two builders and combined exported-GDS validator.')
