"""Guarded one-time migration; approved device-end geometry is untouched."""
from pathlib import Path
import ast
import hashlib

ROOT = Path(__file__).resolve().parent
EXPECTED = {
    '2611_tapered_yue_dual_resist.py': 'a86419b0e02e9e0d1ccd77ad3213a5459c4b399b',
    '2611_tapered_yue_design2.py': 'f0ca088a9dfab43c8cae42ca56645ace0eb9351d',
}

def replace_one(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f'Expected one migration target: {old[:120]!r}; found {text.count(old)}')
    return text.replace(old, new, 1)

for name, expected_sha in EXPECTED.items():
    path = ROOT / name
    raw = path.read_bytes()
    assert hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest() == expected_sha
    text = original = raw.decode()
    text = replace_one(text,
        'Only the middle straight length varies per entry; other settings are shared.',
        'Each entry is the COMMON length of all THREE straights; other settings are shared.\n'
        'Meander centerline length = 3 * L + 2 * actual_180_degree_Euler_bend_length.\n'
        'Complete end assemblies move with the ports; their local geometry is unchanged.')
    text = replace_one(text, '# middle straight lengths; bottom to top',
                       '# EACH of three equal straights; bottom to top')
    text = replace_one(text,
        'S_BEND_HALF_LENGTH = 150.0  # half the horizontal separation of the two ports',
        '# Horizontal half-port separation is derived as L / 2 in make_device().\n'
        '# Do not set a fixed S_BEND_HALF_LENGTH: it would make the outer leads unequal.')
    text = replace_one(text,
        '    Only meander_straight_length varies. Both taper/coupler assemblies,',
        '    All THREE straight segments have length meander_straight_length.\n'
        '    The two port X coordinates are -L/2 and +L/2. Both taper/coupler assemblies,')
    text = replace_one(text,
        '    if not np.isfinite(meander_straight_length) or meander_straight_length < 0:\n'
        '        raise ValueError("meander_straight_length must be finite and non-negative")',
        '    if not np.isfinite(meander_straight_length) or meander_straight_length <= 0:\n'
        '        raise ValueError("meander_straight_length must be finite and positive")\n'
        '    if TAPER_ANGLE != 0.0:\n'
        '        raise ValueError("Equal-straight layout requires TAPER_ANGLE = 0; rotate a complete device instead")\n'
        '    # The native router sets each outer lead to (port_dx + middle_L) / 2.\n'
        '    # Setting port_dx = L makes both outer leads exactly L as well.\n'
        '    meander_half_length = meander_straight_length / 2.0')
    text = replace_one(text, '        s_bend_half_length=S_BEND_HALF_LENGTH,',
                       '        s_bend_half_length=meander_half_length,')
    text = replace_one(text,
        '    total_length = 2 * S_BEND_HALF_LENGTH + 2 * meander_straight_length + 2 * euler_length',
        '    # Centerline between the two meander ports only; excludes end tapers/tethers.\n'
        '    # Euler length is measured from the actual path, not approximated by pi*R.\n'
        '    total_length = 3 * meander_straight_length + 2 * euler_length')
    text = replace_one(text,
        '    component.info["meander_straight_length_um"] = meander_straight_length',
        '    component.info["meander_straight_length_um"] = meander_straight_length\n'
        '    component.info["meander_straight_segment_lengths_um"] = [meander_straight_length] * 3\n'
        '    component.info["meander_straight_total_um"] = 3 * meander_straight_length\n'
        '    component.info["meander_half_port_separation_um"] = meander_half_length\n'
        '    component.info["meander_180_bend_length_um"] = euler_length\n'
        '    component.info["meander_length_formula"] = "3 * L + 2 * Euler_180_length"')
    text = replace_one(text,
        '        if not np.isfinite(float(length)) or float(length) < 0:\n'
        '            raise ValueError("Meander straight lengths must be finite and non-negative")',
        '        if not np.isfinite(float(length)) or float(length) <= 0:\n'
        '            raise ValueError("Meander straight lengths must be finite and positive")')
    a, b = ast.parse(original), ast.parse(text)
    fa = {n.name: ast.dump(n, include_attributes=False) for n in a.body if isinstance(n, ast.FunctionDef)}
    fb = {n.name: ast.dump(n, include_attributes=False) for n in b.body if isinstance(n, ast.FunctionDef)}
    for key in fa:
        if key not in ('make_device', 'validate_settings'):
            assert fa[key] == fb[key], (name, key)
    path.write_text(text)
    print('UPDATED', name)

p = ROOT / 'CAD_EDITING_GUIDE.md'
s = p.read_text()
s = replace_one(s,
    'Both generators now use the same section numbers. The cleanup preserves the\napproved CAD geometry; it does not redesign either device.',
    'Both generators use the same section numbers. The latest meander update uses\n'
    'three equal straight sections joined by two unchanged 180-degree Euler bends.')
s = replace_one(s,
    'approved one-device 800 um setting. Each entry controls the middle straight\n'
    'section, not total optical path length. The outer straight leads also adjust\n'
    'to maintain the two connection ports. Other Section-01 settings are shared by',
    'one-device 800 um setting. Each entry controls ALL THREE straight sections,\n'
    'not the total optical path length. With L=800, the straights are 800/800/800 um.\n'
    'The horizontal port coordinates are -L/2 and +L/2; complete end assemblies\n'
    'follow those ports without any change to their local geometry.\n'
    'Other Section-01 settings are shared by')
s = s.replace('middle-straight sweep control.', 'equal-straight sweep control.')
marker = '### 1E: Tethers and support islands'
addition = '''### Equal-straight length calculation

In both scripts, `make_device(L)` now draws:

```text
straight L -> 180-degree Euler bend -> straight L
           -> opposite 180-degree Euler bend -> straight L
```

The centerline length between the two meander ports is
`L_total = 3 * L + 2 * L_bend180`. The code measures `L_bend180` with
`gf.path.euler(radius=S_BEND_COMPACT_180_RADIUS, angle=180,
p=S_BEND_EULER_P, use_eff=True).length()`, rather than using the circular
approximation `pi*R`. The device ML label and metadata use this same formula.
PNC/optical tapers outside the meander ports and tethers are NOT included.

For the unchanged radius and Euler fraction, the bend contribution is constant;
increasing L by 100 um increases the total meander length by 300 um.
The left bend alignment and Y pitch are retained. End assemblies move along X
as needed; they are not clipped, stretched or copied between devices.
`S_BEND_HALF_LENGTH` is no longer an editable constant: Section 09 calculates
its per-device replacement as L/2. L must be positive. Keep `TAPER_ANGLE=0`;
rotate an entire device instance when a rotated device is needed.

'''
s = replace_one(s, marker, addition + marker)
s = s.replace('## 4. Check changes\n',
    '## 4. Check changes\n\n'
    'Run `python validate_2611_equal_straights.py` for the current equal-straight\n'
    'geometry, complete-end preservation, connected meander/masks, exact length\n'
    'labels, export and placement checks in both designs. This regenerates both\n'
    'production GDS files. The original pre-change comparison is available when\n'
    'the matching baseline is retained in `validation_2611_equal_straights/baseline/`.\n\n')
p.write_text(s)

p = ROOT / 'validate_2611_design2_stages.py'
s = p.read_text()
s = replace_one(s, '    roi=box_region((150.001,-1000,10000,-25),dbu)',
    "    half = float(d.info['meander_half_port_separation_um'])\n"
    '    roi=box_region((half+0.001,-1000,10000,-25),dbu)')
s = replace_one(s,
    '    taper_roi=box_region((195,-205,325,-195),dbu)+box_region((-325,195,-195,205),dbu)',
    '    taper_roi=box_region((half+45,-205,half+175,-195),dbu)+box_region((-half-175,195,-half-45,205),dbu)')
start = s.index("    before=OUT/'before.gds'")
end = s.index('    top=m.build_layout()', start)
s = s[:start] + """    if (OUT/'before.gds').exists():
        raise ValueError('Archive the historical taper-only before.gds; use validate_2611_equal_straights.py for moved-port comparisons')

""" + s[end:]
p.write_text(s)

p = ROOT / 'DESIGN2_README.md'
s = p.read_text()
s = s.replace('800 um middle meander straight', '800 um length for EACH of the three meander straights')
s += '\nEqual-straight meander: total centerline = 3*L + 2*Euler180_length.\nEnd assemblies follow the ports at x=+/-L/2; all taper widths, stages and offset\nare unchanged. See CAD_EDITING_GUIDE.md and validate_2611_equal_straights.py.\n'
p.write_text(s)
