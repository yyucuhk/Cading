from pathlib import Path
import ast
import json
import sys

def edit(path,old,new):
    p=Path(path);s=p.read_text()
    assert s.count(old)==1,(path,old[:80],s.count(old))
    p.write_text(s.replace(old,new))

if sys.argv[1]=='prepare':
    edit('2611_design1_L6000_P10_width_sweep.py','box((x,-204,x+DBU,-196))',
         'box((x,-200-(3 if layer==10 else 2),x+DBU,-200+(3 if layer==10 else 2)))')
    edit('validate_2611_combined_sweep.py','round(-204/dbu),ix+1,round(-196/dbu)',
         'round((-200-(3 if n==10 else 2))/dbu),ix+1,round((-200+(3 if n==10 else 2))/dbu)')
    edit('validate_2611_straight_pitch.py',"and entry['pnc_taper_length_um'] == 45.0)",
         "and abs(entry['pnc_taper_length_um'] - 31.379) < 1e-9\n              and entry['pnc_cutoff_width_um'] == 5.0)")
    p=Path('2611_build_combined_sweep.py');s=p.read_text()
    s=s.replace('full Layer-10 taper after the unchanged 45 um PNC',
                'full Layer-10 taper after the PNC cutoff at 5 um width (31.379 um retained transition)')
    s=s.replace('The taper grid starts at the connection port and its 45 um envelope is unchanged.',
                'The original 45 um reference profile is truncated at 5 um width; its retained length is 31.379 um.')
    p.write_text(s)
    p=Path('validate_2611_pnc5_cutoff.py');s=p.read_text()
    s=s.replace("            c=next(ax._get_lines._cycler_items)['color'] if False else None\n",'')
    s=s.replace("    for n,alpha in ((14,.12),(2,.18),(10,.55),(11,.85),(15,.9)):","    from itertools import cycle\n    colors=cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])\n    for n,alpha in ((14,.12),(2,.18),(10,.55),(11,.85),(15,.9)):")
    s=s.replace("PolyCollection(vertices,alpha=alpha,label=f'Layer {n}')","PolyCollection(vertices,alpha=alpha,facecolors=[next(colors)],label=f'Layer {n}')")
    needle="                prefix=name+'_'+side\n"
    assert s.count(needle)==1
    s=s.replace(needle,needle+'''                if oo is not None:
                    # After the revised first stage, all end structures should
                    # be pure translations, not stretched/clipped/patched.
                    tr=kdb.Trans(round((END-OLD_END)/dbu),0)
                    tail=box((L/2+END+stages[0]+1.,-500,20000,-25),dbu)
                    for n in rr:
                        diff=(rr[n]^oo[n].transformed(tr))&tail
                        check(prefix+f'_L{n}_downstream_shape_preserved',
                              diff.area()*dbu**2<.005,xor_um2=diff.area()*dbu**2)
''')
    p.write_text(s)
    Path('validate_2611_taper_pitch.py').write_text('"""Compatibility entry point for current 5 um PNC-end checks."""\nimport sys\nfrom validate_2611_pnc5_cutoff import capture, run\nif __name__ == "__main__":\n    capture() if len(sys.argv)>1 and sys.argv[1]=="capture" else run()\n')
    for p in Path('.').glob('*.py'):
        ast.parse(p.read_text(),filename=str(p))
else:
    note='''# PNC endpoint at 5 um width - current geometry

Both designs now terminate the PNC-supported transition at a 5 um-wide
central envelope. This does not change the Layer-15 optical ridge to 5 um.
The original 45 um quintic reference profile remains 15 -> 3 um; it is
TRUNCATED where it reaches 5 um, at 31.378626 um, snapped to 31.379 um on
our 1 nm CAD grid. The retained transition is 13.621 um shorter. This is
not a new 15 -> 5 profile stretched over the original full 45 um length.

Five 1 um intervals remain at the narrow PNC endpoint: paired stations
26.379, 27.379, 28.379, 29.379, 30.379, 31.379 um. Earlier regular stations
keep their 5/10 um phase. The last regular-to-dense interval is 1.379 um
for short layouts and 6.379 um for the long sweep; both exceed the 0.626 um
cell extent. No PNC remains beyond this new endpoint except its existing
half-cell longitudinal extent (0.313 um around the station center).

## Corresponding layer changes

Layer 10 starts at this earlier endpoint with width 5 um. Design 1 is now
5 -> 1 -> 0.2 um, with unchanged 40/40 um default lengths. Design 2 is
5 -> 1.65 -> 0.980 -> 0.200 um, with unchanged 20/100/10 um defaults.
The combined 80/100/120 um post-PNC totals and stage ratios are unchanged.
Layer 15 keeps its original widths and stage lengths, including the
Design-2 0.840 um centerline shift. Its constant PNC section ends earlier
and its optical taper follows the new junction.

Layer 2 side bands/release transitions, Layer 14 full opening and Layer 16
ridge protection are rebuilt from the same endpoint. Tethers, support
islands (including Layer 13), Layer-15 alignment tips, long channels and
Layer-51 fiber guides follow the shorter complete end assembly by 13.621 um
toward the meander. Both ends use the same factory and a 180-degree rotation.
No device is clipped or spliced.

Meander straight lengths, bends, straight/bend PNCs, ridge widths, device
counts, 600 um Y pitch and left-bend alignment are unchanged. Existing
label contents and all four production output filenames are unchanged.

## Editing

In Section 01 of either generator:

```python
PNC_TAPER_SUPPORT_CUTOFF_WIDTH = 5.0
PNC_TAPER_TRANSITION_LENGTH = 45.0  # uncut REFERENCE length
PNC_TAPER_PROFILE_END_WIDTH = 3.0  # uncut REFERENCE endpoint width
PNC_CONNECTION_PITCH = 1.0
PNC_CONNECTION_PERIODS = 5
```

`PNC_TAPER_LENGTH` is calculated with `pnc_cutoff_distance`; do not set it
independently. Following stages and masks use the same derived endpoint.
The reusable exported-GDS tests have explicit expectations for this
approved 5 um geometry. Update those expectations for another cutoff.

Run the existing two design generators and two sweep builders to reproduce
the four GDSs. Run `python validate_2611_pnc5_cutoff.py` to check the actual
exported polygons. To compare with a baseline, run the command with
`capture` BEFORE changing and rebuilding geometry. The existing entry
point `validate_2611_taper_pitch.py` now runs the same current checks.

Reproduction uses GDSFactory 9.45.0 with KLayout 0.30.12 for the short
outputs and KLayout 0.30.13 for the 6000 um width sweep, matching their
approved geometry engines. Outputs are measured on the 1 nm grid.

These geometry checks do not establish lower optical loss or mechanical
robustness; optical simulation/measurement and fabrication DRC are separate.
'''
    Path('PNC_5UM_ENDPOINT_README.md').write_text(note)
    for name in ('CAD_EDITING_GUIDE.md','COMBINED_SWEEP_README.md','DESIGN1_L6000_WIDTH_SWEEP_README.md'):
        p=Path(name)
        p.write_text('> CURRENT UPDATE: PNC stops at 5 um width after 31.379 um.\n> See PNC_5UM_ENDPOINT_README.md for current geometry and controls.\n> Older 3 um / full-45 um PNC descriptions below are historical.\n\n'+p.read_text())
    p=Path('gds/2611_pnc5_cutoff_validation.json');r=json.loads(p.read_text())
    r['base_commit']='6d4fcb56e649e80e0270c3af37be16a102799c00'
    r['build_engine_versions']={'gdsfactory':'9.45.0','short_klayout':'0.30.12','long_klayout':'0.30.13'}
    p.write_text(json.dumps(r,indent=2))
    Path('gds/2611_taper_pitch_validation.json').write_text(json.dumps(r,indent=2))
