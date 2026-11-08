"""Restore the approved full-length taper law; run only on the verified base."""
from pathlib import Path
import ast, json, subprocess
ROOT=Path.cwd()
ORIGINAL='6d4fcb56e649e80e0270c3af37be16a102799c00'
WORK=ROOT/'validation_2611_original_profile'
WORK.mkdir(exist_ok=True)
audit={}
for fn in ('2611_tapered_yue_dual_resist.py','2611_tapered_yue_design2.py'):
    old=subprocess.check_output(['git','show',f'{ORIGINAL}:{fn}'],text=True)
    current=(ROOT/fn).read_text()
    assert old.count('PNC_TAPER_SUPPORT_CUTOFF_WIDTH = 3.0')==1
    new=old.replace('PNC_TAPER_SUPPORT_CUTOFF_WIDTH = 3.0','PNC_TAPER_SUPPORT_CUTOFF_WIDTH = 5.0')
    for a,b in [('15 -> 3','15 -> 5'),('15->3','15->5'),('3 -> 1','5 -> 1'),('3->1','5->1'),('3 um','5 um')]:
        new=new.replace(a,b)
    new=new.replace('# PNC taper pitch follows the straight sections; its physical length stays fixed.',
                    '# Full original quintic 15 -> 5 um transition: no truncation or Hermite matching.\n# PNC taper pitch follows the straight sections; its physical length stays fixed.')
    def functions(txt):
        return {f.name:ast.dump(f,include_attributes=False) for f in ast.parse(txt).body if isinstance(f,ast.FunctionDef)}
    def bodies(txt):
        t=ast.parse(txt)
        for n in ast.walk(t):
            if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and n.body and isinstance(n.body[0],ast.Expr) and isinstance(n.body[0].value,ast.Constant) and isinstance(n.body[0].value.value,str):
                n.body=n.body[1:]
        return {f.name:ast.dump(f,include_attributes=False) for f in t.body if isinstance(f,ast.FunctionDef)}
    assert bodies(new)==bodies(old), 'An original drawing function changed unexpectedly'
    a,b=functions(current),functions(new)
    audit[fn]={'original_commit':ORIGINAL,'all_original_function_bodies_preserved':True,
        'removed_helpers':sorted(set(a)-set(b)),
        'functions_restored':[k for k in b if b[k]!=a.get(k)]}
    compile(new,fn,'exec');(ROOT/fn).write_text(new)
for fn in ('2611_build_combined_sweep.py','2611_design1_L6000_P10_width_sweep.py'):
    p=ROOT/fn;s=p.read_text()
    anchor='0 < m.PNC_TAPER_LENGTH < m.PNC_TAPER_TRANSITION_LENGTH'
    assert s.count(anchor)==1
    s=s.replace(anchor,'m.PNC_TAPER_LENGTH == m.PNC_TAPER_TRANSITION_LENGTH == 45.0')
    s=s.replace('list(m.pnc_endpoint_width_derivatives())','[0.0, 0.0]')
    s=s.replace('pnc_junction_continuity="C2",','pnc_junction_continuity="C2",\n                taper_profile="original_quintic_smootherstep",')
    s=s.replace('after the shortened PNC.','after the full 45 um PNC transition.')
    s=s.replace('The PNC ends at the 5 um width point of the original 45 um reference profile.',
                'The PNC uses the original full-length 45 um quintic, now narrowing 15 -> 5 um.')
    compile(s,fn,'exec');p.write_text(s)
p=ROOT/'validate_2611_pnc5um.py';s=p.read_text();a=s.index('def cutoff():');b=s.index('\n\ndef run():',a)
s=s[:a]+'def cutoff():\n    # Full original profile now ends at 5 um; no interior-width root.\n    return 45.0\n'+s[b:];p.write_text(s)
p=ROOT/'validate_2611_smooth_junction.py'
p.write_text('"""Compatibility entry point: original-profile checks supersede Hermite checks."""\nfrom validate_2611_original_profile import run\n\nif __name__ == "__main__":\n    run()\n')
notice='''# Current revision: original full-length taper profiles

The original quintic smootherstep construction is restored. The PNC envelope
is 15 -> 5 um over the full 45 um; no truncation at 31.379 um and no
nonzero-slope Hermite matching remains. Design 1 then uses 5 -> 1 -> 0.2 um.
Design 2 retains 5 -> 1.65 -> 0.980 -> 0.200 um. Existing post-PNC stage
lengths and ratios are unchanged. Both sides of the 5 um join have zero
first and second width derivatives. The final five 1 um reinforcement
periods occupy stations 40,41,42,43,44,45 um. The downstream assemblies
move outward by 13.621373653 um relative to the truncated version.
Meanders, PNC pitches outside the transition, widths, labels and device
origins are otherwise retained. See ORIGINAL_TAPER_PROFILE_README.txt
and gds/2611_original_profile_validation.json for current verification.
Optical loss has not been simulated for these revised widths.

---

'''
for fn in ('CAD_EDITING_GUIDE.md','COMBINED_SWEEP_README.md','DESIGN1_L6000_WIDTH_SWEEP_README.md'):
    p=ROOT/fn
    if p.exists():p.write_text(notice+p.read_text())
for fn in ('PNC5UM_README.txt','SMOOTH_5UM_JUNCTION_README.md'):
    p=ROOT/fn
    if p.exists():p.write_text('SUPERSEDED: see ORIGINAL_TAPER_PROFILE_README.txt for the current full-length profiles.\nThe contents below describe the historical truncated revision.\n\n'+p.read_text())
(ROOT/'ORIGINAL_TAPER_PROFILE_README.txt').write_text('''ORIGINAL TAPER LAW RESTORED -- 15 -> 5 um, then original downstream stages

All dimensions are micrometres. The two design source files again use the
original geometry functions from commit 6d4fcb56, with the PNC endpoint
width changed from 3 to 5 um. The full PNC transition is 45 um long.

For each transition, u=x/L and S(u)=10u^3-15u^4+6u^5:
    w(x) = w_in + (w_out-w_in)*S(u).
The width, slope and curvature are continuous at the 5 um join because
S'(0)=S'(1)=S''(0)=S''(1)=0. The width narrows monotonically. This is the
original taper family, not the recent nonzero-slope Hermite construction.

DESIGN 1: 15 -> 5 over 45; 5 -> 1 over 40; 1 -> 0.2 over 40.
Layer 15 retains 1 -> 0.1 over the first 40 um post-PNC stage.
DESIGN 2: 15 -> 5 over 45; 5 -> 1.65 over 20; 1.65 -> 0.980 over 100;
0.980 -> 0.200 over 10. Layer-15 stages and 0.840 um offset are retained.
COMBINED SWEEP: post-PNC totals 80/100/120, allocated by original ratios.
LONG WIDTH SWEEP: three 6000 um straights, widths 1.62/1.64/1.68 and 10 um PNC.

PNC reinforcement remains five 1 um intervals, with paired supports at
40,41,42,43,44,45 um at both ends. The normal PNC grid is preserved before
this group. Associated L2/L14 openings and L16 protection are regenerated.

RUN from Cading in fresh Python processes:
    python 2611_tapered_yue_dual_resist.py
    python 2611_tapered_yue_design2.py
    python 2611_build_combined_sweep.py
    python 2611_design1_L6000_P10_width_sweep.py
    python validate_2611_original_profile.py

The usual files under gds/ are replaced. Original dependent modules and
tether_variants.gds must remain in this directory. Use GDSFactory 9.45.0,
kfactory 2.6.1; KLayout 0.30.12 for the short layouts and 0.30.13 for the
long-width layout to reproduce prior mask-grid handling.

The analytical envelopes are C2; GDS uses polygons on a 1 nm grid.
Restoring the original profile does NOT independently verify optical
adiabaticity or insertion loss for the new widths. No optical or mechanical
performance certification is implied by the geometric checks.
''')
(WORK/'source_scope_audit.json').write_text(json.dumps(audit,indent=2))
print('Original geometry functions restored; endpoint changed to 5 um.')
