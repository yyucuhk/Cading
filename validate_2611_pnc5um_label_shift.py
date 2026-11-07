"""Independently verify standalone Design-2 label movement versus physical masks.

The standalone label follows its taper tip. The combined layout has labels
outside devices and does not need this correction. Run with --baseline-commit
when the old commit is available to prove the unlabelled core is unchanged.
"""
from pathlib import Path
import argparse
import ast
import json
import subprocess
import sys
import types

ROOT=Path(__file__).resolve().parent
WORK=ROOT/'validation_2611_pnc5um'
SOURCE='2611_tapered_yue_design2.py'

def exact_standalone_text(length, pnc_end, first_stage, dbu):
    import gdsfactory as gf
    import klayout.db as kdb
    gf.gpdk.PDK.activate()
    total=3*length+2*float(gf.path.euler(radius=100,angle=180,p=.2,use_eff=True).length())
    c=gf.components.text(text=f'DESIGN 2_L{length:g}_ML{total:.1f}',size=3,layer=(2,0))
    r=kdb.Region(c.kdb_cell.begin_shapes_rec(gf.get_layer((2,0)))).merged()
    return r.transformed(kdb.Trans(round((-length/2-pnc_end-first_stage+5)/dbu),round(210/dbu)))

def patch_validator():
    p=ROOT/'validate_2611_pnc5um.py';s=p.read_text()
    old='''                    delta=((regs[n]^region(old_ly,old.cell,n))&core).area()*dbu**2
                    check(name+f'_L{n}_meander_preserved',delta<.005,xor_um2=delta)'''
    new='''                    old_region=region(old_ly,old.cell,n)
                    change=(regs[n]^old_region)&core
                    raw_delta=change.area()*dbu**2
                    delta=raw_delta
                    if tag=='design2' and n==2:
                        from validate_2611_pnc5um_label_shift import exact_standalone_text
                        old_text=exact_standalone_text(L,45.,first_stage,dbu)
                        new_text=exact_standalone_text(L,endpoint,first_stage,dbu)
                        # Compare exact glyph regions, not a rectangular exclusion.
                        # Restoring BOTH known labels must recover identical core masks.
                        delta=(((regs[n]|old_text)^(old_region|new_text))&core).area()*dbu**2
                        outside=(change-(old_text|new_text)).area()*dbu**2
                        check(name+'_L2_difference_is_exact_label_translation',outside<.005,
                              outside_glyph_xor_um2=outside,raw_xor_um2=raw_delta)
                    check(name+f'_L{n}_meander_preserved',delta<.005,
                          xor_um2=delta,raw_xor_um2=raw_delta)'''
    if old not in s:raise ValueError('Validator patch anchor not found')
    p.write_text(s.replace(old,new))

def build_isolated(tag, source):
    import gdsfactory as gf
    tree=ast.parse(source,filename=SOURCE)
    count=0
    for f in tree.body:
        if isinstance(f,ast.FunctionDef) and f.name=='make_device':
            for c in ast.walk(f):
                if isinstance(c,ast.Call) and isinstance(c.func,ast.Attribute) and c.func.attr=='text':
                    for k in c.keywords:
                        if k.arg=='layer' and isinstance(k.value,ast.Name) and k.value.id=='MARKER_LAYER':
                            k.value=ast.parse('(998,0)',mode='eval').body;count+=1
    if count!=1:raise ValueError('Expected exactly one final standalone text call')
    ast.fix_missing_locations(tree)
    module=types.ModuleType('label_audit_'+tag);module.__file__=str(ROOT/SOURCE)
    sys.modules[module.__name__]=module
    exec(compile(tree,str(ROOT/SOURCE),'exec'),module.__dict__)
    d=module.make_device(1400.);d.flatten()
    d.write_gds(WORK/f'label_audit_{tag}.gds',with_metadata=False)


def audit(base_commit):
    import klayout.db as kdb
    WORK.mkdir(exist_ok=True)
    old_source=subprocess.check_output(['git','show',f'{base_commit}:{SOURCE}'],text=True,cwd=ROOT)
    for tag,source in [('before',old_source),('after',(ROOT/SOURCE).read_text())]:
        p=WORK/f'label_audit_{tag}.py';p.write_text(source)
        subprocess.run([sys.executable,str(Path(__file__).resolve()),'--build',tag,'--source',str(p)],check=True,cwd=ROOT)
    def read(path):
        ly=kdb.Layout();ly.read(str(path));return ly,ly.top_cell()
    def r(ly,c,n):return kdb.Region(c.begin_shapes_rec(ly.layer(n,0))).merged()
    a,ca=read(WORK/'label_audit_before.gds');b,cb=read(WORK/'label_audit_after.gds');dbu=a.dbu
    def box(v):return kdb.Region(kdb.Box(*(round(x/dbu) for x in v)))
    core=box((-699.5,-240,699.5,240))+box((699.7,-25,920,225))+box((-920,-225,-699.7,25))
    physical_delta=((r(a,ca,2)^r(b,cb,2))&core).area()*dbu**2
    checks=[dict(name='label_free_physical_L2_meander_unchanged',passed=physical_delta<.005,xor_um2=physical_delta)]
    filename='2611_tapered_yue_design2.gds'
    for tag,ly,c,path in [('before',a,ca,WORK/'baseline'/filename),('after',b,cb,ROOT/'gds'/filename)]:
        out,top=read(path)
        devices=[i.cell for i in top.each_inst() if not r(out,i.cell,51).is_empty() and not r(out,i.cell,15).is_empty()]
        if len(devices)!=1:raise ValueError('Expected one exported standalone device')
        expected=(r(ly,c,2)|r(ly,c,998)).merged()
        delta=(expected^r(out,devices[0],2)).area()*dbu**2
        checks.append(dict(name=tag+'_export_equals_physical_mask_union_exact_label',passed=delta<.005,xor_um2=delta))
    result={'status':'PASS' if all(x['passed'] for x in checks) else 'FAIL','checks':checks,'check_count':len(checks),'method':'Two fresh processes; isolate final label on temporary layer 998 before flattening; compare physical masks independently. Production GDS files are read only.'}
    (WORK/'label_shift_validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    if result['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--patch-validator',action='store_true')
    ap.add_argument('--baseline-commit')
    ap.add_argument('--build');ap.add_argument('--source')
    args=ap.parse_args()
    if args.patch_validator:patch_validator()
    elif args.build:build_isolated(args.build,Path(args.source).read_text())
    elif args.baseline_commit:audit(args.baseline_commit)
    else:ap.error('Specify --baseline-commit for the independent physical-mask audit')
