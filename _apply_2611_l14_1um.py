from pathlib import Path
import ast
import hashlib

root = Path.cwd()
p = root / '2611_design1_L6000_straight_output.py'
old = p.read_text()
assert hashlib.sha256(p.read_bytes()).hexdigest() == 'aa9c88a1cdf4e780208d1f7f3b760cd73a3519d3469d2c4caa13423958031746'
new = old.replace('LAYER14_WINDOW_MARGIN = 2.5', 'LAYER14_WINDOW_MARGIN = 1.0', 1)
new = new.replace('L14 starts with L10, at TOTAL width 8 um. Both downstream margins are 2.5 um\nper side.', 'L14 starts with L10, at TOTAL width 8 um. Downstream margins are 2.5 um\nper side on L2 and 1.0 um per side on L14.')
new = new.replace('extra clearances: L2=15 um and L14=8 um. Both downstream margins=2.5 um.', 'extra clearances: L2=15 um and L14=8 um. Downstream margins: L2=2.5 um, L14=1 um.')
# Verify that only one executable setting changes in the production builder.
def normalized(source):
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if node.body and isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Constant) and isinstance(node.body[0].value.value, str):
                node.body.pop(0)
        if isinstance(node, ast.Assign) and any(isinstance(t,ast.Name) and t.id=='LAYER14_WINDOW_MARGIN' for t in node.targets):
            node.value = ast.Constant(1.0)
    return ast.dump(tree, include_attributes=False)
assert normalized(new) == normalized(old)
assert new != old
p.write_text(new)

p = root / 'validate_2611_four_sets.py'
s = p.read_text()
assert hashlib.sha256(p.read_bytes()).hexdigest() == '9406974199b73141418549b9a0a03c35dd8d7dacfb36f5313ae1f1837e9b0fc8'
replacements = {
    "meta['layer2_window_margin_um']==meta['layer14_window_margin_um']==2.5": "meta['layer2_window_margin_um']==2.5 and meta['layer14_window_margin_um']==1.0",
    "            if name=='FRAME_AND_MARKERS':continue\n": "",
    "if name in devices and layer in (2,14):continue": "if name in devices and layer == 14:continue",
    "expected_widths={2:15.,14:8.};offsets={2:.5,14:0.}": "expected_widths={2:15.,14:8.};offsets={2:.5,14:0.};margins={2:2.5,14:1.0}",
    "expected=expected_widths[layer]+(6-expected_widths[layer])*f": "expected=expected_widths[layer]+(1+2*margins[layer]-expected_widths[layer])*f",
    "abs(lower-2.5)<.0025 and abs(upper-2.5)<.0025": "abs(lower-margins[layer])<.0025 and abs(upper-margins[layer])<.0025",
    "downstream_margin_per_side_um=2.5": "downstream_margin_per_side_um={'2':2.5,'14':1.0}",
    "downstream margins 2.5 um')": "downstream margins L2=2.5 um / L14=1 um')",
}
for a,b in replacements.items():
    assert s.count(a)==1, (a,s.count(a))
    s=s.replace(a,b,1)
# Check current twelve-device baseline, including exact L2 and frame preservation.
needle="    for name,pos in devices.items():\n"
extra="""    if old:
        check('all_cell_names_preserved', set(cells)==set(old))
        check('all_cell_references_preserved', all(cells[n]['refs']==old[n]['refs'] for n in cells))
        for name in devices:
            before=union(old[name],14)
            after=geom(name,14)
            check(name+'_L14_margin_really_reduced', before.difference(after).area>1.)
            # Offset polygons are quantized on a 1 nm database grid.
            check(name+'_L14_no_widening_beyond_grid', after.difference(before.buffer(DBU)).area<1e-6)
            single=name.startswith('D1_SINGLE_')
            one=('_1T_' in name) if single else name.endswith('_OUT2000')
            allowed=Polygon()
            for sign,cy in [(-1,0 if single else 200)]+([] if one else [(1,0 if single else -200)]):
                sx=sign*3045
                allowed=allowed.union(box(min(sx,sx+sign*96),cy-30,max(sx,sx+sign*96),cy+30))
            exact(name+'_L14_unchanged_outside_taper_window', after.difference(allowed), before.difference(allowed))
"""
assert s.count(needle)==1
s=s.replace(needle,extra+needle)
ast.parse(s)
p.write_text(s)

p = root / 'DESIGN1_STRAIGHT_OUTPUT_README.txt'
s = p.read_text()
a='Entrance widths are TOTAL widths: L2=15 um, L14=8 um. Both downstream\nexpansions are 2.5 um per side. These margins do not add another 5 um\nto the stated entrance widths.'
b='Entrance widths are TOTAL widths: L2=15 um, L14=8 um. Downstream\nexpansions are 2.5 um per side on L2 and 1.0 um per side on L14.\nThese margins are not added again to the stated entrance widths.'
assert s.count(a)==1
p.write_text(s.replace(a,b,1))
print('Changed only the L14 executable margin; updated independent checks and notes.')
