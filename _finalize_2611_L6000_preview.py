"""Temporary preview-only correction; production GDS remains byte-identical."""
from pathlib import Path
import ast
import hashlib
import importlib.util
import json

p=Path('2611_design1_L6000_P10_width_sweep.py')
s=p.read_text();before=ast.parse(s)
old='(STRAIGHT_LENGTH/2-20,-235,STRAIGHT_LENGTH/2+170,35)'
assert s.count(old)==1
s=s.replace(old,'(STRAIGHT_LENGTH/2-20,-25,STRAIGHT_LENGTH/2+170,225)')
old='fig,ax=plt.subplots(figsize=(15,5))'
assert s.count(old)==1
s=s.replace(old,"fig,ax=plt.subplots(figsize=(7,8) if name=='bend_detail' else (15,5))")
old="            ax.add_collection(LineCollection(lines,linewidths=.6,label=f'Layer {n}'))"
new="""            flat = [point for line in lines for point in line + [(float('nan'),float('nan'))]]
            if flat:
                ax.plot([v[0] for v in flat],[v[1] for v in flat],linewidth=.6,label=f'Layer {n}')"""
assert s.count(old)==1;s=s.replace(old,new)
s=s.replace('    from matplotlib.collections import LineCollection\n','')
old="        ax.set_aspect('equal');ax.legend(loc='upper right');fig.tight_layout()"
new="""        if name == 'overview':
            for entry in result['devices']:
                x,y=entry['label_origin_um']
                ax.text(x,y+20,entry['label'],fontsize=8)
        ax.set_aspect('equal');ax.legend(loc='upper right');fig.tight_layout()"""
assert s.count(old)==1;s=s.replace(old,new)
after=ast.parse(s)
def without_render(tree):
    tree.body=[n for n in tree.body if not(isinstance(n,ast.FunctionDef) and n.name=='render')]
    return ast.dump(tree,include_attributes=False)
assert without_render(before)==without_render(after),'Non-preview code changed'
p.write_text(s)
gds=Path('gds/2611_design1_L6000_P10_width_sweep.gds')
hash_before=hashlib.sha256(gds.read_bytes()).hexdigest()
spec=importlib.util.spec_from_file_location('long_preview',p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
ly,top=m.read_gds(gds)
report=json.loads(gds.with_suffix('.json').read_text())
m.WORK.mkdir(exist_ok=True)
m.render(ly,top,report)
assert hashlib.sha256(gds.read_bytes()).hexdigest()==hash_before==report['gds_sha256']
print('Preview corrected; all non-preview code AST-identical; GDS unchanged:',hash_before)
