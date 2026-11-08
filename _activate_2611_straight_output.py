from pathlib import Path
p=Path('2611_design1_L6000_straight_output.py')
s=p.read_text()
a='    base = load_base()\n    m = base.configured_design(width)'
b='    gf.gpdk.PDK.activate()\n    base = load_base()\n    m = base.configured_design(width)'
assert s.count(a)==1
s=s.replace(a,b)
p.write_text(s)
compile(s,str(p),'exec')
