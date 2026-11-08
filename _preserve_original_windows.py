"""One-use L14-only transplant onto exact approved GDS records."""
from pathlib import Path
import json,struct,math
import klayout.db as kdb
from shapely.geometry import Polygon,box
from shapely.affinity import affine_transform
from shapely.ops import unary_union
from original_design_layer14 import make_terminal_window

FILES=('2611_tapered_yue_dual_resist_L10_L12.gds','2611_tapered_yue_design2.gds','2611_design1_design2_L1400_18device_sweep.gds')
def shapes(r):
 return unary_union([Polygon([(q.x*.001,q.y*.001) for q in p.each_point_hull()],
    [[(q.x*.001,q.y*.001) for q in p.each_point_hole(i)] for i in range(p.holes())]) for p in r.each()])
def records(data):
 pos=0
 while pos<len(data):
  n,t,v=struct.unpack('>HBB',data[pos:pos+4]);yield t,data[pos:pos+n];pos+=n
 assert pos==len(data)
def elements_by_cell(data,only_layer=14):
 name=None;elem=None;layer=None;out={}
 for t,b in records(data):
  if t==6:name=b[4:].rstrip(b'\0').decode();out[name]=[]
  if t in (8,9,10,11,12):elem=[];layer=None
  if elem is not None:elem.append(b)
  if t==13:layer=struct.unpack('>h',b[4:])[0]
  if t==17 and elem is not None:
   if layer==only_layer:out[name].append(b''.join(elem))
   elem=None
 return out
def transplant(old,new14):
 output=[];name=None;elem=None;layer=None
 for t,b in records(old):
  if t==6:name=b[4:].rstrip(b'\0').decode()
  if t in (8,9,10,11,12):elem=[];layer=None
  if elem is not None:
   elem.append(b)
   if t==13:layer=struct.unpack('>h',b[4:])[0]
   if t==17:
    if layer!=14:output.extend(elem)
    elem=None
  else:
   if t==7:output.extend(new14[name])
   output.append(b)
 return b''.join(output)
manifest=json.loads(Path('gds',FILES[2]).with_suffix('.json').read_text())
cases={e['case_id']:e for e in manifest['devices']}
for fn in FILES:
 original=Path('baseline',fn);target=Path('gds',fn)
 ly=kdb.Layout();ly.read(str(original));assert abs(ly.dbu-.001)<1e-12
 l10=ly.layer(10,0);l11=ly.layer(11,0);l14=ly.layer(14,0)
 for cell in ly.each_cell():
  if cell.shapes(l11).is_empty() or cell.shapes(l14).is_empty():continue
  assert not list(cell.each_inst())
  stages=[40.,40.] if fn==FILES[0] else [20.,100.,10.] if fn==FILES[1] else cases[cell.name]['stage_lengths_um']
  original14=kdb.Region(cell.shapes(l14)).merged();result=original14.dup()
  geom=shapes(original14);slab=shapes(kdb.Region(cell.shapes(l10)).merged())
  for sign,cy in [(-1,200.),(1,-200.)]:
   sx=sign*745.;tr=[sign,0,0,sign,-sign*sx,-sign*cy];back=[sign,0,0,sign,sx,cy];roi=box(0,-175,20000,175)
   old=affine_transform(geom,tr).intersection(roi);s=affine_transform(slab,tr).intersection(roi)
   new,meta=make_terminal_window(old,s,stages)
   mapped=affine_transform(new,back);limits=affine_transform(roi,back).bounds
   cut=kdb.Region(kdb.Box(*(round(v*1000) for v in limits)))
   add=kdb.Region(kdb.Polygon([kdb.Point(round(x*1000),round(y*1000)) for x,y in list(mapped.exterior.coords)[:-1]]))
   result=((result-cut)|add).merged()
  cell.shapes(l14).clear();cell.shapes(l14).insert(result)
 temp=target.with_suffix('.candidate.gds');options=kdb.SaveLayoutOptions();options.gds2_max_vertex_count=2000
 ly.write(str(temp),options)
 target.write_bytes(transplant(original.read_bytes(),elements_by_cell(temp.read_bytes())))
 temp.unlink()
 print('Retained every non-L14 GDS record:',fn,flush=True)
