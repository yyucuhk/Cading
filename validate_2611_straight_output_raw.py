"""Independent default-six-device GDSII audit; no generator or KLayout import."""
import json, struct, sys
from pathlib import Path
from shapely.geometry import Polygon, box, LineString
from shapely.ops import unary_union
from shapely import make_valid


def read(path):
    data=Path(path).read_bytes();pos=0;cells={};cell=None;e=None
    while pos<len(data):
        n,t,v=struct.unpack('>HBB',data[pos:pos+4]);b=data[pos+4:pos+n];pos+=n
        if t==6:
            cell=b.rstrip(b'\0').decode();cells[cell]={'polys':{},'refs':[]}
        elif t in (8,9,10,11,12): e={'t':t}
        elif t==13: e['layer']=struct.unpack('>h',b)[0]
        elif t==14: e['datatype']=struct.unpack('>h',b)[0]
        elif t==18: e['ref']=b.rstrip(b'\0').decode()
        elif t==16:
            vals=struct.unpack('>'+str(len(b)//4)+'i',b);e['xy']=list(zip(vals[::2],vals[1::2]))
        elif t==26:
            assert struct.unpack('>H',b)[0]==0, 'Unexpected nontranslation'
        elif t==28:
            assert not any(b), 'Unexpected rotation'
        elif t==17 and e:
            if e['t']==8: cells[cell]['polys'].setdefault(e['layer'],[]).append(e['xy'])
            elif e['t']==10: cells[cell]['refs'].append((e['ref'],e['xy'][0]))
            else: raise ValueError('Unexpected record element '+str(e['t']))
            e=None
    referenced={n for c in cells.values() for n,_ in c['refs']}
    top=[n for n in cells if n not in referenced]
    assert len(top)==1
    return cells,top[0]


def union(c,l):
    return unary_union([make_valid(Polygon([(x*.001,y*.001) for x,y in xy])) for xy in c['polys'].get(l,[])])


def main(root):
    from validate_2611_coupler_release import validate
    validate(Path(root))


if __name__=='__main__':
    main(sys.argv[1] if len(sys.argv)>1 else Path(__file__).resolve().parent)
