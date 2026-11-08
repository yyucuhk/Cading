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
    root=Path(root)
    old,_=read(root/'gds/2611_design1_L6000_P10_width_sweep.gds')
    new,top=read(root/'gds/2611_design1_L6000_P10_straight_output_comparison.gds')
    checks=[]
    def check(name,ok,**kw):
        checks.append(dict(name=name,passed=bool(ok),**kw));assert ok,(name,kw)
    refs=new[top]['refs']
    check('six_devices_six_labels_one_frame',len(refs)==13 and sum(n.startswith('LABEL_') for n,_ in refs)==6)
    devrefs={n:p for n,p in refs if n.startswith('D1_L')}
    check('six_complete_device_cells',len(devrefs)==6)
    check('600um_vertical_pitch',sorted(p[1] for p in devrefs.values())==[-700000,-100000,500000,1100000,1700000,2300000])
    roi=box(3000,-400,20000,-25)
    for w in (1.62,1.64,1.68):
        n=f'D1_L6000_W{round(w*1000)}_P10';v=n+'_OUT2000'
        check(n+'_eight_layers',set(new[n]['polys'])=={2,10,11,13,14,15,16,51})
        check(v+'_eight_layers',set(new[v]['polys'])=={2,10,11,13,14,15,16,51})
        for l in (2,10,11,13,14,15,16,51):
            a=union(old[n],l);b=union(new[n],l);c=union(new[v],l)
            ar=a.symmetric_difference(b).area
            check(n+f'_reference_L{l}_unchanged',ar<1e-7,xor_um2=ar)
            ar=b.symmetric_difference(c).difference(roi).area
            check(v+f'_outside_output_L{l}_unchanged',ar<1e-7,xor_um2=ar)
        optical=union(new[v],15)
        expected=box(3000,-200-w/2,5000,-200+w/2)
        ar=optical.intersection(roi).symmetric_difference(expected).area
        check(v+'_exact_full_2000um_constant_ridge',ar<1e-7,xor_um2=ar)
        for x in (2999.99,3000.01,3050,3500,4000,4500,4999.99):
            seg=optical.intersection(LineString([(x,-203),(x,-197)]))
            check(v+f'_width_{x}',abs(seg.length-w)<1e-7,measured_um=seg.length)
        geoms=list(optical.geoms) if hasattr(optical,'geoms') else [optical]
        check(v+'_one_optical_body_two_input_ticks',len(geoms)==3)
        protection=union(new[v],16)
        check(v+'_protection_connected',protection.geom_type=='Polygon')
        body=max(geoms,key=lambda p:p.area)
        check(v+'_ridge_covered_by_protection',body.difference(protection).area<1e-7)
        check(v+'_ridge_covered_by_release',body.difference(union(new[v],14)).area<1e-7)
        for l in (10,13,51):
            check(v+f'_no_output_taper_hardware_L{l}',union(new[v],l).intersection(roi).area==0)
        pnc=union(new[v],11)
        out=pnc.intersection(box(2999,-230,5001,-170))
        pairs={1:[],-1:[]}
        for p in out.geoms:
            xmin,ymin,xmax,ymax=p.bounds
            pairs[1 if (ymin+ymax)/2>-200 else -1].append((xmin+xmax)/2-3000)
        for side,centers in pairs.items():
            centers.sort()
            check(v+f'_200_PNCs_side{side}',len(centers)==200)
            check(v+f'_exact10um_PNC_side{side}',all(abs(x-i*10)<1e-7 for i,x in enumerate(centers)))
        for c in (n,v):
            x=union(new[c],11).bounds[0]+devrefs[c][0]*.001
            check(c+'_left_bend_minus900',abs(x+900)<1e-7,measured_x_um=x)
    result=dict(status='PASS',check_count=len(checks),checks=checks,
                method='Independent raw GDS record parser and Shapely polygons; no generator or KLayout imported')
    out=root/'validation_2611_straight_output/raw_gds_audit.json';out.write_text(json.dumps(result,indent=2))
    print('PASS',len(checks),'raw GDS checks')


if __name__=='__main__':
    main(sys.argv[1] if len(sys.argv)>1 else Path(__file__).resolve().parent)
