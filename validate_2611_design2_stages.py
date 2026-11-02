"""Polygon checks for the current 20/100/10 um Design 2.

Run python validate_2611_design2.py. An optional before.gds and protected_hashes.json
in validation_2611_design2 enable baseline and unchanged-Design-1 checks.
"""
from __future__ import annotations
from pathlib import Path
from functools import partial
import hashlib
import importlib.util
import json
import traceback

import gdsfactory as gf
import klayout.db as kdb
import numpy as np
from validate_2611_sweep import snapshot, file_snapshot, compare, export, box_region

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'validation_2611_design2'
OUT.mkdir(exist_ok=True)
REPORT = {'checks': [], 'failures': []}
TOL = 0.005


def check(name, passed, **details):
    item = dict(name=name, passed=bool(passed), **details)
    REPORT['checks'].append(item)
    if not passed:
        REPORT['failures'].append(name)
    (OUT / 'validation_report.json').write_text(json.dumps(REPORT, indent=2))
    print(json.dumps(item), flush=True)


def read(path):
    ly = kdb.Layout()
    ly.read(str(path))
    return ly, file_snapshot(ly, ly.top_cell())


def region_polygon(xs, widths, centers, dbu):
    upper = [(x, y + w / 2) for x, w, y in zip(xs, widths, centers)]
    lower = [(x, y - w / 2) for x, w, y in zip(xs[::-1], widths[::-1], centers[::-1])]
    points = [kdb.Point(round(x/dbu), round(y/dbu)) for x, y in upper + lower]
    return kdb.Region(kdb.Polygon(points))


def expected_taper(x0, length, w0, w1, dbu):
    u = np.linspace(0, 1, 201)
    w = w0 + (w1-w0)*(6*u**5-15*u**4+10*u**3)
    return region_polygon(x0 + length*u, w, np.zeros_like(u), dbu)


def section(reg, x, dbu):
    ix = round(x/dbu)
    clip = kdb.Region(kdb.Box(ix, round(-5/dbu), ix+1, round(5/dbu)))
    polys = list((reg & clip).merged().each())
    if not polys:
        return None
    assert len(polys) == 1, ('Disconnected cross-section', x, len(polys))
    b = polys[0].bbox()
    return [(b.bottom+b.top)*dbu/2, (b.top-b.bottom)*dbu]


def actual_device_instances(ly):
    idx = ly.find_layer(51, 0)
    return [i for i in ly.top_cell().each_inst()
            if idx is not None and not kdb.Region(i.cell.begin_shapes_rec(idx)).is_empty()]


def render(regions, dbu, name, bounds, layers, title, origin_x=0, converter=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as MPath
    from matplotlib.patches import PathPatch
    x0, y0, x1, y1 = bounds
    fig, ax = plt.subplots(figsize=(13, 5))
    roi = box_region(bounds, dbu)
    for layer in layers:
        reg = (regions.get((layer,0), kdb.Region()) & roi).merged()
        lines=[]
        for poly in reg.each():
            loops = [list(poly.each_point_hull())] + [list(poly.each_point_hole(i)) for i in range(poly.holes())]
            vertices=[]; codes=[]
            for loop in loops:
                pts=[(p.x*dbu-origin_x,p.y*dbu) for p in loop]
                if not pts:
                    continue
                lines.extend(pts+[pts[0],(np.nan,np.nan)])
                vertices.extend(pts+[pts[0]])
                codes.extend([MPath.MOVETO]+[MPath.LINETO]*(len(pts)-1)+[MPath.CLOSEPOLY])
            if converter and vertices:
                ax.add_patch(PathPatch(MPath(vertices,codes),alpha=0.13))
        if lines:
            p=np.asarray(lines)
            ax.plot(p[:,0],p[:,1],linewidth=1 if converter else 0.6,label=f'Layer {layer}')
    if converter:
        for x in (0,20,120,130):
            ax.axvline(x,linestyle=':',linewidth=0.7)
        ax.text(10,-1.62,'Stage 1: 20 um',ha='center',fontsize=9)
        ax.text(70,-1.62,'Stage 2: 100 um',ha='center',fontsize=9)
        ax.annotate('Layer 10 only\nStage 3: 10 um',xy=(125,0),xytext=(103,-1.35),
                    fontsize=9,arrowprops=dict(arrowstyle='->'))
        ax.annotate('Layer 15 ends\n200 nm at Y = 840 nm',xy=(120,0.84),xytext=(87,1.42),
                    fontsize=9,arrowprops=dict(arrowstyle='->'))
    else:
        ax.set_aspect('equal',adjustable='box')
    ax.set(xlim=(x0-origin_x,x1-origin_x),ylim=(y0,y1),
           xlabel='Distance from the 3 um taper start (um)' if converter else 'Local X (um)',
           ylabel='Local Y (um)',title=title)
    ax.legend(loc='upper right')
    fig.tight_layout()
    fig.savefig(OUT/name,dpi=180)
    plt.close(fig)


def run():
    spec=importlib.util.spec_from_file_location('d2_stage_layout',ROOT/'2611_tapered_yue_design2.py')
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    gf.gpdk.PDK.activate()
    check('single_design2_800um',tuple(m.MEANDER_STRAIGHT_LENGTHS)==(800.0,) and m.DESIGN_ID=='Design 2')
    check('requested_stage_constants',
          (m.POST_PNC_TAPER_LENGTH,m.SLAB_STAGE2_LENGTH,m.SLAB_TAPER_LENGTH)==(20.0,100.0,10.0)
          and (m.RIDGE_TAPER_LENGTH,m.RIDGE_CONVERTER_LENGTH)==(20.0,100.0))
    check('840nm_centerline_shift_retained',m.RIDGE_CONVERTER_Y_SHIFT==0.840)
    check('unchanged_45um_pnc',m.PNC_TAPER_LENGTH==45.0)

    coupler=m.design2_two_layer_tapered_fiber(add_support_tips=False,add_phononic_taper_support=True)
    export(coupler,OUT/'coupler.gds')
    cly,cr=read(OUT/'coupler.gds'); dbu=cly.dbu
    expected10=kdb.Region()
    for x,L,w0,w1 in ((45,20,3,1.65),(65,100,1.65,0.980),(165,10,0.980,0.200)):
        expected10+=expected_taper(x,L,w0,w1,dbu)
    expected15=expected_taper(0,45,1,1,dbu)+expected_taper(45,20,1,0.500,dbu)
    expected15+=region_polygon(np.array([65.,165.]),np.array([.500,.200]),np.array([0.,.840]),dbu)
    for layer,expected in ((10,expected10),(15,expected15)):
        delta=(cr[(layer,0)]^expected).area()*dbu**2
        check(f'complete_layer_{layer}_profile_matches_requested_stages',delta<=TOL,xor_um2=float(delta))
        check(f'coupler_layer_{layer}_all_junctions_connected',len(list(cr[(layer,0)].merged().each()))==1)
    check('layer10_ends_at_130um',abs(cr[(10,0)].bbox().right*dbu-175)<1e-8)
    check('layer15_ends_at_120um',abs(cr[(15,0)].bbox().right*dbu-165)<1e-8)
    no_third=cr[(15,0)] & box_region((165.001,-5,175,5),dbu)
    check('no_layer15_third_stage_or_terminal_extension',no_third.is_empty(),area_um2=float(no_third.area()*dbu**2))
    sampled=[]
    for layer,values in (
        (10,[(45.001,3.,0.),(55.,2.325,0.),(65.001,1.65,0.),(115.,1.315,0.),(164.999,.980,0.),(170.,.590,0.),(174.999,.200,0.)]),
        (15,[(45.001,1.,0.),(55.,.750,0.),(65.001,.500,0.),(115.,.350,.420),(164.999,.200,.840)])):
        for x,width,y in values:
            got=section(cr[(layer,0)],x,dbu)
            passed=got is not None and abs(got[0]-y)<=0.002 and abs(got[1]-width)<=0.002
            check(f'layer{layer}_section_{x:g}',passed,x_from_3um_start=x-45,
                  measured_center_y_um=None if got is None else got[0],
                  measured_width_um=None if got is None else got[1],
                  expected_width_um=width)
            sampled.append(dict(layer=layer,x_um=x-45,center_y_um=got[0],width_um=got[1]))
    v10=section(cr[(10,0)],164.999,dbu);v15=section(cr[(15,0)],164.999,dbu)
    gap=v15[0]-v15[1]/2-v10[1]/2
    check('terminal_vertical_edge_gap_250nm',abs(gap-.250)<=.002,gap_um=float(gap))

    d=m.make_device(800.0);r=snapshot(d);dbu=d.kcl.dbu
    expected_layers={(2,0),(10,0),(11,0),(13,0),(14,0),(15,0),(16,0),(51,0)}
    check('all_device_layers_present',set(r)==expected_layers,layers=[list(l) for l in sorted(r)])
    check('metadata_has_three_slab_two_ridge_stages',
          list(d.info['layer10_stage_lengths_um'])==[20.,100.,10.] and
          list(d.info['layer15_stage_lengths_um'])==[20.,100.])
    ticks=m.alignment_tip_region(d.info['tip_centers'],d.info['connection_orientations'],dbu)
    optical=(r[(15,0)]-ticks).merged()
    check('four_ticks_retained_on_layer15',len(list(ticks.each()))==4 and (ticks-r[(15,0)]).is_empty())
    check('no_ticks_remain_on_layer10',(ticks&r[(10,0)]).is_empty())
    check('optical_ridge_and_both_converters_connected',len(list(optical.each()))==1)
    check('layer16_protection_connected',len(list(r[(16,0)].merged().each()))==1)
    for mask in (14,16):
        missing=(optical-r[(mask,0)]).area()*dbu**2
        check(f'layer{mask}_covers_complete_ridge',missing<=TOL,missing_um2=float(missing))
    label=box_region(d.info['label_bbox_um'],dbu)
    missing=((r[(2,0)]-r[(14,0)])-label).area()*dbu**2
    check('layer14_covers_all_layer2_windows',missing<=TOL,missing_um2=float(missing))
    roi=box_region((150.001,-1000,10000,-25),dbu)
    for layer in sorted(expected_layers):
        reg=r[layer]-label if layer==(2,0) else r[layer]
        diff=((reg^reg.transformed(kdb.Trans(kdb.Trans.R180)))&roi).area()*dbu**2
        check(f'layer{layer[0]}_end_to_end_symmetry',diff<=TOL,rotated_xor_um2=float(diff))
    # Check for residual center holes only in the taper areas, not support pads.
    taper_roi=box_region((195,-205,325,-195),dbu)+box_region((-325,195,-195,205),dbu)
    holes=kdb.Region()
    for poly in r[(14,0)].merged().each():
        for i in range(poly.holes()):
            holes.insert(kdb.Polygon(list(poly.each_point_hole(i))))
    center_holes=holes&taper_roi
    check('no_layer14_center_holes',center_holes.is_empty(),remaining_hole_area_um2=float(center_holes.area()*dbu**2))
    repeated=m.make_device(800.0)
    delta=compare(r,snapshot(repeated),dbu)
    check('factory_repeatable_without_shared_cell_mutation',all(v==0 for v in delta.values()),xor_um2=delta)

    before=OUT/'before.gds'
    if before.exists():
        oldly,_=read(before)
        old_instances=actual_device_instances(oldly)
        assert len(old_instances)==1
        old=file_snapshot(oldly,old_instances[0].cell)
        old_tip_boxes=box_region((-235.25,187.5,-234.75,197.5),dbu)+box_region((-235.25,202.5,-234.75,212.5),dbu)
        old_tip_boxes+=old_tip_boxes.transformed(kdb.Trans(kdb.Trans.R180))
        old_optical=old[(15,0)]-old_tip_boxes
        allowed=box_region((194.999,-205,500,-195),dbu)+box_region((-500,195,-194.999,205),dbu)
        unchanged=((old_optical^optical)-allowed).area()*dbu**2
        check('meander_and_pnc_ridge_unchanged',unchanged<=TOL,unexpected_xor_um2=float(unchanged))
        pd=(old[(11,0)]^r[(11,0)]).area()*dbu**2
        check('entire_layer11_pnc_geometry_unchanged',pd<=TOL,xor_um2=float(pd))

    top=m.build_layout()
    path=ROOT/'gds/2611_tapered_yue_design2.gds'
    export(top,path)
    ly,_=read(path); instances=actual_device_instances(ly)
    check('one_complete_exported_device',len(instances)==1)
    inst=instances[0]
    diff=compare(r,file_snapshot(ly,inst.cell),dbu)
    check('export_preserves_all_device_layers',all(v<=TOL for v in diff.values()),xor_um2=diff)
    edge=inst.trans.disp.x*dbu+r[(11,0)].bbox().left*dbu
    check('left_bend_alignment_retained',abs(edge+900)<1e-8,x_um=edge)
    check('bottom_device_y_retained',abs(inst.trans.disp.y*dbu-m.BOTTOM_DEVICE_Y)<1e-8)
    protected=OUT/'protected_hashes.json'
    if protected.exists():
        for name,sha in json.loads(protected.read_text()).items():
            check('design1_file_unchanged_'+name,hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha)
    render(cr,dbu,'two_layer_coupler.png',(42,-1.8,180,1.85),(10,15),
           'Design 2: Layer 10 = 20/100/10 um; Layer 15 = 20/100 um (transverse scale expanded)',origin_x=45,converter=True)
    render(r,dbu,'device_overview.png',(-1150,-370,1150,370),(2,10,11,13,14,15,16),
           'Design 2: one complete 800 um meander device')
    export(d,OUT/'device.gds')
    REPORT.update(status='PASS' if not REPORT['failures'] else 'FAIL',
                  output_file=str(path.relative_to(ROOT)),
                  layer10_widths_um=[3.,1.65,.980,.200],layer10_lengths_um=[20.,100.,10.],
                  layer15_widths_um=[1.,.500,.200],layer15_lengths_um=[20.,100.],
                  centerline_offset_um=.840,terminal_vertical_edge_gap_um=float(gap),
                  ridge_sections=sampled)
    (OUT/'validation_report.json').write_text(json.dumps(REPORT,indent=2))
    if REPORT['failures']:
        raise AssertionError('Failed checks: '+', '.join(REPORT['failures']))
    print('DESIGN2_REVISED_STAGES_PASS',flush=True)


def main():
    try:
        run()
    except Exception:
        REPORT['status']='FAIL'
        REPORT['exception']=traceback.format_exc()
        (OUT/'validation_report.json').write_text(json.dumps(REPORT,indent=2))
        raise


if __name__=='__main__':
    main()
