import gdsfactory as gf
import numpy as np
import math
gf.gpdk.PDK.activate()
gf.gpdk.PDK.get_layer_views().to_lyp("generic_layers.lyp")


@gf.cell(check_instances=False)
def alignment_markers(
    square_size1: float = 50.0,
    square_size2: float = 50.0,
    cross_arm_length1: float = 300.0,
    cross_arm_width1: float = 10.0,
    cross_arm_length2: float = 300.0,
    cross_arm_width2: float = 10.0,
    square1_pos: tuple = (-200, 300),   # left square position
    square2_pos: tuple = (0, 300), 
    cross1_pos: tuple = (0, 0),   # right square position
    cross2_pos: tuple = (0, 0),          # crosshair position
    layer: tuple = (2, 0),
) -> gf.Component:
    c = gf.Component()

    def add_square(cx, cy,square_size):
        half = square_size / 2
        c.add_polygon([
            (cx - half, cy - half),
            (cx + half, cy - half),
            (cx + half, cy + half),
            (cx - half, cy + half),
        ], layer=layer)

    
    def add_cross(cx, cy,cross_arm_length,cross_arm_width):
        half_l = cross_arm_length / 2
        half_w = cross_arm_width / 2
        c.add_polygon([
            (cx - half_l, cy - half_w),
            (cx + half_l, cy - half_w),
            (cx + half_l, cy + half_w),
            (cx - half_l, cy + half_w),
        ], layer=layer)
        c.add_polygon([
            (cx - half_w, cy - half_l),
            (cx + half_w, cy - half_l),
            (cx + half_w, cy + half_l),
            (cx - half_w, cy + half_l),
        ], layer=layer)

    add_square(*square1_pos,square_size1)
    add_square(*square2_pos,square_size2)
    add_cross(*cross1_pos,cross_arm_length1,cross_arm_width1)
    add_cross(*cross2_pos,cross_arm_length2,cross_arm_width2)

    return c


def snap_angle(angle: float) -> float:
    return float(round(angle / 90) * 90)

def snap_point(pt, grid=0.0001):
    return (
        round(pt[0] / grid) * grid,
        round(pt[1] / grid) * grid,
    )

def snap_w(width: float, grid: float = 0.001) -> float:
    """Snap width to even grid multiple of 0.002 um."""
    dbu = int(round(width / grid))
    if dbu % 2 != 0:
        dbu += 1
    if dbu <= 0:
        dbu = 2
    return float(dbu * grid)


@gf.cell(check_instances=False)   # ← fix 1: added decorator
def half_euler_90(
    radius: float = 80.0,
    wg_width: float = 0.5,
    layer: tuple = (2, 0),
    angle: float = 90,
    euler_ratio: float = 1.0,
    npoints: int = 720,
) -> gf.Component:
    c = gf.Component()
    snapped_width = snap_w(wg_width)
    xs = gf.cross_section.strip(width=snapped_width, layer=layer)
    path_90 = gf.path.euler(radius=radius, angle=angle, p=euler_ratio, npoints=npoints)
    c << gf.path.extrude(path_90, xs)
    c.add_port(
        name="o1",
        center=snap_point(path_90.points[0]),
        width=snapped_width,
        orientation=180.0,
        layer=layer,
    )
    c.add_port(
        name="o2",
        center=snap_point(path_90.points[-1]),
        width=snapped_width,
        orientation=snap_angle(path_90.end_angle),
        layer=layer,
    )
    return c

def make_polygon_component(vertices, layer, port_center, port_width, port_orientation):
    c = gf.Component()
    c.add_polygon(vertices, layer=layer)
    c.add_port(
        name="o1",
        center=port_center,
        width=port_width,
        orientation=port_orientation,
        layer=layer,
    )
    return c



@gf.cell(check_instances=False)
def c_shape_o(  #actually u shape now
    radius_top: float = 80.0,
    length_p: float = 100,
    wgw: float = 10,
    gap: float = 0.25,
    layer: tuple = (2, 0),
    euler_ratio: float = 1.0,

) -> gf.Component:
    c = gf.Component()



    xs = gf.cross_section.strip(width=wgw, layer=layer)
    

    st = c << gf.components.straight(length=length_p/2, cross_section=xs)
    
    bt = c << gf.components.bend_circular(radius=radius_top,angle=180,cross_section=xs,)
    bt.connect("o1", st.ports["o2"])

    stt = c << gf.components.straight(length=length_p/2, cross_section=xs)
    stt.connect("o1", bt.ports["o2"])



    c.add_port("o1", port=st.ports["o1"])
    c.add_port("o2", port=stt.ports["o2"])

    return c


@gf.cell(check_instances=False)
def double_u_shape(  #actually u shape now
    radius_top: float = 80.0,
    radius_bottom: float = 80.0,
    length_p: float = 1000.0,
    wgw: float = 10,
    gap: float = 0.25,
    layer: tuple = (2, 0),
    euler_ratio: float = 1.0,

) -> gf.Component:
    c = gf.Component()

    c1 = c << c_shape_o(
            radius_t = radius_top,
            length_p = length_p,
            wgw = wgw,
            gap = gap,
            layer= layer,
            euler_ratio = euler_ratio,
        )
    c2 = c << c_shape_o(
            radius_t = radius_bottom,
            length_p = length_p,
            wgw = wgw,
            gap = gap,
            layer= layer,
            euler_ratio = euler_ratio,
        )
    
    dx = c1.xmin - c2.xmin
    c2.movex(dx)

    dy = c1.y - c2.y
    c2.movey(dy)


    c.add_port("o1", port=c2.ports["o1"])
    c.add_port("o2", port=c1.ports["o1"])
    c.add_port("o3", port=c1.ports["o2"])
    c.add_port("o4", port=c2.ports["o2"])
    return c














@gf.cell(check_instances=False)
def c_shape(
    radius_top: float = 80.0,
    radius_bottom: float = 80.0,
    length_x: float = 1000.0,
    wgw1: float = 10,
    wgw2: float = 5,
    layer: tuple = (2, 0),
    angle: float = 90,
    euler_ratio: float = 1.0,
) -> gf.Component:
    c = gf.Component()

    # ── snap all widths ──
    w1 = snap_w(wgw1)
    w2 = snap_w(wgw2)

    bt = c << half_euler_90(
        radius=radius_top, wg_width=w1,
        layer=layer, angle=angle, euler_ratio=1.0, npoints=720,
    )
    st = c << gf.components.taper(
        width1=w1, width2=w2,
        length=length_x, layer=layer,
    )
    st.connect("o1", bt.ports["o1"])

    bb = c << half_euler_90(
        radius=radius_bottom, wg_width=w2,
        layer=layer, angle=angle, euler_ratio=1.0, npoints=720,
    )
    bb.connect("o1", st.ports["o2"])

    c.add_port("o1", port=bt.ports["o2"])
    c.add_port("o2", port=bb.ports["o2"])
    return c


@gf.cell(check_instances=False)
def c_shape_e(
    radius_t: float = 100,
    radius_b: float = 100,
    length_x: float = 1000.0,
    wgw: float = 10,
    gap: float = 0.25,
    layer: tuple = (2, 0),
    euler_ratio: float = 1.0,
) -> gf.Component:
    c = gf.Component()
    xs = gf.cross_section.strip(width=snap_w(wgw), layer=layer)

    bt = c << gf.components.bend_circular(radius=radius_t, angle=90, cross_section=xs)
    st = c << gf.components.straight(length=length_x, cross_section=xs)
    st.connect("o1", bt.ports["o1"])
    bb = c << gf.components.bend_circular(radius=radius_b, angle=90, cross_section=xs)
    bb.connect("o1", st.ports["o2"])

    c.add_port("o1", port=bt.ports["o2"])
    c.add_port("o2", port=bb.ports["o2"])
    return c


@gf.cell(check_instances=False)
def meander_SC(  
    n_loops: int = 4,
    o_loops: int = 10,
    length: float = 3000.0,
    radius0: float = 100,
    wgw0: float = 10,
    wgwe: float = 0.2,
    wg_width: float = 0.5,
    delta: float = 0,
    length_after: float = 500,
    length_p: float = 10000,
    length_out: float = 150,
    gap: float = 2.5,
    layer: tuple = (2, 0),
    angle: float = 90,
    euler_ratio: float = 1.0,
    N: float = 2,
) -> gf.Component:
    c = gf.Component()
    segments = []
    length_x = length / n_loops

    # ── snap all widths from start ──
    wgw1 = snap_w(wgw0)
    wgw2 = snap_w(wgw0 - (wgw0-wgwe) * length_x / length)
    radius_top    = snap_w((wgw1 - delta) * N)
    radius_bottom = snap_w((wgw2 - delta) * N)

    for i in range(n_loops):
        if i > 0:
            if i % 2 == 1:
                radius_top    = snap_w(N * (wgw0 - delta - (wgw0-wgwe) * length_x / length))
                wgw1 = snap_w(wgw0 - (wgw0-wgwe) * length_x * (i + 1) / length)
                wgw2 = snap_w(wgw0 - (wgw0-wgwe) * length_x * i / length)
            else:
                radius_bottom = snap_w(N * (wgw0 - delta - (wgw0-wgwe) * length_x  / length))
                wgw1 = snap_w(wgw0 - (wgw0-wgwe) * length_x * i / length)
                wgw2 = snap_w(wgw0 - (wgw0-wgwe) * length_x * (i + 1) / length)

        seg = c << c_shape(
            radius_top=radius_top,
            radius_bottom=radius_bottom,
            length_x=length_x,
            wgw1=wgw1,
            wgw2=wgw2,
            layer=layer,
            angle=angle,
            euler_ratio=euler_ratio,
        )

        if i > 0:
            prev = segments[i - 1]
            if i % 2 == 1:
                seg.mirror_y(seg.xmin)
                seg.connect("o2", prev.ports["o2"])
            else:
                seg.connect("o1", prev.ports["o1"])

        segments.append(seg)

    # ── Input side ──
    r_in = snap_w((wgw0 - delta) * N)
    w_in = snap_w(wgw0)



    bti = c << half_euler_90(
        radius=r_in, wg_width=w_in,
        layer=layer, angle=angle, euler_ratio=1.0, npoints=720,
    )
    bti.connect("o1", segments[0].ports["o1"])

    sti = c << gf.components.taper(
        width1=w_in, width2=w_in,
        length=length_x / 2 - r_in,
        layer=layer,
    )
    sti.connect("o1", bti.ports["o2"])

    bbi = c << half_euler_90(
        radius=r_in, wg_width=w_in,
        layer=layer, angle=angle, euler_ratio=1.0, npoints=720,
    )
    bbi.mirror_x()
    bbi.connect("o1", sti.ports["o2"])

    stii = c << gf.components.taper(
        width1=w_in, width2=w_in,
        length=100, layer=layer,
    )
    stii.connect("o1", bbi.ports["o2"])

    # ── Output side ──
    w_end = snap_w(wgwe)                    # ← fix 4: use wgwe not wgw1

    bto = c << half_euler_90(
        radius=radius_top, wg_width=snap_w(wgw1),
        layer=layer, angle=angle, euler_ratio=1.0, npoints=720,
    )
    bto.mirror_x()
    bto.connect("o1", seg.ports["o1"])

    sto = c << gf.components.taper(
        width1=snap_w(wgw1), width2=w_end,
        length=length_x / 2 - radius_top-50,
        layer=layer,
    )
    sto.connect("o1", bto.ports["o2"])

    bbo = c << half_euler_90(
        radius=radius_top, wg_width=w_end,
        layer=layer, angle=angle, euler_ratio=1.0, npoints=720,
    )
    bbo.connect("o1", sto.ports["o2"])


    xs = gf.cross_section.strip(width=w_end, layer=layer)
    wg_after = c << gf.components.straight(length=length_after, cross_section=xs)
    wg_after.connect("o1", bbo.ports["o2"])

    bend_MZI_seg_top1 = c.add_ref_off_grid(
        gf.components.bend_circular_all_angle(
            radius=radius0,
            angle=angle,
            cross_section=xs,
        )
    )
    bend_MZI_seg_top1.connect("o1", wg_after.ports["o2"])
    
    
    wg_MZI_seg_top1= c << gf.components.straight(length=700, cross_section=xs)
    wg_MZI_seg_top1.connect("o1",bend_MZI_seg_top1.ports["o2"])


    bend_MZI_seg_top2 = c.add_ref_off_grid(
        gf.components.bend_circular_all_angle(
            radius=radius0,
            angle=-angle,
            cross_section=xs,
        )
    )
    bend_MZI_seg_top2.connect("o1", wg_MZI_seg_top1.ports["o2"])
    


    wg_MZI_seg_top2= c << gf.components.straight(length=length_p/2+200, cross_section=xs)
    wg_MZI_seg_top2.connect("o1",bend_MZI_seg_top2.ports["o2"])

    

    c.add_port("o1", port=stii.ports["o2"])   # input pad end
    c.add_port("o2", port=wg_after.ports["o2"])
    c.add_port("o3", port=wg_MZI_seg_top2.ports["o2"])

    return c





@gf.cell(check_instances=False)
def meander_wg(
    n_loops: int = 4,
    o_loops: int = 10,
    length: float = 3000.0,
    radius0: float = 100,
    wgw0: float = 10,
    wgwe: float = 0.2,
    wg_width: float = 0.5,
    delta: float = 0,
    length_after: float = 500,
    length_p: float = 10000,
    length_out: float = 150,
    gap: float = 2.5,
    layer: tuple = (2, 0),
    angle: float = 90,
    euler_ratio: float = 1.0,
    N: float = 2,
) -> gf.Component:
    c = gf.Component()

    input_SC = c << meander_SC(
        n_loops = n_loops,
        o_loops = o_loops,
        length = length,
        radius0 = radius0,
        wgw0 = wgw0,
        wgwe = wgwe,
        wg_width = wg_width,
        delta = delta,
        length_after = length_after,
        length_p = length_p,
        length_out = length_out,
        gap = gap,
        layer = layer,
        angle = angle,
        euler_ratio = 1.0,
        N = N,
    )


    # c.add_port("o1", port=input_SC.ports["o1"])   # input pad end
    # c.add_port("o2", port=input_SC.ports["o2"])


    # ── Optical loops ──
    r_e = radius0 - (gap + wg_width / 2 + (snap_w(wgwe) - delta) / 2)
    segments1 = []

    for i in range(o_loops):
        seg1 = c << c_shape_o(
            radius_top=radius0-gap-(wg_width-delta)/2-wgwe/2,
            length_p=length_p,
            wgw=snap_w(wgwe),
            gap=gap,
            layer=layer,
            euler_ratio=euler_ratio,
        )

        if i > 0:
            prev1 = segments1[i - 1]
            if i % 2 == 1:
                seg1.mirror_x(seg1.xmin)
                seg1.connect("o1", prev1.ports["o2"])
            else:
                seg1.connect("o1", prev1.ports["o2"])

        segments1.append(seg1)

    # ── fix 2: connect after loop, with .ports[] ──
        segments1[0].connect("o1", input_SC.ports["o3"])


    output_SC = c << meander_SC(
        n_loops = n_loops,
        o_loops = o_loops,
        length = length,
        radius0 = radius0,
        wgw0 = wgw0,
        wgwe = wgwe,
        wg_width = wg_width,
        delta = delta,
        length_after = length_after,
        length_p = length_p,
        length_out = length_out,
        gap = gap,
        layer = layer,
        angle = -angle,
        euler_ratio = 1.0,
        N = N,
    )
    #output_SC.mirror_x
    output_SC.connect("o3", segments1[-1].ports["o2"])


    c.add_port("o1", port=input_SC.ports["o1"])   # input pad end
    c.add_port("o2", port=input_SC.ports["o2"])
    c.add_port("o3", port=output_SC.ports["o1"])
    

    for port in c.ports:
        c.add_label(text=port.name, position=port.center, layer=(66, 0))

    #c.flatten()
    return c






# ── Build both layers ──
c1 = meander_wg(
    n_loops=20, o_loops=21, length=30000,radius0=100,
    wgw0=43, wg_width=1, wgwe=0.2, delta=0,
    length_after=200, length_out = 4595, length_p=15517, gap=2.3,
    layer=(1, 0), angle=-90, euler_ratio=1,N=1.5
)

c2 = meander_wg(
    n_loops=20, o_loops=21, length=30000,radius0=100,
    wgw0=43 + 5.6*2, wg_width=1, wgwe=0.2 + 5.6*2, delta=5.6*2,
    length_after=200, length_out = 4595, length_p=15517, gap=2.3,
    layer=(10, 0), angle=-90, euler_ratio=1,N=1.5
)






c = gf.Component("Kp_SC")

ref1    = c << c1
ref2    = c << c2

poly1 = make_polygon_component(
    vertices = [(0, 21.5), (0, -21.5), (-200, -175), (-400,-175),(-400,175),(-200,175)],
    layer=(1, 0),
    port_center=(0, 0),   # edge of polygon where it meets c1
    port_width=ref1.ports["o1"].dwidth,
    port_orientation=0,     # facing right, toward c1
)

poly2 = make_polygon_component(
    vertices = [(0, 5.6+21.5), (0, -5.6-21.5), (-200, -385), (-610,-385),(-610,385),(-200,385)],
    layer=(10, 0),
    port_center=(0, 0),   # edge of polygon where it meets c1
    port_width=ref2.ports["o1"].dwidth,
    port_orientation=0,     # facing right, toward c1
)

poly3 = make_polygon_component(
    vertices = [(0, 21.5), (0, -21.5), (-200, -175), (-400,-175),(-400,175),(-200,175)],
    layer=(1, 0),
    port_center=(0, 0),   # edge of polygon where it meets c1
    port_width=ref1.ports["o3"].dwidth,
    port_orientation=0,     # facing right, toward c1
)

poly4 = make_polygon_component(
    vertices = [(0, 5.6+21.5), (0, -5.6-21.5), (-200, -385), (-610,-385),(-610,385),(-200,385)],
    layer=(10, 0),
    port_center=(0, 0),   # edge of polygon where it meets c1
    port_width=ref2.ports["o3"].dwidth,
    port_orientation=0,     # facing right, toward c1
)



ref_poly1 = c << poly1
ref_poly2 = c << poly2
ref_poly3 = c << poly3
ref_poly4 = c << poly4


ref1.connect("o1", ref_poly1.ports["o1"])
ref2.connect("o1", ref_poly2.ports["o1"])
ref_poly3.connect("o1", ref1.ports["o3"])
ref_poly4.connect("o1", ref2.ports["o3"])

# expose ports
c.add_port("o1", port=ref1.ports["o1"])
c.add_port("o2", port=ref1.ports["o2"])
c.add_port("o3", port=ref1.ports["o3"])


top = gf.Component("top")
ref_all = top << c
ref_all.move(np.array([0, 0]) - np.array(ref_all.ports["o2"].center))


markers1 = alignment_markers(
    square_size1=10,
    square_size2=20,
    cross_arm_length1=200,
    cross_arm_width1=3,
    cross_arm_length2=200,
    cross_arm_width2=2,
    square1_pos=(4000, -1700),
    square2_pos=(4200, -1700),
    cross1_pos=(4950, 2000),
    cross2_pos=(4200, -1200),
    layer=(3, 0),
)

markers2 = alignment_markers(
    square_size1=10,
    square_size2=20,
    cross_arm_length1=200,
    cross_arm_width1=3,
    cross_arm_length2=200,
    cross_arm_width2=2,
    square1_pos=(-800, -1700),
    square2_pos=(-1000, -1700),
    cross1_pos=(-2100, 2000),
    cross2_pos=(-1000, -1200),
    layer=(3, 0),
)


markers3 = alignment_markers(
    square_size1=10,
    square_size2=20,
    cross_arm_length1=200,
    cross_arm_width1=3,
    cross_arm_length2=200,
    cross_arm_width2=2,
    square1_pos=(4000, -14500),
    square2_pos=(4200, -14500),
    cross1_pos=(4950, -17000),
    cross2_pos=(4200, -15000),
    layer=(3, 0),
)

markers4 = alignment_markers(
    square_size1=10,
    square_size2=20,
    cross_arm_length1=200,
    cross_arm_width1=3,
    cross_arm_length2=200,
    cross_arm_width2=2,
    square1_pos=(-800, -14500),
    square2_pos=(-1000, -14500),
    cross1_pos=(-2100, -17000),
    cross2_pos=(-1000, -15000),
    layer=(3, 0),
)

bound = make_polygon_component(
    vertices = [(-3500,1500), (6500, 1500), (6500, -16000), (-3500,-16000)],
    layer=(2, 0),
    port_center=(0, 0),   # edge of polygon where it meets c1
    port_width=10,
    port_orientation=0,     # facing right, toward c1
)

top << markers1
top << markers2
top << markers3
top << markers4
top << bound

print(f"\nPolygons:     {len(top.get_polygons())}")
print(f"Bounding box: {top.bbox()}")
for port in top.ports:
    print(f"  {port.name}: center={port.center}, angle={port.orientation}")
top.flatten()
top.write_gds("Kp_SC.gds")

top.show()
