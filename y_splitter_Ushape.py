import gdsfactory as gf
import numpy as np
gf.gpdk.PDK.activate()

@gf.cell(check_instances=False)
def dice_markers(
    cross_arm_length: float = 300.0,
    cross_arm_width: float = 10.0,
    cross1_pos: tuple = (0, 0),   # right square position       # crosshair position
    layer: tuple = (2, 0),
) -> gf.Component:
    c = gf.Component()

    
    
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

  
    add_cross(*cross1_pos,cross_arm_length,cross_arm_width)
    #add_cross(*cross2_pos,cross_arm_length2,cross_arm_width2)

    return c

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
    #add_cross(*cross2_pos,cross_arm_length2,cross_arm_width2)

    return c

def calc_radius(y_offset: float, angle_deg: float) -> float:
    """Calculate radius needed to achieve a given y_offset at a given angle."""
    a = np.radians(angle_deg)
    return (y_offset / 2) / (1 - np.cos(a))



@gf.cell(check_instances=False)
def half_euler_90(
    radius: float = 80.0,
    wg_width: float = 0.5,
    layer: tuple = (2, 0),
    angle: float = 90,
    euler_ratio: float = 1.0,
    npoints: int = 720,  # <-- Correct argument for modern gdsfactory
) -> gf.Component:
    """True 90 deg Euler bend with explicit high point resolution."""
    c = gf.Component()
    
    snapped_width = snap_width_to_even_grid(wg_width)
    xs = gf.cross_section.strip(width=snapped_width, layer=layer)

    # Use npoints to explicitly force a high density of coordinate vertices
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


@gf.cell(check_instances=False)
def c_shape(  #actually u shape now
    radius_t: float = 100,
    length_p: float = 100,
    wgw: float = 10,
    gap: float = 0.25,
    layer: tuple = (2, 0),
    euler_ratio: float = 1.0,

) -> gf.Component:
    c = gf.Component()



    xs = gf.cross_section.strip(width=wgw, layer=layer)
    

    st = c << gf.components.straight(length=length_p/2, cross_section=xs)
    
    bt = c << gf.components.bend_circular(radius=radius_t,angle=180,cross_section=xs,)
    bt.connect("o1", st.ports["o2"])

    stt = c << gf.components.straight(length=length_p/2, cross_section=xs)
    stt.connect("o1", bt.ports["o2"])



    c.add_port("o1", port=st.ports["o1"])
    c.add_port("o2", port=stt.ports["o2"])

    return c


@gf.cell(check_instances=False)
def double_u_shape(  #actually u shape now
    radius_t: float = 100,
    radius_b: float = 100,
    length_p: float = 1000.0,
    wgw: float = 10,
    gap: float = 0.25,
    layer: tuple = (2, 0),
    euler_ratio: float = 1.0,

) -> gf.Component:
    c = gf.Component()

    c1 = c << c_shape(
            radius_t = radius_t,
            length_p = length_p,
            wgw = wgw,
            gap = gap,
            layer= layer,
            euler_ratio = euler_ratio,
        )
    c2 = c << c_shape(
            radius_t = radius_b,
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


# @gf.cell(check_instances=False)
# def add_holes(
#     n_loops: float = 4,
#     wg_width: float = 1.5,
#     delta: float = 0,
#     radius0: float = 100,
#     layer: tuple = (1, 0),
#     length_in: float = 100,
#     length_after: float = 1000,  # length of waveguides after splitter
#     length_p = 1500,
#     angle: float = 30,
#     angle_o: float = 10,
#     y_offset: float = 300,
#     gap: float =2.5,
#     layer_port: tuple = (66, 0),
#     euler_ratio: float = 1.0
# ) -> gf.Component:
#     c = gf.Component()
#     xs = gf.cross_section.strip(width=wg_width, layer=layer)
#     radius1= calc_radius(y_offset, angle)
#     radius2 = calc_radius(y_offset-gap-(wg_width-delta)/2-0.05, angle)


#     # ── Input waveguide ──
#     wg_in = c << gf.components.straight(length=length_in, cross_section=xs)

  

#     # ── Top arm: single bend up ──
#     bend_top = c << gf.components.bend_euler(
#         radius=radius1, angle=angle, cross_section=xs
#     )
#     bend_top.connect("o1", wg_in.ports["o2"])

#     wg_top = c << gf.components.straight(length=y_str, cross_section=xs)

#     wg_top.connect("o1", bend_top.ports["o2"])

#     # ── Top arm: second bend back to horizontal ──
#     bend_top2 = c << gf.components.bend_euler(
#         radius=radius1, angle=angle, cross_section=xs
#     )
#     bend_top2.mirror_x()
#     bend_top2.connect("o1", wg_top.ports["o2"])

#     # ── Top output straight ──
#     wg_after_top = c << gf.components.straight(length=length_after, cross_section=xs)
#     wg_after_top.connect("o1", bend_top2.ports["o2"])
# # ── Top arm: second bend back to horizontal ──
#     bend_top3 = c << gf.components.bend_euler(
#         radius=radius2, angle=angle, cross_section=xs
#     )
#     bend_top3.mirror_x()
#     bend_top3.connect("o1", wg_after_top.ports["o2"])
# # ── Top arm: single bend up ──
#     bend_top4 = c << gf.components.bend_euler(
#         radius=radius2, angle=angle, cross_section=xs
#     )
#     bend_top4.connect("o1", bend_top3.ports["o2"])








#    # ── Bottom arm: single bend down ──
#     bend_bot = c << gf.components.bend_euler(
#         radius=radius1, angle=angle, cross_section=xs
#     )
#     bend_bot.mirror_x()
#     bend_bot.connect("o1", wg_in.ports["o2"])

#     wg_bot = c << gf.components.straight(length=y_str, cross_section=xs)

#     wg_bot.connect("o1", bend_bot.ports["o2"])

#     # ── Bottom arm: second bend back to horizontal ──
#     bend_bot2 = c << gf.components.bend_euler(
#         radius=radius1, angle=angle, cross_section=xs
#     )
#     bend_bot2.connect("o1", wg_bot.ports["o2"])

#     # ── Bottom output straight ──
#     wg_after_bot = c << gf.components.straight(length=length_after, cross_section=xs)
#     wg_after_bot.connect("o1", bend_bot2.ports["o2"])

#     # ── Top arm: second bend back to horizontal ──
#     bend_bot3 = c << gf.components.bend_euler(
#         radius=radius2, angle=angle, cross_section=xs
#     )

#     bend_bot3.connect("o1", wg_after_bot.ports["o2"])
# # ── Top arm: single bend up ──
#     bend_bot4 = c << gf.components.bend_euler(
#         radius=radius2, angle=angle, cross_section=xs
#     )

#     bend_bot4.mirror_x()
#     bend_bot4.connect("o1", bend_bot3.ports["o2"])


# # ── Expose ports ──
#     c.add_port("o1", port=wg_in.ports["o1"])            # input
#     c.add_port("o2", port=bend_top4.ports["o2"])     # top output
#     c.add_port("o3", port=bend_bot4.ports["o2"])     # bottom output

    
 

   
#     segments1 = []
#     radius = 100

#     for i in range(n_loops):
        
#         seg1 = c << c_shape(
#             radius_t = radius0,
#             radius_b = radius0-gap*2-(wg_width-delta)-0.2,
#             length_x = length_p,
#             wgw = wg_width,
#             gap = gap,
#             layer= layer,
#             euler_ratio = euler_ratio,
#         )

#         if i > 0:
#             prev1 = segments1[i - 1]
#             if i % 2 == 1:
#                 seg1.mirror_y()
#                 seg1.connect("o2", prev1.ports["o2"])
#             else:
#                 seg1.connect("o1", prev1.ports["o1"])

#         segments1.append(seg1)
#         segments1[0].connect("o1", bend_top4.ports["o2"])

#     segments2 = []

#     for i in range(n_loops):
        

#         seg2 = c << c_shape(
#             radius_t = radius0-gap*2-wg_width+delta-0.2,
#             radius_b = radius0,
#             length_x = length_p,
#             wgw = wg_width,
#             gap = gap,
#             layer= layer,
#             euler_ratio = euler_ratio,
#         )


#         if i > 0:
#             prev2 = segments2[i - 1]
#             if i % 2 == 1:
#                 seg2.mirror_y()
#                 seg2.connect("o2", prev2.ports["o2"])
#             else:
#                 seg2.connect("o1", prev2.ports["o1"])

#         segments2.append(seg2)
#         segments2[0].connect("o1", bend_bot4.ports["o2"])






@gf.cell(check_instances=False)
def MZI(
    n_loops: float = 4,
    wg_width: float = 0.5,
    wgwe: float = 0.2,
    wgwe_d: float = 0.5,
    delta: float = 0,
    radius0: float = 100,
    layer: tuple = (1, 0),
    length_in: float = 100,
    length1: float = 100,
    length2: float = 100,
    length3: float = 100,
    length4: float = 100,
    length_mzi_in: float = 1000,  # length of waveguides after splitter
    length_mzi_out: float = 1000,
    length_p: float = 1500,
    angle: float = 30,
    angle_o: float = 10,
    y_offset: float = 300,
    gap: float =2.5,
    layer_port: tuple = (66, 0),
    euler_ratio: float = 1.0,
    y_str1: float = 100,
    y_str2: float = 200
) -> gf.Component:
    c = gf.Component()
    radius1= calc_radius(y_offset, angle)
    radius2 = calc_radius(y_offset-gap-(wg_width-delta)/2-wgwe/2-wgwe_d, angle)
    xs = gf.cross_section.strip(width=wg_width, layer=layer)
    
    wgi = c << gf.components.straight(length=0, cross_section=xs)


    # ── Top arm: single bend up ──
    bend_top = c << gf.components.bend_euler(
        radius=radius1, angle=angle, cross_section=xs
    )
    bend_top.connect("o1", wgi.ports["o2"])

    wg_top = c << gf.components.straight(length=y_str1, cross_section=xs)

    wg_top.connect("o1", bend_top.ports["o2"])

    # ── Top arm: second bend back to horizontal ──
    bend_top2 = c << gf.components.bend_euler(
        radius=radius1, angle=angle, cross_section=xs
    )
    bend_top2.mirror_x()
    bend_top2.connect("o1", wg_top.ports["o2"])

    # ── Top output straight ──
    wg_after_top = c << gf.components.straight(length=length_mzi_in, cross_section=xs)
    wg_after_top.connect("o1", bend_top2.ports["o2"])
# ── Top arm: second bend back to horizontal ──
    bend_top3 = c << gf.components.bend_euler(
        radius=radius2, angle=angle, cross_section=xs
    )
    bend_top3.mirror_x()
    bend_top3.connect("o1", wg_after_top.ports["o2"])

    wg_top2 = c << gf.components.straight(length=y_str1, cross_section=xs)

    wg_top2.connect("o1", bend_top3.ports["o2"])

# ── Top arm: single bend up ──
    bend_top4 = c << gf.components.bend_circular(
        radius=radius2, angle=angle, cross_section=xs
    )
    bend_top4.connect("o1", wg_top2.ports["o2"])


    



   # ── Bottom arm: single bend down ──
    bend_bot = c << gf.components.bend_euler(
        radius=radius1, angle=angle, cross_section=xs
    )
    bend_bot.mirror_x()
    bend_bot.connect("o1", wgi.ports["o2"])

    wg_bot = c << gf.components.straight(length=y_str1, cross_section=xs)

    wg_bot.connect("o1", bend_bot.ports["o2"])


    # ── Bottom arm: second bend back to horizontal ──
    bend_bot2 = c << gf.components.bend_euler(
        radius=radius1, angle=angle, cross_section=xs
    )
    bend_bot2.connect("o1", wg_bot.ports["o2"])

    # ── Bottom output straight ──
    wg_after_bot = c << gf.components.straight(length=length_mzi_in, cross_section=xs)
    wg_after_bot.connect("o1", bend_bot2.ports["o2"])

    # ── Top arm: second bend back to horizontal ──
    bend_bot3 = c << gf.components.bend_euler(
        radius=radius2, angle=angle, cross_section=xs
    )

    bend_bot3.connect("o1", wg_after_bot.ports["o2"])

    wg_bot2 = c << gf.components.straight(length=y_str1, cross_section=xs)

    wg_bot2.connect("o1", bend_bot3.ports["o2"])

# ── Top arm: single bend up ──
    bend_bot4 = c << gf.components.bend_circular(
        radius=radius2, angle=angle, cross_section=xs
    )

    bend_bot4.mirror_x()
    bend_bot4.connect("o1", wg_bot2.ports["o2"])

   

# ── Expose ports ──
    c.add_port("o1", port=bend_top.ports["o1"])            # input
    c.add_port("o2", port=bend_bot.ports["o1"])
    c.add_port("o3", port=bend_top4.ports["o2"])
    c.add_port("o4", port=bend_bot4.ports["o2"])

    return c



@gf.cell(check_instances=False)
def y_splitter_connected(
    n_loops: float = 4,
    wg_width: float = 0.5,
    wgwe: float = 0.2,
    wgwe_d: float = 0.2,
    delta: float = 0,
    radius0: float = 100,
    layer: tuple = (1, 0),
    length_in: float = 100,
    length1: float = 100,
    length2: float = 100,
    length3: float = 100,
    length4: float = 100,
    length_mzi_in: float = 1000,  # length of waveguides after splitter
    length_mzi_out: float = 1000,
    length_p: float = 1500,
    angle: float = 30,
    angle_o: float = 10,
    y_offset: float = 300,
    gap: float =2.5,
    layer_port: tuple = (66, 0),
    euler_ratio: float = 1.0,
    y_str1: float = 100,
    y_str2: float = 200
) -> gf.Component:
    c = gf.Component()
    xs = gf.cross_section.strip(width=wg_width, layer=layer)
    
    radius1= calc_radius(y_offset, angle)
    radius2 = calc_radius(y_offset-gap-(wg_width-delta)/2-wgwe/2-wgwe_d,angle)


    # ── Input waveguide ──
    wg_in1 = c << gf.components.straight(length=length1, cross_section=xs)
    wg_in1.rotate(90)

    bend_in1 = c << gf.components.bend_euler(
        radius=radius0, angle=angle, cross_section=xs
    )
    bend_in1.connect("o1", wg_in1.ports["o2"])

    wg_in2 = c << gf.components.straight(length=length2, cross_section=xs)
    wg_in2.connect("o1", bend_in1.ports["o2"])

    bend_in2 = c << gf.components.bend_euler(
        radius=radius0, angle=-angle, cross_section=xs
    )
    bend_in2.connect("o1", wg_in2.ports["o2"])

    wg_in3 = c << gf.components.straight(length=length3, cross_section=xs)
    wg_in3.connect("o1", bend_in2.ports["o2"])

    bend_in3 = c << gf.components.bend_euler(
        radius=radius0, angle=-angle, cross_section=xs
    )
    
    bend_in3.connect("o1", wg_in3.ports["o2"])

    input_MZI = c << MZI(
        n_loops = n_loops,
        wg_width=wg_width,
        wgwe=wgwe,
        wgwe_d= wgwe_d,
        delta = delta,
        radius0=radius0,
        layer=layer,
        length_in=length_in,
        length1=length1,
        length2=length2,
        length3=length3,
        length4 = length4,
        length_mzi_in=length_mzi_in,
        length_mzi_out=length_mzi_out,
        length_p=length_p,
        angle=angle,
        angle_o=angle_o,
        y_offset = y_offset,
        y_str1=y_str1,
        y_str2=y_str2,
        gap=gap,
        euler_ratio=euler_ratio,
    )
    
    input_MZI.connect("o1",bend_in3.ports["o2"])


    bend_MZI_seg_top1 = c.add_ref_off_grid(
        gf.components.bend_circular_all_angle(
            radius=radius0+gap+wg_width/2-delta/2+wgwe/2+wgwe_d,
            angle=-angle,
            cross_section=xs,
        )
    )
    bend_MZI_seg_top1.connect("o1", input_MZI.ports["o3"])
    
    bend_MZI_seg_bot1 = c.add_ref_off_grid(
        gf.components.bend_circular_all_angle(
            radius=radius0-gap-wg_width/2+delta/2-wgwe/2-wgwe_d,
            angle=-angle,
            cross_section=xs,
        )
    )
    bend_MZI_seg_bot1.connect("o1", input_MZI.ports["o4"])


    wg_MZI_seg_top1= c << gf.components.straight(length=700, cross_section=xs)
    wg_MZI_seg_top1.connect("o1",bend_MZI_seg_top1.ports["o2"])

    wg_MZI_seg_bot1= c << gf.components.straight(length=700, cross_section=xs)
    wg_MZI_seg_bot1.connect("o1",bend_MZI_seg_bot1.ports["o2"])


    bend_MZI_seg_top2 = c << gf.components.bend_circular(
        radius=radius0-gap-wg_width/2+delta/2-wgwe/2, angle=angle, cross_section=xs
    )
    bend_MZI_seg_top2.connect("o1", wg_MZI_seg_top1.ports["o2"])
    
    bend_MZI_seg_bot2 = c << gf.components.bend_circular(
        radius=radius0+gap+wg_width/2-delta/2+wgwe/2, angle=angle, cross_section=xs
    )
    bend_MZI_seg_bot2.connect("o1", wg_MZI_seg_bot1.ports["o2"])


    wg_MZI_seg_top2= c << gf.components.straight(length=length_p/2+200, cross_section=xs)
    wg_MZI_seg_top2.connect("o1",bend_MZI_seg_top2.ports["o2"])

    wg_MZI_seg_bot2= c << gf.components.straight(length=length_p/2+200, cross_section=xs)
    wg_MZI_seg_bot2.connect("o1",bend_MZI_seg_bot2.ports["o2"])


    segments1 = []
   

    for i in range(n_loops):
        
        seg1 = c << double_u_shape(
            radius_t = radius0-gap*2-(wg_width-delta)-wgwe,
            radius_b = radius0,
            length_p = length_p,
            wgw = wg_width,
            gap = gap,
            layer= layer,
            euler_ratio = euler_ratio,
        )

        if i > 0:
            prev1 = segments1[i - 1]
            if i % 2 == 1:
                 
                seg1.mirror_y()
                seg1.connect("o1", prev1.ports["o3"])
            else:
                seg1.connect("o1", prev1.ports["o3"])

        segments1.append(seg1)
        segments1[0].connect("o1", wg_MZI_seg_bot2.ports["o2"])
        segments1[0].connect("o2", wg_MZI_seg_top2.ports["o2"])

   

    
    MZI_wg_seg_top1= c << gf.components.straight(length=length_p/2+200, cross_section=xs)
    MZI_wg_seg_top1.connect("o1",segments1[-1].ports["o4"])

    MZI_wg_seg_bot1= c << gf.components.straight(length=length_p/2+200, cross_section=xs)
    MZI_wg_seg_bot1.connect("o1",segments1[-1].ports["o3"])

    
    MZI_bend_seg_top1 = c << gf.components.bend_circular(
        radius=radius0+gap+wg_width/2-delta/2+wgwe/2, angle=angle, cross_section=xs
    )
    MZI_bend_seg_top1.connect("o1", MZI_wg_seg_top1.ports["o2"])
    
    MZI_bend_seg_bot1 = c << gf.components.bend_circular(
        radius=radius0-gap-wg_width/2+delta/2-wgwe/2, angle=angle, cross_section=xs
    )
    MZI_bend_seg_bot1.connect("o1", MZI_wg_seg_bot1.ports["o2"])


    MZI_wg_seg_top2= c << gf.components.straight(length=700, cross_section=xs)
    MZI_wg_seg_top2.connect("o1",MZI_bend_seg_top1.ports["o2"])

    MZI_wg_seg_bot2= c << gf.components.straight(length=700, cross_section=xs)
    MZI_wg_seg_bot2.connect("o1",MZI_bend_seg_bot1.ports["o2"])


    MZI_bend_seg_top2 = c.add_ref_off_grid(
        gf.components.bend_circular_all_angle(
            radius=radius0-gap-wg_width/2+delta/2-wgwe/2,
            angle=-angle,
            cross_section=xs,
        )
    )
    MZI_bend_seg_top2.connect("o1", MZI_wg_seg_top2.ports["o2"])
    
    MZI_bend_seg_bot2 = c.add_ref_off_grid(
        gf.components.bend_circular_all_angle(
            radius=radius0+gap+wg_width/2-delta/2+wgwe/2,
            angle=-angle,
            cross_section=xs,
        )
    )
    MZI_bend_seg_bot2.connect("o1", MZI_wg_seg_bot2.ports["o2"])


    

    output_MZI= c << MZI(
        n_loops = n_loops,
        wg_width=wg_width,
        wgwe=wgwe,
        wgwe_d= wgwe_d,
        delta = delta,
        radius0=radius0,
        layer=layer,
        length_in=length_in,
        length1=length1,
        length2=length2,
        length3=length3,
        length4 = length4,
        length_mzi_in=length_mzi_in,
        length_mzi_out=length_mzi_out,
        length_p=length_p,
        angle=angle,
        angle_o=angle_o,
        y_offset = y_offset,
        y_str1=y_str1,
        y_str2=y_str2,
        gap=gap,
        euler_ratio=euler_ratio,
    )
    

    output_MZI.connect("o4",MZI_bend_seg_bot2.ports["o2"])




    bend_out3 = c << gf.components.bend_euler(
        radius=radius0, angle=-angle, cross_section=xs
    )
    
    bend_out3.connect("o1", output_MZI.ports["o1"])
    
    wg_out3 = c << gf.components.straight(length=length3, cross_section=xs)
    wg_out3.connect("o1", bend_out3.ports["o2"])
    
    bend_out2 = c << gf.components.bend_euler(
        radius=radius0, angle=-angle, cross_section=xs
    )
    bend_out2.connect("o1", wg_out3.ports["o2"])

    wg_out2 = c << gf.components.straight(length=length4, cross_section=xs)
    wg_out2.connect("o1", bend_out2.ports["o2"])

    bend_out1 = c << gf.components.bend_euler(
        radius=radius0, angle=angle, cross_section=xs
    )
    bend_out1.connect("o1", wg_out2.ports["o2"])

    wg_out1 = c << gf.components.straight(length=length1, cross_section=xs)
    wg_out1.connect("o1", bend_out1.ports["o2"])

    c.add_port("o1", port=input_MZI.ports["o3"])            # input
    #c.add_port("o2", port=input_MZI.ports["o4"])
   
    



    # wg_out = c << gf.components.straight(length=length_in, cross_section=xs)
    # wg_out.connect("o1", bend_bot9.ports["o2"])



    for port in c.ports:
        c.add_label(text=port.name, position=port.center, layer=layer_port)

    
    return c


c1 = y_splitter_connected(
    n_loops = 21,
    wg_width=1.5,
    wgwe=0.2,
    wgwe_d = 0,
    delta = 0,
    radius0=100,
    layer=(10, 0),
    length_in=1500,
    length1=2000,
    length2=23000,
    length3=1000,
    length4 = 23000,
    length_mzi_in=5700,
    length_mzi_out=5700,
    length_p=15517,
    angle=90,
    angle_o=10,
    y_offset = 200,
    y_str1=800,
    y_str2=800,
    gap=2.05,
    euler_ratio=1.0,
)

c2 = y_splitter_connected(
    n_loops = 21,
    wg_width=1.5+2.05*2,
    wgwe=0.2,
    wgwe_d = 0,
    delta = 2.05*2,
    radius0=100,
    layer=(1, 0),
    length_in=1500,
    length1=2000,
    length2=23000,
    length3=1000,
    length4 = 23000,
    length_mzi_in=5700,
    length_mzi_out=5700,
    length_p=15517,
    angle=90,
    angle_o=10,
    y_offset = 200,
    y_str1=800,
    y_str2=800,
    gap=2.05,
    euler_ratio=1.0,
)


c = gf.Component("Kp_SC")
c << c1
c << c2

# expose ports
c.add_port("o1", port=c1.ports["o1"])
c.add_port("o2", port=c2.ports["o1"])

top = gf.Component("top")
ref_all = top << c
ref_all.move(np.array([0, 2.9]) - np.array(ref_all.ports["o1"].center))

markers1 = alignment_markers(
    square_size1=10,
    square_size2=20,
    cross_arm_length1=200,
    cross_arm_width1=3,
    cross_arm_length2=200,
    cross_arm_width2=2,
    square1_pos=(500, 4000),
    square2_pos=(700, 4000),
    cross1_pos=(500, 4200),
    layer=(3, 0),
)

markers2 = alignment_markers(
    square_size1=10,
    square_size2=20,
    cross_arm_length1=200,
    cross_arm_width1=3,
    cross_arm_length2=200,
    cross_arm_width2=2,
    square1_pos=(15000, 4000),
    square2_pos=(14800, 4000),
    cross1_pos=(15000, 4200),
    layer=(3, 0),
)


markers3 = alignment_markers(
    square_size1=10,
    square_size2=20,
    cross_arm_length1=200,
    cross_arm_width1=3,
    cross_arm_length2=200,
    cross_arm_width2=2,
    square1_pos=(500, -1600),
    square2_pos=(700, -1600),
    cross1_pos=(500, -1800),
    layer=(3, 0),
)

markers4 = alignment_markers(
    square_size1=10,
    square_size2=20,
    cross_arm_length1=200,
    cross_arm_width1=3,
    cross_arm_length2=200,
    cross_arm_width2=2,
    square1_pos=(15000, -1600),
    square2_pos=(14800, -1600),
    cross1_pos=(15000, -1800),
    layer=(3, 0),
)

dice_markers1 = dice_markers(
    cross_arm_length=300,
    cross_arm_width=5,
    cross1_pos=(-6250, 5500),
    layer=(3, 0),
)

dice_markers2 = dice_markers(
    cross_arm_length=300,
    cross_arm_width=5,
    cross1_pos=(-6250,-3000),
    layer=(3, 0),
)

dice_markers3 = dice_markers(
    cross_arm_length=300,
    cross_arm_width=5,
    cross1_pos=(16700,5500),
    layer=(3, 0),
)

dice_markers4 = dice_markers(
    cross_arm_length=300,
    cross_arm_width=5,
    cross1_pos=(16700,-3000),
    layer=(3, 0),
)





top << markers1
top << markers2
top << markers3
top << markers4
top << dice_markers1
top << dice_markers2
top << dice_markers3
top << dice_markers4


for port in top.ports:
    print(f"  {port.name}: center={port.center}, angle={port.orientation}")
top.flatten()
top.write_gds("y_splitter.gds")

top.show()
