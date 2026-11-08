"""Layer-14-only terminal windows for the original Design 1/Design 2.

Called by standalone build_layout() and the combined-sweep worker, NOT by
make_device(). This keeps source construction used by the L6000 family intact.
All dimensions are micrometres. Preserve the complete PNC opening, blend from
its existing end width, offset the optical taper by the requested margin, then
use one straight-sided channel to the two unchanged 300 um circles.
"""
from __future__ import annotations
import math


def _cross_section(shape, x):
    from shapely.geometry import LineString
    cut = shape.intersection(LineString([(x, -1000), (x, 1000)]))
    if cut.geom_type == 'MultiLineString':
        from shapely.ops import linemerge
        cut = linemerge(cut)
    if cut.is_empty or cut.geom_type != 'LineString':
        raise ValueError(f'Expected one solid window cross-section at {x}: {cut.geom_type}')
    return cut.bounds[1], cut.bounds[3]


def make_terminal_window(window, slab, stage_lengths, margin=1.0, grid=0.001):
    """Return a normalized replacement: s=0 at the post-PNC L10 entrance.

    The first optical stage inherits the PNC opening width at s=0. Its smooth
    blend avoids narrowing or shortening the PNC support window. Thereafter
    margin is the normal offset from the exported Layer-10 taper polygon.
    """
    from shapely.geometry import Polygon, box
    from shapely.ops import unary_union
    if margin <= 0 or not math.isfinite(margin) or min(stage_lengths) <= 0:
        raise ValueError('Margin and stage lengths must be positive')
    if window.geom_type != 'Polygon':
        raise ValueError('Expected one connected terminal window')
    first, total = float(stage_lengths[0]), float(sum(stage_lengths))
    radius = 150.0
    upper = sorted({round(x, 6) for x,y in window.exterior.coords if abs(y-radius)<grid/10})
    lower = sorted({round(x, 6) for x,y in window.exterior.coords if abs(y+radius)<grid/10})
    if len(upper) != 2 or upper != lower:
        raise ValueError('Cannot identify the two existing 300 um circles')
    far = window.bounds[2]
    if upper[0]-radius <= total or far <= upper[-1]+radius:
        raise ValueError('Unexpected circle/channel positions')
    low0, high0 = _cross_section(window, 0.0)
    farlow, farhigh = _cross_section(window, far)
    optical = slab.intersection(box(0,-30,total,30))
    if optical.geom_type != 'Polygon' or optical.is_empty:
        raise ValueError('Expected a connected Layer-10 optical taper')
    expanded = optical.buffer(margin, join_style=2)
    low1, high1 = _cross_section(expanded, first)
    low2, high2 = _cross_section(expanded, total)
    # Use the same 201-point quintic sampling as the existing taper profile.
    xs = [first*i/200 for i in range(201)]
    f = lambda u: u*u*u*(10+u*(-15+6*u))
    lo = [(x,low0+(low1-low0)*f(x/first)) for x in xs]
    hi = [(x,high0+(high1-high0)*f(x/first)) for x in xs]
    blend = Polygon(hi + lo[::-1])
    near_taper = expanded.intersection(box(first,-100,total,100))
    channel = Polygon([(total,low2),(far,farlow),(far,farhigh),(total,high2)])
    circles=[]
    snap=lambda v: round(v/grid)*grid
    for cx in upper:
        circles.append(Polygon([(snap(cx+radius*math.cos(math.radians(i*2.5))),
                                 snap(radius*math.sin(math.radians(i*2.5))))
                                for i in range(144)]))
    replacement=unary_union([blend,near_taper,channel,*circles])
    if not replacement.is_valid or replacement.geom_type!='Polygon' or replacement.interiors:
        raise ValueError('Replacement must be one solid, hole-free terminal window')
    return replacement, dict(stage_lengths_um=list(stage_lengths),
        post_pnc_margin_um=margin, pnc_end_width_um=high0-low0,
        first_stage_end_width_um=high1-low1, linear_start_from_l10_um=total,
        linear_start_width_um=high2-low2, far_end_from_l10_um=far,
        far_edges_um=[farlow,farhigh], circle_centers_from_l10_um=upper,
        circle_diameter_um=300.0, intermediate_frame_removed=True)


def apply_post_pnc_layer14(component, pnc_length, stage_lengths, margin=1.0):
    """Edit only owned Layer-14 polygons; no shared cell or other layer changes."""
    import gdsfactory as gf
    import klayout.db as kdb
    from shapely.geometry import Polygon, box
    from shapely.affinity import affine_transform
    from shapely.ops import unary_union
    component.flatten()
    dbu=float(component.kcl.dbu)
    index=gf.get_layer((14,0)); slab_index=gf.get_layer((10,0))
    def to_shape(region):
        shapes=[]
        for p in region.each():
            hull=[(q.x*dbu,q.y*dbu) for q in p.each_point_hull()]
            holes=[[(q.x*dbu,q.y*dbu) for q in p.each_point_hole(i)] for i in range(p.holes())]
            shapes.append(Polygon(hull,holes))
        return unary_union(shapes)
    original=kdb.Region(component.kdb_cell.shapes(index)).merged()
    slab=to_shape(kdb.Region(component.kdb_cell.shapes(slab_index)).merged())
    window=to_shape(original)
    modified=original.dup(); records=[]
    for center,orientation in zip(component.info['connection_centers'],
                                  component.info['connection_orientations'],strict=True):
        angle=math.radians(float(orientation))
        sign=int(round(-math.cos(angle)))
        if abs(math.sin(angle))>1e-9 or abs(sign)!=1:
            raise ValueError('Expected horizontal, opposing taper assemblies')
        sx=float(center[0])+sign*pnc_length; cy=float(center[1])
        tr=[sign,0,0,sign,-sign*sx,-sign*cy]
        back=[sign,0,0,sign,sx,cy]
        roi=box(0,-175,20000,175)
        local_window=affine_transform(window,tr).intersection(roi)
        local_slab=affine_transform(slab,tr).intersection(roi)
        new,record=make_terminal_window(local_window,local_slab,stage_lengths,margin,dbu)
        # Cut on the exact PNC endpoint, leaving every upstream coordinate intact.
        rect=affine_transform(roi,back).bounds
        removal=kdb.Region(kdb.Box(*(round(v/dbu) for v in rect)))
        mapped=affine_transform(new,back)
        added=kdb.Region(kdb.Polygon([kdb.Point(round(x/dbu),round(y/dbu))
                                      for x,y in list(mapped.exterior.coords)[:-1]]))
        modified=((modified-removal)|added).merged()
        record.update(side='input' if sign<0 else 'output',
                      l10_start_local_x_um=sx,center_y_um=cy,
                      pnc_length_um=pnc_length,outward_sign=sign)
        records.append(record)
    component.kdb_cell.shapes(index).clear()
    component.kdb_cell.shapes(index).insert(modified)
    component.info['layer14_post_pnc_windows']=records
    return component
