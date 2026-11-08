"""Layer-14-only linear continuation; all dimensions are in micrometres.

Keep the upstream release window unchanged. Replace all downstream window
fragments with one straight-sided channel and the two existing glue circles.
This changes exposure-mask geometry only, never optical tapers or Layer 2.
"""
from __future__ import annotations
import math


def linearize_layer14_windows(layout, entries, policies, start_distance_um=58.847,
                              circle_diameter_um=300.0):
    import klayout.db as kdb
    dbu = float(layout.dbu)
    if not math.isfinite(start_distance_um) or start_distance_um <= 40.0:
        raise ValueError('Linear continuation must start after the first 40 um stage')
    if abs(circle_diameter_um - 300.0) > 1e-9:
        raise ValueError('This revision preserves the existing 300 um glue circles')
    unit = lambda v: round(v / dbu)
    join = unit(start_distance_um)
    radius = unit(circle_diameter_um / 2)
    li = layout.layer(14, 0)
    records = []

    def outer_edges(region, x):
        ys = []
        for poly in region.each():
            pts = list(poly.each_point_hull())
            if poly.holes():
                raise ValueError('Unexpected hole in the terminal release window')
            for a, b in zip(pts, pts[1:] + pts[:1]):
                if a.x == b.x:
                    if a.x == x:
                        ys.extend((a.y, b.y))
                elif min(a.x, b.x) <= x <= max(a.x, b.x):
                    ys.append(a.y + (b.y-a.y) * (x-a.x) / (b.x-a.x))
        if len(ys) < 2:
            raise ValueError('Cannot measure the existing channel cross-section')
        low, high = round(min(ys)), round(max(ys))
        if high <= low:
            raise ValueError('Release window has zero width')
        return low, high

    for entry in entries:
        cell = layout.cell(entry['case_id'])
        if cell is None or list(cell.each_inst()):
            raise ValueError('Expected independently owned flat device cells')
        original = kdb.Region(cell.shapes(li)).merged()
        modified = original.dup()
        for end in policies[entry['case_id']]['ends']:
            sx, cy = end['l10_start_local_x_um'], end['center_y_um']
            sign = -1 if end['side'] == 'input' else 1
            to_local = kdb.Trans(2 if sign < 0 else 0, False,
                                 -sign*unit(sx), -sign*unit(cy))
            back = to_local.inverted()
            roi = kdb.Region(kdb.Box(0, -radius-unit(25), unit(20000), radius+unit(25)))
            window = (original.transformed(to_local) & roi).merged()
            if window.count() != 1:
                raise ValueError('Expected one connected terminal window')
            points = [p for poly in window.each() for p in poly.each_point_hull()]
            upper = sorted({p.x for p in points if p.y == radius})
            lower = sorted({p.x for p in points if p.y == -radius})
            if len(upper) != 2 or upper != lower:
                raise ValueError('Cannot identify the two unchanged 300 um circle centers')
            if upper[0]-radius <= join:
                raise ValueError('Linear start must precede the first glue circle')
            far_x = window.bbox().right
            low, high = outer_edges(window, join)
            far_low, far_high = outer_edges(window, far_x)
            if far_x <= upper[-1]+radius:
                raise ValueError('Unexpected outer channel endpoint')
            channel = kdb.Region(kdb.Polygon([
                kdb.Point(join, low), kdb.Point(far_x, far_low),
                kdb.Point(far_x, far_high), kdb.Point(join, high)]))
            circles = kdb.Region()
            for cx in upper:
                vertices = [kdb.Point(cx+round(radius*math.cos(math.radians(i*2.5))),
                                     round(radius*math.sin(math.radians(i*2.5))))
                            for i in range(144)]
                circles.insert(kdb.Polygon(vertices))
            prefix = window & kdb.Region(kdb.Box(0, -unit(1000), join, unit(1000)))
            replacement = (prefix | channel | circles).merged()
            if replacement.count() != 1 or any(p.holes() for p in replacement.each()):
                raise ValueError('Linear channel and both circles must form one solid opening')
            modified = ((modified-window.transformed(back)) | replacement.transformed(back)).merged()
            records.append(dict(case_id=entry['case_id'], end=end['side'],
                start_distance_from_l10_um=join*dbu,
                start_local_x_um=sx+sign*join*dbu,
                start_global_x_um=entry['origin_um'][0]+sx+sign*join*dbu,
                start_width_um=(high-low)*dbu,
                start_edges_local_normalized_um=[low*dbu,high*dbu],
                outer_endpoint_distance_um=far_x*dbu,
                outer_edges_local_normalized_um=[far_low*dbu,far_high*dbu],
                circle_diameter_um=circle_diameter_um,
                circle_center_distances_from_l10_um=[x*dbu for x in upper],
                circle_centers_global_um=[[entry['origin_um'][0]+sx+sign*x*dbu,
                                           entry['origin_um'][1]+cy] for x in upper],
                middle_frame_removed=True))
        cell.shapes(li).clear()
        cell.shapes(li).insert(modified)
        entry['layer14_linear_start_distance_um'] = start_distance_um
    if len(records) != sum(len(p['ends']) for p in policies.values()):
        raise ValueError('Wrong number of linearized taper windows')
    return dict(enabled=True,window_count=len(records),circle_count=2*len(records),
                preserved_upstream_margin_um=1.0,windows=records)
