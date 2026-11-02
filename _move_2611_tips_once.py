"""One-time, guarded migration of final alignment-tip layer assignment."""
from pathlib import Path
import hashlib

ROOT = Path(__file__).resolve().parent
p = ROOT / '2611_tapered_yue_dual_resist.py'
data = p.read_bytes()
sha = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
assert sha == '8c7692ebcfd4f413db17898b6faf0ff15c9f8f6f', sha
source = data.decode()
helper = '''def alignment_tip_region(tip_centers, orientations, dbu: float):
    """Exact 0.5 x 10 um alignment ticks: two at each taper end.

    Matches tapered_local.tapered_fiber_coupler's local tick definition.
    The returned region is in integer database units, with no size bias.
    """
    region = gf.kdb.Region()
    for tip, orientation in zip(tip_centers, orientations, strict=True):
        angle = radians(float(orientation) + 180.0)
        tangent = np.array([cos(angle), sin(angle)])
        normal = np.array([-sin(angle), cos(angle)])
        for side in (-1.0, 1.0):
            center = np.asarray(tip, dtype=float) + side * 7.5 * normal
            corners = [center + dx * tangent + dy * normal
                       for dx, dy in ((-0.25, -5.0), (0.25, -5.0),
                                      (0.25, 5.0), (-0.25, 5.0))]
            region.insert(gf.kdb.Polygon([
                gf.kdb.Point(int(round(point[0] / dbu)),
                             int(round(point[1] / dbu)))
                for point in corners
            ]))
    return region.merged()


def move_alignment_tips_to_layer15(component, tip_centers, orientations):
    """Move only the four existing ticks after all exposure masks are built.

    Late reassignment deliberately leaves the L2/L14 windows and L16
    protection mask unchanged. The input must be independently flattened;
    no cached child cells or other Layer-10/Layer-15 polygons are edited.
    """
    component.flatten()
    expected = alignment_tip_region(tip_centers, orientations,
                                    float(component.kcl.dbu))
    expected_polygons = list(expected.each())
    source_shapes = component.kdb_cell.shapes(gf.get_layer(NEGATIVE_RESIST_LAYER))
    destination = component.kdb_cell.shapes(gf.get_layer(RIDGE_LAYER))
    selected = []
    selected_region = gf.kdb.Region()
    for shape in source_shapes.each():
        if not (shape.is_polygon() or shape.is_box()):
            continue
        polygon = gf.kdb.Polygon(shape.box) if shape.is_box() else shape.polygon
        for target in expected_polygons:
            if polygon.bbox() == target.bbox() and (
                gf.kdb.Region(polygon) ^ gf.kdb.Region(target)
            ).is_empty():
                selected.append(shape)
                selected_region.insert(polygon)
                break
    if len(selected) != 2 * len(tip_centers) or not (
        selected_region ^ expected
    ).is_empty():
        raise ValueError("Expected exactly two complete Layer-10 ticks per taper end")
    # Copy each original shape exactly, then erase only that source shape.
    for shape in selected:
        destination.insert(shape)
        source_shapes.erase(shape)
    component.info["alignment_tip_layer"] = list(RIDGE_LAYER)
    component.info["alignment_tip_count"] = len(selected)


'''
anchor = 'def make_device(meander_straight_length: float) -> gf.Component:'
assert source.count(anchor) == 1
source = source.replace(anchor, helper + anchor)
old = '    component.flatten()\n    component.info["builder"] = "make_device"'
new = '    component.flatten()\n    move_alignment_tips_to_layer15(component, tips, connection_orientations)\n    component.info["builder"] = "make_device"'
assert source.count(old) == 1
p.write_text(source.replace(old, new))

# The central ridge must remain connected; the relocated alignment marks
# are deliberately separate polygons and must not count as broken ridges.
p = ROOT / 'validate_2611_sweep.py'
s = p.read_text()
old = '            components = list(regions[(layer, 0)].merged().each())'
new = '''            connected_region = regions[(layer, 0)].dup()
            if layer == 15:
                connected_region -= m.alignment_tip_region(
                    device.info['tip_centers'], device.info['connection_orientations'], dbu)
            components = list(connected_region.merged().each())'''
assert s.count(old) == 1
p.write_text(s.replace(old, new))
print('Applied tip-only L10 -> L15 reassignment and ridge connectivity expectation')
