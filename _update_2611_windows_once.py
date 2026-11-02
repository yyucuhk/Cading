"""One-time source migration; removed after validated generation."""
from pathlib import Path
import hashlib

path = Path('2611_tapered_yue_dual_resist.py')
raw = path.read_bytes()
blob_sha = hashlib.sha1(f'blob {len(raw)}\0'.encode() + raw).hexdigest()
assert blob_sha == '0884b64ded623a6ddb2b23ca68b879dca59795c5', blob_sha
source = raw.decode()

start = source.index('    # Keep the Layer-2 companion windows around the two Layer-10 tip markers,')
stop = source.index('    layer14_central_pattern = gf.boolean(', start)
source = source[:start] + '''    # Remove the companion-window blocks at BOTH taper tips.  Work in
    # device-local coordinates so every straight-length variant receives
    # exactly the same two-ended cleanup. Layer-10 alignment tips are kept.
    layer14_tip_exclusions = gf.Component()
    tick_half_separation = 15.0 / 2 - ALIGNMENT_TIP_WINDOW_INWARD_SHIFT
    tick_window_x_offset = 0.25 + ALIGNMENT_TIP_WINDOW_WIDTH / 2
    tip_orientations = window_source.info["s_bend_connection_orientations"]
    for tip_start, tip_orientation in zip(
        layer14_start_points, tip_orientations, strict=True
    ):
        tip_angle = radians(tip_orientation)
        outward = np.array([-cos(tip_angle), -sin(tip_angle)])
        normal = np.array([-outward[1], outward[0]])
        for transverse_side in (-1.0, 1.0):
            for longitudinal_side in (-1.0, 1.0):
                window_center = (
                    np.asarray(tip_start, dtype=float)
                    + longitudinal_side * tick_window_x_offset * outward
                    + transverse_side * tick_half_separation * normal
                )
                exclusion = layer14_tip_exclusions.add_ref(
                    gf.components.rectangle(
                        size=(ALIGNMENT_TIP_WINDOW_WIDTH,
                              ALIGNMENT_TIP_WINDOW_LENGTH),
                        centered=True,
                        layer=Layers.via,
                    )
                )
                exclusion.rotate(tip_orientation + 180)
                exclusion.move(tuple(window_center))

''' + source[stop:]

old = '''    bridge = smooth_tip_bridge.add_ref(
        gf.components.rectangle(
            size=(bridge_length, bridge_width),
            centered=True,
            layer=positive_layer,
        )
    )
    bridge.rotate(orientation)
    bridge.move(tuple(np.asarray(start_point, dtype=float)))'''
new = '''    for tip_start, tip_orientation in zip(
        layer14_start_points, tip_orientations, strict=True
    ):
        bridge = smooth_tip_bridge.add_ref(
            gf.components.rectangle(
                size=(bridge_length, bridge_width),
                centered=True,
                layer=positive_layer,
            )
        )
        bridge.rotate(tip_orientation)
        bridge.move(tuple(np.asarray(tip_start, dtype=float)))'''
assert source.count(old) == 1
source = source.replace(old, new)
source = source.replace('Remove the annoying companion blocks near the larger-x taper tip',
                        'Remove the companion-window blocks at both taper tips')

helper = '''def fill_layer14_taper_seams(
    layer14: gf.Component,
    connection_centers,
    connection_orientations,
) -> gf.Component:
    """Fill enclosed L14 holes only inside the two taper center regions.

    Detect actual hole polygons, rather than drawing a guessed rectangle.
    The allowed region follows each port's local outward direction, from
    its meander connection to its slab tip. Larger support-pad openings
    farther out are explicitly outside this region and remain unchanged.
    Input cells are only read; the returned polygons are independently owned.
    """
    layer_index = gf.get_layer(Layers.via)
    dbu = float(layer14.kcl.dbu)
    region = gf.kdb.Region(
        layer14.kdb_cell.begin_shapes_rec(layer_index)
    ).merged()
    fills = gf.kdb.Region()
    filled_count = 0
    half_width = SLAB_EXTENSION_END_WIDTH / 2
    tolerance = dbu / 2
    frames = []
    for center, orientation in zip(
        connection_centers, connection_orientations, strict=True
    ):
        angle = radians(float(orientation))
        outward = np.array([-cos(angle), -sin(angle)])
        normal = np.array([-outward[1], outward[0]])
        frames.append((np.asarray(center, dtype=float), outward, normal))

    for polygon in region.each():
        for hole_index in range(polygon.holes()):
            points = list(polygon.each_point_hole(hole_index))
            coords = np.asarray([(p.x * dbu, p.y * dbu) for p in points])
            for center, outward, normal in frames:
                relative = coords - center
                longitudinal = relative @ outward
                transverse = relative @ normal
                inside_taper = (
                    np.all(longitudinal >= -tolerance)
                    and np.all(longitudinal <= SLAB_EXTENSION_LENGTH + tolerance)
                    and np.all(np.abs(transverse) <= half_width + tolerance)
                )
                if inside_taper:
                    fills.insert(gf.kdb.Polygon(points))
                    filled_count += 1
                    break

    repaired = gf.Component()
    repaired.kdb_cell.shapes(layer_index).insert((region | fills).merged())
    repaired.info["taper_seam_holes_filled"] = filled_count
    repaired.info["taper_seam_fill_area_um2"] = float(fills.area() * dbu**2)
    return repaired


'''
anchor = 'def add_shared_mla_window(\n'
assert source.count(anchor) == 1
source = source.replace(anchor, helper + anchor)

anchor = '''    # Replace only Layer 14 with the merged open-window geometry.
    result.flatten()'''
replacement = '''    # Fill the tiny enclosed taper-junction holes on Layer 14 only.
    # Do this after the union so mask fragments cannot recreate the seams.
    merged_layer14 = fill_layer14_taper_seams(
        merged_layer14,
        window_source.info["s_bend_connection_centers"],
        window_source.info["s_bend_connection_orientations"],
    )
    result.info["layer14_taper_holes_filled"] = int(
        merged_layer14.info["taper_seam_holes_filled"]
    )
    result.info["layer14_taper_fill_area_um2"] = float(
        merged_layer14.info["taper_seam_fill_area_um2"]
    )
    result.info["tip_window_cleanup"] = "both_ends_with_continuous_central_opening"

    # Replace only Layer 14 with independently owned merged polygons.
    result.flatten()'''
assert source.count(anchor) == 1
source = source.replace(anchor, replacement)
source = source.replace('shared MLA exposure. No Boolean geometry merge is applied.',
                        'shared MLA exposure. Layer 14 merges the cleaned release masks.')
assert 'tip_exclusion_index' not in source
compile(source, str(path), 'exec')
path.write_text(source)
print('Applied two-ended window cleanup and bounded Layer14 hole fill.')
