"""One-time migration. Removed after the candidate has been verified."""
from pathlib import Path
import ast
import importlib.util
import sys
import textwrap

ROOT = Path(__file__).resolve().parent
TARGET = ROOT / '2611_tapered_yue_dual_resist.py'
OUT = ROOT / 'validation_2611'
OUT.mkdir(exist_ok=True)


def capture_baseline():
    import gdsfactory as gf
    import numpy as np
    spec = importlib.util.spec_from_file_location('baseline2611', TARGET)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    gf.gpdk.PDK.activate()
    original = m.add_shared_mla_window

    def capture(device, window_source, *args, **kwargs):
        before = {
            tuple(layer): gf.kdb.Region(device.kdb_cell.begin_shapes_rec(gf.get_layer(tuple(layer)))).merged()
            for layer in device.layers
        }
        result = original(device, window_source, *args, **kwargs)
        changed = []
        for layer, region in before.items():
            after = gf.kdb.Region(device.kdb_cell.begin_shapes_rec(gf.get_layer(layer))).merged()
            area = (region ^ after).area() * device.kcl.dbu**2
            if area:
                changed.append({'layer': list(layer), 'changed_area_um2': area})
        print('BASELINE_SHARED_CHILD_MUTATIONS', changed, flush=True)
        frozen = result.copy()
        frozen.flatten()
        euler_length = float(gf.path.euler(radius=m.S_BEND_COMPACT_180_RADIUS, angle=180,
                            p=m.S_BEND_EULER_P, use_eff=True).length())
        length = float(m.MEANDER_STRAIGHT_LENGTHS[0])
        total = 2 * m.S_BEND_HALF_LENGTH + 2 * length + 2 * euler_length
        endpoints = []
        for center, orientation in zip(window_source.info['s_bend_connection_centers'],
                                       window_source.info['s_bend_connection_orientations'], strict=True):
            angle = np.radians(orientation)
            endpoints.append(np.asarray(center) + (m.PNC_TAPER_LENGTH + m.RIDGE_TAPER_LENGTH)
                             * np.array([-np.cos(angle), -np.sin(angle)]))
        endpoint = min(endpoints, key=lambda p: p[0])
        label = frozen.add_ref(gf.components.text(
            text=f'ML{total:.1f}_TL{m.POST_PNC_TAPER_LENGTH + m.SLAB_TAPER_LENGTH:g}',
            size=m.DETAIL_LABEL_SIZE, layer=m.MARKER_LAYER))
        label.move((float(endpoint[0]) + m.DETAIL_LABEL_TIP_X_OFFSET,
                    float(endpoint[1]) + m.DETAIL_LABEL_TIP_Y_OFFSET))
        frozen.flatten()
        options = gf.kdb.SaveLayoutOptions()
        options.gds2_max_vertex_count = 2000
        frozen.write_gds(OUT / 'baseline_bottom.gds', save_options=options, with_metadata=False)
        return result

    m.add_shared_mla_window = capture
    m.MEANDER_STRAIGHT_LENGTHS = (200.0,)
    m.build_layout()
    assert (OUT / 'baseline_bottom.gds').stat().st_size > 0


def apply():
    source = TARGET.read_text()
    tree = ast.parse(source)
    functions = {node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)}
    assert 'make_device' not in functions
    assert 'splice_reference_with_swept_meander' in functions
    build = functions['build_layout']
    lines = source.splitlines(keepends=True)
    build_source = ''.join(lines[build.lineno - 1:build.end_lineno])
    start = build_source.index('        connection_width = REVIEW_TETHER_CONNECTION_WIDTH')
    stop = build_source.index('        # Build each sweep device as ONE intact component.')
    construction = textwrap.indent(textwrap.dedent(build_source[start:stop]), '    ')

    # Keep the approved first-device mask construction, but never edit its
    # shared child cells. Flatten into the NEW result before each replacement.
    prefix = source[:source.index('def repair_layer14_local_artifacts(')]
    for target in ('positive_layer', 'Layers.via', 'MEANDER_CONSTRUCTION_LAYER'):
        old = f'    result.remove_layers([{target}], recursive=True)'
        assert prefix.count(old) == 1, old
        prefix = prefix.replace(old,
            f'    result.flatten()\n    result.remove_layers([{target}], recursive=False)')
    end = '    return result\n'
    assert prefix.count(end) == 1
    prefix = prefix.replace(end, '    result.flatten()\n    return result\n')
    prefix = prefix.replace('shared MLA exposure. No Boolean geometry merge is applied.',
                            'shared MLA exposure. Each completed device owns all its mask polygons.')
    prefix = prefix.replace('# Align the left meander connection for every swept straight length. Three\n# 600 um EBL subfield columns can then share the same x boundaries.',
                            '# Use only the bottom device to establish the common X origin.\n# All three fixed end assemblies then remain vertically aligned.')

    new = '''def make_device(meander_straight_length: float) -> gf.Component:
    """Build one complete, independent all-layer device in local coordinates.

    Only meander_straight_length varies. Both taper/coupler assemblies,
    butterfly tethers, support pads, alignment tips, exposure windows, and
    all device layers are constructed here through the same helper calls.
    No geometry is copied from another sweep device or cut at a global X.
    """
    meander_straight_length = float(meander_straight_length)
    if not np.isfinite(meander_straight_length) or meander_straight_length < 0:
        raise ValueError("meander_straight_length must be finite and non-negative")
    gf.gpdk.PDK.activate()
'''
    new += construction
    new += '''
    # Snapshot placement information BEFORE assembling the independent masks.
    left_pnc_x = float(device.bbox(layer=gf.get_layer(PNC_EBL_LAYER)).left)
    connection_centers = [list(map(float, p)) for p in
                          window_source.info["s_bend_connection_centers"]]
    connection_orientations = list(map(float,
        window_source.info["s_bend_connection_orientations"]))
    component = add_shared_mla_window(device, window_source)
    component.flatten()

    euler_length = float(gf.path.euler(
        radius=S_BEND_COMPACT_180_RADIUS, angle=180,
        p=S_BEND_EULER_P, use_eff=True).length())
    total_length = (2 * S_BEND_HALF_LENGTH + 2 * meander_straight_length
                    + 2 * euler_length)
    tips = []
    for center, orientation in zip(connection_centers,
                                   connection_orientations, strict=True):
        angle = radians(orientation)
        tips.append(np.asarray(center) + SLAB_EXTENSION_LENGTH
                    * np.array([-cos(angle), -sin(angle)]))
    left_tip = min(tips, key=lambda point: float(point[0]))
    label = component.add_ref(gf.components.text(
        text=f"ML{total_length:.1f}_TL{POST_PNC_TAPER_LENGTH + SLAB_TAPER_LENGTH:g}",
        size=DETAIL_LABEL_SIZE, layer=MARKER_LAYER))
    label.move((float(left_tip[0]) + DETAIL_LABEL_TIP_X_OFFSET,
                float(left_tip[1]) + DETAIL_LABEL_TIP_Y_OFFSET))
    label_box = label.bbox()
    label_bounds = [float(label_box.left), float(label_box.bottom),
                    float(label_box.right), float(label_box.top)]
    component.flatten()
    component.info["builder"] = "make_device"
    component.info["meander_straight_length_um"] = meander_straight_length
    component.info["total_meander_length_um"] = total_length
    component.info["left_pnc_bend_x"] = left_pnc_x
    component.info["connection_centers"] = connection_centers
    component.info["connection_orientations"] = connection_orientations
    component.info["tip_centers"] = [list(map(float, p)) for p in tips]
    component.info["label_bbox_um"] = label_bounds
    component.info["mask_layers"] = [list(map(int, layer)) for layer in component.layers]
    return component


def build_layout() -> gf.Component:
    """Call the same all-layer factory three times; only L and Y vary."""
    gf.gpdk.PDK.activate()
    fin = gf.Component(f"{Path(__file__).stem}_L10_L12")
    lengths = tuple(float(length) for length in MEANDER_STRAIGHT_LENGTHS)
    devices = [make_device(length) for length in lengths]
    x_shift = LEFT_BEND_ALIGNMENT_X - float(devices[0].info["left_pnc_bend_x"])
    y_positions = [BOTTOM_DEVICE_Y + i * DEVICE_Y_PITCH for i in range(len(devices))]
    for device, y_position in zip(devices, y_positions, strict=True):
        fin.add_ref(device).move((x_shift, y_position))

    # Chip outline and global overlay markers are not part of the length sweep.
    chip_shape = (5000.0, max(10000.0, (len(devices) - 1) * DEVICE_Y_PITCH + 1000))
    fin.add_ref(gf.components.rectangle(size=chip_shape, centered=True, layer=(0, 0)))
    for marker_chip_shape, arm_length in (
        (chip_shape, 50), ((chip_shape[0] - 1000, chip_shape[1]), 80)
    ):
        fin.add_ref(calipermarker_array(
            chip_shape=marker_chip_shape, arm_length=arm_length,
            arm_width=LAYER10_MARKER_ARM_WIDTH,
            center_square_size=LAYER10_MARKER_CENTER_SQUARE_SIZE,
            layer1=NEGATIVE_RESIST_LAYER, layer2=None))
        fin.add_ref(calipermarker_array(
            chip_shape=marker_chip_shape, arm_length=arm_length + MLA_MARKER_SIZE_INCREASE,
            arm_width=MLA_MARKER_ARM_WIDTH, center_square_size=MLA_MARKER_CENTER_SQUARE_SIZE,
            layer1=MARKER_LAYER, layer2=None))
    fin.info["device_count"] = len(devices)
    fin.info["device_cells"] = [device.name for device in devices]
    fin.info["meander_straight_lengths"] = list(lengths)
    fin.info["total_meander_lengths"] = [float(d.info["total_meander_length_um"]) for d in devices]
    fin.info["device_x_shifts"] = [float(x_shift)] * len(devices)
    fin.info["device_y_positions"] = y_positions
    fin.info["device_y_pitch"] = DEVICE_Y_PITCH
    fin.info["bottom_device_y"] = BOTTOM_DEVICE_Y
    fin.info["placement_mode"] = "independent_all_layer_factory_fixed_endpoints"
    return fin


'''
    main_guard = source[source.index('if __name__ == "__main__":'):]
    updated = prefix + new + main_guard
    ast.parse(updated)
    assert 'splice_reference_with_swept_meander' not in updated
    assert 'tip_patch_size' not in updated
    assert 'recursive=True' not in updated
    TARGET.write_text(updated)
    print('REFACTOR_APPLIED', len(source.splitlines()), '->', len(updated.splitlines()), flush=True)


if __name__ == '__main__':
    if sys.argv[1:] == ['baseline']:
        capture_baseline()
    elif sys.argv[1:] == ['apply']:
        apply()
    else:
        raise SystemExit('Use baseline or apply')
