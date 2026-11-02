"""One-time placement migration; no device geometry changes."""
from pathlib import Path
import ast

ROOT = Path(__file__).resolve().parent

def replace_once(text, old, new):
    if text.count(old) != 1:
        raise RuntimeError(f'Expected exactly one patch target: {old[:100]!r}')
    return text.replace(old, new, 1)

path = ROOT / '2611_tapered_yue_dual_resist.py'
original = path.read_text()
source = replace_once(original,
    '    """Call the same all-layer factory three times; only L and Y vary."""',
    '    """Place complete devices with aligned left bends and 600 um Y pitch."""')
source = replace_once(source,
    '    x_shift = LEFT_BEND_ALIGNMENT_X - float(devices[0].info["left_pnc_bend_x"])\n'
    '    y_positions = [BOTTOM_DEVICE_Y + i * DEVICE_Y_PITCH for i in range(len(devices))]\n'
    '    for device, y_position in zip(devices, y_positions, strict=True):\n'
    '        fin.add_ref(device).move((x_shift, y_position))',
    '    # Translate each COMPLETE device. Do not edit, clip, or remap any layer.\n'
    '    # The reference is the leftmost Layer-11 PNC edge at the left Euler bend.\n'
    '    # Device 1 stays in place; the longer devices move right to align here.\n'
    '    x_shifts = [\n'
    '        LEFT_BEND_ALIGNMENT_X - float(device.info["left_pnc_bend_x"])\n'
    '        for device in devices\n'
    '    ]\n'
    '    y_positions = [BOTTOM_DEVICE_Y + i * DEVICE_Y_PITCH for i in range(len(devices))]\n'
    '    for device, x_shift, y_position in zip(devices, x_shifts, y_positions, strict=True):\n'
    '        fin.add_ref(device).move((x_shift, y_position))')
source = replace_once(source,
    '    fin.info["device_x_shifts"] = [float(x_shift)] * len(devices)',
    '    fin.info["device_x_shifts"] = [float(value) for value in x_shifts]\n'
    '    fin.info["left_bend_alignment_x"] = LEFT_BEND_ALIGNMENT_X\n'
    '    fin.info["alignment_reference_layer"] = list(PNC_EBL_LAYER)\n'
    '    fin.info["alignment_reference"] = "leftmost_pnc_edge_of_left_euler_bend"')
source = replace_once(source,
    '    fin.info["placement_mode"] = "independent_all_layer_factory_fixed_endpoints"',
    '    fin.info["placement_mode"] = "independent_all_layer_factory_left_bend_aligned"')
# Prove that the complete all-layer factory and every geometry helper are untouched.
def nonplacement_tree(text):
    tree = ast.parse(text)
    tree.body = [node for node in tree.body if not
                 (isinstance(node, ast.FunctionDef) and node.name == 'build_layout')]
    return ast.dump(tree, include_attributes=False)
assert nonplacement_tree(original) == nonplacement_tree(source)
path.write_text(source)

path = ROOT / 'validate_2611_sweep.py'
source = path.read_text()
source = replace_once(source,
    "        check('export_actual_600um_y_pitch_and_same_x',\n"
    '              np.allclose(np.diff(np.asarray(actual_positions)[:, 1]), 600.0, rtol=0, atol=1e-9)\n'
    '              and np.ptp(np.asarray(actual_positions)[:, 0]) == 0,\n'
    '              device_origins_um=actual_positions)',
    '        left_bend_edges = [\n'
    '            kdb.Region(final.cell(name).begin_shapes_rec(final.layer(11, 0))).bbox().left\n'
    '            * final.dbu + position[0]\n'
    '            for name, position in zip(names, actual_positions, strict=True)\n'
    '        ]\n'
    "        check('export_actual_600um_y_pitch_and_aligned_left_bends',\n"
    '              np.allclose(np.diff(np.asarray(actual_positions)[:, 1]), 600.0, rtol=0, atol=1e-9)\n'
    '              and np.allclose(left_bend_edges, m.LEFT_BEND_ALIGNMENT_X, rtol=0, atol=1e-9),\n'
    '              device_origins_um=actual_positions, left_bend_edges_um=left_bend_edges)')
path.write_text(source)

path = ROOT / 'validate_2611_windows.py'
source = path.read_text()
source = replace_once(source,
    "    check('export_keeps_fixed_x_and_600um_pitch',\n"
    '          len(positions)==3 and np.ptp(np.asarray(positions)[:,0])==0 and\n'
    '          np.allclose(np.diff(np.asarray(positions)[:,1]),600,atol=1e-9,rtol=0), positions_um=positions)',
    '    left_bend_edges = [\n'
    '        kdb.Region(inst.cell.begin_shapes_rec(final.layer(11, 0))).bbox().left\n'
    '        * final.dbu + inst.trans.disp.x * final.dbu\n'
    '        for inst in instances\n'
    '    ]\n'
    "    check('export_keeps_aligned_left_bends_and_600um_pitch',\n"
    '          len(positions)==3 and\n'
    '          np.allclose(left_bend_edges,module.LEFT_BEND_ALIGNMENT_X,atol=1e-9,rtol=0) and\n'
    '          np.allclose(np.diff(np.asarray(positions)[:,1]),600,atol=1e-9,rtol=0),\n'
    '          positions_um=positions, left_bend_edges_um=left_bend_edges)')
path.write_text(source)
print('PLACEMENT_ONLY_PATCH_APPLIED; all geometry function ASTs unchanged')
