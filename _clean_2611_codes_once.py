"""One-time guarded source cleanup. Geometry is checked separately before merge."""
from pathlib import Path
import ast
import hashlib
import sys

EXPECTED = {
    '2611_tapered_yue_dual_resist.py': '03b8c0ee5b99ee954d3bbcd4965141c83330e536',
    '2611_tapered_yue_design2.py': '8ce5330b3fb0155c0e9cc16072d9f69934e26c47',
}


def section(number, title, description=''):
    return f'\n\n# %% [{number:02d}] {title}\n# ' + '=' * 77 + '\n' + (description + '\n' if description else '')


def replace_once(s, a, b):
    if s.count(a) != 1:
        raise ValueError(f'Expected one match, found {s.count(a)}: {a[:100]}')
    return s.replace(a, b, 1)


def clean(path):
    raw = path.read_bytes()
    sha = hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest()
    if sha != EXPECTED[path.name]:
        raise ValueError(f'Upstream changed; do not overwrite {path.name}: {sha}')
    old = raw.decode()
    tree = ast.parse(old)
    lines = old.splitlines()
    d2 = 'design2' in path.name
    prefix = 'design2_' if d2 else ''
    assigns = {n.targets[0].id: ast.get_source_segment(old, n) for n in tree.body
               if isinstance(n, ast.Assign) and isinstance(n.targets[0], ast.Name)}
    functions = {}
    for n in tree.body:
        if isinstance(n, ast.FunctionDef):
            start = min([n.lineno] + [d.lineno for d in n.decorator_list])
            functions[n.name] = '\n'.join(lines[start-1:n.end_lineno])
    used = {'REPOSITORY_ROOT', 'GDS', 'EBL_SUBFIELD_SIZE'}

    def param(name, comment=''):
        used.add(name)
        return assigns[name] + ('  # ' + comment if comment else '') + '\n'

    def params(items):
        return ''.join(param(*item) if isinstance(item, tuple) else param(item) for item in items)

    title = 'Design 2 - offset two-stage ridge' if d2 else 'Design 1 - aligned ridge taper'
    header = f'''"""{title}: editable, all-layer GDS generator.

QUICK START
    1. Edit Section 01 (all lengths/widths in um; angles in degrees).
    2. Run the WHOLE file in a fresh Python process after changing settings.
    3. The output path is printed at the end. Existing GDS is overwritten.

SECTION MAP
    00 imports/paths              06 coupler packaging and long channels
    01 USER PARAMETERS           07 Layer 2 / 14 / 16 exposure masks
    02 layers/derived values     08 move alignment tips to Layer 15
    03 taper polygon primitives  09 make_device(): one COMPLETE device
    04 PNC/meander support masks  10 chip placement and global markers
    05 two-layer taper assembly  11 settings checks and GDS export

DRAWING FLOW
    build_layout() -> make_device(L) -> taper/meander/support geometry
                   -> release/protection masks -> labels/tip-layer move
                   -> whole-device placement -> export_gds()

MEANDER_STRAIGHT_LENGTHS has one entry per device, ordered bottom to top.
Only the middle straight length varies per entry; other settings are shared.
Both ends are constructed together, not copied or clipped from another device.
Use CAD_EDITING_GUIDE.md for the parameter and layer reference.
"""
'''
    out = header + section(0, 'IMPORTS AND REPOSITORY PATHS') + '''from functools import partial
from math import cos, radians, sin
from pathlib import Path
import sys

import gdsfactory as gf
import numpy as np

REPOSITORY_ROOT = Path(__file__).resolve().parent
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))
GDS = REPOSITORY_ROOT / "gds"

from layers import Layers
from markers import calipermarker_array
from phononics import phononic_tether
from tapered_local import (
    edge_coupler,
    tapered_fiber_3d,
    tapered_fiber_coupler,
    tapered_fiber_inout,
)
from tether_variants import TETHER_FACTORIES
'''
    out += section(1, 'USER PARAMETERS - EDIT HERE', '''# Units: um for dimensions, degrees for angles, dimensionless for scales.
# Change values here, then rerun the entire file; decorated-cell defaults and
# the derived values in Section 02 are evaluated when the module is loaded.''')
    out += '\n# --- 1A. Device count, straight-length sweep, and whole-device placement ---\n'
    out += '# One entry = one device. Example: (200.0, 800.0, 1400.0).\n'
    out += params([
        ('MEANDER_STRAIGHT_LENGTHS', 'middle straight lengths; bottom to top'),
        ('DEVICE_Y_PITCH', 'origin-to-origin vertical pitch'),
        ('BOTTOM_DEVICE_Y', 'first device origin, NOT its lower taper tip'),
        ('LEFT_BEND_ALIGNMENT_X', 'common leftmost Layer-11 PNC bend edge'),
    ])
    out += '\n# --- 1B. Optical taper stages (outward from the 3 um PNC endpoint) ---\n'
    if d2:
        out += '''#                   Stage 1          Stage 2           Stage 3
# Layer 10:       3 -> 1.65         1.65 -> 0.980      0.980 -> 0.200
# Layer 15:       1 -> 0.500        0.500 -> 0.200     NO third stage
# Lengths:          20                 100                 10
# L10 and L15 stage 1 share the first length; stage 2 shares the second.
'''
        out += params([
            'TAPER_STAGE1_LENGTH', 'TAPER_STAGE2_LENGTH',
            ('TAPER_STAGE3_LENGTH', 'Layer 10 only'),
            'SLAB_STAGE1_END_WIDTH', 'SLAB_STAGE2_END_WIDTH',
            ('RIDGE_MEANDER_WIDTH', 'Layer-15 meander and stage-1 input width'),
            ('RIDGE_EXTENSION_START_WIDTH', 'Layer-15 stage-1 END width'),
            ('RIDGE_CONVERTER_END_WIDTH', 'Layer-15 terminal width'),
            ('RIDGE_CONVERTER_Y_SHIFT', 'final CENTERLINE offset, not edge gap'),
            ('REVIEW_TETHER_CONNECTION_WIDTH', 'Layer-10 terminal width AND tether neck'),
        ])
    else:
        out += '''# Layer 10: 3 -> 1 over RIDGE_TAPER_LENGTH, then 1 -> 0.2 over
# SLAB_TAPER_LENGTH. Layer 15: 1 -> 0.1 over RIDGE_TAPER_LENGTH, then ends.
# These are post-PNC lengths; the separate 15 -> 3 PNC region precedes them.
'''
        out += params([
            ('RIDGE_TAPER_LENGTH', 'common first-stage length for Layers 10/15'),
            ('SLAB_TAPER_LENGTH', 'second/final Layer-10 stage length'),
            ('SLAB_EXTENSION_START_WIDTH', 'Layer-10 intermediate width'),
            ('RIDGE_MEANDER_WIDTH', 'Layer-15 meander and taper input width'),
            ('RIDGE_EXTENSION_START_WIDTH', 'Layer-15 terminal width'),
            ('REVIEW_TETHER_CONNECTION_WIDTH', 'Layer-10 terminal width AND tether neck'),
        ])
    out += params([
        ('PNC_TAPER_SUPPORT_CUTOFF_WIDTH', 'Layer-10 input width after the PNC region'),
        ('TAPER_ANGLE', 'paired taper placement angle; validated design uses 0'),
        ('ADIABATIC_TAPER_POINTS', 'polygon samples per quintic transition'),
    ])
    out += '\n# --- 1C. Euler meander ---\n'
    out += params([
        ('S_BEND_COMPACT_180_RADIUS', 'effective radius of each 180-degree bend'),
        ('S_BEND_EULER_P', 'Euler fraction: 0 circular, 1 full Euler'),
        ('S_BEND_HALF_LENGTH', 'half the horizontal separation of the two ports'),
        ('S_BEND_WIDTH', 'construction/support envelope; NOT optical ridge width'),
        ('SLAB_EXTENSION_END_WIDTH', 'width at the start of the 15 -> 3 PNC transition'),
        ('S_BEND_TAPER_LENGTH', 'legacy noncompact-route setting; compact path ignores it'),
    ])
    out += '# Optical width is RIDGE_MEANDER_WIDTH in 1B.\n'
    out += '\n# --- 1D. Phononic-crystal (PNC) support pattern ---\n'
    out += params([
        ('REVIEW_SUPPORT_ANCHOR_PERIOD', 'requested PNC pitch on straight regions'),
        ('REVIEW_SUPPORT_BEND_PERIOD', 'requested PNC pitch on bends'),
        ('PHONONIC_SUPPORT_SEGMENTS', 'six square-mass sizes in the PNC library cell'),
        ('PHONONIC_SUPPORT_GAPS', 'gap between PNC masses'),
        ('PNC_NARROW_END_EXTENSION', 'narrow neck extension at both radial ends'),
    ])
    out += '# The PNC taper length is three straight periods; see Section 02.\n'
    out += '\n# --- 1E. Butterfly tethers and support islands ---\n'
    out += params([
        ('TETHER_VARIANT', '2 = split/forked butterfly in tether_variants.py'),
        ('TETHER_PAD_SCALE', 'library pad scale'),
        ('TETHER_CONNECTION_BLEND_LENGTH', 'neck blending length'),
        ('SUPPORT_RECTANGLE_SIZE', '(longitudinal length, transverse width)'),
        ('EXTRA_RECTANGLE_GAP', 'distance between opposing support islands'),
        ('SUPPORT_SURROUND_SLOT_WIDTH', 'Layer-2 slot around each support island'),
    ])
    out += '''# Optional extra side-support tips are disabled in make_device(). They
# are NOT the fiber-alignment ticks (which are kept on Layer 15).
'''
    out += params(['SUPPORT_TIP_POSITION', 'SUPPORT_TIP_WIDTH', 'SUPPORT_TIP_LENGTH', 'SUPPORT_TIP_OVERLAP'])
    out += '\n# --- 1F. Layer 2 / 14 release windows and Layer 16 protection ---\n'
    out += params([
        ('POSITIVE_WINDOW_MARGIN', 'Layer-2 taper-window expansion'),
        ('POSITIVE_WINDOW_END_EXTENSION', 'must be >= POSITIVE_WINDOW_MARGIN'),
        ('LAYER2_LAYER10_COPY_DELAY', 'inset of L10 copy from the PNC boundary'),
        ('RIDGE_PROTECTION_MARGIN', 'Layer-16 offset around ridge; excludes alignment ticks'),
        ('MLA_CHANNEL_ANGLE', 'full opening angle on Layer 2'),
        ('LAYER14_CHANNEL_ANGLE', 'full opening angle on Layer 14'),
        'MLA_CHANNEL_START_WIDTH', 'LAYER14_CHANNEL_START_WIDTH',
        ('LAYER14_CHANNEL_START_OFFSET', 'from stage-1 end, NOT final optical endpoint'),
        ('POSITIVE_CHANNEL_LENGTH', 'Layer-2 long-channel length'),
    ])
    out += 'LAYER14_CHANNEL_LENGTH = 2000.0  # Layer-14 long-channel length\n'
    out += '''# Dimensions below define REMOVED companion windows, not alignment ticks.
# Actual alignment ticks are fixed at 0.5 x 10 um in tapered_local.py and
# matched by alignment_tip_region(); edit both together to change tick shape.
'''
    out += params(['ALIGNMENT_TIP_WINDOW_WIDTH', 'ALIGNMENT_TIP_WINDOW_LENGTH', 'ALIGNMENT_TIP_WINDOW_INWARD_SHIFT'])
    out += '\n# --- 1G. Labels, chip outline, global overlay marks, and output ---\n'
    out += params(['TETHER_LABEL_SIZE', 'DETAIL_LABEL_SIZE', 'DETAIL_LABEL_TIP_X_OFFSET', 'DETAIL_LABEL_TIP_Y_OFFSET',
                   'LAYER10_MARKER_ARM_WIDTH', 'LAYER10_MARKER_CENTER_SQUARE_SIZE', 'MLA_MARKER_SIZE_INCREASE',
                   'MLA_MARKER_ARM_WIDTH', 'MLA_MARKER_CENTER_SQUARE_SIZE'])
    out += '''CHIP_WIDTH = 5000.0
CHIP_MIN_HEIGHT = 10000.0
CHIP_Y_ALLOWANCE = 1000.0
INNER_MARKER_CHIP_WIDTH_REDUCTION = 1000.0
GLOBAL_MARKER_ARM_LENGTHS = (50.0, 80.0)  # outer row, inner row
GDS_MAX_VERTEX_COUNT = 2000
'''
    out += 'TOP_CELL_NAME = "Design2_L800"\n' if d2 else 'TOP_CELL_NAME = f"{Path(__file__).stem}_L10_L12"\n'
    out += 'OUTPUT_GDS_FILENAME = "2611_tapered_yue_design2.gds"\n' if d2 else 'OUTPUT_GDS_FILENAME = f"{TOP_CELL_NAME}.gds"\n'
    out += '# TOP_CELL_NAME is a naming label, not a device-count or geometry control.\n'
    out += section(2, 'MASK LAYERS AND DERIVED VALUES - NORMALLY DO NOT EDIT')
    out += params([
        ('NEGATIVE_RESIST_LAYER', '10/0: slab/tethers; alignment ticks move off this layer'),
        ('FULL_ETCH_LAYER', '2/0: continuous positive-resist MLA mask'),
        ('POSITIVE_RESIST_LAYER', 'alias for the shared Layer-2 mask'),
        ('MARKER_LAYER', '2/0: labels and coarse alignment markers'),
        ('PNC_EBL_LAYER', '11/0: solid PNC structures'),
        ('RIDGE_LAYER', '15/0: optical ridge AND isolated alignment ticks'),
        ('RIDGE_PROTECTION_LAYER', '16/0: optical-ridge protection'),
        ('MEANDER_CONSTRUCTION_LAYER', 'temporary; never exported'),
    ])
    out += '# Other inherited layers: Layers.via = 14/0 (release), Layers.ring =\n# 13/0 (island fill), Layers.fiber = 51/0 (fiber guide), 0/0 = chip outline.\n\n'
    if d2:
        out += params(['DESIGN_ID', 'RIDGE_SEPARATION_DEFINITION', 'RIDGE_TAPER_LENGTH', 'SLAB_STAGE2_LENGTH', 'SLAB_TAPER_LENGTH',
                       'SLAB_EXTENSION_START_WIDTH', 'RIDGE_CONVERTER_LENGTH'])
    out += params(['S_BEND_VERTICAL_HALF_SPAN', 'PNC_TAPER_LENGTH', 'POST_PNC_TAPER_LENGTH',
                   ('SLAB_EXTENSION_LENGTH', 'stage-1 endpoint / alignment landmark, not total taper length'),
                   'PNC_TAPER_END_WIDTH', 'RIDGE_EXTENSION_END_WIDTH', 'PHONONIC_SUPPORT_OVERSHOOT', 'MLA_CHANNEL_START_OFFSET'])
    missing = set(assigns) - used
    if missing:
        raise ValueError(f'Unclassified settings: {sorted(missing)}')

    blocks = dict(functions)
    for name, s in list(blocks.items()):
        s = s.replace('# Make one continuous band around the complete meander. PNC cells are\n    # subtracted only after they have been placed as complete library cells,\n    # so no cell is clipped at a straight/Euler section boundary.',
                      '# Draw continuous Layer-2 side bands and separate solid Layer-11\n    # PNC cells. Keep whole cells at straight/Euler boundaries; no PNC-shaped\n    # holes are subtracted from the positive mask.')
        s = s.replace('# Clip the former Layer-12 shape at the Layer-15 endpoint before assigning\n    # it to the shared Layer-2 MLA mask.',
                      '# Start the expanded Layer-2 taper window at the END OF STAGE 1.\n    # For Design 2 this is not the final ridge endpoint. There is no Layer 12.')
        s = s.replace('# Copy only the device\'s meander/PNC release geometry.  Deliberately do\n    # not copy ``positive_window`` because that is the expanded taper-tip\n    # region on Layer 2 and is not wanted on Layer 14.',
                      '# First copy the central PNC/taper bands. The COMPLETE cleaned\n    # Layer-2 window is merged into Layer 14 later in this function.')
        s = s.replace('# exactly the same two-ended cleanup. Layer-10 alignment tips are kept.',
                      '# exactly the same two-ended cleanup. Alignment tips are retained\n    # here and moved from Layer 10 to Layer 15 at the end of make_device().')
        s = s.replace('# Ignore the old top/bottom tip blocks entirely and rebuild only the',
                      '# Rebuild the central opening after companion-window removal. Keep the')
        s = s.replace('# (optical taper width + 2*MLA margin), so it does not recreate the',
                      '# (stage-1 slab width + 2*MLA margin), so it does not recreate the')
        blocks[name] = s

    coupled = prefix + 'two_layer_tapered_fiber'
    for anchor, note in [
        ('    post_pnc_slab_taper =', '5A. Layer 10: first post-PNC slab stage'),
        ('    ridge_pnc_section =', '5B. Layer 15: constant ridge through PNC, then first taper'),
        ('    outer_slab_taper =', '5D. Layer 10: final slab taper, ending at the tether neck'),
        ('    release_tether =', '5E. Layer 2: transition windows and delayed slab copy'),
        ('    if add_phononic_taper_support:', '5F. Layers 11/14: PNC taper supports and full release opening'),
        ('    support_u =', '5G. Optional extra side-supports (disabled by current factory)'),
        ('    component.add_port(\n        "o1",', '5H. Ports and metadata used by coupler packaging'),
    ]:
        if anchor not in blocks[coupled]: raise ValueError(anchor)
        blocks[coupled] = blocks[coupled].replace(anchor, '    # --- ' + note + ' ---\n' + anchor, 1)
    if d2:
        blocks[coupled] = blocks[coupled].replace('    ridge_converter =', '    # --- 5C. Layer 15: shifted second stage; no third ridge stage ---\n    ridge_converter =', 1)
        blocks[coupled] = blocks[coupled].replace('    middle_slab_taper =', '    # --- 5C continued. Layer 10: independent second stage ---\n    middle_slab_taper =', 1)
    blocks['make_tapered_coupler'] = replace_once(blocks['make_tapered_coupler'], 'oxide_open_duplicate_length=2000.0,', 'oxide_open_duplicate_length=LAYER14_CHANNEL_LENGTH,')
    for anchor, note in [
        ('    positive_window = gf.Component()', '7A. Layer 2: expanded optical-taper exposure window'),
        ('    ridge_protection = gf.Component()', '7B. Layer 16: optical-ridge offset (before moving alignment tips)'),
        ('    layer14_start_points = []', '7C. Layer 14: central opening and local start coordinates'),
        ('    layer14_tip_exclusions = gf.Component()', '7D. Both ends: remove only companion-window blocks'),
        ('    for center, orientation, taper_width in zip(', '7E. Layer 2: extend the outer window to meet its long channel'),
        ('    layer2_tip_exclusions =', '7F. Clean Layer 2 and restore its continuous central openings'),
        ('    layer2_full_copy =', '7G. Layer 14: merge the ENTIRE cleaned Layer-2 mask'),
        ('    merged_layer14 = fill_layer14_taper_seams(', '7H. Layer 14 only: fill bounded taper-junction seams'),
        ('    result.info["mask_mode"] =', '7I. Record final mask metadata; construction layer already removed'),
    ]:
        if anchor not in blocks['add_shared_mla_window']: raise ValueError(anchor)
        blocks['add_shared_mla_window'] = blocks['add_shared_mla_window'].replace(anchor, '    # --- ' + note + ' ---\n' + anchor, 1)
    for anchor, note in [
        ('    tether_spec =', '9A. Configure butterfly tether and common meander geometry'),
        ('    device = tapered_fiber_inout(', '9B. Draw BOTH ends, meander, PNCs, tethers, and support pads'),
        ('    window_source = tapered_fiber_inout(', '9C. Separate source for exposure-window construction'),
        ('    component = add_shared_mla_window(', '9D. Complete Layer 2 / 14 / 16 masks on owned polygons'),
        ('    euler_length =', '9E. Length label and alignment-tip coordinates'),
        ('    move_alignment_tips_to_layer15(', '9F. Move alignment ticks AFTER all exposure/protection masks'),
        ('    component.info["builder"] =', '9G. Return complete independent geometry plus metadata'),
    ]:
        blocks['make_device'] = blocks['make_device'].replace(anchor, '    # --- ' + note + ' ---\n' + anchor, 1)

    for number, title, names, description in [
        (3, 'TAPER POLYGON PRIMITIVES', [prefix+'tapered_fiber_adiabatic'] + ([prefix+'offset_ridge_converter'] if d2 else []),
         '# Edit profiles here only when changing the taper law, not dimensions.'),
        (4, 'PNC SUPPORT AND MEANDER MASK GEOMETRY', [prefix+'phononic_meander_shield', prefix+'phononic_taper_shield'],
         '# Writes Layer 11 PNCs, Layer 2 bands, and optional solid Layer 14 fill.\n# The connected optical meander path itself comes from tapered_local.py.'),
        (5, 'TWO-LAYER TAPER ASSEMBLY', [coupled], '# Local +x points outward from the meander. Widths are along local y.\n# The second device end is placed by rotation; never hand-patch one end.'),
        (6, 'COUPLER PACKAGING AND LONG CHANNELS', ['make_tapered_coupler'],
         '# Uses tapered_local.py for fiber guide, glue spots, ticks and channels.'),
        (7, 'EXPOSURE AND PROTECTION MASKS: LAYERS 2 / 14 / 16', ['fill_layer14_taper_seams','add_shared_mla_window'],
         '# Flatten independent copies before editing. Do NOT recursively remove\n# layers from shared/cached child cells: that can erase other devices.'),
        (8, 'ALIGNMENT TICKS: MOVE LAYER 10 TO LAYER 15', ['alignment_tip_region','move_alignment_tips_to_layer15'],
         '# Two 0.5 x 10 um ticks per taper end; shapes/positions are unchanged.\n# This late move prevents them from changing the Layer-16 protection mask.'),
        (9, 'ONE COMPLETE DEVICE: make_device(meander_straight_length)', ['make_device'],
         '# This is the sole all-layer device factory. Every sweep entry calls it.'),
    ]:
        out += section(number, title, description)
        out += '\n\n\n'.join(blocks[name] for name in names) + '\n'

    layout = blocks['build_layout']
    a = layout.index('    # Chip outline and global overlay markers')
    b = layout.index('    fin.info[', a)
    chip_body = layout[a:b].replace('len(devices)', 'device_count')
    chip_body = chip_body.replace('(5000.0, max(10000.0, (device_count - 1) * DEVICE_Y_PITCH + 1000))',
                                  '(CHIP_WIDTH, max(CHIP_MIN_HEIGHT, (device_count - 1) * DEVICE_Y_PITCH + CHIP_Y_ALLOWANCE))')
    chip_body = chip_body.replace('(chip_shape, 50), ((chip_shape[0] - 1000, chip_shape[1]), 80)',
                                  '(chip_shape, GLOBAL_MARKER_ARM_LENGTHS[0]),\n        ((chip_shape[0] - INNER_MARKER_CHIP_WIDTH_REDUCTION, chip_shape[1]), GLOBAL_MARKER_ARM_LENGTHS[1])')
    chip_body = chip_body.replace('fin.', 'component.')
    chip_func = '''def add_chip_outline_and_markers(component: gf.Component, device_count: int) -> None:
    """Add only the chip outline and global overlay markers to a layout.

    All devices are already placed. No device polygons are edited here.
    """
''' + chip_body
    layout = layout[:a] + '    add_chip_outline_and_markers(fin, len(devices))\n' + layout[b:]
    layout = layout.replace('"""Place complete devices with aligned left bends and 600 um Y pitch."""',
                            '"""Call make_device for each length; align entire devices, then add chip marks."""\n    validate_settings()')
    layout = layout.replace('fin = gf.Component("Design2_L800")' if d2 else 'fin = gf.Component(f"{Path(__file__).stem}_L10_L12")',
                            'fin = gf.Component(TOP_CELL_NAME)')
    out += section(10, 'CHIP PLACEMENT AND GLOBAL MARKERS', '# Only whole-device translations occur here. All layers move together.')
    out += chip_func + '\n\n' + layout + '\n'
    out += section(11, 'SETTINGS CHECKS AND GDS EXPORT')
    out += '''def validate_settings() -> None:
    """Catch empty sweeps and invalid common dimensions before writing a GDS.

    These are basic input checks, NOT fabrication DRC or optical validation.
    """
    if not MEANDER_STRAIGHT_LENGTHS:
        raise ValueError("MEANDER_STRAIGHT_LENGTHS must contain at least one device")
    for length in MEANDER_STRAIGHT_LENGTHS:
        if not np.isfinite(float(length)) or float(length) < 0:
            raise ValueError("Meander straight lengths must be finite and non-negative")
    dimensions = {
        "DEVICE_Y_PITCH": DEVICE_Y_PITCH,
        "S_BEND_COMPACT_180_RADIUS": S_BEND_COMPACT_180_RADIUS,
        "RIDGE_MEANDER_WIDTH": RIDGE_MEANDER_WIDTH,
        "RIDGE_TAPER_LENGTH": RIDGE_TAPER_LENGTH,
        "SLAB_TAPER_LENGTH": SLAB_TAPER_LENGTH,
        "REVIEW_TETHER_CONNECTION_WIDTH": REVIEW_TETHER_CONNECTION_WIDTH,
        "POSITIVE_WINDOW_MARGIN": POSITIVE_WINDOW_MARGIN,
        "CHIP_WIDTH": CHIP_WIDTH,
        "CHIP_MIN_HEIGHT": CHIP_MIN_HEIGHT,
    }
'''
    if d2:
        out += '''    dimensions.update({
        "TAPER_STAGE2_LENGTH": TAPER_STAGE2_LENGTH,
        "SLAB_STAGE1_END_WIDTH": SLAB_STAGE1_END_WIDTH,
        "SLAB_STAGE2_END_WIDTH": SLAB_STAGE2_END_WIDTH,
        "RIDGE_CONVERTER_END_WIDTH": RIDGE_CONVERTER_END_WIDTH,
    })
'''
    out += '''    for name, value in dimensions.items():
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if not 0 <= S_BEND_EULER_P <= 1:
        raise ValueError("S_BEND_EULER_P must be between 0 and 1")
    if POSITIVE_WINDOW_END_EXTENSION < POSITIVE_WINDOW_MARGIN:
        raise ValueError("POSITIVE_WINDOW_END_EXTENSION must be >= POSITIVE_WINDOW_MARGIN")
    if not all(np.isfinite(v) for v in (BOTTOM_DEVICE_Y, LEFT_BEND_ALIGNMENT_X)):
        raise ValueError("Device placement coordinates must be finite")


def export_gds(output_path: str | Path | None = None) -> Path:
    """Build and write the complete layout; return the actual saved path.

    Default: GDS / OUTPUT_GDS_FILENAME. Supplying a path does not change
    geometry. An existing file at that path is overwritten.
    """
    path = Path(output_path) if output_path is not None else GDS / OUTPUT_GDS_FILENAME
    if path.suffix.lower() != ".gds":
        raise ValueError("The output filename must end in .gds")
    component = build_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    options = gf.kdb.SaveLayoutOptions()
    options.gds2_max_vertex_count = GDS_MAX_VERTEX_COUNT
    component.write_gds(path, save_options=options, with_metadata=False)
    return path


if __name__ == "__main__":
    saved_path = export_gds()
    print(f"Wrote {saved_path}")
'''
    ast.parse(out)
    path.write_text(out)
    print('CLEANED', path.name, len(lines), '->', len(out.splitlines()))


if __name__ == '__main__':
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    for name in EXPECTED:
        clean(root / name)
