"""Create Design 2 without changing the approved Design 1 source."""
from pathlib import Path
import hashlib
import re
import ast

ROOT = Path(__file__).resolve().parent
base = ROOT / '2611_tapered_yue_dual_resist.py'
raw = base.read_bytes()
sha = hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest()
if sha != '03b8c0ee5b99ee954d3bbcd4965141c83330e536':
    raise RuntimeError(f'Design 1 changed; review before applying: {sha}')
s = raw.decode('utf-8')

def replace_one(old, new):
    global s
    if s.count(old) != 1:
        raise RuntimeError(f'Expected one patch anchor: {old[:90]!r}')
    s = s.replace(old, new, 1)

replace_one('"""Three-device two-EBL Euler-meander review.', '''"""Design 2: one 800 um Euler meander with two-layer offset converters.

Design 1 remains in 2611_tapered_yue_dual_resist.py and is not modified.
Layer 10 retains the Design-1 profile; its FINAL 1.0 -> 0.2 um taper is 100 um.
Layer 15 remains 1.0 um across the PNC, narrows 1.0 -> 0.6 um over the
inherited 40 um post-PNC transition, then 0.6 -> 0.2 um over 60 um while
its centerline shifts laterally by 0.840 um. Both ends use the same factory.
IMPORTANT: 0.840 um is interpreted as centerline-to-centerline offset,
NOT edge-to-edge gap; confirm this convention before fabrication.''')
replace_one('MEANDER_STRAIGHT_LENGTHS = (800,)', 'MEANDER_STRAIGHT_LENGTHS = (800.0,)')
replace_one('SLAB_TAPER_LENGTH = 40.0', 'SLAB_TAPER_LENGTH = 100.0  # Final Layer-10 1.0 -> 0.2 um taper only')
replace_one('RIDGE_EXTENSION_START_WIDTH = 0.1', 'RIDGE_EXTENSION_START_WIDTH = 0.6  # End of the 40 um first ridge transition')
replace_one('RIDGE_MEANDER_WIDTH = 1.0', '''RIDGE_MEANDER_WIDTH = 1.0
# DESIGN-2 CONVERTER SETTINGS (um). Widths are measured along local Y.
DESIGN_ID = "Design 2"
RIDGE_CONVERTER_LENGTH = 60.0
RIDGE_CONVERTER_END_WIDTH = 0.2
RIDGE_CONVERTER_Y_SHIFT = 0.840
RIDGE_SEPARATION_DEFINITION = "centerline_to_centerline"
# The image specifies the 60 um linear section. The earlier 1.0 -> 0.6 um
# transition retains the Design-1 40 um quintic profile (RIDGE_TAPER_LENGTH).
'''.rstrip())
anchor = '@gf.cell\ndef two_layer_tapered_fiber('
helper = '''@gf.cell
def design2_offset_ridge_converter(
    length: float = RIDGE_CONVERTER_LENGTH,
    width1: float = RIDGE_EXTENSION_START_WIDTH,
    width2: float = RIDGE_CONVERTER_END_WIDTH,
    y_shift: float = RIDGE_CONVERTER_Y_SHIFT,
    layer=RIDGE_LAYER,
) -> gf.Component:
    """Linear width taper and linear lateral shift, matching the sketch.

    The input cross-section is at x=0, y=0. The output is at x=length,
    y=y_shift. Widths are vertical cross-section widths, not perpendicular
    to the oblique centerline. The far end is intentionally a taper tip.
    """
    if not all(np.isfinite(v) for v in (length, width1, width2, y_shift)):
        raise ValueError("Converter dimensions must be finite")
    if min(length, width1, width2) <= 0:
        raise ValueError("Converter length and widths must be positive")
    component = gf.Component()
    component.add_polygon([
        (0.0, -width1 / 2),
        (length, y_shift - width2 / 2),
        (length, y_shift + width2 / 2),
        (0.0, width1 / 2),
    ], layer=layer)
    component.add_port("o1", center=(0.0, 0.0), width=width1,
                       orientation=180, layer=layer)
    component.add_port("o2", center=(length, y_shift), width=width2,
                       orientation=0, layer=layer)
    component.info["length_um"] = float(length)
    component.info["widths_um"] = [float(width1), float(width2)]
    component.info["centerline_offset_um"] = float(y_shift)
    component.info["profile"] = "linear_width_and_centerline"
    return component


'''
replace_one(anchor, helper + anchor)
replace_one('    ridge_extension.connect("o1", ridge_pnc_section.ports["o2"])', '''    ridge_extension.connect("o1", ridge_pnc_section.ports["o2"])

    # The 600 -> 200 nm displaced converter starts exactly at the first
    # ridge taper endpoint. It is not spliced into or clipped out of a mask.
    ridge_converter = component.add_ref(design2_offset_ridge_converter())
    ridge_converter.connect("o1", ridge_extension.ports["o2"])
    component.info["ridge_converter_end_um"] = [
        float(value) for value in ridge_converter.ports["o2"].center
    ]''')
replace_one('    """Piecewise-adiabatic slab and ridge transition after the meander.\n\n    Along +x, the 15->3 um PNC transition is excluded from Layer 10. Layer 10\n    starts again at the 3 um endpoint, then tapers 3->1 um while Layer 15\n    begins its aligned 1->0.1 um taper. Layer 10 finally continues from\n    1 um to the tether width.\n    """', '''    """Design-2 two-layer coupler after the unchanged PNC region.

    Layer 10: 3 -> 1 um over 40 um, then 1 -> 0.2 um over 100 um.
    Layer 15: 1 -> 0.6 um over 40 um, then a laterally shifted
    0.6 -> 0.2 um linear converter over 60 um.
    """''')
replace_one('        text=f"ML{total_length:.1f}_TL{POST_PNC_TAPER_LENGTH + SLAB_TAPER_LENGTH:g}",', '        text=f"DESIGN 2_L{meander_straight_length:g}_ML{total_length:.1f}",')
replace_one('    component.info["builder"] = "make_device"', '''    component.info["builder"] = "make_device"
    component.info["design"] = DESIGN_ID
    component.info["layer10_final_taper_length_um"] = SLAB_TAPER_LENGTH
    component.info["ridge_stage1_length_um"] = RIDGE_TAPER_LENGTH
    component.info["ridge_stage2_length_um"] = RIDGE_CONVERTER_LENGTH
    component.info["ridge_width_sequence_um"] = [1.0, 0.6, 0.2]
    component.info["ridge_tip_centerline_offset_um"] = RIDGE_CONVERTER_Y_SHIFT
    component.info["ridge_separation_definition"] = RIDGE_SEPARATION_DEFINITION
    converter_ends = []
    for center, orientation in zip(connection_centers, connection_orientations, strict=True):
        angle = radians(float(orientation) + 180.0)
        tangent = np.array([cos(angle), sin(angle)])
        normal = np.array([-sin(angle), cos(angle)])
        converter_ends.append(list(map(float, np.asarray(center)
            + (PNC_TAPER_LENGTH + RIDGE_TAPER_LENGTH + RIDGE_CONVERTER_LENGTH) * tangent
            + RIDGE_CONVERTER_Y_SHIFT * normal)))
    component.info["ridge_converter_end_centers_um"] = converter_ends''')
replace_one('    fin = gf.Component(f"{Path(__file__).stem}_L10_L12")', '    fin = gf.Component("Design2_L800")')
replace_one('    fin.info["device_count"] = len(devices)', '    fin.info["design"] = DESIGN_ID\n    fin.info["device_count"] = len(devices)')
replace_one('        os.path.join(GDS, f"{fin.name}.gds"),', '        os.path.join(GDS, "2611_tapered_yue_design2.gds"),')
# Avoid collisions with Design 1 cached component names in one Python session.
for name in ('phononic_meander_support', 'phononic_meander_shield',
             'phononic_taper_shield', 'tapered_fiber_adiabatic', 'two_layer_tapered_fiber'):
    s = re.sub(r'\b' + name + r'\b', 'design2_' + name, s)
ast.parse(s)
(ROOT / '2611_tapered_yue_design2.py').write_text(s)
assert base.read_bytes() == raw
print('DESIGN_2_CREATED; DESIGN_1_SOURCE_UNCHANGED')
