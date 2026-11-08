from pathlib import Path
import ast
import hashlib

ROOT = Path.cwd()
BUILDER = ROOT / '2611_design1_L6000_straight_output.py'
VALIDATOR = ROOT / 'validate_2611_release_offset.py'

def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('Expected one exact match: '+old[:100])
    return text.replace(old, new, 1)

source = BUILDER.read_text()
blob = hashlib.sha1(b'blob '+str(len(BUILDER.read_bytes())).encode()+b'\0'+BUILDER.read_bytes()).hexdigest()
assert blob == 'dd46e5ef80ce58326c33f51f844bafe3c18052ca', blob
s = replace_once(source,
    'RELEASE_WINDOW_START_OFFSET = 10.0    # um OUTWARD after L10 starts; applies to BOTH masks',
    'RELEASE_WINDOW_START_OFFSET = 0.0     # start BOTH masks exactly at the Layer 10 entrance')
s = replace_once(s,
    'a configurable offset after Layer-10 entrances: two windows on two-coupler devices and one',
    'the Layer-10 entrances: two windows on two-coupler devices and one')
s = replace_once(s,
    'controls. Their retained masks begin 10 um outward from each Layer-10\nentrance, not at a 10 um waveguide width or at a global X coordinate.',
    'controls. Both masks begin exactly at each Layer-10 entrance (zero setback).\nThe optical taper and all downstream mask shapes remain unchanged.')
s = replace_once(s,
    '"""Narrow L2/L14 and start both masks after the L10 entrance.',
    '"""Narrow L2/L14 and align both masks with the L10 entrance by default.')
s = replace_once(s,
    'narrow both masks and apply a 10 um outward start offset',
    'narrow both masks and align their starts with the L10 entrances')
s = replace_once(s,
    'Input: unchanged optical taper; both masks 1 um expansion, start offset 10 um',
    'Input: unchanged optical taper; both masks 1 um expansion, starts aligned to L10')
ast.parse(s)
class Normalize(ast.NodeTransformer):
    def visit_Constant(self, node):
        if isinstance(node.value, str):
            return ast.copy_location(ast.Constant(value='TEXT'), node)
        return node
    def visit_Assign(self, node):
        if any(isinstance(t, ast.Name) and t.id == 'RELEASE_WINDOW_START_OFFSET' for t in node.targets):
            node.value = ast.Constant(value=0.0)
        return self.generic_visit(node)
assert ast.dump(Normalize().visit(ast.parse(source))) == ast.dump(Normalize().visit(ast.parse(s)))
BUILDER.write_text(s)

v = VALIDATOR.read_text()
v = replace_once(v, 'Measure both 1 um release margins and the 10 um longitudinal mask setback.',
                 'Measure both 1 um release margins and exact L2/L14/L10 start alignment.')
v = replace_once(v, 'non-mask layers, labels, markers, transforms, and the remaining Layer14 shape.',
                 'non-mask layers, labels, markers, transforms, and mask shapes beyond the former setback.')
v = replace_once(v, 'OFFSET = 10.0', 'OFFSET = 0.0\nPREVIOUS_OFFSET = 10.0  # baseline from the superseded setback revision')
v = replace_once(v, '        keep = Polygon()\n', '        keep = Polygon()\n        previous_keep = Polygon()\n')
v = replace_once(v,
    '            end_boxes.append(end)\n            keep = keep.union(end)',
    '''            end_boxes.append(end)
            keep = keep.union(end)
            previous_start = sx + sign*PREVIOUS_OFFSET
            previous_end = (box(-20000,25,previous_start,20000) if sign<0
                            else box(previous_start,-20000,20000,-25))
            previous_keep = previous_keep.union(previous_end)''')
v = replace_once(v, "_exact_10um_start'", "_starts_exactly_at_L10'")
v = replace_once(v,
    '''                lead = box(min(sx,start)+.001,y-20,max(sx,start)-.001,y+20)
                check(key+f'_L{layer}_no_mask_in_first10um', rs[layer].intersection(lead).area<1e-9)
                for distance in (10.1,15,20,30,35,39):''',
    '''                ahead = sx-sign*10
                upstream = box(min(sx,ahead)+.001,y-20,max(sx,ahead)-.001,y+20)
                check(key+f'_L{layer}_no_mask_before_L10', rs[layer].intersection(upstream).area<1e-9)
                for distance in (.001,.1,1,5,9.9,10.1,15,20,30,35,39):''')
v = replace_once(v,
    '''                check(name+f'_L{layer}_only_reduced',rs[layer].difference(before).area<1e-6)
                if layer==14:
                    exact(name+'_L14_only10um_setback_no_other_change',rs[14],before.intersection(keep))
                else:
                    check(name+'_L2_expansion_reduced_not_just_cropped',
                          before.intersection(keep).area-rs[2].area>100,
                          removed_after_crop_um2=before.intersection(keep).area-rs[2].area)''',
    '''                exact(name+f'_L{layer}_existing_mask_not_removed',before.difference(rs[layer]),Polygon())
                exact(name+f'_L{layer}_downstream_of_former10um_start_unchanged',
                      rs[layer].intersection(previous_keep),before)
                added = rs[layer].difference(before)
                restored_zone = keep.difference(previous_keep)
                exact(name+f'_L{layer}_additions_only_in_restored_first10um',
                      added.difference(restored_zone),Polygon())
                check(name+f'_L{layer}_former_setback_restored',added.area>1,
                      restored_area_um2=added.area)''')
v = replace_once(v, "offset_reference='Outward along taper from actual Layer10 entrance'",
                 "offset_reference='Zero: both masks start exactly at actual Layer10 entrance'")
v = replace_once(v, "print('PASS',len(checks),'independent mask-margin and setback checks')",
                 "print('PASS',len(checks),'independent mask-margin and common-start checks')")
v = replace_once(v,
    '''        ax.axvline(sx+sign*10,linestyle=':',linewidth=.6)
        ax.annotate('10 um mask setback',xy=(sx+sign*5,y+8),xytext=(sx+sign*45,y+8),
                    arrowprops={'arrowstyle':'->'},fontsize=10)''',
    '''        ax.annotate('L2 / L14 / L10 common start',xy=(sx,y+6),xytext=(sx+sign*45,y+8),
                    arrowprops={'arrowstyle':'->'},fontsize=10)''')
v = replace_once(v, "title='Both release masks: 1 um expansion; start 10 um after Layer 10 begins'",
                 "title='Both release masks: 1 um expansion; start aligned with Layer 10'")
ast.parse(v)
VALIDATOR.write_text(v)

p = ROOT / 'DESIGN1_STRAIGHT_OUTPUT_README.txt'
r = p.read_text()
r = replace_once(r, 'RELEASE_WINDOW_START_OFFSET = 10.0', 'RELEASE_WINDOW_START_OFFSET = 0.0')
r = replace_once(r,
    '''The 10 um offset is measured OUTWARD along the taper from its Layer10 entrance.
It does not move or shorten Layer10 or Layer15. The first 10 um has no L2/L14 mask.
Local starts: -3055 um (input), +3055 um (output).
Global starts: -812.847 um and +5297.153 um, respectively.''',
    '''Both release masks now start exactly at the Layer10 entrance (zero setback).
The former 10 um setback is removed. This restores only that first 10 um of
release-mask geometry; no optical geometry is moved or shortened, and the
existing downstream masks and both 1 um expansion settings are unchanged.
Local starts: -3045 um (input), +3045 um (output), identical to Layer10.
Global starts: -802.847 um and +5287.153 um, respectively.''')
r = replace_once(r, 'Audit: python validate_2611_release_offset.py --baseline PREVIOUS.gds --render',
    'Audit: python validate_2611_release_offset.py --render\nOptional --baseline PREVIOUS.gds compares against the previous 10 um-setback revision.')
p.write_text(r)
print('Patched zero offset; AST confirms only the offset changes executable builder behavior.')
