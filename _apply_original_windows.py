from pathlib import Path
import ast
for design,filename in ((1,'2611_tapered_yue_dual_resist.py'),(2,'2611_tapered_yue_design2.py')):
    p=Path(filename);text=p.read_text();old=ast.parse(text)
    anchor='POSITIVE_WINDOW_MARGIN = 3.0  # Layer-2 taper-window expansion'
    assert text.count(anchor)==1
    text=text.replace(anchor,anchor+'\nLAYER14_POST_PNC_MARGIN = 1.0  # per-side optical-taper clearance; PNC unchanged\nLAYER14_SIMPLE_TERMINAL_WINDOW = True  # remove terminal frame on L14 only',1)
    marker='def build_layout() -> gf.Component:'
    assert text.count(marker)==1
    stages='(RIDGE_TAPER_LENGTH, SLAB_TAPER_LENGTH)' if design==1 else '(TAPER_STAGE1_LENGTH, TAPER_STAGE2_LENGTH, TAPER_STAGE3_LENGTH)'
    wrapper='''def apply_layer14_post_pnc_window(component: gf.Component) -> gf.Component:
    """Final L14 policy for original layouts; raw make_device stays unchanged.

    Keep the full PNC support window, narrow the optical-taper window only,
    and replace the downstream rectangular frame with a linear channel.
    The L6000 family uses make_device and its own independent final policy.
    """
    if LAYER14_SIMPLE_TERMINAL_WINDOW:
        from original_design_layer14 import apply_post_pnc_layer14
        apply_post_pnc_layer14(component, PNC_TAPER_LENGTH,
                              STAGES, LAYER14_POST_PNC_MARGIN)
    return component


'''.replace('STAGES',stages)
    text=text.replace(marker,wrapper+marker,1)
    anchor='    devices = [make_device(length) for length in lengths]'
    assert text.count(anchor)==1
    text=text.replace(anchor,anchor+'\n    for device in devices:\n        apply_layer14_post_pnc_window(device)',1)
    new=ast.parse(text);a={n.name:ast.dump(n,include_attributes=False) for n in old.body if isinstance(n,ast.FunctionDef)}
    b={n.name:ast.dump(n,include_attributes=False) for n in new.body if isinstance(n,ast.FunctionDef)}
    assert all(a[n]==b[n] for n in a if n!='build_layout')
    p.write_text(text)
p=Path('2611_build_combined_sweep.py');text=p.read_text()
anchor='    removed_label = remove_original_length_label(d, m)'
assert text.count(anchor)==1
text=text.replace(anchor,anchor+'\n    m.apply_layer14_post_pnc_window(d)',1)
anchor="                tip_centers=d.info['tip_centers'], orientations=d.info['connection_orientations'],"
assert text.count(anchor)==1
text=text.replace(anchor,"                layer14_post_pnc_windows=d.info['layer14_post_pnc_windows'],\n"+anchor,1)
ast.parse(text);p.write_text(text)
print('Only original standalone exports and combined sweep final L14 policy changed.')
