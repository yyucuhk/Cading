L6000 FOUR-SET COMPARISON

Run: python 2611_design1_L6000_straight_output.py
Output: gds/2611_design1_L6000_P10_straight_output_comparison.gds

Bottom to top: three meander/two-coupler devices, three meander/one-coupler
devices, three single-straight/two-coupler devices, three single-straight/
one-coupler devices. Each set sweeps widths 1.62, 1.64, 1.68 um.
600 um row pitch; original six origins unchanged. New straight center
sections are 6000 um, NOT the approximately 18.742 mm total meander path.
One-coupler sets retain a 2000 um bare extension, giving 8000 um in the
new bend-free route, excluding the retained input transition/taper.

Layer 2 starts 0.5 um outward from Layer 10: global input X=-803.347 um
and output X=5287.653 um. Layer 14 starts with L10: -802.847/5287.153 um.
The 0.5 um original delayed-copy interface is retained, not doubled.
Entrance widths are TOTAL widths: L2=15 um, L14=8 um. Both downstream
expansions are 2.5 um per side. These margins do not add another 5 um
to the stated entrance widths. The smooth entrance blend joins the
original downstream mask generator; distant channels and pads stay unchanged.

Section 01 controls widths, start offsets and margins. The default width
sweep uses the approved source GDS; new widths require regenerating that
source and manifest first. Four-set default requires reference devices.
No L11, L16, L51 anywhere. Labels remain L2; global markers remain L15.
Old P10 text in the existing six labels is only a legacy source identifier.
Frame height increases for the additional rows; optical couplers and the
first six device placements remain unchanged.

Independent check: python validate_2611_four_sets.py --render
Geometry checks are not optical, mechanical or full fabrication DRC.
