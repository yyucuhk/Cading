CURRENT SIX-DEVICE L6000 COMPARISON

Run from the existing Cading folder:
    python 2611_design1_L6000_straight_output.py

Section 01 controls:
    LAYER2_WINDOW_MARGIN = 1.0
    LAYER14_WINDOW_MARGIN = 1.0
    RELEASE_WINDOW_START_OFFSET = 0.0

Margin settings are in um per side. Both use the original smooth first-stage
mask blend and 1 um expansion on the downstream narrow taper. This does not
replace the entire first-stage opening with a constant offset.
Both release masks now start exactly at the Layer10 entrance (zero setback).
The former 10 um setback is removed. This restores only that first 10 um of
release-mask geometry; no optical geometry is moved or shortened, and the
existing downstream masks and both 1 um expansion settings are unchanged.
Local starts: -3045 um (input), +3045 um (output), identical to Layer10.
Global starts: -802.847 um and +5287.153 um, respectively.

Applies to both sets: six devices and nine coupler windows. Two windows are kept
on each two-coupler device; only the input window on each straight-output device.
The 2000 um straight outputs remain unchanged and have no release windows.
Long channels and pads keep their original definitions. Optical profiles, widths,
alignment tips, all six labels, device origins and 600 um Y pitch are unchanged.
Layers 11/16/51 remain absent; global markers remain on L15; labels stay on L2.
P10 in labels is a historical identifier, not an exported PNC pattern.

Output: gds/2611_design1_L6000_P10_straight_output_comparison.gds
Audit: python validate_2611_release_offset.py --render
Optional --baseline PREVIOUS.gds compares against the previous 10 um-setback revision.
Current report: gds/2611_release_offset_validation.json

The source GDS, standalone designs and 18-device sweep are not modified.
Keep supporting modules, source GDS and its manifest in the existing Cading folder.
Engine: gdsfactory 9.45.0, KLayout 0.30.13, kfactory 2.6.1, NumPy 2.4.6.
Geometry checks are not optical/mechanical simulation or full fabrication DRC.
