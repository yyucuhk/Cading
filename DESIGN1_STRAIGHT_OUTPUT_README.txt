L6000 DESIGN 1: FOUR-SET COMPARISON

Run: python 2611_design1_L6000_straight_output.py
Output: gds/2611_design1_L6000_P10_straight_output_comparison.gds
Keep linear_layer14_window.py alongside the main builder.

DEVICE SETS (bottom to top)
1. Three meanders with two couplers.
2. Three meanders with one coupler and a 2000 um bare output extension.
3. Three single 6000 um straights with two couplers, no bends.
4. Three single 6000 um straights with one coupler and a 2000 um bare output.
All sets sweep ridge widths 1.62, 1.64, 1.68 um; vertical pitch is 600 um.
The 6000 um value is each original meander straight, not total meander length.
No optical taper or waveguide routing changes in this revision.

RELEASE ENTRANCES
Layer 2: total entrance width 15 um, 2.5 um downstream margins, and the
original 0.5 um outward interface offset from Layer 10. Global input start
is -803.347 um; two-coupler output start is 5287.653 um.
Layer 14: total entrance width 8 um, starting with Layer 10 at global
-802.847/5287.153 um. Its upstream taper profile retains 1 um clearance.
Entrance widths are total widths: do not add the margins a second time.
The taper-free 2000 um outputs have no release windows.

LAYER 14 LINEAR CONTINUATION
Keep existing L14 through 58.847 um from the Layer-10 entrance. This is
output X=5346 um, or input X=-861.694 um. Its measured opening width at
that location is 2.644 um on the 1 nm layout grid.
Beyond that location, draw ONE linear trapezoid to the original outer
channel endpoint, joined to the two 300 um circles at their existing centers.
Remove the old intermediate rectangular frame and its cutouts entirely.
Circles use the original 2.5 degree sampling. The 1 um clearance applies to
the preserved upstream profile, not the downstream widening linear channel.
Layer 2, including its own original frame/channel shapes, is NOT changed.

LOCAL SUPPORTS AND LAYERS
Each of the 18 tapered ends retains three paired plain beams on Layer 11,
with 0.5 um longitudinal width and 1 um pitch (stations 0, 1, 2 um).
Their 5.096 um anchor span and 0.5 um overlap at each end are unchanged.
No PNC masses, meander supports or bare-output supports are restored.
Layers 16/51 stay absent. Labels remain on Layer 2; global markers on Layer 15.

EDITING AND CHECKS
Section 01 exposes LAYER14_LINEAR_START_FROM_L10=58.847 and
LAYER14_GLUE_DIAMETER=300.0. The main builder imports the isolated L14 helper.
Use the existing source width-sweep GDS/manifest and supporting Python modules.
python validate_2611_four_sets.py --render
python validate_2611_linear_window.py --baseline <previous-comparison.gds> --render
The package includes baseline/comparison.gds for the latter comparison.
These are geometry checks, not optical/mechanical qualification or full DRC.
The straight-sided opening is width-continuous; no C2 continuity is claimed.
