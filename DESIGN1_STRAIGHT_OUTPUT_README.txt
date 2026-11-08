DESIGN 1 L6000: CURRENT SIX-DEVICE COMPARISON
=============================================
Run: python 2611_design1_L6000_straight_output.py
Output: gds/2611_design1_L6000_P10_straight_output_comparison.gds

Bottom three: two couplers, ridge widths 1.62 / 1.64 / 1.68 um.
Top three: tapered input and 2000 um constant-width straight output.
Each meander straight is 6000 um. Y pitch stays 600 um. All placements,
optical tapers, waveguides, support islands and alignment tips are retained.

CURRENT FINAL-EXPORT POLICY (BOTH SETS)
Layers 11, 16 and 51 are removed, every datatype in every cell.
Both Layer 2 and Layer 14 release masks begin at the actual Layer-10
entrance and extend outward through the existing terminal assemblies.
Two-coupler devices retain two windows; straight-output devices retain
only their tapered input window. No mask remains along the meander,
preceding transition, bends, or taper-free 2000 um output extension.
The straight output ridge itself is unchanged on Layer 15.

Layer 14 ONLY uses LAYER14_WINDOW_MARGIN = 1.0 um, instead of 3.0.
It is rebuilt using the original smooth mask construction with an
isolated margin setting. The broad first-stage entrance still joins
the inherited band and narrows smoothly into the 1 um-expanded taper.
Long channel and pad definitions are unchanged. Intersection with the
approved mask prevents any widening. Layer 2 is only cropped, not
resized. Its surviving terminal polygons are unchanged.

Labels remain on Layer 2, global alignment markers on Layer 15,
and Layer-10 markers and the frame remain unchanged. P10 in the
preserved labels is a legacy identifier; there are no exported PNCs.
Controls are in Section 01. Intermediate worker masks still contain
construction layers for baseline checks; only the final GDS is the
deliverable. Existing standalone designs, 18-device sweep and the
original three-device source GDS are NOT modified.

Current independent audit:
python validate_2611_coupler_release.py --render
Add --baseline PATH for before/after preservation checks.
Report: gds/2611_coupler_release_validation.json
Old raw/two-coupler validator entrypoints delegate to this audit.
Geometry checks do not establish optical loss, mechanical stability,
or full fabrication-rule compliance.

Engine: Python 3.11, gdsfactory 9.45.0, KLayout 0.30.13, NumPy 2.4.6.
Keep the source width-sweep GDS/manifest, two source builders, local
support modules and tether_variants.gds in the Cading directory.
