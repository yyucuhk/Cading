PNC CUTOFF AT 5 UM - CURRENT GEOMETRY
====================================
The original 45 um quintic 15-to-3 um reference profile is truncated at
its 5 um width point: active PNC length = 31.378626347 um. The 45 um value
is a reference-curve scale, NOT the active PNC length.

Five 1 um reinforcement periods end at this new cutoff. Six paired centers
are at 26.378626347 through 31.378626347 um. Normal straight, bend, and earlier
taper PNC positions are preserved. Both input/output ends use this rule.

Layer 10 starts at width 5 um. Design 1 retains 5 -> 1 -> 0.2 um stages;
Design 2 retains 5 -> 1.65 -> 0.980 -> 0.200 um. Downstream stage lengths,
ridge widths, the 840 nm Design-2 shift and all sweep values are unchanged.
Layer 15's constant PNC section shortens. Layers 2/14 windows and Layer 16
protection follow the new endpoint. Tapers, tips, tethers, pads and channels
move 13.621373653 um toward the meander. No device is clipped or spliced.

Run the complete scripts in a fresh Python process from Cading:
  python 2611_tapered_yue_dual_resist.py
  python 2611_tapered_yue_design2.py
  python 2611_build_combined_sweep.py
  python 2611_design1_L6000_P10_width_sweep.py
  python validate_2611_combined_sweep.py
  python validate_2611_pnc5um.py

Use GDSFactory 9.45.0 / kfactory 2.6.1 / NumPy 2.4.6; KLayout 0.30.12
for the first three outputs and 0.30.13 for the separate 6000 um output.
Each generator writes its previous gds/ filename. Back up before rebuilding.
Old validate_2611_taper_pitch.py targets the retired 3 um cutoff and is not
the validator for this revision. Use validate_2611_pnc5um.py, which measures
the saved polygons. Report: gds/2611_pnc5um_validation.json.
Validation concerns geometry only; it does not demonstrate optical loss
or mechanical performance and does not replace fabrication DRC.
