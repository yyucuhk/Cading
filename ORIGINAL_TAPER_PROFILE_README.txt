ORIGINAL TAPER LAW RESTORED -- 15 -> 5 um, then original downstream stages

All dimensions are micrometres. The two design source files again use the
original geometry functions from commit 6d4fcb56, with the PNC endpoint
width changed from 3 to 5 um. The full PNC transition is 45 um long.

For each transition, u=x/L and S(u)=10u^3-15u^4+6u^5:
    w(x) = w_in + (w_out-w_in)*S(u).
The width, slope and curvature are continuous at the 5 um join because
S'(0)=S'(1)=S''(0)=S''(1)=0. The width narrows monotonically. This is the
original taper family, not the recent nonzero-slope Hermite construction.

DESIGN 1: 15 -> 5 over 45; 5 -> 1 over 40; 1 -> 0.2 over 40.
Layer 15 retains 1 -> 0.1 over the first 40 um post-PNC stage.
DESIGN 2: 15 -> 5 over 45; 5 -> 1.65 over 20; 1.65 -> 0.980 over 100;
0.980 -> 0.200 over 10. Layer-15 stages and 0.840 um offset are retained.
COMBINED SWEEP: post-PNC totals 80/100/120, allocated by original ratios.
LONG WIDTH SWEEP: three 6000 um straights, widths 1.62/1.64/1.68 and 10 um PNC.

PNC reinforcement remains five 1 um intervals, with paired supports at
40,41,42,43,44,45 um at both ends. The normal PNC grid is preserved before
this group. Associated L2/L14 openings and L16 protection are regenerated.

RUN from Cading in fresh Python processes:
    python 2611_tapered_yue_dual_resist.py
    python 2611_tapered_yue_design2.py
    python 2611_build_combined_sweep.py
    python 2611_design1_L6000_P10_width_sweep.py
    python validate_2611_original_profile.py

The usual files under gds/ are replaced. Original dependent modules and
tether_variants.gds must remain in this directory. Use GDSFactory 9.45.0,
kfactory 2.6.1; KLayout 0.30.12 for the short layouts and 0.30.13 for the
long-width layout to reproduce prior mask-grid handling.

The analytical envelopes are C2; GDS uses polygons on a 1 nm grid.
Restoring the original profile does NOT independently verify optical
adiabaticity or insertion loss for the new widths. No optical or mechanical
performance certification is implied by the geometric checks.
