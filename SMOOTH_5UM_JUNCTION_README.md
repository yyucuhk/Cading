# Smooth 5 um PNC-to-optical-taper junction

The truncated reference PNC envelope arrives at width 5 um with nonzero
slope and curvature. The first Layer-10 stage now uses a quintic Hermite
polynomial matching that width, slope AND curvature, and ends with zero
slope and curvature to join the next stage. No lengths or endpoint widths
are changed. All earlier PNC supports, the dense five 1 um periods and the
5 um cutoff position are unchanged.

Layer-2 transition edges and its delayed Layer-10 copy use the same
matched profile; Layer 14 is rebuilt from the resulting masks. Layer 15
already joins its constant-width PNC section to a zero-slope quintic and
is unchanged, as is Layer 16. The approved laterally shifted converter in
Design 2 is not redesigned by this local junction edit.

The shared profile and analytic monotonicity guard are in Section 03 of
both generators. The first-stage call and the matching positive-mask edges
are in Section 05. The coefficients are computed from endpoint constraints,
not hand-adjusted. A future length/width edit that creates widening or an
undershoot raises an error rather than exporting it.

The mathematical profile is C2. GDS polygons approximate it on a 1 nm grid;
the exported-geometry report separately measures that approximation.
Smooth geometry is not proof of optical adiabaticity: wavelength, film and
etch thicknesses, material indices, polarization, local mode separation,
and the PNC perturbations are needed for EME/FDTD validation. No insertion
loss, reflection or higher-order-mode conversion has been simulated here.

Outputs keep their four existing names under gds/. Short layouts use
GDSFactory 9.45.0 and KLayout 0.30.12; the L6000 width sweep uses 0.30.13.
Run validate_2611_smooth_junction.py after regeneration. Its before/after
comparisons require saved old GDSs in validation_2611_smooth_junction/baseline.
