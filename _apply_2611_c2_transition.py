"""Guarded one-time migration: match the 5 um junction without moving geometry."""
from pathlib import Path
import ast

HELPER = '''def pnc_endpoint_width_derivatives():
    """Derivatives of the FULL slab width at the active PNC cutoff.

    The 15 -> 3 um reference quintic is truncated at width 5 um. Its
    derivative is not zero there. Units: dw/dx is dimensionless; d2w/dx2
    is 1/um. The two sidewall derivatives are +/- one half these values.
    """
    u = PNC_TAPER_LENGTH / PNC_TAPER_TRANSITION_LENGTH
    delta = PNC_TAPER_PROFILE_END_WIDTH - SLAB_EXTENSION_END_WIDTH
    slope = delta * 30 * u**2 * (1 - u)**2 / PNC_TAPER_TRANSITION_LENGTH
    curvature = delta * 60 * u * (1 - u) * (1 - 2*u) / PNC_TAPER_TRANSITION_LENGTH**2
    return float(slope), float(curvature)


def quintic_taper_coefficients(length, width1, width2, start_slope=0.0, start_curvature=0.0):
    """C2 Hermite polynomial in t=x/length, with a flat final cross-section.

    Match width, slope AND curvature at the entrance; set first and second
    derivatives to zero at the exit. Check all derivative extrema, not only
    sampled points, so unsupported parameter edits fail rather than create
    a neck/overshoot. This is a geometric check, not an optical-loss model.
    """
    values = (length, width1, width2, start_slope, start_curvature)
    if not all(np.isfinite(v) for v in values) or length <= 0 or min(width1, width2) < 0:
        raise ValueError("Taper dimensions/derivatives must be finite and widths nonnegative")
    delta = width2 - width1
    a = length * start_slope
    b = 0.5 * length**2 * start_curvature
    coeff = np.array([width1, a, b, 10*delta - 6*a - 3*b,
                      -15*delta + 8*a + 3*b, 6*delta - 3*a - b], dtype=float)
    d1 = np.polynomial.polynomial.polyder(coeff)
    roots = np.polynomial.polynomial.polyroots(np.polynomial.polynomial.polyder(d1))
    test_t = [0.0, 1.0] + [float(r.real) for r in roots
                           if abs(r.imag) < 1e-8 and 0 < r.real < 1]
    slopes = np.polynomial.polynomial.polyval(test_t, d1) / length
    if delta < 0 and max(slopes) > 1e-10:
        raise ValueError("Matched taper is not monotonic; revise its length or endpoint widths")
    if delta > 0 and min(slopes) < -1e-10:
        raise ValueError("Matched taper is not monotonic; revise its length or endpoint widths")
    if delta == 0 and max(abs(slopes)) > 1e-10:
        raise ValueError("Equal-width section cannot have nonzero matched derivatives")
    return coeff


def quintic_taper_widths(x, length, width1, width2, start_slope=0.0, start_curvature=0.0):
    """Shared profile for Layer 10 and the corresponding positive-mask edges."""
    t = np.asarray(x, dtype=float) / length
    if np.any(t < -1e-12) or np.any(t > 1 + 1e-12):
        raise ValueError("Sample coordinates must stay inside the taper")
    if start_slope == 0.0 and start_curvature == 0.0:
        # Preserve the exact approved arithmetic for all unaffected tapers.
        weight = 6*t**5 - 15*t**4 + 10*t**3
        return width1 + (width2 - width1) * weight
    coeff = quintic_taper_coefficients(length, width1, width2, start_slope, start_curvature)
    values = np.polynomial.polynomial.polyval(t, coeff)
    values = np.where(t == 0, width1, values)
    return np.where(t == 1, width2, values)


'''

def change(text, old, new, count=1):
    if text.count(old) != count:
        raise RuntimeError(f'Expected {count} anchors; found {text.count(old)}: {old[:90]}')
    return text.replace(old,new)

for fn in ['2611_tapered_yue_dual_resist.py','2611_tapered_yue_design2.py']:
    p=Path(fn); s=p.read_text()
    s=change(s,'# Edit profiles here only when changing the taper law, not dimensions.\n',
             '# Edit profiles here only when changing the taper law, not dimensions.\n'+HELPER)
    s=change(s,'    layer=NEGATIVE_RESIST_LAYER,\n) -> gf.Component:\n    """Quintic adiabatic taper with zero slope/curvature at both ends."""',
        '    layer=NEGATIVE_RESIST_LAYER,\n    start_slope: float = 0.0,\n    start_curvature: float = 0.0,\n) -> gf.Component:\n    """Smooth geometric taper; optional C2 entrance matching, flat exit.\n\n    Optical adiabaticity must be checked with local modes/EME or FDTD.\n    """')
    s=change(s,'    straight_length = length - taper_length\n',
        '    straight_length = length - taper_length\n    if straight_length > 0 and (start_slope != 0 or start_curvature != 0):\n        raise ValueError("Matched entrance cannot follow an internal flat lead")\n')
    s=change(s,'    # Quintic smootherstep: first and second derivatives vanish at u=0,1.\n    weight = 6 * u**5 - 15 * u**4 + 10 * u**3\n    transition_width = width1 + (width2 - width1) * weight',
        '    transition_width = quintic_taper_widths(\n        transition_x, taper_length, width1, width2, start_slope, start_curvature\n    )')
    s=change(s,'    component.info["profile"] = "quintic_adiabatic"\n    component.info["npoints"] = npoints',
        '    component.info["profile"] = ("quintic_C2_matched" if start_slope or start_curvature\n                                 else "quintic_adiabatic")\n    component.info["start_width_derivatives"] = [float(start_slope), float(start_curvature)]\n    component.info["npoints"] = npoints')
    anchor='    # --- 5A. Layer 10: first post-PNC slab stage ---\n'
    s=change(s,anchor,anchor+'    pnc_slope, pnc_curvature = pnc_endpoint_width_derivatives()\n')
    s=change(s,'            width2=SLAB_EXTENSION_START_WIDTH,\n            taper_length=POST_PNC_TAPER_LENGTH,\n            layer=slab_layer,',
        '            width2=SLAB_EXTENSION_START_WIDTH,\n            taper_length=POST_PNC_TAPER_LENGTH,\n            layer=slab_layer,\n            start_slope=pnc_slope,\n            start_curvature=pnc_curvature,')
    s=change(s,'    release_inner = start_inner * (1 - transition_weight)\n    release_outer = start_outer + (release_target_half_y - start_outer) * transition_weight',
        '    # Match the existing PNC-band edges, then flatten at the stage exit.\n    release_inner = quintic_taper_widths(transition_x, RIDGE_TAPER_LENGTH,\n        start_inner, 0.0, pnc_slope/2, pnc_curvature/2)\n    release_outer = quintic_taper_widths(transition_x, RIDGE_TAPER_LENGTH,\n        start_outer, release_target_half_y, pnc_slope/2, pnc_curvature/2)')
    s=change(s,'    copy_width = (\n        PNC_TAPER_END_WIDTH + (SLAB_EXTENSION_START_WIDTH - PNC_TAPER_END_WIDTH) * copy_weight\n    )',
        '    copy_width = quintic_taper_widths(\n        copy_x_local, POST_PNC_TAPER_LENGTH, PNC_TAPER_END_WIDTH,\n        SLAB_EXTENSION_START_WIDTH, pnc_slope, pnc_curvature\n    )')
    s=change(s,'    component.info["pnc_cutoff_width_um"] = PNC_TAPER_SUPPORT_CUTOFF_WIDTH\n',
        '    component.info["pnc_cutoff_width_um"] = PNC_TAPER_SUPPORT_CUTOFF_WIDTH\n    component.info["pnc_junction_continuity"] = "C2: matched width, slope and curvature"\n    component.info["pnc_junction_width_derivatives"] = [pnc_slope, pnc_curvature]\n    component.info["optical_adiabaticity_simulated"] = False\n')
    s=s.replace('    transition_u = transition_x / RIDGE_TAPER_LENGTH\n    transition_weight = 6 * transition_u**5 - 15 * transition_u**4 + 10 * transition_u**3\n','')
    s=s.replace('    copy_u = copy_x_local / POST_PNC_TAPER_LENGTH\n    copy_weight = 6 * copy_u**5 - 15 * copy_u**4 + 10 * copy_u**3\n','')
    s=s.replace('    # width is greater than 3 um.', '    # width is greater than 5 um.')
    s=s.replace('    # exact point where the narrow ridge taper begins.', '    # exact point where the narrow ridge taper begins (zero slope on L15).')
    ast.parse(s);p.write_text(s)

for fn in ['2611_build_combined_sweep.py','2611_design1_L6000_P10_width_sweep.py']:
    p=Path(fn);s=p.read_text()
    anchor='                pnc_cutoff_width_um=float(m.PNC_TAPER_SUPPORT_CUTOFF_WIDTH),\n'
    s=change(s,anchor,anchor+'                pnc_junction_continuity="C2",\n                pnc_junction_width_derivatives=list(m.pnc_endpoint_width_derivatives()),\n                optical_adiabaticity_simulated=False,\n')
    ast.parse(s);p.write_text(s)

Path('SMOOTH_5UM_JUNCTION_README.md').write_text('''# Smooth 5 um PNC-to-optical-taper junction

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
''')
print('Patched both generators and both sweep manifests.')
