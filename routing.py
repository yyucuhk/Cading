"""
Routers connecting gdsfactory ports with waveguides.

* ``connect_intersection`` -- single euler bend filleting the port-ray
  intersection, with a ``connect_ata`` fallback (draws into a component;
  ``route0`` is a deprecated alias).
* ``connect_spline`` -- smooth spline with an adaptive corner radius
  (``route`` is a deprecated alias).
* ``connect_circular`` -- two equal-radius arcs (ported from pGDS CIRCULAR).
* ``connect_ata`` -- arc-tangent-arc at a fixed radius; shortest of the four
  left/right combinations (ported from pGDS::connectThickenShortestDistance);
  falls back to ``connect_aaa`` for close-quarters wrap-arounds.
* ``connect_aaa`` -- three arcs (Dubins CCC) at a fixed radius; the loop-free
  alternative where ``connect_ata`` self-intersects.
* ``connect_manhattan`` -- axis-aligned straights filleted with 90 deg bends
  at a fixed radius; requires both port orientations on the 90 deg grid.
  Turns cluster at the *port1* end, so a fan of parallel routes never
  crosses (shape policy in ``_solve_manhattan_waypoints``).
* ``connect`` -- unified entry point dispatching to the above (plus
  ``"linear"``, a straight fill between the ports) and layering a three-width
  policy on top: each port is tapered to ``bend_width`` for the bends and
  long straights are widened to ``routing_width``, every step joined by a
  ``taper_length`` linear taper.

All routers share the gdsfactory port convention: orientations are *outward*
normals, so a route leaves ``port1`` along ``port1.orientation`` and arrives
at ``port2`` traveling against ``port2.orientation`` (head-on into the port,
as if the two ports faced each other).

Arc geometry is solved by ``_solve_circular`` / ``_solve_ata`` /
``_solve_aaa`` / ``_solve_manhattan_waypoints`` (pure NumPy); every solver
result is verified by forward simulation (``_simulate_segments``) before
being returned, and is then realized by ``_place_segments``; unequal widths
are tapered along the path with ``gf.path.extrude_transition``.

The port-level helpers (``extend``, ``bend_towards`` and its wrappers,
``get_direction``, ``get_intersection``, ``flipped``, ...) live in
:mod:`layout.core.ports` -- also available as methods on every ``gf.Port`` --
and are re-exported here for backwards compatibility.
"""
import numbers
import warnings
from dataclasses import dataclass
from functools import partial
from typing import Literal

import numpy as np

import gdsfactory as gf
from gdsfactory import logger
from gdsfactory.port import Port
from gdsfactory.typings import CrossSectionSpec

from layout.core.ports import (  # noqa: F401  (re-exported public API)
    extend, extend_to, bend_towards, bend_west, bend_east,
    bend_north, bend_south, bend_by, get_direction, get_intersection, flipped,
    translated, mirrored, measure,
)
from layout.core.ports import (
    _DEFAULT_BEND_RADIUS, _component_length, _copy_info, _draw_straight,
    _split_angle, get_pathlength,
)

bend_euler = partial(gf.components.bend_euler_all_angle)
bend_circular = partial(gf.components.bend_circular_all_angle, angular_step=1)
bend_default = bend_euler


# Anti-parallel and collinear tolerances: past the first a "straight" would kink at each
# port, and a residual kink past the second is drawn but warned about.
_COLLINEAR_TOL_COS = float(np.cos(np.deg2rad(0.5)))
_COLLINEAR_WARN_SIN = float(np.sin(np.deg2rad(0.05)))
# Arcs below this angle are skipped as numerically zero.
_ZERO_ANGLE_DEG = 1e-6
#: How far a placed route may end from port2 before the routers warn (um).
MISS_TOLERANCE = 0.005
# Safety margin applied to the largest fillet radius a control polygon accepts.
_FILLET_MARGIN = 0.95


class _FallbackToAta(Exception):
    """Internal: the strict single-fillet route cannot serve this pair.

    Raised by `_connect_intersection_strict` and caught only by
    `connect_intersection`, which then falls back to `connect_ata` -- a typed
    signal instead of exception-message matching.
    """


def _classify_collinear(p1, d1, p2, d2, caller, stacklevel=3):
    """Classify a port pair against the 0.5 deg collinearity policy.

    Returns ``"facing"`` (a straight fill is acceptable; warns through
    *caller* when the residual kink exceeds 0.05 deg), ``"back_to_back"``
    (both ports point away from each other along one line), or ``None``.
    *p1*/*p2* are positions, *d1*/*d2* unit direction vectors; the pair must
    not be coincident.
    """
    D = p2 - p1
    dhat = D / float(np.linalg.norm(D))
    a1 = float(dhat @ d1)
    a2 = float(dhat @ d2)
    if a1 > _COLLINEAR_TOL_COS and a2 < -_COLLINEAR_TOL_COS:
        # Either port may carry the residual kink; report the worse one.
        kink_sin = max(abs(dhat[0] * d1[1] - dhat[1] * d1[0]),
                       abs(dhat[0] * d2[1] - dhat[1] * d2[0]))
        if kink_sin > _COLLINEAR_WARN_SIN:
            warnings.warn(
                f"{caller}: ports are anti-parallel but laterally offset; the "
                f"straight fill kinks {np.degrees(np.arcsin(kink_sin)):.3f} deg "
                "at each port",
                stacklevel=stacklevel,
            )
        return "facing"
    if a1 < -_COLLINEAR_TOL_COS and a2 > _COLLINEAR_TOL_COS:
        return "back_to_back"
    return None


def connect_intersection(
    component: gf.Component,
    port1: Port,
    port2: Port,
    taper_length: float = 100.0,
    routing_width: "float | None" = None,
    *,
    radius: "float | None" = None,
    **kwargs,
) -> float:
    """Route *port1* to *port2* by filleting their ray intersection, drawing
    into *component*.

    Collinear opposing ports get a straight fill; ports whose rays cross
    forward of both get a single euler bend at the corner; everything else
    (same orientation, anti-parallel, no forward intersection, bend or taper
    does not fit) falls back to `connect_ata` with a warning, so the route is
    always drawn.  Unequal port widths are joined with a *taper_length* taper;
    *routing_width* widens a long straight.

    Parameters
    ----------
    component : gf.Component
        Component to draw all routing geometry into.
    port1, port2 : Port
        Endpoints (module-level outward-normal convention).
    taper_length : float
        Length of each linear width taper.
    routing_width : float, optional
        Widened width for the straight sections.
    radius : float, optional
        Fillet radius, overriding the radius baked into *bend*; it also seeds
        the ata fallback.  ``None`` (default) leaves the bend's own radius.
    **kwargs
        Forwarded to `bend_towards` (``bend=`` sets the fillet, and its
        radius also seeds the ata fallback; default ``_DEFAULT_BEND_RADIUS``).
        ``cross_section=`` (optional) extrudes every straight/taper/bend with
        that cross-section instead of a bare quad on ``port1.layer``.

    Returns
    -------
    float
        Total route path length.  The primitives stamp ``pathlength`` on the
        ports they draw; the route-level ``o2`` stamp comes from `connect`.
    """
    if "add_to_component" in kwargs:
        raise TypeError(
            "connect_intersection: geometry is drawn into the *component* "
            "argument; add_to_component is not accepted")

    try:
        return _connect_intersection_strict(
            component, port1, port2, taper_length=taper_length,
            routing_width=routing_width, radius=radius, **kwargs)
    except _FallbackToAta as e:
        reason = str(e)

    # Fallback bend factory and radius from the bend= keyword.
    bend = kwargs.get("bend")
    ata_bend = bend_default
    ata_radius = _DEFAULT_BEND_RADIUS
    if isinstance(bend, numbers.Real):
        ata_radius = float(bend)
    elif bend is not None:
        ata_bend = bend
        baked = getattr(bend, "keywords", {}).get("radius")
        if baked is not None:
            ata_radius = float(baked)
    # An explicit radius= (forwarded to bend_towards) overrides the bend's own.
    if radius is not None:
        ata_radius = float(radius)

    warnings.warn(
        f"connect_intersection: falling back to connect_ata "
        f"(radius {ata_radius:g}) -- {reason}",
        stacklevel=2,
    )
    c = connect_ata(
        port1, port2, radius=ata_radius, bend=ata_bend,
        cross_section=kwargs.get("cross_section"),
        routing_width=routing_width, taper_length=taper_length)
    component.add_ref(c)
    return float(c.info["length"])


def _connect_intersection_strict(
    component: gf.Component,
    port1: Port,
    port2: Port,
    taper_length: float = 100.0,
    routing_width: "float | None" = None,
    *,
    radius: "float | None" = None,
    **kwargs,
) -> float:
    """Single-fillet core of `connect_intersection`; raises `_FallbackToAta`
    where the public wrapper falls back.  Nothing is drawn into *component*
    before any raise point, so the fallback never leaves partial geometry."""
    # Popped so the bend probe and fillet are stamped at the right width.
    cross_section = kwargs.pop("cross_section", None)
    dport = np.array(port2.center, float) - np.array(port1.center, float)
    dist = float(np.linalg.norm(dport))
    if dist < 1e-9:
        raise ValueError("connect_intersection: port1 and port2 are coincident")

    d1 = get_direction(port1)
    d2 = get_direction(port2)
    topology = _classify_collinear(np.array(port1.center, float), d1,
                                   np.array(port2.center, float), d2,
                                   caller="connect_intersection",
                                   stacklevel=4)  # helper -> strict -> public -> user

    if topology == "facing":
        if routing_width is not None:
            n_tapers = (int(not np.isclose(routing_width, port1.width))
                        + int(not np.isclose(routing_width, port2.width)))
            if dist < n_tapers * taper_length:
                raise _FallbackToAta(
                    f"cannot fit {n_tapers} routing-width taper(s) of "
                    f"{taper_length} um in a {dist:.3f} um span")
            if not np.isclose(routing_width, port1.width):
                port1 = extend(port1, distance=taper_length, new_width=routing_width,
                               add_to_component=component, cross_section=cross_section)
            if not np.isclose(routing_width, port2.width):
                port2 = extend(port2, distance=taper_length, new_width=routing_width,
                               add_to_component=component, cross_section=cross_section)
        elif not np.isclose(port1.width, port2.width):
            # Taper at the thinner port and run the wider width.
            if dist < taper_length:
                raise _FallbackToAta(
                    f"cannot fit a {taper_length} um width taper in a "
                    f"{dist:.3f} um span")
            if port1.width < port2.width:
                port1 = extend(port1, distance=taper_length, new_width=port2.width,
                               add_to_component=component, cross_section=cross_section)
            else:
                port2 = extend(port2, distance=taper_length, new_width=port1.width,
                               add_to_component=component, cross_section=cross_section)

        _draw_straight(component, port1, port2, cross_section)

        return float(dist)

    if topology == "back_to_back":
        raise _FallbackToAta(
            "back-to-back collinear ports need a loop, which a single bend "
            "cannot make")
    if float(d1 @ d2) > _COLLINEAR_TOL_COS:
        raise _FallbackToAta(
            "ports share the same orientation; a single bend cannot route "
            "them")

    # Dry probe for where the bend lands; unbind it, else it draws a phantom bend.
    probe = port1.copy()
    probe.info["bind_component"] = ""
    port12 = bend_towards(probe, orientation=180 + port2.orientation,
                          radius=radius, **kwargs)
    bend_intersection = get_intersection(port1, port12)

    if bend_intersection is None:
        raise _FallbackToAta(
            "anti-parallel ports with a lateral offset need an S-bend, which "
            "a single bend cannot make")

    L_bend = float(d1 @ (bend_intersection - port1.center))
    if L_bend < 0:
        raise ValueError(
            "connect_intersection: internal error: the probe bend lands "
            "behind port1")

    intersection = get_intersection(port1, port2)
    if intersection is None:
        raise _FallbackToAta(
            "the port rays are parallel and never intersect, so there is no "
            "corner to fillet")

    L1 = float(d1 @ (intersection - port1.center))
    L2 = float(d2 @ (intersection - port2.center))

    if min(L1, L2) <= 0:
        raise _FallbackToAta(
            f"the port rays intersect behind a port (L1={L1:.1f}, "
            f"L2={L2:.1f}); there is no forward corner to fillet")
    if min(L1, L2) < L_bend:
        raise _FallbackToAta(
            f"bend does not fit between the ports (L1={L1:.1f}, L2={L2:.1f}, "
            f"L_bend={L_bend:.1f})")

    # port1_bend/port2_bend are the quad endpoints; bend_port seeds bend_towards.
    port1_bend = port1.copy()
    port1_bend.center = port1.center + (L1 - L_bend) * d1
    bend_port = port1_bend.copy()
    port2_bend = port2.copy()
    port2_bend.center = port2.center + (L2 - L_bend) * d2

    if not np.isclose(port1.width, port2.width):
        if max(L1, L2) < taper_length + L_bend:
            raise _FallbackToAta(
                f"no room for a {taper_length} um width taper "
                f"(L1={L1:.1f}, L2={L2:.1f}, L_bend={L_bend:.1f})")

        w_lo = min(port1.width, port2.width)   # o1 is the small end
        w_hi = max(port1.width, port2.width)   # o2 is the large end
        if cross_section is None:
            t = gf.components.taper(length=taper_length, width1=w_lo,
                                    width2=w_hi, layer=port1.layer)
        else:
            xs = gf.get_cross_section(cross_section)
            t = gf.path.extrude_transition(gf.path.straight(taper_length),
                gf.path.transition(xs.copy(width=w_lo), xs.copy(width=w_hi),
                                   width_type="linear"))
        t1 = component.add_ref(t)

        # Place the taper at the thin port, or against the bend when only the long leg has room.
        thin1 = port1.width < port2.width
        at_bend = (min(L1, L2) <= taper_length + L_bend
                   and (L1 > L2) != thin1)
        if at_bend:
            bend_w = min(port1.width, port2.width)
            if L1 > L2:
                t1.connect("o1", port1_bend.copy(width=bend_w).flipped(),
                           allow_type_mismatch=True)
                port1_bend = t1.ports["o2"].copy()
            else:
                t1.connect("o1", port2_bend.copy(width=bend_w).flipped(),
                           allow_type_mismatch=True)
                port2_bend = t1.ports["o2"].copy()
        else:
            bend_w = max(port1.width, port2.width)
            if thin1:
                t1.connect("o1", port1, allow_type_mismatch=True)
                port1 = t1.ports["o2"].copy()
            else:
                t1.connect("o1", port2, allow_type_mismatch=True)
                port2 = t1.ports["o2"].copy()
            port1_bend = port1_bend.copy(width=bend_w)
            port2_bend = port2_bend.copy(width=bend_w)
        bend_port = bend_port.copy(width=bend_w)

    bent = bend_towards(
        bend_port,
        orientation=180 + port2.orientation,
        add_to_component=component,
        cross_section=cross_section,
        radius=radius,
        **kwargs
    )
    # Arc length of the fillet, read off the port pathlength accounting.
    bend_arc = get_pathlength(bent) - get_pathlength(bend_port)
    distance = bend_arc + (L1 + L2 - 2 * L_bend)

    if routing_width is not None:
        # Room checks use the *current* port positions; a width taper may have moved one.
        len1 = float((np.asarray(port1_bend.center) - np.asarray(port1.center)) @ d1)
        if len1 > 2 * taper_length:
            port1 = extend(port1, distance=taper_length, new_width=routing_width,
                           add_to_component=component, cross_section=cross_section)
            # A width taper against the bend leaves its port facing back down the leg.
            back1 = -taper_length * np.sign(np.asarray(get_direction(port1_bend)) @ d1)
            port1_bend = extend(port1_bend, distance=back1, new_width=routing_width,
                                add_to_component=component, cross_section=cross_section)
        len2 = float((np.asarray(port2_bend.center) - np.asarray(port2.center)) @ d2)
        if len2 > 2 * taper_length:
            port2 = extend(port2, distance=taper_length, new_width=routing_width,
                           add_to_component=component, cross_section=cross_section)
            back2 = -taper_length * np.sign(np.asarray(get_direction(port2_bend)) @ d2)
            port2_bend = extend(port2_bend, distance=back2, new_width=routing_width,
                                add_to_component=component, cross_section=cross_section)

    _draw_straight(component, port1, port1_bend, cross_section)
    _draw_straight(component, port2_bend, port2, cross_section)

    return float(distance)


# Deprecated alias kept for existing tapeout scripts.
route0 = connect_intersection


def _intersection_segments(port1, port2, radius) -> "list | None":
    """Single-fillet centerline segments anchored at *port1*.

    The corner at the port-ray crossing is filleted with an arc of *radius*,
    which eats ``radius * tan(theta/2)`` of each leg.  ``None`` for the pairs
    `_connect_intersection_strict` hands to arc-tangent-arc.
    """
    p0, d0 = np.array(port1.center, float), get_direction(port1)
    p1, d1 = np.array(port2.center, float), get_direction(port2)
    topology = _classify_collinear(p0, d0, p1, d1, caller="connect_intersection")
    if topology == "facing":
        return [("straight", float(np.linalg.norm(p1 - p0)))]
    if topology == "back_to_back" or float(d0 @ d1) > _COLLINEAR_TOL_COS:
        return None
    x = get_intersection(port1, port2)
    if x is None:
        return None
    L1 = float(d0 @ (x - p0))
    L2 = float(d1 @ (x - p1))
    t = -d1  # arrival heading after the fillet
    angle = float(np.degrees(np.arctan2(d0[0] * t[1] - d0[1] * t[0], d0 @ t)))
    L_bend = radius * np.tan(np.radians(abs(angle)) / 2.0)
    if min(L1, L2) <= 0 or min(L1, L2) < L_bend:
        return None
    return [("straight", L1 - L_bend), ("arc", angle, radius),
            ("straight", L2 - L_bend)]


def _max_fillet_radius(points) -> float:
    """Largest corner radius ``gf.path.smooth`` can fit on polyline *points*.

    smooth() replaces each interior corner of turn angle theta with an arc
    that consumes ``R*tan(theta/2)`` of tangent length on *both* adjacent
    segments, so each segment must accommodate the demands of its two ends.
    Returns 0.0 for an exact reversal (theta = 180 deg, no radius fits) and
    ``inf`` when the points are collinear (no corner consumes anything).
    """
    pts = np.asarray(points, float)
    vecs = np.diff(pts, axis=0)
    lens = np.linalg.norm(vecs, axis=1)
    if np.any(lens < 1e-12):
        return 0.0
    dirs = vecs / lens[:, None]
    tans = []
    for i in range(len(dirs) - 1):
        theta = np.arccos(float(np.clip(dirs[i] @ dirs[i + 1], -1.0, 1.0)))
        if theta > np.pi - 1e-6:
            return 0.0  # exact reversal: the fillet tangent length diverges
        tans.append(np.tan(theta / 2.0))
    r_max = np.inf
    for i, seg_len in enumerate(lens):
        demand = (tans[i - 1] if i > 0 else 0.0) + (tans[i] if i < len(tans) else 0.0)
        if demand > 1e-12:
            r_max = min(r_max, seg_len / demand)
    return float(r_max)


def _spline_path(port1: Port, port2: Port, *, radius: float = 100.0,
                 smoothing=gf.path.smooth) -> gf.Path:
    """Smoothed centerline from *port1* to *port2* (see `connect_spline`)."""
    p1 = np.array(port1.center, dtype=float)
    p2 = np.array(port2.center, dtype=float)
    d1 = get_direction(port1)
    d2 = get_direction(port2)

    if float(np.linalg.norm(p2 - p1)) < 1e-9:
        raise ValueError("connect_spline: port1 and port2 are coincident")

    topology = _classify_collinear(p1, d1, p2, d2, caller="connect_spline")
    if topology == "back_to_back":
        raise ValueError(
            "connect_spline: back-to-back collinear ports cannot be routed "
            "with a spline; use connect(..., method='ata')"
        )

    if topology == "facing":
        return smoothing(points=[p1, p2], radius=radius)
    else:
        # Control arms on the *outward* side of each port, so the spline meets the port
        # normals; capped at half the separation, else the polygon collapses into a reversal.
        w_max = max(port1.width, port2.width)
        dist = float(np.linalg.norm(p2 - p1))
        arm = min(dist / 2.0, max(2.0 * radius, 4.0 * w_max))
        points = [p1, p1 + arm * d1, p2 + arm * d2, p2]
        # *radius* is a maximum: shrink it to what the control polygon accepts.
        r_fit = _FILLET_MARGIN * _max_fillet_radius(points)
        r_used = min(radius, r_fit)
        if r_used + 1e-12 < min(radius, w_max):
            raise ValueError(
                f"connect_spline: the control polygon needs a hairpin "
                f"reversal the spline cannot fillet (max feasible radius "
                f"{r_fit:.3f} um < waveguide width {w_max:g} um); use "
                "connect(..., method='ata')"
            )
        if r_used < 0.5 * radius:
            warnings.warn(
                f"connect_spline: smoothing radius reduced from {radius:g} "
                f"to {r_used:.3f} um to fit the control polygon",
                stacklevel=2,
            )
        # Only the default smoother takes a bend= keyword; others get (points, radius).
        if smoothing is gf.path.smooth:
            return smoothing(points=points, radius=r_used, bend=gf.path.arc)
        return smoothing(points=points, radius=r_used)


def connect_spline(
    port1: Port,
    port2: Port,
    *,
    smoothing=gf.path.smooth,
    radius: float = 100.0,
    cross_section: "CrossSectionSpec | None" = None,
) -> gf.Component:
    """Create a smooth spline waveguide component from *port1* to *port2*.

    Builds a four-point control polygon whose end tangents match the ports
    (module-level outward-normal convention: head-on arrival at *port2*),
    then fillets it with *smoothing*.  Collinear opposing ports get a
    straight path; unequal port widths are tapered along the spline (warns).

    Parameters
    ----------
    port1, port2 : Port
        Endpoints.
    smoothing : callable
        Path-smoothing function ``smoothing(points, radius) -> gf.Path``.
        Defaults to ``gf.path.smooth``.
    radius : float
        *Maximum* corner radius: automatically reduced to the largest fillet
        the control polygon accepts (a severe reduction warns; a true hairpin
        reversal, where no usable radius exists, raises ``ValueError``).
        Also scales the control-arm length (see the arm comment in the
        source), so larger radii sweep wider curves.
    cross_section : CrossSectionSpec, optional
        Cross-section for the spline; ``None`` inherits width/layer from
        *port1*.  Its width is replaced by the port widths.

    Returns
    -------
    gf.Component
        Extruded waveguide with ``o1``/``o2`` ports and ``info["length"]``;
        ``o2`` continues *port1*'s ``pathlength`` chain.
    """
    path = _spline_path(port1, port2, radius=radius, smoothing=smoothing)

    resolved = gf.get_cross_section(cross_section) if cross_section is not None else None

    def make_cs(w):
        return (resolved.copy(width=w) if resolved is not None
                else gf.cross_section.cross_section(width=w, layer=port1.layer))

    tapering, _ = _decide_widths(port1.width, port2.width, None,
                                 caller="connect_spline")
    if not tapering:
        extruded = gf.path.extrude(path, cross_section=make_cs(port1.width))
    else:
        transition = gf.path.transition(
            make_cs(port1.width), make_cs(port2.width), width_type="linear")
        extruded = gf.path.extrude_transition(path, transition)

    c = gf.Component()
    ref = c.add_ref(extruded)
    c.add_ports(ref.ports)
    # Both ports continue port1's chain; o2 adds the route length.
    for name in ("o1", "o2"):
        _copy_info(port1.info, c.ports[name].info)
    c.ports["o2"].info["pathlength"] = get_pathlength(port1) + float(path.length())
    c.info["length"] = float(path.length())
    return c


# Deprecated alias kept for existing tapeout scripts.
route = connect_spline


# Path planners: pure-geometry solvers and drawing helpers, from pGDS/primitives.cpp.

@dataclass
class _Arc:
    """Parameters for a single circular arc segment."""
    angle_deg: float  # signed: positive = CCW, negative = CW
    radius: float     # always positive


# Pure-geometry helpers (no gdsfactory)

def _perp_ccw(v: np.ndarray) -> np.ndarray:
    """90 deg CCW rotation: [x, y] -> [-y, x]."""
    return np.array([-v[1], v[0]])


def _rotate2d(v: np.ndarray, angle: float) -> np.ndarray:
    """Rotate a 2D vector by *angle* radians CCW."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([v[0]*c - v[1]*s, v[0]*s + v[1]*c])


def _directed_angle_deg(u: np.ndarray, v: np.ndarray, ccw: bool) -> float:
    """Directed angle in degrees from vector *u* to vector *v* in a fixed sense.

    Returns a value in [0, 360) when *ccw* and in (-360, 0] otherwise;
    0 when the vectors coincide.  Unlike an arccos + cross-product-side
    heuristic this is unambiguous for every configuration.
    """
    a0 = np.arctan2(u[1], u[0])
    a1 = np.arctan2(v[1], v[0])
    if ccw:
        return float(np.degrees(np.mod(a1 - a0, 2.0 * np.pi)))
    return -float(np.degrees(np.mod(a0 - a1, 2.0 * np.pi)))


def _sample_segments(p, d, segments, step_deg: float = 5.0) -> "tuple[np.ndarray, np.ndarray]":
    """Polyline samples along ``("arc", angle_deg, radius)`` /
    ``("straight", length)`` segments from point *p*, unit direction *d*.

    Returns ``(points, end_direction)``.  Arc steps are exact rotations, so
    the samples (and in particular the end point) carry no discretization
    error regardless of *step_deg*.
    """
    pts = [np.asarray(p, float).copy()]
    p = np.asarray(p, float).copy()
    d = np.asarray(d, float).copy()
    for seg in segments:
        if seg[0] == "arc":
            _, angle_deg, radius = seg
            if abs(angle_deg) < 1e-12:
                continue
            s = 1.0 if angle_deg > 0 else -1.0
            c = p + s * radius * _perp_ccw(d)
            n = max(2, int(np.ceil(abs(angle_deg) / step_deg)))
            for k in range(1, n + 1):
                pts.append(c + _rotate2d(p - c, np.deg2rad(angle_deg * k / n)))
            p = pts[-1].copy()
            d = _rotate2d(d, np.deg2rad(angle_deg))
        else:
            p = p + d * seg[1]
            pts.append(p.copy())
    return np.array(pts), d


def _simulate_segments(p: np.ndarray, d: np.ndarray, segments) -> "tuple[np.ndarray, np.ndarray]":
    """End point and end direction of arc/straight *segments* from (p, d).

    This is the exact geometry `_place_arcs` realizes, so a solver candidate
    that simulates onto the target is guaranteed to place correctly (for
    endpoint-faithful bends such as arcs or eulers with
    ``with_arc_floorplan=True``).
    """
    pts, d_end = _sample_segments(p, d, segments, step_deg=180.0)
    return pts[-1], d_end


def _crosses_ray(points, origin, direction, eps: float = 1e-3) -> bool:
    """True when polyline *points* crosses the ray from *origin* along
    *direction* at a distance > *eps* from the origin.

    Sign-change walk on the lateral coordinate that skips exactly-on-ray
    samples: a path endpoint sitting at the ray origin never counts, while a
    genuine crossing that happens to land on a sample point still does.
    """
    u = np.asarray(direction, float)
    rel = np.asarray(points, float) - np.asarray(origin, float)
    lat = rel[:, 0] * u[1] - rel[:, 1] * u[0]
    along = rel @ u
    last_sign, last_i = 0.0, -1
    for i in range(len(rel)):
        sign = np.sign(lat[i])
        if sign == 0.0:
            continue
        if last_sign != 0.0 and sign != last_sign:
            t = lat[last_i] / (lat[last_i] - lat[i])
            if along[last_i] + t * (along[i] - along[last_i]) > eps:
                return True
        last_sign, last_i = sign, i
    return False


def _self_intersects(points, eps: float = 1e-9) -> bool:
    """True when two non-adjacent segments of polyline *points* cross."""
    pts = np.asarray(points, float)
    n = len(pts) - 1
    for i in range(n - 2):
        a = pts[i]
        ab = pts[i + 1] - a
        c = pts[i + 2:n]
        cd = pts[i + 3:n + 1] - c
        denom = ab[0] * cd[:, 1] - ab[1] * cd[:, 0]
        ac = c - a
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (ac[:, 0] * cd[:, 1] - ac[:, 1] * cd[:, 0]) / denom
            u = (ac[:, 0] * ab[1] - ac[:, 1] * ab[0]) / denom
        if ((np.abs(denom) > 1e-18) & (t > eps) & (t < 1 - eps)
                & (u > eps) & (u < 1 - eps)).any():
            return True
    return False


def _verified(p0, d0, segments, p1, e1, scale) -> bool:
    """True if *segments* run from (p0, d0) to (p1, e1) within tolerance."""
    pe, de = _simulate_segments(p0, d0, segments)
    return (np.linalg.norm(pe - p1) <= 1e-6 * max(scale, 1.0)
            and np.linalg.norm(de - e1) <= 1e-9)


def _segments_length(segments) -> float:
    """Centerline length of arc/straight *segments* (arcs as ``|angle| * radius``)."""
    return sum(abs(np.deg2rad(s[1])) * s[2] if s[0] == "arc" else s[1] for s in segments)


def _is_facing_collinear(p0, d0, p1, d1) -> bool:
    """True when the poses face each other along one line, or coincide.

    The exact route is then a plain straight, which the arc solvers cannot
    express; they return ``None`` and `_no_solution` draws it instead.
    """
    D = p1 - p0
    dist = float(np.linalg.norm(D))
    if dist < 1e-9:
        return True
    return bool(np.dot(d0 + d1, d0 + d1) < 1e-6 and 1.0 - np.dot(D / dist, d0) < 1e-6)


def _solve_circular(
    p0: np.ndarray,
    d0: np.ndarray,
    p1: np.ndarray,
    d1: np.ndarray,
    min_radius: float = 0.0,
) -> "list[_Arc] | None":
    """Solve for the two equal-radius arc path between two directed points.

    The path starts at *p0* traveling along *d0* and ends at *p1* traveling
    along ``-d1`` (module-level outward-normal convention).  A smooth
    equal-radius biarc joins two externally tangent circles, so the two arcs
    turn in opposite senses (an S); a same-sense pair only occurs as the
    degenerate single-circle case (second arc 0 deg).  For each first-arc
    sense the radius follows from ``|c0 - c1|^2 = 4 r^2``, a quadratic in r.
    Every candidate is verified by forward simulation before being
    considered; the shortest verified path wins.

    Returns a list of two ``_Arc`` specs, or ``None`` when the geometry
    degenerates to a straight line (collinear facing ports) or no solution
    exists above *min_radius* (e.g. back-to-back collinear ports).
    """
    e1 = -np.asarray(d1, float)  # arrival travel direction
    dist = float(np.linalg.norm(p1 - p0))
    if _is_facing_collinear(p0, d0, p1, d1):  # includes coincident (dist < 1e-9)
        return None

    candidates = []
    dp = p0 - p1
    for s0 in (1.0, -1.0):
        # Opposite-sense biarc: |c0 - c1| = 2r with c0/c1 offset by the senses.
        s1 = -s0
        u = s0 * _perp_ccw(d0) - s1 * _perp_ccw(e1)
        a_c = float(u @ u) - 4.0
        b_c = 2.0 * float(u @ dp)
        c_c = float(dp @ dp)
        roots: "list[float]" = []
        if abs(a_c) < 1e-12:
            if abs(b_c) > 1e-12:
                roots = [-c_c / b_c]
        else:
            disc = b_c * b_c - 4.0 * a_c * c_c
            if disc >= 0.0:
                sq = np.sqrt(disc)
                roots = [(-b_c + sq) / (2 * a_c), (-b_c - sq) / (2 * a_c)]
        for r in roots:
            if r <= 1e-9 or (min_radius > 0.0 and r < min_radius):
                continue
            c0 = p0 + s0 * r * _perp_ccw(d0)
            c1 = p1 + s1 * r * _perp_ccw(e1)
            m = 0.5 * (c0 + c1)
            a1 = _directed_angle_deg(p0 - c0, m - c0, s0 > 0)
            a2 = _directed_angle_deg(m - c1, p1 - c1, s1 > 0)
            candidates.append([("arc", a1, r), ("arc", a2, r)])

        # Same-sense degenerate: both ports on one circle -> a single arc.
        u1 = s0 * (_perp_ccw(d0) - _perp_ccw(e1))
        uu = float(u1 @ u1)
        if uu > 1e-12:
            r = -float(u1 @ dp) / uu
            if (r > 1e-9 and np.linalg.norm(dp + r * u1) < 1e-6 * max(dist, 1.0)
                    and not (min_radius > 0.0 and r < min_radius)):
                c0 = p0 + s0 * r * _perp_ccw(d0)
                a1 = _directed_angle_deg(p0 - c0, p1 - c0, s0 > 0)
                candidates.append([("arc", a1, r), ("arc", 0.0, r)])

    best, best_len = None, np.inf
    for segs in candidates:
        if not _verified(p0, d0, segs, p1, e1, dist):
            continue
        length = _segments_length(segs)
        if length < best_len:
            best_len = length
            best = [_Arc(a, r) for _, a, r in segs]
    return best


def _solve_ata(
    p0: np.ndarray,
    d0: np.ndarray,
    p1: np.ndarray,
    d1: np.ndarray,
    radius: float,
) -> "tuple[_Arc, float, _Arc] | None":
    """Solve for the shortest arc-straight-arc (Dubins CSC) path at a fixed radius.

    The path starts at *p0* traveling along *d0* and ends at *p1* traveling
    along ``-d1`` (module-level outward-normal convention).  All four
    turn-sense combinations are evaluated: same-sense pairs (LSL/RSR) connect
    along the outer tangent (length ``|c1 - c0|``); opposite-sense pairs
    (LSR/RSL) along the inner tangent (length ``sqrt(|c1 - c0|^2 - 4 r^2)``,
    which requires ``|c1 - c0| >= 2 r``).  Every candidate is verified by
    forward simulation before being considered; the shortest verified path is
    returned as ``(arc1, straight_length, arc2)``, or ``None`` if nothing
    verifies (degenerate inputs such as coincident ports).
    """
    r = float(radius)
    e1 = -np.asarray(d1, float)  # arrival travel direction
    scale = max(float(np.linalg.norm(p1 - p0)), r)

    candidates = []
    for s0 in (1.0, -1.0):
        for s1 in (1.0, -1.0):
            c0 = p0 + s0 * r * _perp_ccw(d0)
            c1 = p1 + s1 * r * _perp_ccw(e1)
            L = c1 - c0
            L_len = float(np.linalg.norm(L))
            if L_len < 1e-9:
                if s0 == s1:  # same circle -> single arc, no straight
                    a1 = _directed_angle_deg(p0 - c0, p1 - c0, s0 > 0)
                    candidates.append([("arc", a1, r), ("straight", 0.0), ("arc", 0.0, r)])
                continue
            L_unit = L / L_len
            if s0 == s1:
                # Outer tangent: straight runs parallel to the center line.
                tangent_dirs = [L_unit]
            else:
                # Inner tangent: crosses between the circles.
                if L_len < 2.0 * r - 1e-9:
                    continue
                phi = float(np.arcsin(np.clip(2.0 * r / L_len, -1.0, 1.0)))
                tangent_dirs = [_rotate2d(L_unit, phi), _rotate2d(L_unit, -phi)]
            for t in tangent_dirs:
                tp0 = c0 - s0 * r * _perp_ccw(t)   # leave circle 0 traveling t
                tp1 = c1 - s1 * r * _perp_ccw(t)   # enter circle 1 traveling t
                seg = tp1 - tp0
                seg_len = float(np.linalg.norm(seg))
                if seg_len > 1e-9:
                    if float(seg @ t) < 0:
                        continue  # tangent points ordered against travel
                    if abs(seg[0] * t[1] - seg[1] * t[0]) > 1e-6 * seg_len:
                        continue  # wrong inner-tangent branch
                a1 = _directed_angle_deg(p0 - c0, tp0 - c0, s0 > 0)
                a2 = _directed_angle_deg(tp1 - c1, p1 - c1, s1 > 0)
                candidates.append([("arc", a1, r), ("straight", seg_len), ("arc", a2, r)])

    best, best_len = None, np.inf
    for segs in candidates:
        if not _verified(p0, d0, segs, p1, e1, scale):
            continue
        length = _segments_length(segs)
        if length < best_len:
            best_len = length
            best = (_Arc(segs[0][1], r), segs[1][1], _Arc(segs[2][1], r))
    return best


def _solve_aaa(
    p0: np.ndarray,
    d0: np.ndarray,
    p1: np.ndarray,
    d1: np.ndarray,
    radius: float,
) -> "list[_Arc] | None":
    """Solve for a three-arc (Dubins CCC) path at a fixed radius.

    The path starts at *p0* traveling along *d0* and ends at *p1* traveling
    along ``-d1`` (module-level outward-normal convention).  Both outer arcs
    are tangent to the ports and turn with the same sense; the middle arc
    opposes it, its center at distance ``2r`` from both tangent-circle
    centers -- which requires ``|c1 - c0| <= 4r``.  Candidates are tried outer
    sense first (the one whose tangent circles sit farthest apart), then
    shortest; a self-intersecting candidate is taken only when *every*
    verified candidate of either sense crosses itself.

    Returns three ``_Arc`` specs, or ``None`` when the geometry degenerates
    to a straight line (collinear facing ports) or no middle circle exists
    for either sense (ports too far apart at this radius -- use
    `_solve_ata`).
    """
    r = float(radius)
    e1 = -np.asarray(d1, float)  # arrival travel direction
    dist = float(np.linalg.norm(p1 - p0))
    if _is_facing_collinear(p0, d0, p1, d1):  # includes coincident (dist < 1e-9)
        return None
    scale = max(dist, r)

    # Outer senses ordered by the pushed-apart policy: farthest tangent circles first.
    options = []
    for s in (1.0, -1.0):
        c0 = p0 + s * r * _perp_ccw(d0)
        c1 = p1 + s * r * _perp_ccw(e1)
        options.append((float(np.linalg.norm(c1 - c0)), s, c0, c1))
    options.sort(reverse=True, key=lambda o: o[0])

    fallback = None
    for dist_c, s, c0, c1 in options:
        candidates = []
        if dist_c < 1e-9:
            # Ports share one tangent circle: a single arc.
            a0 = _directed_angle_deg(p0 - c0, p1 - c0, s > 0)
            candidates.append([("arc", a0, r), ("arc", 0.0, r), ("arc", 0.0, r)])
        elif dist_c <= 4.0 * r:
            u = (c1 - c0) / dist_c
            h = float(np.sqrt(max(4.0 * r * r - 0.25 * dist_c * dist_c, 0.0)))
            for side in (1.0, -1.0):
                cm = c0 + 0.5 * dist_c * u + side * h * _perp_ccw(u)
                t0 = c0 + 0.5 * (cm - c0)   # tangency points: |cm - c0,1| = 2r
                t1 = c1 + 0.5 * (cm - c1)
                a0 = _directed_angle_deg(p0 - c0, t0 - c0, s > 0)
                am = _directed_angle_deg(t0 - cm, t1 - cm, s < 0)  # opposite sense
                a1 = _directed_angle_deg(t1 - c1, p1 - c1, s > 0)
                candidates.append([("arc", a0, r), ("arc", am, r), ("arc", a1, r)])
                if h < 1e-9:
                    break  # circles exactly 4r apart: a single middle circle

        verified = sorted(
            (segs for segs in candidates if _verified(p0, d0, segs, p1, e1, scale)),
            key=_segments_length)
        for segs in verified:
            if not _self_intersects(_sample_segments(p0, d0, segs)[0]):
                return [_Arc(a, rr) for _, a, rr in segs]
        if fallback is None and verified:
            fallback = [_Arc(a, rr) for _, a, rr in verified[0]]
    # Every verified candidate self-crosses: return the preferred one anyway.
    return fallback


# Manhattan solver (pure geometry)

# How far off the 90 deg grid an orientation may sit and still snap to a cardinal.
_MANHATTAN_ANGLE_TOL = 1e-3
# Lateral offsets below this (one database unit) are treated as exactly on-axis (um).
_MANHATTAN_POS_TOL = 1e-3


def _cardinal(orientation: float) -> "int | None":
    """Index of the nearest cardinal direction (0=E, 1=N, 2=W, 3=S), or
    ``None`` when *orientation* is more than `_MANHATTAN_ANGLE_TOL` degrees
    off the 90 deg grid."""
    k = int(np.rint(orientation / 90.0))
    if abs(orientation - 90.0 * k) > _MANHATTAN_ANGLE_TOL:
        return None
    return k % 4


def _rot90_xy(x: float, y: float, k: int) -> "tuple[float, float]":
    """Rotate (x, y) by ``k * 90`` deg CCW (exact, no trigonometry)."""
    for _ in range(k % 4):
        x, y = -y, x
    return x, y


def _solve_manhattan_waypoints(p0, k0, p1, k1, radius, lead) -> np.ndarray:
    """Rectilinear waypoints from pose (*p0*, cardinal *k0*) to (*p1*, *k1*).

    Returns an (n, 2) array of waypoints (*p0* first, *p1* last) whose legs
    are axis-aligned, leaving *p0* along *k0* and arriving at *p1* traveling
    against *k1* (module-level outward-normal convention).

    Shape policy, which keeps a fan of pairs with varying pitch crossing-free:

    * the first turn heads into the half-plane containing port2 (the
      "quadrant rule"); ties break to port1's left (+90 deg);
    * jogs sit a fixed *lead* from port1, so parallel routes all turn at the
      same station and may share a jog channel;
    * the leg arriving at port2 is kept as long as the pair allows.

    *radius* only sets the corner clearances (``2*radius`` between adjacent
    corners, ``radius`` at the ports); the caller fits the actual fillet and
    shrinks it for the cramped shapes, which violate the clearances rather
    than fail.
    """
    r = float(radius)
    L = max(float(lead), r)
    tol = _MANHATTAN_POS_TOL

    # Canonical frame: port1 at the origin traveling +x, port2 in the upper half-plane.
    X, Y = _rot90_xy(float(p1[0] - p0[0]), float(p1[1] - p0[1]), -k0 % 4)
    a = (k1 + 2 - k0) % 4  # arrival travel direction (0=E, 1=N, 2=W, 3=S)
    my = 1.0
    if Y < -tol:
        my, Y, a = -1.0, -Y, (-a) % 4
    if Y <= tol:
        Y = 0.0

    pts = [(0.0, 0.0)]
    if a == 0:  # arrive traveling +x, head-on into a west-facing port
        if Y == 0.0 and X > tol:
            pass  # facing collinear: pure straight
        elif X >= 2.0 * r:
            # S-jog a lead from port1, kept in the corridor so parallel routes never cross.
            j = max(r, min(L, X - r))
            pts += [(j, 0.0), (j, Y)]
        else:
            # Port2 beside/behind: wrap over the top and come back east.
            s = max(L, r)
            back = X - s
            first = max(L, back + 2.0 * r)
            ym = Y + 2.0 * r
            pts += [(first, 0.0), (first, ym), (back, ym), (back, Y)]
    elif a == 2:  # arrive traveling -x, head-on into an east-facing port
        if Y >= 2.0 * r:
            # C-shape around port2's east side.
            xt = max(L, X + L)
            pts += [(xt, 0.0), (xt, Y)]
        else:
            s = max(L, r)
            ym = Y + 2.0 * r
            if X > tol:
                # Port2 dead ahead: stop short, dodge over it, come back west into its east face.
                first = min(L, X - 2.0 * r)
                if first < r:
                    first = X / 3.0  # cramped: fillet fit shrinks the radius
                pts += [(first, 0.0), (first, ym), (X + s, ym), (X + s, Y)]
            else:
                # Beside/behind: go east past the port, over the top, drop.
                if Y == 0.0:
                    # Port2 on port1's axis: the west leg stops short of the east leg.
                    s = min(max(L, r), -X - r)
                    if s < r:
                        s = max(-X / 2.0, tol)  # cramped: fillet shrinks
                first = max(L, X + s + 2.0 * r)
                pts += [(first, 0.0), (first, ym), (X + s, ym), (X + s, Y)]
    elif a == 1:  # arrive traveling +y, head-on into a south-facing port
        if X >= r and Y >= r:
            pts += [(X, 0.0)]  # the natural single-corner L
        elif Y >= 3.0 * r:
            # Behind with headroom: hook up and back, staying in port2's half-plane.
            first = max(L, X + 2.0 * r)
            pts += [(first, 0.0), (first, 2.0 * r), (X, 2.0 * r)]
        elif X < -tol or X - 2.0 * r >= r:
            # Full-radius dip below the axis: port2 behind, or the approach column clears port1.
            yc = min(-2.0 * r, Y - L)
            first = (min(L, X - 2.0 * r) if X - 2.0 * r >= r
                     else max(L, X + 2.0 * r))
            pts += [(first, 0.0), (first, yc), (X, yc)]
        elif Y >= X:
            # Cramped hook: the riser shrinks with Y (fillet fit follows).
            yc = 2.0 * Y / 3.0
            first = max(L, X + 2.0 * r)
            pts += [(first, 0.0), (first, yc), (X, yc)]
        else:
            # Cramped dip tucked in front of the column: fillet shrinks.
            yc = min(-2.0 * r, Y - L)
            first = min(L, max(X / 3.0, X - 2.0 * r))
            pts += [(first, 0.0), (first, yc), (X, yc)]
    else:  # a == 3: arrive traveling -y, head-on into a north-facing port
        ym = max(Y + L, 2.0 * r)
        if X - 2.0 * r >= r:
            first = min(L, X - 2.0 * r)
        elif X > tol and Y < r:
            first = X / 3.0  # cramped: fillet fit shrinks the radius
        else:
            first = max(L, X + 2.0 * r)
        pts += [(first, 0.0), (first, ym), (X, ym)]
    pts.append((X, Y))

    return np.array([
        np.asarray(p0, float) + _rot90_xy(x, my * y, k0) for x, y in pts])


def _waypoints_to_segments(points, radius) -> list:
    """Fillet the rectilinear polyline *points* into ``("straight", L)`` /
    ``("arc", +-90, radius)`` segments.

    Each corner fillet consumes *radius* of tangent length on both adjacent
    legs (90 deg corners: ``tan(45 deg) = 1``); the manhattan waypoint
    builder plus the caller's fillet fit guarantee that this fits.
    """
    pts = np.asarray(points, float)
    vecs = np.diff(pts, axis=0)
    lens = np.linalg.norm(vecs, axis=1)
    dirs = vecs / lens[:, None]
    turns = [90.0 * np.sign(dirs[i][0] * dirs[i + 1][1]
                            - dirs[i][1] * dirs[i + 1][0])
             for i in range(len(dirs) - 1)]
    segments: list = []
    for i, length in enumerate(lens):
        cut = ((radius if i > 0 else 0.0)
               + (radius if i < len(turns) else 0.0))
        if length - cut > 1e-9:
            segments.append(("straight", float(length - cut)))
        if i < len(turns):
            segments.append(("arc", float(turns[i]), float(radius)))
    return segments


# gdsfactory drawing helpers

def _warn_if_route_bad(component, port2, *, caller):
    """Warn when the placed route misses *port2*, or when it crosses itself.

    Guards against bend factories whose endpoints deviate from the solved
    arc floorplan (e.g. ``bend_default(..., with_arc_floorplan=False)``) and
    against solver regressions.  A self-crossing ribbon shows up as polygon holes.
    """
    if any(p.holes() > 0 for polys in component.get_polygons().values() for p in polys):
        warnings.warn(f"{caller}: the drawn waveguide crosses itself (optical short)",
                      stacklevel=3)
    err = float(np.linalg.norm(np.array(component.ports["o2"].center) - np.array(port2.center)))
    if err > 0.001:
        logger.debug("{}: route miss distance={:.6f} tol={:.6f}", caller, err, MISS_TOLERANCE)
    if err > MISS_TOLERANCE:
        warnings.warn(
            f"{caller}: route ends {err:.4f} um away from port2 "
            f"{tuple(np.round(port2.center, 4))}; the bend factory may not "
            "match the solved arc floorplan",
            stacklevel=3,
        )


def _decide_widths(width1, width2, routing_width, *, caller, stacklevel=3):
    """Shared width/taper decision used by the arc and spline routers.

    Returns ``(tapering, routing_width)``.  When the end widths differ the
    route is tapered along its whole length: *tapering* is ``True``, a warning
    is emitted, and *routing_width* is dropped (a port-to-port taper and a
    wide-straight section cannot be combined).
    """
    tapering = not np.isclose(width1, width2)
    if tapering:
        warnings.warn(
            f"{caller}: port widths differ ({width1:g} != {width2:g}); "
            "tapering the width along the path (bend is ignored -- the path is extruded "
            "as one transition)"
            + ("; routing_width dropped" if routing_width is not None else ""),
            stacklevel=stacklevel,
        )
        routing_width = None
    return tapering, routing_width


def _arcs_to_segments(arcs: "list[_Arc]", straight_len: float) -> list:
    """The arc-solver result as a generic segment list.

    Arc solvers return ``_Arc`` specs plus at most one straight, which by
    their shared convention sits after the first arc; with no arcs the
    straight is the whole route.
    """
    segments: list = []
    for i, arc in enumerate(arcs):
        segments.append(("arc", arc.angle_deg, arc.radius))
        if i == 0:
            segments.append(("straight", straight_len))
    if not arcs:
        segments.append(("straight", straight_len))
    return segments


def _segments_to_path(segments, *, euler: bool = True) -> gf.Path:
    """Concatenate arc/straight *segments* into one centerline Path.

    Euler arcs use ``use_eff=True`` so their effective radius matches the
    arc-floorplan bend the solvers assume; windings beyond 180 deg are split
    into equal sub-arcs (see :func:`layout.core.ports._split_angle`).
    ``euler=False`` uses circular arcs instead, matching ``bend_circular``.
    """
    segs: list = []
    for seg in segments:
        if seg[0] == "arc":
            _, angle_deg, radius = seg
            if abs(angle_deg) > _ZERO_ANGLE_DEG:
                n_sub, sub_angle = _split_angle(angle_deg)
                for _ in range(n_sub):
                    segs.append(gf.path.euler(radius=radius, angle=sub_angle,
                                              use_eff=True) if euler
                                else gf.path.arc(radius=radius, angle=sub_angle))
        elif seg[1] > 1e-9:
            segs.append(gf.path.straight(seg[1]))
    if not segs:
        total = sum(s[1] for s in segments if s[0] == "straight")
        return gf.path.straight(max(total, 1e-6))
    path = segs[0]
    for s in segs[1:]:
        path = path + s
    return path


def _path_at_port(path: gf.Path, port) -> gf.Path:
    """Anchor a local centerline (start at the origin, heading +x) at *port*."""
    path.rotate(port.orientation)
    path.move((0.0, 0.0), tuple(float(x) for x in port.center))
    return path


def _place_arcs(
    arcs: "list[_Arc]",
    straight_len: float,
    port1: Port,
    bend,
    **kwargs,
) -> gf.Component:
    """Realize a solved arc/straight centerline (see `_place_segments`).

    Arc-solver calling convention: the single straight is placed after the
    first arc.
    """
    return _place_segments(_arcs_to_segments(arcs, straight_len), port1,
                           bend, **kwargs)


def _round_segments(segments):
    """*segments* with solver float noise rounded off.

    A fan of congruent routes differs in the last bits of its solved angles and lengths,
    which would key a body per route.
    """
    return tuple(tuple(round(v, 9) if isinstance(v, float) else v for v in seg)
                 for seg in segments)


@gf.cell(basename="route")
def _route_body(
    segments,
    orientation: float,
    bend,
    *,
    width1: float,
    width2: float,
    routing_width: "float | None" = None,
    taper_length: "float | None" = 100.0,
    cross_section: "CrossSectionSpec | None" = None,
    layer=None,
) -> gf.Component:
    """The drawn route, leaving the origin along *orientation*, as a shared cell.

    Two width modes:

    * **uniform** (``width1 == width2``) -- place each arc from the *bend* factory;
      every straight is widened to *routing_width* (up-taper / flat / down-taper)
      when it is longer than ``2*taper_length``.
    * **tapering** -- extrude the whole centerline as one ``gf.path.transition``
      (*routing_width* does not apply).

    Keyed on the shape alone, so routes congruent up to a shift share one cell.  Exposes
    ``o1``/``o2`` and the centerline length in ``info["length"]``.
    """
    # Layer 0 placeholder: the connects below are position-only (see allow).
    start = gf.Port(name="o1", center=(0.0, 0.0), width=width1,
                    orientation=orientation, layer=0)
    resolved = gf.get_cross_section(cross_section) if cross_section is not None else None

    def make_cs(w):
        return resolved.copy(width=w) if resolved is not None else gf.cross_section.cross_section(width=w, layer=layer)

    # Position-only connects: route pieces carry optical ports even on RF routes.
    allow = {"allow_type_mismatch": True, "allow_width_mismatch": True,
             "allow_layer_mismatch": True}

    c = gf.Component()

    if not np.isclose(width1, width2):
        # Circular arcs: what the solvers laid out, so the drawn radius is the requested one.
        path = _segments_to_path(segments, euler=False)
        transition = gf.path.transition(make_cs(width1), make_cs(width2), width_type="linear")
        ref = c.add_ref(gf.path.extrude_transition(path, transition))
        ref.connect("o1", start, **allow)
        end_port = ref.ports["o2"].copy()
        total = float(path.length())
    else:
        wb = width1
        cursor = start
        total = 0.0

        def chain(comp, off_grid=False):
            nonlocal cursor
            ref = c.add_ref_off_grid(comp) if off_grid else c.add_ref(comp)
            ref.connect("o1", cursor, **allow)
            cursor = ref.ports["o2"]
            return ref

        def place_bend(angle_deg, radius):
            nonlocal total
            if abs(angle_deg) < _ZERO_ANGLE_DEG:
                return
            n_sub, sub_angle = _split_angle(angle_deg)
            comp = gf.get_component(bend, angle=sub_angle, radius=radius,
                                    cross_section=make_cs(wb))
            for _ in range(n_sub):
                chain(comp, off_grid=True)
            total += n_sub * _component_length(
                comp, abs(np.deg2rad(sub_angle)) * radius)

        def taper(w_a, w_b, length):
            return gf.path.extrude_transition(
                gf.path.straight(length), gf.path.transition(make_cs(w_a), make_cs(w_b), width_type="linear"))

        def straight(length, w):
            return gf.path.extrude(gf.path.straight(length), cross_section=make_cs(w))

        def place_straight(length):
            nonlocal total
            if length <= 1e-9:
                return
            if (routing_width is not None and taper_length is not None
                    and not np.isclose(routing_width, wb) and length > 2 * taper_length):
                chain(taper(wb, routing_width, taper_length))
                chain(straight(length - 2 * taper_length, routing_width))
                chain(taper(routing_width, wb, taper_length))
            else:
                chain(straight(length, wb))
            total += length

        for seg in segments:
            if seg[0] == "arc":
                place_bend(seg[1], seg[2])
            else:
                place_straight(seg[1])

        end_port = cursor.copy()

    # Flatten realizes the off-grid virtual instances all-angle bends place.
    c.flatten()
    c.add_ports([start.copy(name="o1"), end_port.copy(name="o2")])
    c.info["length"] = float(total)
    return c


def _place_segments(
    segments,
    port1: Port,
    bend,
    *,
    width1: float,
    width2: float,
    routing_width: "float | None" = None,
    taper_length: "float | None" = 100.0,
    cross_section: "CrossSectionSpec | None" = None,
    layer=None,
    caller: str,
) -> gf.Component:
    """Place the `_route_body` for a solved centerline at *port1*.

    Exposes ``o1`` (= *port1*) and ``o2`` (route end); the centerline length goes to
    ``info["length"]``.
    """
    _, routing_width = _decide_widths(
        width1, width2, routing_width, caller=caller,
        stacklevel=5)  # helper -> _place_segments -> _place_arcs -> router -> user

    body = _route_body(
        _round_segments(segments), port1.orientation, bend,
        width1=width1, width2=width2, routing_width=routing_width,
        taper_length=taper_length,
        layer=None if layer is None else int(gf.get_layer(layer)),
        cross_section=None if cross_section is None else gf.get_cross_section(cross_section),
    )

    c = gf.Component()
    ref = c.add_ref(body).dmove(port1.center)
    length = float(body.info["length"])
    # o2 continues port1's chain: info plus the accumulated pathlength.
    end_port = ref.ports["o2"].copy()
    _copy_info(port1.info, end_port.info)
    end_port.info["pathlength"] = get_pathlength(port1) + length
    c.add_ports([port1.copy(name="o1"), end_port.copy(name="o2")])
    c.info["length"] = length
    c.info["shared_body"] = True
    return c


def _no_solution(caller, port1, port2, detail):
    """Straight length for facing-collinear ports; raise otherwise.

    The solvers return ``None`` both for the degenerate straight-line case
    and for genuinely unsolvable pairs; only the former may be drawn.
    """
    p0, d0 = np.array(port1.center, float), get_direction(port1)
    p1, d1 = np.array(port2.center, float), get_direction(port2)
    if float(np.linalg.norm(p1 - p0)) < 1e-9:
        raise ValueError(f"{caller}: port1 and port2 are coincident")
    if not _is_facing_collinear(p0, d0, p1, d1):
        hint = "" if caller == "connect_ata" else "; try connect(..., method='ata')"
        raise ValueError(
            f"{caller}: no {detail} for these ports "
            f"(port1 at {tuple(np.round(p0, 3))} facing {port1.orientation:g} deg, "
            f"port2 at {tuple(np.round(p1, 3))} facing {port2.orientation:g} deg)"
            f"{hint}"
        )
    return float(np.linalg.norm(p1 - p0))


# Public path-planner API

def connect_circular(
    port1: Port,
    port2: Port,
    *,
    cross_section: "CrossSectionSpec | None" = None,
    expected_radius: "float | None" = None,
    min_radius: float = 0.0,
    bend=bend_default,
) -> gf.Component:
    """Route *port1* to *port2* using two equal-radius arcs of *bend*.

    `_solve_circular` computes the shared radius so the path departs *port1*
    and arrives at *port2* head-on (module-level outward-normal convention) --
    the arc radius is therefore *set by the port geometry*, not chosen;
    collinear opposing ports get a straight.  Unequal port widths are tapered
    along the arcs (warns) -- for a wide routing-width straight use
    :func:`connect_ata` or :func:`connect_intersection` instead.

    Parameters
    ----------
    port1, port2 : Port
        Endpoints.
    cross_section : CrossSectionSpec, optional
        Cross-section for the arcs; ``None`` inherits width/layer from *port1*.
    expected_radius : float, optional
        Bend radius to check against: warns when the solved radius comes out
        smaller (the route is still drawn -- pass *min_radius* to reject
        instead).
    min_radius : float
        Reject arc solutions below this radius.
    bend : callable
        Bend factory called as ``bend(angle=, radius=, ...)``.

    Returns
    -------
    gf.Component
        Route anchored at *port1*, with ``o1``/``o2`` ports and
        ``info["length"]``; ``o2`` continues *port1*'s ``pathlength`` chain.
    """
    arcs = _solve_circular(np.array(port1.center, float), get_direction(port1),
                           np.array(port2.center, float), get_direction(port2),
                           min_radius)
    straight_len = 0.0 if arcs else _no_solution(
        "connect_circular", port1, port2, "two-equal-arc solution")

    if expected_radius is not None and arcs:
        solved = min(a.radius for a in arcs)
        if solved < expected_radius * (1.0 - 1e-9):
            warnings.warn(
                f"connect_circular: solved radius {solved:.3f} um is smaller "
                f"than the expected {expected_radius:g} um (the biarc radius is "
                "set by the port geometry; pass min_radius= to reject instead)",
                stacklevel=2,
            )

    c = _place_arcs(
        arcs or [], straight_len, port1, bend,
        width1=port1.width, width2=port2.width,
        cross_section=cross_section, layer=port1.layer, caller="connect_circular",
    )
    _warn_if_route_bad(c, port2, caller="connect_circular")
    return c


def _ata_prefers_aaa(p0, d0, p1, d1, radius, segments) -> bool:
    """True when the arc-tangent-arc route should defer to three arcs.

    Close ports (``< 2*radius``) whose path sweeps behind both are a
    wrap-around; the three-arc route avoids it -- unless that replacement
    self-crosses and this path does not (never trade for a worse route).
    """
    if float(np.linalg.norm(p1 - p0)) >= 2.0 * radius:
        return False
    pts, _ = _sample_segments(p0, d0, segments)
    if not (_crosses_ray(pts, p0, -d0) and _crosses_ray(pts, p1, -d1)):
        return False
    aaa_arcs = _solve_aaa(p0, d0, p1, d1, radius) or []
    if not aaa_arcs:
        return False
    aaa_pts, _ = _sample_segments(
        p0, d0, [("arc", a.angle_deg, a.radius) for a in aaa_arcs])
    return (not _self_intersects(aaa_pts)) or _self_intersects(pts)


def connect_ata(
    port1: Port,
    port2: Port,
    *,
    radius: float,
    cross_section: "CrossSectionSpec | None" = None,
    bend=bend_default,
    routing_width: "float | None" = None,
    taper_length: "float | None" = 100.0,
) -> gf.Component:
    """Route *port1* to *port2* via arc-tangent-arc at a fixed bend *radius*.

    `_solve_ata` evaluates all four left/right combinations and the shortest
    verified path is realized (module-level outward-normal convention).  With
    equal port widths, *routing_width* widens the middle straight
    (taper / flat / taper); unequal port widths are tapered along the path
    instead (warns).

    Ports within ``2*radius`` whose path sweeps behind both ports fall back
    to :func:`connect_aaa` (warns), unless that would trade a clean path for
    a self-crossing one.

    Parameters
    ----------
    port1, port2 : Port
        Endpoints.
    radius : float
        Bend radius for both arcs.
    cross_section : CrossSectionSpec, optional
        Cross-section for the arcs/straight; ``None`` inherits from *port1*.
    bend : callable
        Bend factory called as ``bend(angle=, radius=, ...)``.
    routing_width : float, optional
        Width of the widened middle straight (requires equal port widths).
    taper_length : float
        Length of each width taper.

    Returns
    -------
    gf.Component
        Route anchored at *port1*, with ``o1``/``o2`` ports and
        ``info["length"]``; ``o2`` continues *port1*'s ``pathlength`` chain.
    """
    p0, d0 = np.array(port1.center, float), get_direction(port1)
    p1, d1 = np.array(port2.center, float), get_direction(port2)
    result = _solve_ata(p0, d0, p1, d1, radius)
    if result is None:
        arcs, straight_len = [], _no_solution(
            "connect_ata", port1, port2,
            f"arc-straight-arc solution at radius {radius:g}")
    else:
        arc1, straight_len, arc2 = result
        arcs = [arc1, arc2]
        if _ata_prefers_aaa(p0, d0, p1, d1, radius, _arcs_to_segments(arcs, straight_len)):
            warnings.warn(
                "connect_ata: falling back to connect_aaa -- the "
                "arc-tangent-arc path sweeps behind both ports "
                f"(ports {np.linalg.norm(p1 - p0):.1f} um apart < 2*radius = "
                f"{2 * radius:g})"
                + ("; routing_width ignored" if routing_width is not None else ""),
                stacklevel=2)
            return connect_aaa(port1, port2, radius=radius,
                               cross_section=cross_section, bend=bend)

    c = _place_arcs(
        arcs, straight_len, port1, bend,
        width1=port1.width, width2=port2.width,
        routing_width=routing_width, taper_length=taper_length,
        cross_section=cross_section, layer=port1.layer, caller="connect_ata",
    )
    _warn_if_route_bad(c, port2, caller="connect_ata")
    return c


def connect_aaa(
    port1: Port,
    port2: Port,
    *,
    radius: float,
    cross_section: "CrossSectionSpec | None" = None,
    bend=bend_default,
) -> gf.Component:
    """Route *port1* to *port2* via three arcs (Dubins CCC) at a fixed *radius*.

    Two port-tangent outer arcs joined by an opposite-sense middle arc (see
    `_solve_aaa` for the sense policy) -- the loop-free alternative to
    `connect_ata` where its big-winding solutions self-intersect.  Feasible
    while the tangent circles are within ``4*radius``; collinear opposing
    ports get a straight; unequal port widths taper along the arcs (warns).

    Parameters
    ----------
    port1, port2 : Port
        Endpoints.
    radius : float
        Radius of all three arcs.
    cross_section : CrossSectionSpec, optional
        Cross-section for the arcs; ``None`` inherits width/layer from *port1*.
    bend : callable
        Bend factory called as ``bend(angle=, radius=, ...)``.

    Returns
    -------
    gf.Component
        Route anchored at *port1*, with ``o1``/``o2`` ports and
        ``info["length"]``; ``o2`` continues *port1*'s ``pathlength`` chain.
    """
    arcs = _solve_aaa(np.array(port1.center, float), get_direction(port1),
                      np.array(port2.center, float), get_direction(port2),
                      radius)
    straight_len = 0.0 if arcs else _no_solution(
        "connect_aaa", port1, port2,
        f"three-arc solution at radius {radius:g}")

    c = _place_arcs(
        arcs or [], straight_len, port1, bend,
        width1=port1.width, width2=port2.width,
        cross_section=cross_section, layer=port1.layer, caller="connect_aaa",
    )
    _warn_if_route_bad(c, port2, caller="connect_aaa")
    return c


def _manhattan_segments(port1, port2, radius, lead) -> list:
    """Filleted rectilinear centerline segments (see `connect_manhattan`)."""
    k0 = _cardinal(port1.orientation)
    k1 = _cardinal(port2.orientation)
    if k0 is None or k1 is None:
        raise ValueError(
            "connect_manhattan: manhattan routing requires port orientations "
            f"on the 90 deg grid (got {port1.orientation:g} and "
            f"{port2.orientation:g} deg)")
    p0 = np.array(port1.center, float)
    p1 = np.array(port2.center, float)
    if float(np.linalg.norm(p1 - p0)) < 1e-9:
        raise ValueError("connect_manhattan: port1 and port2 are coincident")

    r = float(radius)
    lead_resolved = r if lead is None else float(lead)
    wps = _solve_manhattan_waypoints(p0, k0, p1, k1, r, lead_resolved)

    r_fit = _max_fillet_radius(wps)
    r_used = r if r_fit >= r else _FILLET_MARGIN * r_fit
    w_max = max(port1.width, port2.width)
    if len(wps) > 2 and r_used < w_max:
        raise ValueError(
            "connect_manhattan: ports are too cramped for a manhattan route "
            f"(max feasible fillet radius {r_fit:.3f} um < waveguide width "
            f"{w_max:g} um)")
    if r_used < 0.5 * r:
        warnings.warn(
            f"connect_manhattan: fillet radius reduced from {r:g} to "
            f"{r_used:.3f} um to fit the corners",
            stacklevel=3,
        )

    segments = _waypoints_to_segments(wps, r_used)
    d0 = np.array(_rot90_xy(1.0, 0.0, k0), float)
    e1 = -np.array(_rot90_xy(1.0, 0.0, k1), float)
    # The last waypoint, not p1: a sub-grid offset was clamped (_MANHATTAN_POS_TOL).
    if not _verified(p0, d0, segments, wps[-1], e1, float(np.linalg.norm(p1 - p0))):
        raise ValueError(
            "connect_manhattan: internal error: the filleted rectilinear "
            "path does not land on port2")
    return segments


def connect_manhattan(
    port1: Port,
    port2: Port,
    *,
    radius: float,
    cross_section: "CrossSectionSpec | None" = None,
    bend=bend_default,
    routing_width: "float | None" = None,
    taper_length: "float | None" = 100.0,
    lead: "float | None" = None,
) -> gf.Component:
    """Route *port1* to *port2* with axis-aligned straights filleted by
    90 deg bends of *bend* at *radius*.

    Both port orientations must lie on the 90 deg grid (within
    `_MANHATTAN_ANGLE_TOL`); anything else raises ``ValueError`` -- this
    router is deliberately orientation-specific.  The rectilinear centerline
    comes from `_solve_manhattan_waypoints`, whose shape policy makes the
    route deterministic across sweeps: the first turn heads into the
    half-plane containing port2 (ties break to port1's left), jogs sit a
    fixed *lead* from port1, and port2 keeps the longest straight approach
    the pair allows -- port1 does the bending, and a fan of routes with
    varying pitch never crosses (parallel routes may share a jog channel).

    When the pair is too cramped for *radius* the fillet is shrunk to the
    largest radius the waypoints accept (a severe reduction warns; below the
    waveguide width, where no usable fillet exists, raises ``ValueError``).
    With equal port widths, *routing_width* widens every straight longer
    than ``2*taper_length``; unequal port widths are tapered along the path
    instead (warns).

    Parameters
    ----------
    port1, port2 : Port
        Endpoints (module-level outward-normal convention).
    radius : float
        Fillet bend radius.
    cross_section : CrossSectionSpec, optional
        Cross-section for the route; ``None`` inherits from *port1*.
    bend : callable
        Bend factory called as ``bend(angle=, radius=, ...)``.
    routing_width : float, optional
        Width of the widened straights (requires equal port widths).
    taper_length : float
        Length of each width taper.
    lead : float, optional
        Distance port1 travels before the jog it takes when the pair needs one
        (clamped to at least *radius*; default *radius*); a pair that meets in a
        single corner ignores it.  Sweeping it across a bundle staggers the jog
        channels of parallel routes.

    Returns
    -------
    gf.Component
        Route anchored at *port1*, with ``o1``/``o2`` ports and
        ``info["length"]``; ``o2`` continues *port1*'s ``pathlength`` chain.
    """
    segments = _manhattan_segments(port1, port2, radius, lead)

    c = _place_segments(
        segments, port1, bend,
        width1=port1.width, width2=port2.width,
        routing_width=routing_width, taper_length=taper_length,
        cross_section=cross_section, layer=port1.layer,
        caller="connect_manhattan",
    )
    _warn_if_route_bad(c, port2, caller="connect_manhattan")
    return c


def _taper_to_width(port, target_width, taper_length, component, cross_section=None):
    """Taper *port* to *target_width* over *taper_length*, returning the inner port.

    No-op when *target_width* is ``None`` or already equal; otherwise draws
    the taper into *component* via `extend` (with *cross_section*, so the
    entry tapers match the route body).
    """
    if target_width is None or np.isclose(target_width, port.width):
        return port
    return extend(port, distance=taper_length, new_width=target_width,
                  add_to_component=component, cross_section=cross_section)


def _taper_to_common_width(port1, port2, taper_length, component, cross_section=None):
    """Bring two unequal-width ports to the smaller width with a straight taper.

    The wider port is tapered down over *taper_length* so the route between
    the returned inner ports runs at one uniform width.  No-op when the
    widths match or *taper_length* is ``None`` (the caller then tapers along
    the curve instead).
    """
    if taper_length is None or np.isclose(port1.width, port2.width):
        return port1, port2
    separation = float(np.linalg.norm(np.asarray(port2.center) - np.asarray(port1.center)))
    if separation < taper_length:
        warnings.warn(
            f"connect: taper_length={taper_length:g} does not fit in the {separation:g} um port "
            "separation; the entry taper overshoots and the route may cross it",
            stacklevel=3)
    w = min(port1.width, port2.width)
    return (_taper_to_width(port1, w, taper_length, component, cross_section),
            _taper_to_width(port2, w, taper_length, component, cross_section))


def _centerline(method, port1, port2, *, radius, bend=bend_default,
                lead=None, smoothing=gf.path.smooth, min_radius=0.0) -> gf.Path:
    """Route centerline as a `gf.Path` (width/layer-free).

    Backs `connect(..., add_to_component=False)`, running the routers' own
    pure solvers and drawing their arcs with *bend*'s own profile, euler
    unless *bend* resolves to ``bend_circular_all_angle``.  Accepts the same
    router options `connect` forwards (``lead``, ``smoothing``,
    ``min_radius``).
    """
    if method == "linear":
        return gf.path.Path([tuple(float(x) for x in port1.center),
                             tuple(float(x) for x in port2.center)])
    if method == "spline":
        return _spline_path(port1, port2, radius=radius, smoothing=smoothing)

    euler = getattr(bend, "func", bend) is not gf.components.bend_circular_all_angle
    p0, d0 = np.array(port1.center, float), get_direction(port1)
    p1, d1 = np.array(port2.center, float), get_direction(port2)

    if method == "intersection":
        segments = _intersection_segments(port1, port2, radius)
        if segments is not None:
            return _path_at_port(_segments_to_path(segments, euler=euler), port1)
        method = "ata"  # same fallback the drawn route takes

    if method == "circular":
        arcs = _solve_circular(p0, d0, p1, d1, min_radius)
        segments = _arcs_to_segments(arcs or [], 0.0 if arcs else _no_solution(
            "connect_circular", port1, port2, "two-equal-arc solution"))
    elif method == "aaa":
        arcs = _solve_aaa(p0, d0, p1, d1, radius)
        segments = _arcs_to_segments(arcs or [], 0.0 if arcs else _no_solution(
            "connect_aaa", port1, port2,
            f"three-arc solution at radius {radius:g}"))
    elif method == "manhattan":
        segments = _manhattan_segments(port1, port2, radius, lead)
    else:  # ata
        result = _solve_ata(p0, d0, p1, d1, radius)
        if result is None:
            segments = _arcs_to_segments([], _no_solution(
                "connect_ata", port1, port2,
                f"arc-straight-arc solution at radius {radius:g}"))
        else:
            arc1, straight_len, arc2 = result
            segments = _arcs_to_segments([arc1, arc2], straight_len)
            if _ata_prefers_aaa(p0, d0, p1, d1, radius, segments):
                segments = _arcs_to_segments(
                    _solve_aaa(p0, d0, p1, d1, radius) or [], 0.0)
    return _path_at_port(_segments_to_path(segments, euler=euler), port1)


def connect(
    port1: Port,
    port2: Port,
    method: str = "intersection",
    bend=bend_default,
    bend_width: "float | None" = None,
    routing_width: "float | None" = None,
    taper_length: "float | None" = 100.0,
    radius: "float | None" = None,
    cross_section: "CrossSectionSpec | None" = None,
    add_to_component: "gf.Component | Literal[False] | None" = None,
    **kwargs,
) -> "gf.Component | gf.Path":
    """Unified router: connect two ports with up to three waveguide widths.

    A thin dispatcher over `connect_intersection`, `connect_ata`,
    `connect_manhattan`, `connect_circular` and `connect_spline` adding a
    uniform width policy: each port is tapered to *bend_width* over
    *taper_length* before routing, and *routing_width* widens long straights
    (``"intersection"``, ``"ata"`` and ``"manhattan"`` only -- the other
    methods have no straight and warn).  When
    *bend_width* is ``None`` and the port widths differ, they are matched
    with a straight taper down to the smaller width; ``taper_length=None``
    instead tapers the width along the curve (warns).

    Parameters
    ----------
    port1, port2 : Port
        Endpoints.
    method : str
        ``"intersection"`` | ``"ata"`` | ``"aaa"`` | ``"circular"`` |
        ``"spline"`` | ``"manhattan"`` | ``"linear"`` (a straight fill between the
        ports: a ``route_quad`` unless a ``cross_section`` is given, which is
        extruded instead).  ``"route0"`` is a deprecated alias for
        ``"intersection"``.
    bend : callable
        Bend factory forwarded to the chosen router (ignored, with a
        warning, for ``"spline"`` and ``"linear"``).  A named function or
        ``functools.partial``; a lambda cannot key the route cache.
    bend_width : float, optional
        Single-mode width carried through every bend.
    routing_width : float, optional
        Wide width on long straights.
    taper_length : float, optional
        Straight-taper length; ``None`` tapers along the curve instead.
    radius : float, optional
        Bend radius: fixed for ``"ata"``, ``"aaa"`` and ``"manhattan"``, the
        fillet radius for ``"intersection"``, a *maximum* for ``"spline"``,
        and only an expectation to warn against for ``"circular"``, whose
        radius is solved from the port geometry.  ``None`` (the default)
        means 100 for the radius-taking methods and no expectation for
        ``"circular"``.
    cross_section : CrossSectionSpec, optional
        Base cross-section for the body (every drawn method; its width is
        replaced by *bend_width* / the port widths).
    add_to_component : gf.Component or False, optional
        A ``gf.Component`` references the finished route into it; ``None``
        (default) returns a fresh component; ``False`` returns the route
        centerline as a ``gf.Path`` (width params ignored, with a warning).
    **kwargs
        Router-specific options forwarded to the selected router: ``lead``
        (manhattan), ``smoothing`` (spline), ``min_radius`` (circular).  An
        option the chosen router does not accept raises (``intersection``
        forwards extras to ``bend_towards``).

    Returns
    -------
    gf.Component OR gf.Path
        A component (anchored at *port1*, exposing ``o1``/``o2`` with the
        total length in ``info["length"]`` and ``o2`` continuing *port1*'s
        ``pathlength`` chain), or -- when ``add_to_component is False`` -- the
        bare centerline ``gf.Path``.
    """
    if method == "route0":  # deprecated alias
        method = "intersection"
    if method not in ("intersection", "ata", "aaa", "circular", "spline",
                      "manhattan", "linear"):
        raise ValueError(
            f"connect: unknown method {method!r}; expected "
            "intersection/ata/aaa/circular/spline/manhattan/linear")
    if add_to_component is False:
        # A gf.Path is a bare centerline: width/extrusion params don't apply.
        ignored = [n for n, v in (("cross_section", cross_section),
                                   ("bend_width", bend_width),
                                   ("routing_width", routing_width)) if v is not None]
        if ignored:
            warnings.warn(
                f"connect: {', '.join(ignored)} ignored for a path return "
                "(add_to_component=False yields a width-free centerline)",
                stacklevel=2)
        return _centerline(method, port1, port2,
                           radius=100.0 if radius is None else radius, bend=bend,
                           **kwargs)
    if taper_length is None and (bend_width is not None or routing_width is not None):
        raise ValueError(
            "connect: taper_length=None (taper along the curve) cannot be "
            "combined with bend_width or routing_width, which both need "
            "straight tapers"
        )
    if bend_width is not None and bend_width <= 0:
        raise ValueError(f"connect: bend_width must be positive, got {bend_width!r}")
    if routing_width is not None and routing_width <= 0:
        raise ValueError(
            f"connect: routing_width must be positive, got {routing_width!r}")

    # radius=None means "no expectation" for circular, and 100 elsewhere.
    circular_expectation = radius
    if radius is None:
        radius = 100.0

    c = gf.Component()
    pl: list = []

    if bend_width is not None:
        sep = float(np.linalg.norm(np.array(port2.center, float)
                                   - np.array(port1.center, float)))
        if sep < 2 * taper_length:
            warnings.warn(
                f"connect: the two bend-width port tapers ({taper_length} um "
                f"each) exceed the port separation ({sep:.1f} um) and may "
                "overlap",
                stacklevel=2,
            )
    p1 = _taper_to_width(port1, bend_width, taper_length, c, cross_section)
    p2 = _taper_to_width(port2, bend_width, taper_length, c, cross_section)
    if bend_width is None:
        p1, p2 = _taper_to_common_width(p1, p2, taper_length, c, cross_section)
    # Straight-taper lengths, read off the port pathlength accounting.
    pl.append(get_pathlength(p1) - get_pathlength(port1))
    pl.append(get_pathlength(p2) - get_pathlength(port2))

    if routing_width is not None and method in ("aaa", "circular", "spline", "linear"):
        warnings.warn(
            f"connect: routing_width has no straight to widen for "
            f"method={method!r}; ignored",
            stacklevel=2,
        )
    if method in ("spline", "linear") and bend is not bend_default:
        why = ("the spline is smoothed, not built from bend cells"
               if method == "spline" else "a linear route has no bend")
        warnings.warn(f"connect: bend is ignored for method={method!r} ({why})",
                      stacklevel=2)

    bend_cs = cross_section
    if cross_section is not None and bend_width is not None:
        # Resolve specs (e.g. "strip") to a CrossSection before overriding width.
        bend_cs = gf.get_cross_section(cross_section).copy(width=bend_width)

    if method == "intersection":
        # connect_intersection takes only straight tapers; default their length here.
        pl.append(connect_intersection(
            c, p1, p2, routing_width=routing_width,
            taper_length=100.0 if taper_length is None else taper_length,
            bend=partial(bend, radius=radius), cross_section=bend_cs, **kwargs))
    elif method == "linear":
        if kwargs:
            raise TypeError(
                f"connect: method='linear' takes no router options, got {sorted(kwargs)}")
        _draw_straight(c, p1, p2, bend_cs)
        pl.append(float(np.linalg.norm(np.array(p2.center, float)
                                       - np.array(p1.center, float))))
    else:
        if method == "ata":
            sub = connect_ata(p1, p2, radius=radius, cross_section=bend_cs,
                              bend=bend, routing_width=routing_width,
                              taper_length=taper_length, **kwargs)
        elif method == "manhattan":
            sub = connect_manhattan(
                p1, p2, radius=radius, cross_section=bend_cs, bend=bend,
                routing_width=routing_width, taper_length=taper_length, **kwargs)
        elif method == "aaa":
            sub = connect_aaa(p1, p2, radius=radius, cross_section=bend_cs,
                              bend=bend, **kwargs)
        elif method == "circular":
            sub = connect_circular(p1, p2, cross_section=bend_cs, bend=bend,
                                   expected_radius=circular_expectation, **kwargs)
        else:
            sub = connect_spline(p1, p2, radius=radius, cross_section=bend_cs, **kwargs)
        c.add_ref(sub)
        pl.append(float(sub.info["length"]))

    # A shared body must survive unflattened; every other route must not.
    if c.vinsts or not c.insts or not all(i.cell.info.get("shared_body") for i in c.insts):
        c.flatten()
    c.add_port("o1", port=port1)
    c.add_port("o2", port=port2)
    # o2 continues port1's chain: pathlength accumulated through the route.
    c.ports["o2"].info["pathlength"] = get_pathlength(port1) + float(sum(pl))
    c.info["length"] = float(sum(pl))
    if add_to_component is not None:
        add_to_component.add_ref(c)
    return c
