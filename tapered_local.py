"""
Fiber-facing couplers and the chip-scale assemblies built from them.

`tapered_fiber_linear` draws the inverse taper that meets a lensed/tapered fiber,
`tapered_fiber_3d` adds its double-etch conversion, and `tapered_fiber_coupler` wraps
that in the alignment furniture (fiber footprint, ticks, glue spots, oxide opening,
logos, labels). `edge_coupler` is the facet-side alternative, and `make_ec` runs
every port of a device out to the facets and caps each one with it. `tapered_fiber_test`,
`tapered_fiber_inout` and `tapered_fiber_device` compose them into whole chips.

Every composite takes its parts as ``*_spec`` factories (`taper_spec`,
`tapered_fiber_spec`, `edge_coupler_spec`, `device_spec`) called with no arguments, so a
caller customises with `functools.partial` — `tapered_fiber_780` and
`tapered_fiber_1550` are the two in-house presets.
"""
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import numpy as np
import gdsfactory as gf
from functools import partial

from layout.core.layers import Layers
from layout.core.logos import get_cal, get_ionq
from layout.core.text import text
from layout.core.routing import connect

@gf.cell
def tapered_fiber_linear(
    length=50,
    width1=.5,
    width2=.2,
    taper_length=None,
    smooth=False,
):
    """The part that couples with the fiber, tapered from *width1* to *width2*.

    *taper_length* defaults to the whole *length*; anything shorter puts the balance into
    a straight at *width1* behind ``o1``. Set *smooth* to use gdsfactory's sine
    cross-section taper, whose zero-slope endpoints join straight sections and
    the reshaped tether neck without a boundary kink.
    """
    c = gf.Component()

    if taper_length is None:
        taper_length = length

    if smooth:
        cross_section1 = gf.cross_section.cross_section(
            width=width1,
            layer=Layers.wg,
        )
        cross_section2 = gf.cross_section.cross_section(
            width=width2,
            layer=Layers.wg,
        )
        taper_component = gf.components.taper_cross_section_sine(
            cross_section1=cross_section1,
            cross_section2=cross_section2,
            length=taper_length,
        )
    else:
        taper_component = gf.components.taper(
            length=taper_length,
            width1=width1,
            width2=width2,
            layer=Layers.wg,
        )

    taper = c.add_ref(taper_component)
    p = taper.ports["o1"]
    q = taper.ports["o2"]

    if taper_length < length:
        dx = length - taper_length
        taper.move((dx, 0))
        p = taper.ports["o1"]
        q = taper.ports["o2"]
        p = p.extend(dx, add_to_component=c)

    c.add_ports(
        [p, q]
    )

    return c

@gf.cell
def tapered_fiber_3d(
    taper_spec=tapered_fiber_linear,
    add_doubleetch_strip=True,
    input_extension_length=20.0,
):
    """
    Part that couples to the fiber, plus the converter to the waveguide layer,
    in this case the 3D design.

    ``taper_spec`` is registered at the origin. ``input_extension_length``
    sets the straight section before the taper; use zero when the taper must
    begin directly at the preceding S-bend. With ``add_doubleetch_strip``
    enabled, the taper is covered by a strip overhanging both ends.
    """
    c = gf.Component()
    if input_extension_length < 0:
        raise ValueError("input_extension_length must be non-negative")

    taper = c.add_ref(
        taper_spec()
    )
    taper.move((-taper.ports["o1"].center[0], 0))

    p1 = taper.ports["o1"]
    if input_extension_length > 0:
        p1 = p1.extend(input_extension_length, add_to_component=c)
    p2 = taper.ports["o2"]

    x1 = p1.center[0]
    x2 = p2.center[0] + 50
    w = 5.

    if add_doubleetch_strip:
        c.add_polygon(
            [
                (x1, -w/2),
                (x1, w/2),
                (x2, w/2),
                (x2, -w/2),
            ],
            layer=Layers.wg_doubleetch,
        )

    c.add_ports(
        [p1, p2]
    )

    return c

@gf.cell
def tapered_fiber_coupler(
    tapered_fiber_spec=tapered_fiber_3d,
    tapered_fiber_angle=1.5,
    tapered_fiber_length=5e3,
    bend_radius=40,
    bend_side=1,
    output_bend=True,
    center_taper_length=20,
    center_straight_length=10,
    full_etch=False,
    label="1550",
    label2="",
    waveguide_label_size=5,
    doubleetch_label_size=20,
    waveguide_label_layer=Layers.wg,
    doubleetch_label_layer=Layers.wg_doubleetch,
    alignment_tick_layer=Layers.wg,
    alignment_tick_window_layer=None,
    alignment_tick_window_width=2.0,
    alignment_tick_window_length=8.0,
    alignment_tick_window_inward_shift=0.0,
    alignment_tick_x_offset=0.0,
    glue_layers=(Layers.wg_doubleetch, Layers.via),
    oxide_open_layer=Layers.via,
    oxide_open_angle=None,
    oxide_open_length=2e3,
    oxide_open_start_width=2.0,
    oxide_open_start_at_taper_end=False,
    oxide_open_start_offset=0.0,
    oxide_open_duplicate_layer=None,
    oxide_open_duplicate_angle=None,
    oxide_open_duplicate_length=None,
    oxide_open_duplicate_start_width=None,
    oxide_open_duplicate_start_at_taper_end=None,
    oxide_open_duplicate_start_offset=None,
):
    """
    Tapered fiber coupler plus mode converter, with glue spots, alignment ticks
    and logos.

    The fiber is drawn as a footprint running east from the coupler tip at the origin,
    with two glue-spot windows under it and an oxide opening around the taper. ``o1``
    leaves west at the taper's back; ``o2`` leaves east through a 30 deg bend of
    *bend_radius*, turned up or down by the signs of *bend_side* and *bend_radius*.

    Parameters
    ----------
    tapered_fiber_angle : float
        Drawn taper angle of the fiber, degrees.
    tapered_fiber_length : float
        Distance from the tip to where the fiber reaches its full 250 um width.
    center_taper_length : float
        Length of the width transition on the center/S-bend side. Set to zero
        to expose the underlying taper endpoint directly.
    center_straight_length : float
        Straight extension after the center-side transition. Set to zero for
        no additional straight section before the S-bend/tether junction.
    full_etch : bool
        Replaces the double-etch geometry with its complement inside a window over the
        taper, for a fully etched release.
    label, label2 : str
        Text written above and below the coupler, twice: small on the waveguide layer at
        the taper, large on the double-etch layer between the glue spots.
    waveguide_label_size, doubleetch_label_size : float
        Text sizes in um for the waveguide- and double-etch-layer copies.
    waveguide_label_layer, doubleetch_label_layer, alignment_tick_layer : LayerSpec
        Layers for the small taper label, large packaging label, and fiber-tip
        alignment ticks.
    alignment_tick_window_layer : LayerSpec or None
        Optional layer for paired side windows along each 0.5 x 10 um
        fiber-tip alignment tick.
    alignment_tick_window_width, alignment_tick_window_length : float
        Width and vertical length of each side window.
    alignment_tick_window_inward_shift : float
        Moves both side windows toward the optical centerline so they can
        overlap the main release opening without moving the alignment ticks.
    alignment_tick_x_offset : float
        Longitudinal offset of the alignment ticks and their side windows from
        the center-side taper port. Use this to retain a process landmark at
        an internal taper junction.
    glue_layers : sequence of LayerSpec
        Layers receiving the two 300 um-diameter glue circles. Use one layer
        when the circles belong only to a coarse-write mask.
    oxide_open_layer : LayerSpec or None
        Layer for the tapered oxide-opening polygon. Set to ``None`` to omit
        the primary oxide-opening pattern.
    oxide_open_angle : float or None
        Independent full taper angle in degrees for the oxide opening. ``None``
        preserves the historical behavior of using ``tapered_fiber_angle``.
    oxide_open_length : float
        Longitudinal length in um of the Layer 14 tapered opening.
    oxide_open_start_width : float
        Transverse width in um at the primary tapered-slot start point.
    oxide_open_start_at_taper_end : bool
        When true, start the opening at the outer ``o2`` taper endpoint rather
        than the historical center/S-bend-side ``o1`` endpoint.
    oxide_open_start_offset : float
        Additional outward distance from the selected primary start point.
    oxide_open_duplicate_layer : LayerSpec or None
        Optional layer receiving an independently configurable duplicate.
    oxide_open_duplicate_angle, oxide_open_duplicate_length : float or None
        Duplicate geometry. ``None`` inherits the primary angle or length.
    oxide_open_duplicate_start_width : float or None
        Duplicate start width. ``None`` inherits the primary start width.
    oxide_open_duplicate_start_at_taper_end : bool or None
        Duplicate start-point rule. ``None`` inherits the primary rule.
    oxide_open_duplicate_start_offset : float or None
        Duplicate outward start offset. ``None`` inherits the primary offset.
    """
    c = gf.Component()
    if alignment_tick_window_width <= 0:
        raise ValueError("alignment_tick_window_width must be positive")
    if alignment_tick_window_length <= 0:
        raise ValueError("alignment_tick_window_length must be positive")
    if alignment_tick_window_inward_shift < 0:
        raise ValueError("alignment_tick_window_inward_shift cannot be negative")

    # Fiber taper.
    taper = c.add_ref(
        tapered_fiber_spec()
    )

    # Optional output bend and width transition
    if output_bend:
        radius = abs(bend_radius)
        angle = -bend_side * np.sign(bend_radius) * 30

        p = taper.ports["o2"].bend_by(
            angle,
            radius=radius,
            add_to_component=c,
        )

        p = p.extend(
            20,
            add_to_component=c,
            new_width=0.2,
        )
    else:
        # Expose the original straight taper endpoint directly.
        p = taper.ports["o2"]


    # Center-side endpoint. Optional transition/straight lengths make the
    # coupler compact without changing the historical defaults.
    q2 = taper.ports["o1"]
    if center_taper_length > 0:
        q2 = q2.extend(
            center_taper_length,
            add_to_component=c,
            new_width=0.8,
        )
    if center_straight_length > 0:
        q2 = q2.extend(
            center_straight_length,
            add_to_component=c,
            new_width=0.8,
        )

    c.add_ports(
        [p, q2]
    )

    # Triangle representation of the fiber.
    tapered_fiber_width = tapered_fiber_length * np.tan(np.radians(tapered_fiber_angle))
    tapered_fiber_extension = 5e3
    c.add_polygon(
        [
            (tapered_fiber_length, -tapered_fiber_width/2),
            (tapered_fiber_length, -125),
            (tapered_fiber_length+tapered_fiber_extension, -125),
            (tapered_fiber_length+tapered_fiber_extension, 125),
            (tapered_fiber_length, 125),
            (tapered_fiber_length, tapered_fiber_width/2),
            (0, 0),
        ],
        layer=Layers.fiber,
    )

    # Ticks for aligning the tip of the fiber.
    for i in [-1, 1]:
        x1 = alignment_tick_x_offset
        w = 15
        d1 = .25
        d2 = 5
        c.add_polygon(
            [
                (x1+d1, i*w/2+d2),
                (x1+d1, i*w/2-d2),
                (x1-d1, i*w/2-d2),
                (x1-d1, i*w/2+d2),
            ],
            layer=alignment_tick_layer,
        )
        if alignment_tick_window_layer is not None:
            x_offset = d1 + alignment_tick_window_width / 2
            for side in (-1, 1):
                window = c.add_ref(
                    gf.components.rectangle(
                        size=(
                            alignment_tick_window_width,
                            alignment_tick_window_length,
                        ),
                        centered=True,
                        layer=alignment_tick_window_layer,
                    )
                )
                window.move(
                    (
                        x1 + side * x_offset,
                        i * (w / 2 - alignment_tick_window_inward_shift),
                    )
                )

    # Glue spots
    glue_spot_width = 300
    glue_spot_offset = 450 + taper.ports["o2"].center[0]

    X = [glue_spot_offset, glue_spot_offset+500]

    for x in X:
        for glue_layer in glue_layers:
            c.add_ref(
                gf.components.circle(
                    radius=glue_spot_width/2,
                    layer=glue_layer,
                )
            ).move((x, 0))

    # Open oxide: primary pattern plus an optional independent duplicate.
    x0 = (
        p.center[0] if oxide_open_start_at_taper_end else q2.center[0]
    ) + oxide_open_start_offset
    if oxide_open_length <= 0:
        raise ValueError("oxide_open_length must be positive")
    if oxide_open_start_width <= 0:
        raise ValueError("oxide_open_start_width must be positive")
    if oxide_open_angle is None:
        oxide_open_angle = tapered_fiber_angle
    if not 0 <= oxide_open_angle < 90:
        raise ValueError("oxide_open_angle must satisfy 0 <= angle < 90 degrees")

    def add_oxide_opening(layer, angle, length, start_x, start_width):
        if length <= 0:
            raise ValueError("oxide opening length must be positive")
        if not 0 <= angle < 90:
            raise ValueError(
                "oxide opening angle must satisfy 0 <= angle < 90 degrees"
            )
        if start_width <= 0:
            raise ValueError("oxide opening start width must be positive")
        start_half_width = start_width / 2
        open_width = length * np.tan(np.radians(angle))
        c.add_polygon(
            [
                (start_x + length, 2 - open_width / 2),
                (start_x, -start_half_width),
                (start_x, start_half_width),
                (start_x + length, 2 + open_width / 2),
            ],
            layer=layer,
        )

    if oxide_open_layer is not None:
        add_oxide_opening(
            oxide_open_layer,
            float(oxide_open_angle),
            float(oxide_open_length),
            float(x0),
            float(oxide_open_start_width),
        )
    if oxide_open_duplicate_layer is not None:
        duplicate_angle = (
            oxide_open_angle
            if oxide_open_duplicate_angle is None
            else oxide_open_duplicate_angle
        )
        duplicate_length = (
            oxide_open_length
            if oxide_open_duplicate_length is None
            else oxide_open_duplicate_length
        )
        duplicate_start_width = (
            oxide_open_start_width
            if oxide_open_duplicate_start_width is None
            else oxide_open_duplicate_start_width
        )
        duplicate_starts_at_taper_end = (
            oxide_open_start_at_taper_end
            if oxide_open_duplicate_start_at_taper_end is None
            else oxide_open_duplicate_start_at_taper_end
        )
        duplicate_x0 = (
            p.center[0]
            if duplicate_starts_at_taper_end
            else q2.center[0]
        )
        duplicate_x0 += (
            oxide_open_start_offset
            if oxide_open_duplicate_start_offset is None
            else oxide_open_duplicate_start_offset
        )
        add_oxide_opening(
            oxide_open_duplicate_layer,
            float(duplicate_angle),
            float(duplicate_length),
            float(duplicate_x0),
            float(duplicate_start_width),
        )

    # Logos
    x = taper.ports["o2"].center[0] / 2
    #c.add_ref(get_cal(width=20, layer=Layers.wg)).move((x, 25))
    #c.add_ref(get_ionq(width=16, layer=Layers.wg)).move((x, -25))

    # Text
    for x, layer, size in zip(
        [-20, glue_spot_offset+250],
        [waveguide_label_layer, doubleetch_label_layer],
        [waveguide_label_size, doubleetch_label_size],
    ):
    # for x, layer, size in zip(
    #     [taper.ports["o1"].center[0]],
    #     [Layers.wg],
    #     [5],
    # ):
        for s, label_ in zip([1, -1], [label, label2]):
            c.add_ref(
                text(
                    text=label_,
                    layer=layer,
                    size=size,
                    justify="center",
                )
            ).move((x, (1*s - .5)*size + s*15))

    if full_etch:
        cc = gf.Component(c.name + "_full_etch")
        full_etch_shift = 250
        cc.add_ref(
            gf.components.rectangle(
                size=(open_length-full_etch_shift, 100),
                centered=True,
                layer=Layers.label,
            )
        ).move((full_etch_shift + (open_length-full_etch_shift)/2, 0))

        ccc = gf.boolean(
            cc,
            c,
            operation="not",
            layer=Layers.label,
            layer1=Layers.label,
            layer2=Layers.wg_doubleetch,
        )
        c.add_ref(ccc)
        c.remove_layers([Layers.wg_doubleetch])

    return c

@gf.cell
def edge_coupler(
    length=100,
    width=3,
    taper_length=100,
    taper_width=.8,
):
    """Facet edge coupler: a straight of *width* over *length*, tapered down to
    *taper_width*.

    ``o1`` is the facet end, ``o2`` the chip end.
    """
    c = gf.Component()

    p = gf.Port(
        name="o1",
        center=(-length/2, 0),
        width=width,
        orientation=0,
        layer=Layers.wg,
    )

    c.add_port("o1", port=p)

    p = p.extend(length, add_to_component=c)
    p = p.extend(taper_length, new_width=taper_width, add_to_component=c)

    c.add_port("o2", port=p)

    return c

@gf.cell
def make_ec(
    component,
    position,
    chip_length,
    horizontal=True,
    label="",
    dy=None,
    stub_length=20,
    routing_width=.8,
    radius=200,
    taper_length=500,
    edge_length=1e3,
    edge_width=3,
    label_size=20,
):
    """Run every port of *component* out to a chip facet and cap it with an `edge_coupler`.

    *component* is placed at *position* and each of its ports is bent toward a facet at
    ``+-chip_length/2``, routed out with an ``ata`` route and terminated in an
    `edge_coupler`. A port already on the *horizontal* axis runs straight out to its own
    side; a port on the other axis takes a quarter turn, putting the two ends of that
    axis into a Z: with *horizontal* unset the west-facing ports go north and the
    east-facing ports go south, with it set the north-facing ports go west and the
    south-facing ports go east. A port off the 90 deg grid raises.

    The returned component carries *component*'s own ports plus one ``ec_<name>`` facet
    port per coupler, continuing that port's ``pathlength`` chain.

    Parameters
    ----------
    position : (float, float)
        Placement of *component*.
    chip_length : float
        Facet-to-facet span; each coupler's facet end lands at ``+-chip_length/2``.
    horizontal : bool
        Axis the couplers run along.
    label : str
        Written beside each coupler, suffixed with the port name.
    dy : float OR None
        Coupler pitch at the facet, so a fan of ports lands on a regular spacing centred
        on the fan. ``None`` keeps each port on its own line.
    stub_length : float
        Narrowing stub drawn at each port ahead of the route.
    routing_width : float
        Width of the stub, of the routed straights and of the coupler's chip end.
    taper_length, edge_length, edge_width : float
        `edge_coupler` geometry: widening taper, facet run, facet width.
    """
    c = gf.Component()
    ref = c.add_ref(component).move(position)
    c.add_ports(list(ref.ports))

    # A 2x2 device names each port twice (o-number and corner); route each guide once,
    # under its o-number where it has one.
    physical = {}
    for port in ref.ports:
        key = (round(float(port.center[0]), 3), round(float(port.center[1]), 3),
               round(float(port.orientation), 3))
        numbered = port.name[:1] == "o" and port.name[1:].isdigit()
        kept = physical.get(key)
        if kept is None or (numbered and not (kept.name[:1] == "o" and kept.name[1:].isdigit())):
            physical[key] = port

    # Which facet each port leaves through, as a quarter-turn index east/north/west/south.
    routes = []
    for port in physical.values():
        turns = float(port.orientation) / 90
        if abs(turns - round(turns)) > 1e-6:
            raise ValueError(
                f"make_ec: port {port.name!r} at {port.orientation:g} deg is off the "
                "90 deg grid")
        axis = int(round(turns)) % 4
        out = axis if (axis % 2 == 0) == bool(horizontal) else 3 - axis
        p = port.extend(stub_length, new_width=routing_width, add_to_component=c)
        routes.append(
            (out, port.name, p.bend_towards(90*out, radius=radius, add_to_component=c))
        )

    for out in sorted({r[0] for r in routes}):
        lat = 1 - out % 2
        sign = 1 if out < 2 else -1
        group = sorted((r for r in routes if r[0] == out), key=lambda r: r[2].center[lat])

        offsets = [r[2].center[lat] for r in group]
        if dy is not None:
            offsets = np.mean(offsets) + dy*(np.arange(len(group)) - (len(group) - 1)/2)

        for (_, name, p), offset in zip(group, offsets):
            facet = np.zeros(2)
            facet[out % 2] = sign * chip_length/2
            facet[lat] = offset

            ec = c.add_ref(
                edge_coupler(
                    length=edge_length,
                    width=edge_width,
                    taper_length=taper_length,
                    taper_width=routing_width,
                )
            ).rotate((90*out + 180) % 360)
            ec.move(tuple(float(v) for v in facet - ec.ports["o1"].position))

            q = p.route_to(
                ec.ports["o2"],
                method="ata",
                add_to_component=c,
                routing_width=routing_width,
                radius=radius,
            )
            c.add_ports([
                ec.ports["o1"].copy(
                    name=f"ec_{name}",
                    pathlength=q.pathlength + edge_length + taper_length,
                )
            ])

            text_position = facet.copy()
            text_position[out % 2] -= sign * edge_length/2
            text_position[lat] += 2*label_size
            c.add_ref(
                text(text=label + str(name), size=label_size, justify="center")
            ).move(tuple(float(v) for v in text_position))

    return c

@gf.cell
def tapered_fiber_test(
    tapered_fiber_spec,
    spacing_y=300,
    tapered_fiber_angle=10,
):
    """
    Same-side, U-shaped tapered fiber coupler test.

    Two couplers *spacing_y* apart, each splayed out by *tapered_fiber_angle*, with
    their ``o1`` ends joined by a `connect` circular route.
    """
    c = gf.Component()

    ports = []

    for s in [-1, 1]:
        dev = c.add_ref(
            tapered_fiber_spec()
        ).rotate(tapered_fiber_angle*s).move((0, s*spacing_y/2))
        ports.append(dev.ports["o1"])

    connect(
        ports[0],
        ports[1],
        method="circular",
        add_to_component=c,
    )

    c.flatten()

    return c

# 780 nm tapered fiber coupler
tapered_fiber_780 = partial(
    tapered_fiber_coupler,
    tapered_fiber_spec=partial(
        tapered_fiber_3d,
        taper_spec=partial(
            tapered_fiber_linear,
            width1=.5,
            width2=.2,
            length=20,
        )
    ),
)

# 1550 nm tapered fiber coupler
tapered_fiber_1550 = partial(
    tapered_fiber_coupler,
    tapered_fiber_spec=partial(
        tapered_fiber_3d,
        taper_spec=partial(
            tapered_fiber_linear,
            width1=.7,
            width2=.2,
            length=50,
        )
    ),
)

@gf.cell
def tapered_fiber_device(
    tapered_fiber_spec1=None,
    tapered_fiber_spec2=tapered_fiber_1550,
    edge_coupler_spec=edge_coupler,
    chip_length=21e3,
    tapered_fiber_angle=10,
    taper_same_side=False,
    device_spec=None,
    **kwargs,
):
    """
    Full device with tapered fiber couplers and long waveguide inbetween,
    potentially with a device like a mux or MZI in the middle for more
    complicated DFG or SQZ experiments.

    Fiber 2 sits near the east facet; fiber 1 goes to the west facet, or beside fiber 2
    when *taper_same_side* is set. Each fiber's ``o2`` is routed out to its own
    *edge_coupler_spec* at the facet, and the free ``o1`` ends are the arms the middle
    device is spliced between. With no fiber 1 and no device, the west edge coupler is
    wired straight across to fiber 2's arm.

    Parameters
    ----------
    tapered_fiber_spec1, tapered_fiber_spec2 : ComponentFactory OR None
        Fiber-coupler factories, called with ``bend_side`` and *kwargs*. None omits one.
    chip_length : float
        Facet-to-facet width; the edge couplers land at ``+-chip_length/2``.
    tapered_fiber_angle : float
        Rotation applied to each fiber coupler, degrees.
    taper_same_side : bool OR float
        Puts fiber 1 on the same (east) side as fiber 2, offset this many um in y.
        ``True`` means 300; ``False`` sends fiber 1 to the west facet instead.
    device_spec : ComponentFactory OR None
        4-port device to splice between the two fiber arms. Its two unused ports are
        routed to their own edge couplers as diagnostics.
    **kwargs
        Forwarded to both fiber-coupler specs.

    Raises
    ------
    ValueError
        For a combination whose arms are never drawn: *device_spec* needs both
        fibers, and omitting *tapered_fiber_spec1* needs *tapered_fiber_spec2* to
        wire the west edge coupler across to, and rules out *taper_same_side*.
    """
    if device_spec is not None:
        missing = [n for n, s in (("tapered_fiber_spec1", tapered_fiber_spec1),
                                  ("tapered_fiber_spec2", tapered_fiber_spec2)) if s is None]
        if missing:
            raise ValueError(
                "tapered_fiber_device: device_spec is spliced between the two fiber "
                f"arms, so it requires {' and '.join(missing)}")
    if tapered_fiber_spec1 is None and tapered_fiber_spec2 is None:
        raise ValueError(
            "tapered_fiber_device: with no tapered_fiber_spec1 the west edge coupler is "
            "wired across to fiber 2, so it requires tapered_fiber_spec2")
    if tapered_fiber_spec1 is None and taper_same_side:
        raise ValueError(
            "tapered_fiber_device: taper_same_side puts fiber 1 beside fiber 2, so it "
            "requires tapered_fiber_spec1")

    c = gf.Component()

    taper_position = chip_length/2 - 1.5e3

    if tapered_fiber_spec2 is not None:
        tf2 = c.add_ref(
            tapered_fiber_spec2(bend_side=+1, **kwargs)
        ).rotate(tapered_fiber_angle).move((taper_position, 75))

        p = tf2.ports["o1"].bend_west()
        tf2.move((0, -p.center[1]))
        p = tf2.ports["o1"].bend_west(add_to_component=c, radius=500)

        ec2 = c.add_ref(
            edge_coupler_spec()
        ).rotate(180).move((chip_length/2, -40))

        connect(
            ec2.ports["o2"],
            tf2.ports["o2"],
            method="intersection",
            add_to_component=c,
            radius=500,
            # bend_
        )

    if tapered_fiber_spec1 is not None:
        tf1 = c.add_ref(
            tapered_fiber_spec1(bend_side=(1 if taper_same_side else -1), **kwargs)
        ).rotate(tapered_fiber_angle)

        q = tf1.ports["o1"].bend_west()
        tf1.move((0, -q.center[1]))

        if taper_same_side:
            if taper_same_side is True:
                taper_same_side = 300
            tf1.move((taper_position, taper_same_side))
        else:
            tf1.move((-taper_position-100, 100))

        q = tf1.ports["o1"].bend_west(add_to_component=c)

        if taper_same_side:
            ec1 = c.add_ref(
                edge_coupler_spec()
            ).rotate(180).move((chip_length/2, taper_same_side-40))
            connect(
                ec1.ports["o2"],
                tf1.ports["o2"],
                method="intersection",
                add_to_component=c,
                radius=500,
                # bend_
            )

    if not taper_same_side:
        ec1 = c.add_ref(
            edge_coupler_spec()
        )

        if tapered_fiber_spec1 is not None:
            qq = tf1.ports["o2"].bend_west(add_to_component=c, radius=150)
            ec1.move((-chip_length/2, qq.center[1]))
            connect(
                ec1.ports["o2"],
                qq,
                method="intersection",
                add_to_component=c,
                radius=150,
                # bend_
            )
        else:
            ec1.move((-chip_length/2, 0))



    kwargs_connect = dict(
        method="intersection",
        add_to_component=c,
        routing_width=2,
        bend_width=.8,
        radius=150,
        taper_length=50,
    )

    if device_spec is not None:
        device_postion = chip_length/2 - 1e3
        dev = c.add_ref(
            device_spec()
        )
        if taper_same_side:
            dev.move((
                -dev.ports["o1"].center[0]-device_postion,
                -dev.ports["o1"].center[1]+taper_same_side,
            ))

            connect(
                dev.ports["o1"],
                p,
                **kwargs_connect,
            )
            connect(
                dev.ports["o4"],
                q,
                **kwargs_connect,
            )
        else:
            dev.move((
                -dev.ports["o2"].center[0]-device_postion,
                -dev.ports["o2"].center[1],
            ))

            connect(
                dev.ports["o4"],
                p,
                **kwargs_connect,
            )

            connect(
                dev.ports["o2"],
                q,
                **kwargs_connect,
            )

        # Diagnostic ports
        if taper_same_side:
            ports = ["o2", "o3"]
        else:
            ports = ["o1", "o3"]
        yy = [-90, 400]
        for i, port in enumerate(ports):
            pp = dev.ports[port]
            if i == 0:
                y = pp.center[1]
            else:
                y = yy[taper_same_side != 0]

            j = -(2*i - 1)
            ec = c.add_ref(
                edge_coupler_spec()
            )
            ec.rotate(i*180)
            ec.move((-j * chip_length/2, y))

            if i == 1:
                pp = pp.copy(center=(pp.center[0]+300, y))
                connect(
                    dev.ports[port],
                    pp.flipped(),
                    method="circular",
                    add_to_component=c,
                    radius=150,
                    # bend_
                )

            connect(
                ec.ports["o2"],
                pp,
                method="intersection",
                add_to_component=c,
                radius=150,
                # bend_
            )



    if tapered_fiber_spec1 is None: # Simple case
        connect(
            ec1.ports["o2"],
            p,
            **kwargs_connect,
        )
    else:
        pass




    c.flatten()

    return c

@gf.cell
def tapered_fiber_inout(
    tapered_fiber_spec,
    edge_coupler_spec=edge_coupler,
    tether_spec=None,
    include_edge_couplers=True,
    extra_rectangle_size=None,
    extra_rectangle_gap=100,
    extra_rectangle_layer=Layers.wg_doubleetch,
    extra_rectangle_fill_layer=None,
    extra_rectangle_margin=0.0,
    extra_rectangle_margin_layer=None,
    extra_rectangle_surround_slot_width=None,
    tether_rectangle_slot_width=None,
    tether_rectangle_slot_layer=Layers.wg_doubleetch,
    tether_rectangle_slot_overlap=1.0,
    tether_rectangle_slot_clearance=0.0,
    tether_rectangle_slot_start_offset=0.0,
    tether_rectangle_slot_tether_margin=None,
    tether_rectangle_slot_transition_length=20.0,
    tether_rectangle_slot_transition_points=32,
    s_bend_width=0.8,
    s_bend_taper_length=10,
    s_bend_half_length=200,
    s_bend_compact_180_radius=None,
    s_bend_compact_euler_p=0.5,
    s_bend_compact_straight_length=0.0,
    s_bend_vertical_half_span=100.0,
    s_bend_layer=None,
    s_bend_endpoint_spec=None,
    ridge_s_bend_layer=None,
    ridge_s_bend_width=None,
    s_bend_support_shield_spec=None,
    s_bend_support_spec=None,
    s_bend_support_layer=None,
    s_bend_support_pitch=5.0,
    s_bend_support_neck_width=0.2,
    s_bend_support_outer_width=1.0,
    s_bend_support_neck_length=1.5,
    s_bend_support_taper_length=1.5,
    s_bend_support_overlap=0.2,
    s_bend_support_end_clearance=5.0,
    chip_length=5e3,
    tapered_fiber_angle=5,
    flatten=True,
    **kwargs,
):
    """
    Different-side tapered fiber in/out coupler test, connected by an S-bend waveguide.

    The two couplers face opposite facets, each ``o2`` bending out to its own
    `edge_coupler` at ``+-chip_length/2`` while their ``o1`` ends meet in the middle.

    If ``include_edge_couplers`` is false, the long outer bends, diagonal routes,
    and facet edge couplers are omitted. In that mode, ``tether_spec`` connects
    directly to each tapered-fiber coupler's exposed ``o2`` endpoint. The center
    S-bend is retained unchanged.

    If ``extra_rectangle_size`` is provided as ``(length, width)``, a rectangle
    is placed beyond each taper ``o2`` port. ``extra_rectangle_gap`` is the
    edge-to-edge distance from the taper endpoint to the nearest rectangle edge.
    A positive ``extra_rectangle_margin`` adds a second, independent rectangle
    enlarged by that amount on all four sides. It is placed on
    ``extra_rectangle_margin_layer`` (or ``extra_rectangle_layer`` when none is
    specified) and is not boolean-merged with the nominal rectangle.

    If ``tether_rectangle_slot_width`` is provided, a slot on
    ``tether_rectangle_slot_layer`` joins the far pad of each tether to the
    nearest edge of its extra rectangle. ``tether_rectangle_slot_overlap``
    extends the slot under both endpoints to avoid a lithographic seam, while
    ``tether_rectangle_slot_clearance`` shortens it at both ends to enforce a
    non-overlapping gap. ``tether_rectangle_slot_start_offset`` moves only the
    tether-side endpoint farther toward the rectangle, which can make the slot
    merge with an offset release halo. If
    ``tether_rectangle_slot_tether_margin`` is provided, only the section next
    to the tether is smoothly flared to the tether's actual transverse extent
    plus that margin. ``tether_rectangle_slot_transition_length`` controls the
    flare length; the remainder stays at ``tether_rectangle_slot_width``.
    ``tether_rectangle_slot_transition_points`` controls curve sampling.
    Overlap and clearance may be combined; their net endpoint extension is
    ``overlap - clearance``.

    ``s_bend_width`` fixes the center S-bend body width. Short transitions of
    ``s_bend_taper_length`` connect the coupler ports to that width.
    ``s_bend_half_length`` controls the horizontal placement of each taper
    relative to the layout center. If ``ridge_s_bend_layer`` is provided, a
    second constant-width copy of the complete center S-bend is drawn on that
    layer using ``ridge_s_bend_width``.
    If ``s_bend_compact_180_radius`` is provided, the broad circular S-route
    is replaced by straight leads and two mirrored 180-degree Euler bends.
    ``s_bend_compact_euler_p`` is the Euler fraction: 0 is circular, 0.5 is
    the gdsfactory default partial Euler bend, and 1 is a full Euler bend.
    The effective radius must match one quarter of the port y-separation.
    ``s_bend_compact_straight_length`` inserts a straight meander arm between
    the two 180-degree bends; the two outer leads are adjusted symmetrically
    so the taper endpoint coordinates remain unchanged.
    ``s_bend_vertical_half_span`` places the two taper centers at equal and
    opposite y coordinates and therefore sets half their vertical separation.
    ``s_bend_layer`` optionally draws the meander body on a layer independent
    of the taper ports. If ``s_bend_endpoint_spec`` is provided, one component
    is inserted between each taper and the meander. It must expose anchor
    ports ``a1`` (taper side) and ``a2`` (meander side).
    If ``s_bend_support_shield_spec`` is provided, it is called with the
    complete centerline path and should return the continuous outer-slab PNC
    mask. Otherwise, if ``s_bend_support_spec`` is provided, that component is sampled along
    both sides of the complete S-bend with its ``a1`` port touching the
    waveguide and its ``a2`` port pointing outward. Otherwise, when
    ``s_bend_support_layer`` is provided, symmetric tapered anchors are used.
    relative to the layout center; reducing it shortens the S-bend.

    Set ``flatten`` false to retain the support rectangles and other child
    geometry as separate GDS hierarchy objects. This is useful when multiple
    raw mask sources share one layer but must not be geometrically merged.
    """
    c = gf.Component()

    if s_bend_half_length <= 0:
        raise ValueError("s_bend_half_length must be positive")
    taper_position = s_bend_half_length

    ports = []
    outer_ports = []
    taper_end_ports = []
    tether_slot_lengths = []
    tether_slot_transverse_extents = []

    for i in range(2):
        j = 2*i - 1
        ang = i*180

        tf = c.add_ref(
            tapered_fiber_spec(bend_side=+1, **kwargs)
        ).rotate(180+ang+tapered_fiber_angle).move(
            (j * taper_position, -j * s_bend_vertical_half_span)
        )

        # Center connection
        ports.append(tf.ports["o1"])
        taper_end_ports.append(tf.ports["o2"])

        if include_edge_couplers:
            # Optional long route to a facet edge coupler.
            p = tf.ports["o2"].bend_towards(
                180 + ang,
                add_to_component=c,
                radius=1000,
            )
            ec = c.add_ref(
                edge_coupler_spec()
            ).rotate(ang).move((j * chip_length / 2, p.center[1]))
            outer_ports.append(ec.ports["o1"])
            connect(
                p,
                ec.ports["o2"],
                add_to_component=c,
            )
        else:
            # Direct tether connection at the tapered-fiber coupler endpoint.
            outer_ports.append(tf.ports["o2"])


    route_ports = ports
    if s_bend_endpoint_spec is not None:
        endpoint_component = s_bend_endpoint_spec()
        try:
            endpoint_component.ports["a1"]
            endpoint_component.ports["a2"]
        except KeyError as error:
            raise ValueError(
                "s_bend_endpoint_spec must provide anchor ports 'a1' and 'a2'"
            ) from error
        route_ports = []
        for port in ports:
            endpoint = c.add_ref(endpoint_component)
            # A phononic neck intentionally meets the center of a much wider
            # slab/taper face, as in the source PNC geometry.
            endpoint.connect(
                "a1",
                port,
                allow_width_mismatch=True,
                allow_layer_mismatch=True,
                allow_type_mismatch=True,
            )
            route_ports.append(endpoint.ports["a2"])

    meander_layer = (
        s_bend_layer if s_bend_layer is not None else route_ports[0].layer
    )

    if s_bend_compact_180_radius is None:
        s_bend_centerline = connect(
            route_ports[0],
            route_ports[1],
            method="circular",
            add_to_component=False,
        )
        connect(
            route_ports[0],
            route_ports[1],
            method="circular",
            bend_width=s_bend_width,
            taper_length=s_bend_taper_length,
            cross_section=partial(
                gf.cross_section.cross_section,
                width=s_bend_width,
                layer=meander_layer,
            ),
            add_to_component=c,
        )
    else:
        if s_bend_compact_180_radius <= 0:
            raise ValueError("s_bend_compact_180_radius must be positive")
        if not 0 <= s_bend_compact_euler_p <= 1:
            raise ValueError("s_bend_compact_euler_p must be between 0 and 1")
        if s_bend_compact_straight_length < 0:
            raise ValueError(
                "s_bend_compact_straight_length must be non-negative"
            )
        theta = np.radians(float(route_ports[0].orientation))
        rotation = np.array(
            [[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]]
        )
        local_delta = (
            np.asarray(route_ports[1].center, dtype=float)
            - np.asarray(route_ports[0].center, dtype=float)
        ) @ rotation
        expected_dy = 4 * s_bend_compact_180_radius
        if not np.isclose(abs(local_delta[1]), expected_dy, atol=1e-3):
            raise ValueError(
                "Two 180-degree Euler bends require abs(port dy) = "
                f"4 * radius; got {abs(local_delta[1]):g} and "
                f"radius {s_bend_compact_180_radius:g}"
            )
        turn_sign = 1 if local_delta[1] > 0 else -1
        s_bend_centerline = gf.Path()
        # The middle straight is traversed in the reverse x direction after
        # the first 180-degree bend. Add that length back to the two outer
        # leads, split evenly, so the input/output coordinates stay fixed.
        half_lead = (
            local_delta[0] + s_bend_compact_straight_length
        ) / 2
        if half_lead < 0:
            raise ValueError(
                "Meander straight length and port separation produce a "
                "negative outer lead"
            )
        s_bend_centerline.append(gf.path.straight(length=half_lead))
        s_bend_centerline.append(
            gf.path.euler(
                radius=s_bend_compact_180_radius,
                angle=turn_sign * 180,
                p=s_bend_compact_euler_p,
                use_eff=True,
            )
        )
        if s_bend_compact_straight_length > 0:
            s_bend_centerline.append(
                gf.path.straight(
                    length=s_bend_compact_straight_length
                )
            )
        s_bend_centerline.append(
            gf.path.euler(
                radius=s_bend_compact_180_radius,
                angle=-turn_sign * 180,
                p=s_bend_compact_euler_p,
                use_eff=True,
            )
        )
        s_bend_centerline.append(gf.path.straight(length=half_lead))
        s_bend_centerline.rotate(float(route_ports[0].orientation))
        s_bend_centerline.move(route_ports[0].center)
        c.add_ref(
            gf.path.extrude(
                s_bend_centerline,
                cross_section=partial(
                    gf.cross_section.cross_section,
                    width=s_bend_width,
                    layer=meander_layer,
                ),
            )
        )
    if s_bend_support_shield_spec is not None:
        c.add_ref(s_bend_support_shield_spec(path=s_bend_centerline))
    elif s_bend_support_spec is not None or s_bend_support_layer is not None:
        if s_bend_support_pitch <= 0:
            raise ValueError("s_bend_support_pitch must be positive")
        if s_bend_support_spec is None and min(
                s_bend_support_neck_width,
                s_bend_support_outer_width,
                s_bend_support_neck_length,
                s_bend_support_taper_length,
        ) <= 0:
            raise ValueError("S-bend support dimensions must be positive")

        support_component = (
            s_bend_support_spec()
            if s_bend_support_spec is not None
            else None
        )
        if support_component is not None:
            try:
                support_component.ports["a1"]
                support_component.ports["a2"]
            except KeyError as error:
                raise ValueError(
                    "s_bend_support_spec must provide anchor ports 'a1' and 'a2'"
                ) from error

        points = np.asarray(s_bend_centerline.points, dtype=float)
        segment_vectors = np.diff(points, axis=0)
        segment_lengths = np.linalg.norm(segment_vectors, axis=1)
        cumulative = np.concatenate(([0.0], np.cumsum(segment_lengths)))
        total_length = float(cumulative[-1])
        sample_distances = np.arange(
            s_bend_support_end_clearance,
            total_length - s_bend_support_end_clearance + 1e-9,
            s_bend_support_pitch,
        )
        for distance in sample_distances:
            segment_index = min(
                int(np.searchsorted(cumulative, distance, side="right") - 1),
                len(segment_lengths) - 1,
            )
            segment_length = segment_lengths[segment_index]
            if segment_length <= 0:
                continue
            fraction = (
                distance - cumulative[segment_index]
            ) / segment_length
            center = points[segment_index] + fraction * segment_vectors[
                segment_index
            ]
            tangent = segment_vectors[segment_index] / segment_length
            normal = np.array([-tangent[1], tangent[0]])
            for side in (-1.0, 1.0):
                inner = center + side * normal * (
                    s_bend_width / 2 - s_bend_support_overlap
                )
                if support_component is not None:
                    outward = side * normal
                    angle = float(np.arctan2(outward[1], outward[0]))
                    rotation = np.array(
                        [
                            [np.cos(angle), -np.sin(angle)],
                            [np.sin(angle), np.cos(angle)],
                        ]
                    )
                    anchor = np.asarray(
                        support_component.ports["a1"].center,
                        dtype=float,
                    )
                    # Add transformed polygons directly. This preserves the
                    # arbitrary tangent angle without leaving off-grid GDS
                    # instances in the returned cell hierarchy.
                    for layer, polygons in support_component.get_polygons_points(
                        by="tuple"
                    ).items():
                        for polygon in polygons:
                            transformed = (
                                (np.asarray(polygon) - anchor) @ rotation.T
                                + np.asarray(inner, dtype=float)
                            )
                            c.add_polygon(transformed, layer=layer)
                    continue
                outer = center + side * normal * (
                    s_bend_width / 2
                    + s_bend_support_neck_length
                    + s_bend_support_taper_length
                )
                shoulder = center + side * normal * (
                    s_bend_width / 2 + s_bend_support_neck_length
                )
                polygon = [
                    inner - tangent * s_bend_support_neck_width / 2,
                    inner + tangent * s_bend_support_neck_width / 2,
                    shoulder + tangent * s_bend_support_neck_width / 2,
                    outer + tangent * s_bend_support_outer_width / 2,
                    outer - tangent * s_bend_support_outer_width / 2,
                    shoulder - tangent * s_bend_support_neck_width / 2,
                ]
                c.add_polygon(polygon, layer=s_bend_support_layer)
    if ridge_s_bend_layer is not None:
        if ridge_s_bend_width is None or ridge_s_bend_width <= 0:
            raise ValueError(
                "ridge_s_bend_width must be positive when "
                "ridge_s_bend_layer is provided"
            )
        if s_bend_compact_180_radius is not None:
            c.add_ref(
                gf.path.extrude(
                    s_bend_centerline,
                    cross_section=partial(
                        gf.cross_section.cross_section,
                        width=ridge_s_bend_width,
                        layer=ridge_s_bend_layer,
                    ),
                )
            )
        else:
            ridge_ports = [
                port.copy(
                    width=ridge_s_bend_width,
                    layer=ridge_s_bend_layer,
                )
                for port in route_ports
            ]
            connect(
                ridge_ports[0],
                ridge_ports[1],
                method="circular",
                bend_width=ridge_s_bend_width,
                taper_length=s_bend_taper_length,
                cross_section=partial(
                    gf.cross_section.cross_section,
                    width=ridge_s_bend_width,
                    layer=ridge_s_bend_layer,
                ),
                add_to_component=c,
            )

    # Preserve the exact taper/S-bend interfaces for downstream mask clipping.
    # Each port orientation points from the taper into the center S-bend, so
    # the taper-side half-plane is the opposite direction.
    c.info["s_bend_connection_centers"] = [
        [float(value) for value in port.center]
        for port in ports
    ]
    c.info["s_bend_connection_orientations"] = [
        float(port.orientation)
        for port in ports
    ]
    c.info["taper_end_centers"] = [
        [float(value) for value in port.center]
        for port in taper_end_ports
    ]
    c.info["taper_end_orientations"] = [
        float(port.orientation)
        for port in taper_end_ports
    ]
    c.info["taper_end_widths"] = [
        float(port.width)
        for port in taper_end_ports
    ]

    # Put the mechanical tethers at the two outer extremes of the pattern,
    # directly on the facet-side ends of the left and right edge couplers.
    if tether_spec is not None:
        for port in outer_ports:
            tether_component = tether_spec()
            tether = c.add_ref(tether_component)
            # Keep this strict: the swept taper endpoint and the reshaped
            # tether face must have the same physical width.
            tether.connect("o1", port)
            tether_length = float(
                tether_component.info.get(
                    "tether_length",
                    tether_component.xsize,
                )
            )
            tether_slot_lengths.append(tether_length)
            tether_slot_transverse_extents.append(
                (float(tether_component.ymin), float(tether_component.ymax))
            )

    if extra_rectangle_size is not None:
        rectangle_length, rectangle_width = extra_rectangle_size
        if rectangle_length <= 0 or rectangle_width <= 0:
            raise ValueError("extra_rectangle_size values must be positive")
        if extra_rectangle_gap < 0:
            raise ValueError("extra_rectangle_gap must be non-negative")
        if extra_rectangle_margin < 0:
            raise ValueError("extra_rectangle_margin must be non-negative")
        if (
            extra_rectangle_surround_slot_width is not None
            and extra_rectangle_surround_slot_width <= 0
        ):
            raise ValueError(
                "extra_rectangle_surround_slot_width must be positive"
            )
        if tether_rectangle_slot_width is not None:
            if tether_spec is None:
                raise ValueError(
                    "tether_rectangle_slot_width requires tether_spec"
                )
            if tether_rectangle_slot_width <= 0:
                raise ValueError(
                    "tether_rectangle_slot_width must be positive"
                )
            if tether_rectangle_slot_overlap < 0:
                raise ValueError(
                    "tether_rectangle_slot_overlap must be non-negative"
                )
            if tether_rectangle_slot_clearance < 0:
                raise ValueError(
                    "tether_rectangle_slot_clearance must be non-negative"
                )
            if tether_rectangle_slot_start_offset < 0:
                raise ValueError(
                    "tether_rectangle_slot_start_offset must be non-negative"
                )
            if (
                tether_rectangle_slot_tether_margin is not None
                and tether_rectangle_slot_tether_margin < 0
            ):
                raise ValueError(
                    "tether_rectangle_slot_tether_margin must be non-negative"
                )
            if tether_rectangle_slot_transition_length <= 0:
                raise ValueError(
                    "tether_rectangle_slot_transition_length must be positive"
                )
            if tether_rectangle_slot_transition_points < 4:
                raise ValueError(
                    "tether_rectangle_slot_transition_points must be at least 4"
                )

        distance_to_center = extra_rectangle_gap + rectangle_length / 2
        for index, port in enumerate(taper_end_ports):
            angle = np.radians(port.orientation)
            outward = np.array((np.cos(angle), np.sin(angle)))
            center = (
                port.center[0] + distance_to_center * outward[0],
                port.center[1] + distance_to_center * outward[1],
            )
            if extra_rectangle_surround_slot_width is None:
                rectangle = c.add_ref(
                    gf.components.rectangle(
                        size=(rectangle_length, rectangle_width),
                        centered=True,
                        layer=extra_rectangle_layer,
                    )
                )
                rectangle.rotate(port.orientation)
                rectangle.move(center)
            else:
                slot_width = extra_rectangle_surround_slot_width
                local_bars = (
                    (
                        (rectangle_length + 2 * slot_width, slot_width),
                        (0.0, rectangle_width / 2 + slot_width / 2),
                    ),
                    (
                        (rectangle_length + 2 * slot_width, slot_width),
                        (0.0, -rectangle_width / 2 - slot_width / 2),
                    ),
                    (
                        (slot_width, rectangle_width),
                        (rectangle_length / 2 + slot_width / 2, 0.0),
                    ),
                    (
                        (slot_width, rectangle_width),
                        (-rectangle_length / 2 - slot_width / 2, 0.0),
                    ),
                )
                angle_degrees = float(port.orientation)
                angle_radians = np.radians(angle_degrees)
                rotation = np.array(
                    [
                        [np.cos(angle_radians), -np.sin(angle_radians)],
                        [np.sin(angle_radians), np.cos(angle_radians)],
                    ]
                )
                for bar_size, local_center in local_bars:
                    bar_center = (
                        np.asarray(center, dtype=float)
                        + rotation @ np.asarray(local_center, dtype=float)
                    )
                    bar = c.add_ref(
                        gf.components.rectangle(
                            size=bar_size,
                            centered=True,
                            layer=extra_rectangle_layer,
                        )
                    )
                    bar.rotate(angle_degrees)
                    bar.move(tuple(float(value) for value in bar_center))

            if extra_rectangle_fill_layer is not None:
                fill = c.add_ref(
                    gf.components.rectangle(
                        size=(rectangle_length, rectangle_width),
                        centered=True,
                        layer=extra_rectangle_fill_layer,
                    )
                )
                fill.rotate(port.orientation)
                fill.move(center)

            if extra_rectangle_margin > 0:
                margin_layer = (
                    extra_rectangle_layer
                    if extra_rectangle_margin_layer is None
                    else extra_rectangle_margin_layer
                )
                surrounding_rectangle = c.add_ref(
                    gf.components.rectangle(
                        size=(
                            rectangle_length + 2 * extra_rectangle_margin,
                            rectangle_width + 2 * extra_rectangle_margin,
                        ),
                        centered=True,
                        layer=margin_layer,
                    )
                )
                surrounding_rectangle.rotate(port.orientation)
                surrounding_rectangle.move(center)

            if tether_rectangle_slot_width is not None:
                tether_length = tether_slot_lengths[index]
                slot_start = (
                    tether_length
                    - tether_rectangle_slot_overlap
                    + tether_rectangle_slot_clearance
                    + tether_rectangle_slot_start_offset
                )
                slot_end = (
                    extra_rectangle_gap
                    + tether_rectangle_slot_overlap
                    - tether_rectangle_slot_clearance
                )
                slot_length = slot_end - slot_start
                if slot_length <= 0:
                    raise ValueError(
                        "tether overlaps the extra rectangle; no positive slot "
                        "length remains"
                    )

                port_center = np.asarray(port.center, dtype=float)
                if tether_rectangle_slot_tether_margin is None:
                    # Legacy constant-width slot.
                    slot_center = (
                        port_center
                        + (slot_start + slot_end) / 2 * outward
                    )
                    slot = c.add_ref(
                        gf.components.rectangle(
                            size=(slot_length, tether_rectangle_slot_width),
                            centered=True,
                            layer=tether_rectangle_slot_layer,
                        )
                    )
                    slot.rotate(port.orientation)
                    slot.move(tuple(float(value) for value in slot_center))
                else:
                    # Local tether +y after reference.connect(). This preserves
                    # the offset of asymmetric tether variants.
                    transverse = np.array(
                        (np.sin(angle), -np.cos(angle)),
                        dtype=float,
                    )
                    tether_ymin, tether_ymax = (
                        tether_slot_transverse_extents[index]
                    )
                    tether_ymin -= tether_rectangle_slot_tether_margin
                    tether_ymax += tether_rectangle_slot_tether_margin
                    half_width = tether_rectangle_slot_width / 2

                    # Flare only beside the tether window. Smoothstep makes
                    # the sidewalls tangent to both the broad tether opening
                    # and the constant-width slot, removing the circled kink.
                    transition_length = min(
                        tether_rectangle_slot_transition_length,
                        slot_length,
                    )
                    lower_curve = []
                    upper_curve = []
                    for u in np.linspace(
                        0.0,
                        1.0,
                        tether_rectangle_slot_transition_points,
                    ):
                        blend = u * u * (3.0 - 2.0 * u)
                        distance = slot_start + u * transition_length
                        lower = tether_ymin + blend * (
                            -half_width - tether_ymin
                        )
                        upper = tether_ymax + blend * (
                            half_width - tether_ymax
                        )
                        lower_curve.append(
                            port_center
                            + distance * outward
                            + lower * transverse
                        )
                        upper_curve.append(
                            port_center
                            + distance * outward
                            + upper * transverse
                        )

                    # Continue horizontally at the selected slot width until
                    # the support rectangle. Duplicate endpoints are harmless
                    # and are removed when the component is flattened.
                    end_lower = (
                        port_center
                        + slot_end * outward
                        - half_width * transverse
                    )
                    end_upper = (
                        port_center
                        + slot_end * outward
                        + half_width * transverse
                    )
                    slot_points = (
                        lower_curve
                        + [end_lower, end_upper]
                        + list(reversed(upper_curve))
                    )
                    c.add_polygon(
                        [
                            tuple(float(value) for value in point)
                            for point in slot_points
                        ],
                        layer=tether_rectangle_slot_layer,
                    )

    if flatten:
        c.flatten()

    return c
