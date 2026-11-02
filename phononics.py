"""
Phononic crystal shielding for a suspended nanowire.

:func:`phononic_tether` is one unit cell of the crystal -- a beam carrying a chain of
square masses -- and :func:`phononic_shield` tiles it along a straight run or extrudes
it along a curve, leaving a band of tethers either side of whatever it suspends.
:func:`phononic_cap` is that band's terminator, disabled and drawing nothing, and
:func:`phononic_test` is the lithography test structure for the band. Everything here
is drawn on ``wg_doubleetch``.
"""
import inspect
import numpy as np

import gdsfactory as gf
from gdsfactory.cross_section import ComponentAlongPath

from layers import Layers



@gf.cell
def phononic_tether(
    width=0.25,
    segments=(0.626, 0.6, 0.56, 0.54, 0.516, 0.504),
    gaps=0.25,
    overshoot=1
):
    """One unit cell of the phononic crystal: a chain of square masses on a beam.

    Parameters
    ----------
    width : float
        Beam width, in microns.
    segments, gaps : float or list of float
        Mass side lengths and the beam runs between them, in microns; ``gaps`` is one longer
        than ``segments``, opening and closing the cell. Either may be a scalar and is
        broadcast against the other -- both scalar raises -- and either ``None`` leaves
        the bare beam.
    overshoot : float
        Beam drawn past each anchor, in microns, so neighboring cells overlap.

    Returns
    -------
    gf.Component
        Anchor ports ``a1`` and ``a2`` at the two ends of the cell.
    """
    c = gf.Component()


    if segments is None or gaps is None:
        x = width
    else:
        if np.isscalar(gaps) and np.isscalar(segments):
            raise ValueError("Segments and gaps cannot both be scalars.")

        if np.isscalar(gaps):
            gaps = [gaps]*(len(segments)+1)

        if np.isscalar(segments):
            segments = [segments]*(len(gaps)-1)

        x = gaps[0]
        for segment, gap in zip(segments, gaps[1:]):    # Currently hard-code square
            rect = gf.components.rectangle(size=(segment, segment), centered=True, layer=Layers.wg_doubleetch)
            c.add_ref(rect).move((x + segment/2, 0))
            x += segment + gap

    c.add_ref(gf.components.rectangle(size=(overshoot*2 + x, width), centered=True, layer=Layers.wg_doubleetch)).move((x/2, 0))
    c.add_port(name="a1", port_type="anchor", center=(0, 0), width=width, orientation=180, layer=Layers.wg_doubleetch)
    c.add_port(name="a2", port_type="anchor", center=(x, 0), width=width, orientation=0, layer=Layers.wg_doubleetch)

    c.flatten()

    return c

def get_phononic_tether_length(
    tether_spec=phononic_tether,
) -> float:
    """Anchor-to-anchor length of one tether, in microns."""
    tether = tether_spec()
    return tether.ports["a2"].center[0] - tether.ports["a1"].center[0]


@gf.cell
def phononic_shield(
    path=None,
    tether_spec=phononic_tether,
    pitch=1.,
    slab_width=15.,
    side="both",
):
    """Bands of phononic tethers either side of a path, suspending what runs between them.

    ``path`` picks the geometry: a number tiles the unit cell along a straight run of
    that length starting at the origin, a ``gf.Path`` extrudes the bands along it as a
    ``ComponentAlongPath``, and ``None`` is a 100 micron straight.

    Parameters
    ----------
    path : float or gf.Path, optional
        Straight length, in microns, or the curve to follow.
    pitch : float
        Tether-to-tether spacing along the path's centerline, in microns; opened out just
        enough that a whole number of tethers spans the run.
    slab_width : float
        Clear opening between the two bands, in microns; whatever runs down the middle is
        suspended across it.
    side : str
        ``"both"``, ``"left"`` or ``"right"``. Sides are named against the direction of
        travel, so ``"right"`` is starboard -- the way a positive ``gf.Section`` offset
        already points.
    """
    c = gf.Component()

    try:
        sides = {"both": [-1, 1], "right": [1], "left": [-1]}[side]
    except KeyError:
        raise ValueError(f"side must be 'both', 'left' or 'right', got {side!r}") from None

    tether = tether_spec()
    tether_length = get_phononic_tether_length(tether_spec=tether_spec)

    if path is None:
        path_is_linear = True
        path_length = 100
    elif isinstance(path, (int, float)):
        path_is_linear = True
        path_length = path
    else:
        path_is_linear = False
        path_length = path.length()

    if path_length < pitch:
        raise ValueError(
            f"phononic_shield: path_length={path_length:g} um is shorter than "
            f"pitch={pitch:g} um, so not one tether fits along it"
        )

    if path_is_linear:
        n = int(np.floor(path_length / pitch))
        pitch_ = path_length / n

        subtract_from = gf.components.rectangle(size=(pitch_, tether_length), centered=False, layer=Layers.wg_doubleetch)

        subtract_with = gf.Component(c.name + "_subtract_with")
        subtract_with.add_ref(tether).rotate(90)
        subtract_with.add_ref(tether).rotate(90).move((pitch_, 0))
    else:
        offset = (slab_width + tether_length) / 2
        layer = Layers.wg_doubleetch
        s0 = gf.Section(width=1, offset=0, layer=layer, port_names=("o1", "o2"))

        cs_from = gf.CrossSection(
            sections=tuple(
                gf.Section(width=tether_length, offset=s*offset, layer=layer,
                           **({"port_names": ("o1", "o2")} if i == 0 else {}))
                for i, s in enumerate(sides)
            ),
        )
        subtract_from = gf.path.extrude(path, cross_section=cs_from)

        # Round the pitch down to fit a whole number of tethers, the last one included.
        pitch_ = path_length / np.floor(path_length / pitch) - 1e-6

        tethers = gf.Component(c.name + "_tethers")
        for s in sides:
            # Negated: the s section above sits at -s in y, else the boolean leaves a solid bar.
            tethers.add_ref(tether).rotate(-90*s).move((0, -s*slab_width/2))

        cs_with = gf.CrossSection(
            sections=(s0,),
            components_along_path=(ComponentAlongPath(component=tethers, spacing=pitch_),),
        )
        subtract_with = gf.path.extrude(path, cross_section=cs_with)

    # Cut the tethers out of the slab.
    void = gf.boolean(
        A=subtract_from,
        B=subtract_with,
        operation="not",
        layer=Layers.wg_doubleetch
    )

    if path_is_linear:
        unit_cell = gf.Component(c.name + "_unit_cell")
        if 1 in sides:
            unit_cell.add_ref(void).move((0, slab_width/2)).mirror_y()
        if -1 in sides:
            unit_cell.add_ref(void).move((0, slab_width/2))

        for i in range(n):
            c.add_ref(unit_cell).move((i*pitch_, 0))
    else:
        c.add_ref(void)

    return c

@gf.cell
def phononic_cap(
    tether_spec=phononic_tether,
    slab_width=15.,
    inset=-0.5,
    outset=1,
    offset=2,
):
    """Bracket closing off the end of a :func:`phononic_shield`.

    A slab from ``inset`` to ``outset + offset`` across both tether bands, joined over
    the ``slab_width`` between them by a spine from ``outset`` to ``outset + offset``,
    so it caps the bands clear of whatever that opening suspends. A negative ``inset``
    reaches back over the last tether and merges with the shield.

    Warnings
    --------
    The polygon is disabled, so this draws nothing, and ``inset`` is overridden before use.
    """
    tether_length = get_phononic_tether_length(tether_spec=tether_spec)

    c = gf.Component()

    inset = -0.5

    if False:
        c.add_polygon(
            [
                (outset+offset, -tether_length-slab_width/2),
                (inset, -tether_length-slab_width/2),
                (inset, -slab_width/2),
                (outset, -slab_width/2),
                (outset, slab_width/2),
                (inset, slab_width/2),
                (inset, tether_length+slab_width/2),
                (outset+offset, tether_length+slab_width/2),
            ],
            layer=Layers.wg_doubleetch
        )

    return c


@gf.cell
def phononic_test(
    length=1000,
    phononic_shield_spec=phononic_shield,
):
    """Lithography test structure: a straight phononic shield with a cap at each end.

    :func:`phononic_cap` draws nothing, so what comes out is the bare shield.
    """
    c = gf.Component()

    c.add_ref(phononic_shield_spec(length)).move((-length/2, 0))

    tether_spec = inspect.signature(phononic_shield_spec).parameters["tether_spec"].default

    c.add_ref(phononic_cap(tether_spec=tether_spec)).rotate(180).move((-length/2, 0))
    c.add_ref(phononic_cap(tether_spec=tether_spec)).move((length/2, 0))

    return c