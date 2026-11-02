"""
The group's GDS layer map and helpers built on it.

:class:`Layers` names every process layer as a ``(gdslayer, gdspurpose)`` pair and joins the
generic PDK's layer map at import time; a member whose GDS number a generic layer already
carries joins as an alias of it, so the generic name is what ``gf.get_layer_name`` reports.
:func:`layer` resolves a name, member or tuple to a plain tuple, and :func:`generate_lyp`
writes a KLayout ``.lyp`` file of layer colors and names.
"""
from aenum import extend_enum

import gdsfactory as gf
from gdsfactory.typings import Layer
from gdsfactory.technology.layer_views import LayerViews, LayerView
from gdsfactory.technology.layer_map import LayerMap

class Layers(LayerMap):
    """The group's process layers, as ``(gdslayer, gdspurpose)`` pairs.

    Six families, by leading GDS number: chip outline, labels and alignment
    markers (0-9); waveguide-etch features (10-19); metal and heater
    electrodes (20-29); poling electrodes and their windows (31-39); the
    adaptive-poling set (41-49); and packaging features (50-59). Every
    member has purpose 0; non-zero purposes on the same number carry text
    labels.

    Resolve a name, member or tuple to a plain tuple with :func:`layer`.
    """
    # Chip outlines, labels, and markers.
    outline:        Layer = (0, 0)
    label:          Layer = (1, 0) # Fast e-beam layer
    globalmarker:   Layer = (2, 0) # Pre-written markers
    patternmarker:  Layer = (3, 0)
    aligncheck:     Layer = (4, 0)
    dieetch:        Layer = (5, 0)
    title:          Layer = (6, 0)

    # Waveguide related.
    wg:             Layer = (10, 0) # Slow e-beam layer
    wg_doubleetch:  Layer = (11, 0) # Double-etched region
    grating:        Layer = (12, 0)
    ring:           Layer = (13, 0)
    via:            Layer = (14, 0) # Oxide etch
    mc:             Layer = (15, 0)
    mc2:            Layer = (16, 0)
    shield:         Layer = (17, 0)

    # Metal related.
    metal1:         Layer = (20, 0)
    metal2:         Layer = (21, 0)
    heater:         Layer = (22, 0)
    metal_remove:   Layer = (23, 0)

    # Poling related.
    polemarker:     Layer = (31, 0)
    poleelectrode:  Layer = (32, 0)
    polewindow:     Layer = (33, 0)
    poleetchwindow: Layer = (34, 0)
    poleelectrodetop:   Layer = (35, 0)

    # Adaptive poling variants.
    apolemarker:    Layer = (41, 0)
    apoleelectrode: Layer = (42, 0)
    apolewindow:    Layer = (43, 0)

    # Packaging-related.
    packaging:      Layer = (50, 0)
    fiber:          Layer = (51, 0)


# PDK activation prunes kcl layers the PDK does not own, so join the generic layer map.
for _member in Layers:
    if _member.name not in gf.gpdk.PDK.layers.__members__:
        extend_enum(gf.gpdk.PDK.layers, _member.name, _member.layer, _member.datatype)


def layer(name):
    """Resolve a semantic layer name to its ``(gdslayer, gdspurpose)`` tuple.

    A tuple or list is passed through as a tuple, so callers may also supply an explicit layer.

    Parameters
    ----------
    name : str or Layers or tuple of int
        A ``Layers`` attribute name (e.g. ``"wg"``), a member, or a ``(layer, purpose)`` tuple.

    Returns
    -------
    tuple of int

    Raises
    ------
    KeyError
        If ``name`` names no attribute of :class:`Layers`.
    AttributeError
        If it names one that is not a layer, such as an inherited ``Enum`` method.
    """
    # Always a tuple, never the member: a member's int is a klayout layer *index*.
    if isinstance(name, Layers):
        return (int(name.layer), int(name.datatype))
    if isinstance(name, (tuple, list)):
        return tuple(name)
    try:
        member = getattr(Layers, name)
    except (AttributeError, TypeError):
        raise KeyError("Unknown layer name {!r}; not an attribute of Layers.".format(name)) from None
    return (int(member.layer), int(member.datatype))


def generate_lyp(filename="layer_properties.lyp"):
    """Write a KLayout .lyp file of layer colors and names (File > Load Layer Properties)."""
    layer_views = LayerViews()

    def add_layer(lv, n, lay, c, d):
        """Add one named, colored, described view of ``lay`` to the ``lv`` LayerViews."""
        lv.add_layer_view(
            name=n,
            layer=lay,
            layer_view=LayerView(
                layer=lay,
                color=c,
                visible=True,
                name=n,
                description=d
            )
        )

    # Label layers
    add_layer(layer_views, "outline", Layers.outline, "White", "Chip outline")
    add_layer(layer_views, "label", Layers.label, "Pink", "General label")
    add_layer(layer_views, "globalmarker", Layers.globalmarker, "LightPink", "Global marker")
    add_layer(layer_views, "patternmarker", Layers.patternmarker, "HotPink", "EBL marker")
    add_layer(layer_views, "aligncheck", Layers.aligncheck, "DeepPink", "alignment")
    add_layer(layer_views, "dieetch", Layers.dieetch, "PaleVioletRed", "deep etch")
    add_layer(layer_views, "title", Layers.title, "MediumVioletRed", "General label")
    # WG layers
    add_layer(layer_views, "wg", Layers.wg, "CadetBlue", "LN wg")
    add_layer(layer_views, "wg_doubleetch", Layers.wg_doubleetch, "SteelBlue", "LN wg deepetched")
    add_layer(layer_views, "grating", Layers.grating, "LightSteelBlue", "bragg_grating")
    add_layer(layer_views, "ring", Layers.ring, "LightBlue", "LN ring")
    add_layer(layer_views, "via", Layers.via, "PowderBlue", "hole")
    add_layer(layer_views, "mc", Layers.mc, "LightSkyBlue", "mode converter")
    add_layer(layer_views, "mc2", Layers.mc2, "SkyBlue", "mode converter")
    add_layer(layer_views, "shield", Layers.shield, "CornflowerBlue", "shield")
    # Metal layers
    add_layer(layer_views, "metal1", Layers.metal1, "gold", "Gold electrode")
    add_layer(layer_views, "metal2", Layers.metal2, "goldenrod", "2nd Gold electrode")
    add_layer(layer_views, "heater", Layers.heater, "#667FCC", "heater")
    add_layer(layer_views, "metal_remove", Layers.metal_remove, "khaki", "Gold removal")
    # Poling Layers
    add_layer(layer_views, "polemarker", Layers.polemarker, "Cornsilk", "polemarker")
    add_layer(layer_views, "poleelectrode", Layers.poleelectrode, "BlanchedAlmond", "Poling electrode")
    add_layer(layer_views, "polewindow", Layers.polewindow, "Bisque", "Windows for Poling electrode")
    add_layer(layer_views, "poleetchwindow", Layers.poleetchwindow, "NavajoWhite", "Etching windows for Poling electrode")
    add_layer(layer_views, "poleelectrodetop", Layers.poleelectrodetop, "Wheat", "Top poling electrode (2-layer)")
    add_layer(layer_views, "apolemarker", Layers.apolemarker, "Moccasin", "Adaptive-poling marker")
    add_layer(layer_views, "apoleelectrode", Layers.apoleelectrode, "PeachPuff", "Adaptive-poling electrode")
    add_layer(layer_views, "apolewindow", Layers.apolewindow, "PapayaWhip", "Windows for adaptive-poling electrode")

    layer_views.to_lyp(filename)
