"""
Collaborator logos, vendored as ``.gds`` files beside this module.

One ``get_*`` wrapper per file (`get_cal`, `get_mj`, `get_alp`, `get_harq`,
`get_ionq`); each recolors the imported polygons onto a chosen ``layer`` and
optionally rescales the logo to a target ``width``. A multi-layer logo takes a
``{input_layer: output_layer}`` dict instead of a single layer.
"""
from pathlib import Path

import gdsfactory as gf
import numpy as np
from kfactory import kdb

from layout.core.layers import Layers, layer as resolve_layer

_LOGO_DIR = Path(__file__).parent


def _import_logo(filename, width=None, layer=Layers.wg):
    """Import a logo from a GDS file onto ``layer``, optionally rescaling it.

    Parameters
    ----------
    filename : str or pathlib.Path
        GDS file (in this directory) to import.
    width : float, optional
        If given, scale the logo to this total width (microns).  The logo is centered
        on the origin either way.
    layer : LayerSpec or dict[LayerSpec, LayerSpec]
        Output layer for the placed polygons; a dict maps input layers to output
        layers, anything else maps every input layer to the one given.

    Returns
    -------
    gdsfactory.Component
        Component holding the logo polygons on ``layer``.
    """
    c = gf.Component()

    # A private layout, so the file's own cells never reach the layout being built.
    logo = kdb.Layout()
    logo.read(str(_LOGO_DIR / filename))
    top = logo.top_cell()

    def points(indexes):
        """Vertices, in microns, of every logo polygon on the given layer indexes."""
        return [
            np.array([(p.x, p.y) for p in
                      poly.to_simple_polygon().to_dtype(logo.dbu).each_point()])
            for index in indexes
            for poly in kdb.Region(top.begin_shapes_rec(index)).each()
        ]

    # Bounding-box center and scale factor for the whole logo.
    all_polygons = points(logo.layer_indexes())
    xs = np.concatenate([p[:, 0] for p in all_polygons])
    ys = np.concatenate([p[:, 1] for p in all_polygons])

    center = [(xs.min() + xs.max()) / 2, (ys.min() + ys.max()) / 2]

    if width is not None:
        scale = width / (xs.max() - xs.min())
    else:
        scale = 1.0

    if not isinstance(layer, dict):
        layer = {None: layer}

    for input_layer, output_layer in layer.items():
        if input_layer is None:
            indexes = logo.layer_indexes()
        else:
            specs = input_layer if isinstance(input_layer, list) else [input_layer]
            indexes = [logo.layer(kdb.LayerInfo(*resolve_layer(spec))) for spec in specs]

        polygons = [(p - center) * scale for p in points(indexes)]

        for p in polygons:
            c.add_polygon(p, layer=output_layer)

    c.add_label("__" + Path(filename).stem + "__", position=(0, 0), layer=(Layers.label[0], 1))

    return c


@gf.cell
def get_cal(width=100, layer=Layers.metal1):
    return _import_logo("cal.gds", layer=layer, width=width)

@gf.cell
def get_mj(width=None, layer=Layers.wg):
    return _import_logo("mj.gds", layer=layer, width=width)

@gf.cell
def get_alp(width=None, layer=Layers.wg):
    return _import_logo("alp.gds", layer=layer, width=width)

@gf.cell
def get_harq(width=None, layer=None):
    """HarQ logo; ``layer=None`` maps its two input layers to metal1 and wg."""
    if layer is None:
        layer = {(1,0):Layers.metal1, (2,0):Layers.wg}
    return _import_logo("harq.gds", layer=layer, width=width)

@gf.cell
def get_ionq(width=None, layer=None):
    """IonQ logo, single-layer; ``layer=None`` maps it onto wg."""
    if layer is None:
        layer = Layers.wg
    layer = {(1,0):layer}
    return _import_logo("ionq.gds", layer=layer, width=width)
