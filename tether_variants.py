"""Reusable gdsfactory components for the four tethers in tether_variants.gds.

Every source polygon is translated so its native 200 nm endface is centered at
(0, 0). The end of the imported tether itself is smoothly reshaped to the
requested connection width; no extra taper is appended outside the source
geometry. Port ``o1`` points right, allowing direct attachment with
``reference.connect``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import gdsfactory as gf

REPOSITORY_ROOT = Path(__file__).resolve().parent
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from extract_layer2_gdsfactory import read_layer_polygons


DEFAULT_SOURCE_GDS = REPOSITORY_ROOT / "tether_variants.gds"
DEFAULT_LAYER = (2, 0)
PORT_WIDTH = 0.2
DEFAULT_CONNECTION_TAPER_LENGTH = 10.0
MAX_CONNECTION_BLEND_FRACTION = 0.45


def _tether_variant(
    variant: int,
    layer: tuple[int, int] = DEFAULT_LAYER,
    source_gds: str | Path = DEFAULT_SOURCE_GDS,
    length_scale: float = 1.0,
    width_scale: float = 1.0,
    pad_scale: float | None = None,
    connection_width: float = PORT_WIDTH,
    connection_taper_length: float = DEFAULT_CONNECTION_TAPER_LENGTH,
) -> gf.Component:
    """Build one exact tether polygon from the source GDS.

    Args:
        variant: Tether number in the source file, from 1 through 4.
        layer: Output ``(layer, datatype)``.
        source_gds: Source file containing the four Layer-2 polygons.
        length_scale: Legacy pad scale. Must equal ``width_scale``.
        width_scale: Legacy pad scale. Must equal ``length_scale``.
        pad_scale: Uniform scale applied rigidly to the far pad in both axes.
        connection_width: Width presented to the adjoining optical taper.
        connection_taper_length: Maximum distance, measured inward from the
            connection face, over which a narrower connection is blended back
            to the native tether shape. The actual distance is capped at 45
            percent of the tether length so compact variants retain their pad
            geometry. Wider connections instead widen the native tether so it
            never falls below ``connection_width``.
    """
    if variant not in (1, 2, 3, 4):
        raise ValueError(f"variant must be 1, 2, 3, or 4; got {variant}")
    if pad_scale is None:
        if abs(length_scale - width_scale) > 1e-12:
            raise ValueError(
                "Use one uniform pad_scale; independent length/width scaling "
                "distorts tilted tether pads"
            )
        pad_scale = length_scale
    if pad_scale <= 0:
        raise ValueError("pad_scale must be positive")
    if connection_width <= 0:
        raise ValueError("connection_width must be positive")
    if connection_taper_length < 0:
        raise ValueError("connection_taper_length must be non-negative")
    if connection_width < PORT_WIDTH and connection_taper_length == 0:
        raise ValueError(
            "connection_taper_length must be positive for a narrower connection"
        )

    _cell_name, polygons = read_layer_polygons(source_gds, layer_number=2)
    if len(polygons) != 4:
        raise ValueError(
            f"Expected four Layer-2 tether polygons in {source_gds}, found {len(polygons)}"
        )

    _datatype, source_points = polygons[variant - 1]

    # The connection endface is the rightmost vertical edge. In every source
    # variant it is exactly 0.2 um wide.
    x_port = max(x for x, _y in source_points)
    endface_y = [y for x, y in source_points if x == x_port]
    if len(endface_y) != 2:
        raise ValueError(
            f"Variant {variant} does not have the expected two-point right endface"
        )
    y_port = (min(endface_y) + max(endface_y)) / 2
    measured_width = max(endface_y) - min(endface_y)
    if abs(measured_width - PORT_WIDTH) > 1e-6:
        raise ValueError(
            f"Variant {variant} endface width is {measured_width:g} um, "
            f"expected {PORT_WIDTH:g} um"
        )

    # Scale the far pad rigidly about the pad's own center so tilted features
    # retain their angle and aspect ratio. Only the pad transition and the
    # short connection-side neck are interpolated.
    tether_length = x_port - min(x for x, _y in source_points)
    if tether_length <= 0:
        raise ValueError(f"Variant {variant} has zero or negative tether length")

    transition_start = 0.45
    pad_start = 0.7

    def normalized_distance(x: float) -> float:
        return (x_port - x) / tether_length

    pad_points = [
        (x, y)
        for x, y in source_points
        if normalized_distance(x) >= pad_start
    ]
    if not pad_points:
        raise ValueError(f"Could not identify the far pad of variant {variant}")

    # Bounding-box center is stable for asymmetric and forked pad polygons.
    pad_center_x = (
        min(x for x, _y in pad_points) + max(x for x, _y in pad_points)
    ) / 2
    pad_center_y = (
        min(y for _x, y in pad_points) + max(y for _x, y in pad_points)
    ) / 2

    blend_length = min(
        connection_taper_length,
        MAX_CONNECTION_BLEND_FRACTION * tether_length,
    )
    connection_scale = connection_width / PORT_WIDTH

    def reshape_point(x: float, y: float) -> tuple[float, float]:
        # Reshape the source polygon itself instead of placing a second taper
        # beyond it. At the port, both endface vertices therefore have the
        # exact requested width. Smoothstep gives zero slope at both ends of
        # the blend and rejoins the untouched source outline tangentially.
        distance_from_port = x_port - x
        if (
            connection_width < PORT_WIDTH
            and blend_length > 0
            and distance_from_port < blend_length
        ):
            u = max(0.0, min(1.0, distance_from_port / blend_length))
            neck_weight = u * u * (3.0 - 2.0 * u)
            local_scale = connection_scale + (1.0 - connection_scale) * neck_weight
            y = y_port + local_scale * (y - y_port)

        u = normalized_distance(x)
        t = (u - transition_start) / (pad_start - transition_start)
        t = max(0.0, min(1.0, t))
        weight = t * t * (3.0 - 2.0 * t)  # smoothstep

        # For weight=1 this is a true uniform affine scale about the pad
        # center. For weight=0 the original stem and port remain exact.
        x_scaled = x + weight * (pad_scale - 1.0) * (x - pad_center_x)
        y_scaled = y + weight * (pad_scale - 1.0) * (y - pad_center_y)
        return (
            x_scaled - x_port,
            y_scaled - y_port,
        )

    points = [reshape_point(x, y) for x, y in source_points]

    component = gf.Component()
    component.add_polygon(points, layer=layer)

    if connection_width > PORT_WIDTH:
        # A requested connection wider than the native 0.2 um neck must not
        # taper back down to 0.2 um inside the tether. Dilating the imported
        # tether by half the width difference raises every native minimum to
        # the requested width while preserving its original centerline. Crop
        # the right side back to x=0 so the connection face stays fixed and
        # remains exactly ``connection_width`` wide.
        offset = (connection_width - PORT_WIDTH) / 2
        component.offset(layer=layer, distance=offset, corner_mode=2)

        x_min = min(x for x, _y in points) - offset - 1.0
        y_min = min(y for _x, y in points) - offset - 1.0
        y_max = max(y for _x, y in points) + offset + 1.0
        clip = gf.Component()
        clip.add_polygon(
            [(x_min, y_min), (0.0, y_min), (0.0, y_max), (x_min, y_max)],
            layer=layer,
        )
        cropped = gf.boolean(
            component,
            clip,
            operation="and",
            layer=layer,
            layer1=layer,
            layer2=layer,
        )
        component = gf.Component()
        component.add_ref(cropped)
        component.flatten()

    component.add_port(
        name="o1",
        center=(0.0, 0.0),
        width=connection_width,
        orientation=0,
        layer=layer,
        port_type="optical",
    )
    component.info["native_connection_width"] = PORT_WIDTH
    component.info["connection_width"] = float(connection_width)
    component.info["connection_blend_length"] = float(
        blend_length if connection_width < PORT_WIDTH else 0.0
    )
    # Use the final drawn extent, including any dilation for connection widths
    # above the native 0.2 um. Slots can then terminate at the exact tether edge.
    component.info["tether_length"] = float(component.xsize)
    return component


@gf.cell
def tether_1(
    layer: tuple[int, int] = DEFAULT_LAYER,
    source_gds: str | Path = DEFAULT_SOURCE_GDS,
    length_scale: float = 1.0,
    width_scale: float = 1.0,
    pad_scale: float | None = None,
    connection_width: float = PORT_WIDTH,
    connection_taper_length: float = DEFAULT_CONNECTION_TAPER_LENGTH,
) -> gf.Component:
    """Variant 1: compact, vertically symmetric tether."""
    return _tether_variant(
        1, layer, source_gds, length_scale, width_scale, pad_scale,
        connection_width, connection_taper_length,
    )


@gf.cell
def tether_2(
    layer: tuple[int, int] = DEFAULT_LAYER,
    source_gds: str | Path = DEFAULT_SOURCE_GDS,
    length_scale: float = 1.0,
    width_scale: float = 1.0,
    pad_scale: float | None = None,
    connection_width: float = PORT_WIDTH,
    connection_taper_length: float = DEFAULT_CONNECTION_TAPER_LENGTH,
) -> gf.Component:
    """Variant 2: split/forked tether."""
    return _tether_variant(
        2, layer, source_gds, length_scale, width_scale, pad_scale,
        connection_width, connection_taper_length,
    )


@gf.cell
def tether_3(
    layer: tuple[int, int] = DEFAULT_LAYER,
    source_gds: str | Path = DEFAULT_SOURCE_GDS,
    length_scale: float = 1.0,
    width_scale: float = 1.0,
    pad_scale: float | None = None,
    connection_width: float = PORT_WIDTH,
    connection_taper_length: float = DEFAULT_CONNECTION_TAPER_LENGTH,
) -> gf.Component:
    """Variant 3: broad curved single-arm tether."""
    return _tether_variant(
        3, layer, source_gds, length_scale, width_scale, pad_scale,
        connection_width, connection_taper_length,
    )


@gf.cell
def tether_4(
    layer: tuple[int, int] = DEFAULT_LAYER,
    source_gds: str | Path = DEFAULT_SOURCE_GDS,
    length_scale: float = 1.0,
    width_scale: float = 1.0,
    pad_scale: float | None = None,
    connection_width: float = PORT_WIDTH,
    connection_taper_length: float = DEFAULT_CONNECTION_TAPER_LENGTH,
) -> gf.Component:
    """Variant 4: narrow curved single-arm tether."""
    return _tether_variant(
        4, layer, source_gds, length_scale, width_scale, pad_scale,
        connection_width, connection_taper_length,
    )


TETHER_FACTORIES = {
    1: tether_1,
    2: tether_2,
    3: tether_3,
    4: tether_4,
}


def add_tether_connected(
    parent: gf.Component,
    current_port: gf.Port,
    variant: int = 1,
    layer: tuple[int, int] = DEFAULT_LAYER,
    source_gds: str | Path = DEFAULT_SOURCE_GDS,
    length_scale: float = 1.0,
    width_scale: float = 1.0,
    pad_scale: float | None = None,
    connection_width: float = PORT_WIDTH,
    connection_taper_length: float = DEFAULT_CONNECTION_TAPER_LENGTH,
):
    """Add a selected tether to ``parent`` and connect its ``o1`` port.

    Returns the new ComponentReference so the caller can move, inspect, or
    retain it. The connection automatically handles translation and rotation.
    """
    try:
        factory = TETHER_FACTORIES[variant]
    except KeyError as error:
        raise ValueError(f"variant must be 1, 2, 3, or 4; got {variant}") from error

    reference = parent.add_ref(
        factory(
            layer=layer,
            source_gds=source_gds,
            length_scale=length_scale,
            width_scale=width_scale,
            pad_scale=pad_scale,
            connection_width=connection_width,
            connection_taper_length=connection_taper_length,
        )
    )
    reference.connect("o1", current_port)
    return reference


__all__ = [
    "tether_1",
    "tether_2",
    "tether_3",
    "tether_4",
    "add_tether_connected",
    "TETHER_FACTORIES",
]
