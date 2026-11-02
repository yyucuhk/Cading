r"""
Custom geometric text font translated from old C++ code.

Uppercase and lowercase letters share the same glyph. Subscript and superscript
are written inside backtick-delimited LaTeX regions (see ``text``); outside them
``_ ^ { }`` are literal. Supports a-z, A-Z, 0-9, and common punctuation/symbols:
  + - = * / \ ^ ~ _ . , : ; ' " ! ? | ( ) [ ] { } < > @ # $ %
The @ # $ % characters are rendered as directional arrows (right, left, up,
down triangles). The backtick toggles the LaTeX region and is consumed there,
so it draws only where it lands inside a sub/superscript target.
"""

from __future__ import annotations

__all__ = ["text"]

import math
import warnings

import gdsfactory as gf
from gdsfactory.component import Component
from gdsfactory.typings import Coordinate, LayerSpec
from kfactory import kdb
from shapely.geometry import LineString, MultiPolygon, Polygon
from shapely.ops import unary_union

from layout.core.layers import Layers

_WIDTH_RATIO = 0.8   # character width / height
_THICK_RATIO = _WIDTH_RATIO / 5   # stroke thickness / height
_SCRIPT_SCALE = 0.5   # sub/superscript size relative to full glyph height


def _rect(x1: float, y1: float, x2: float, y2: float) -> Polygon:
    return Polygon([(x1, y1), (x2, y1), (x2, y2), (x1, y2)])


def _thicken(pts: list[tuple[float, float]], thick: float) -> Polygon:
    """Stroke a polyline to produce a filled polygon."""
    if len(pts) < 2:
        return Polygon()
    return LineString(pts).buffer(
        thick / 2, cap_style="flat", join_style="mitre", mitre_limit=10.0
    )


def _add_geom(component: Component, geom, x_off: float, y_off: float, layer_idx: int, dbu: float, x_scale :float=1, y_scale :float=1) -> None:
    """Add a Shapely geometry to a kfactory component, preserving interior holes."""
    if geom is None or geom.is_empty:
        return
    if isinstance(geom, MultiPolygon) or hasattr(geom, "geoms"):
        for g in geom.geoms:
            _add_geom(component, g, x_off, y_off, layer_idx, dbu, x_scale=x_scale, y_scale=y_scale)
        return
    if not isinstance(geom, Polygon):
        return

    def to_ipts(coords):
        return [kdb.Point(round((x * x_scale + x_off) / dbu), round((y * y_scale + y_off) / dbu)) for x, y in coords]

    kpoly = kdb.Polygon(to_ipts(geom.exterior.coords[:-1]))
    for interior in geom.interiors:
        kpoly.insert_hole(to_ipts(interior.coords[:-1]))
    component.shapes(layer_idx).insert(kpoly)


def _char_polys(
    c: str, thick: float, height: float, width: float, cut: float | bool = True,
) -> list[Polygon]:
    """Shapely polygons for a single character at unit position."""
    bot = thick / 2
    mid = height / 2
    top = height - thick / 2
    left = thick / 2
    right = width - thick / 2
    chamfer = width / 6

    if cut is True:
        cut = thick / 2

    geoms: list[Polygon] = []
    cu = c.upper()

    if cu == "A":
        p = [(left, 0), (left, top - chamfer), (left + chamfer, top), (right, top), (right, 0)]
        geoms.append(_thicken(p, thick))
        geoms.append(_rect(thick + cut, mid - bot, width - thick, mid + bot))

    elif cu == "B":
        p = [
            (left, 0), (left, top),
            (right - 1.5 * chamfer, top), (right - 0.5 * chamfer, top - chamfer),
            (right - 0.5 * chamfer, mid), (right, mid),
            (right, bot), (thick+cut, bot),
        ]
        geoms.append(_thicken(p, thick))
        geoms.append(_rect(thick+cut, mid - bot, width - thick - 0.5 * chamfer, mid + bot))

    elif cu == "C":
        p = [
            (width, bot), (left + chamfer, bot), (left, bot + chamfer),
            (left, top - chamfer), (left + chamfer, top), (width, top),
        ]
        geoms.append(_thicken(p, thick))

    elif cu == "D":
        p = [
            (left, mid-cut-thick), (left, bot), (right - chamfer, bot), (right, bot + chamfer),
            (right, top - chamfer), (right - chamfer, top), (left, top), (left, mid-thick),
        ]
        geoms.append(_thicken(p, thick))

    elif cu == "E":
        p = [
            (width, bot), (left + chamfer, bot), (left, bot + chamfer),
            (left, top - chamfer), (left + chamfer, top), (width, top),
        ]
        geoms.append(_thicken(p, thick))
        geoms.append(_rect(thick, mid - bot, width - thick, mid + bot))

    elif cu == "F":
        p = [(left, 0), (left, top - chamfer), (left + chamfer, top), (width, top)]
        geoms.append(_thicken(p, thick))
        geoms.append(_rect(thick, mid - bot, width - thick, mid + bot))

    elif cu == "G":
        p = [
            (right - 1.5 * thick, mid), (right, mid), (right, bot),
            (left + chamfer, bot), (left, bot + chamfer),
            (left, top - chamfer), (left + chamfer, top), (width, top),
        ]
        geoms.append(_thicken(p, thick))

    elif cu == "H":
        geoms.append(_rect(0, 0, thick, height))
        geoms.append(_rect(width - thick, 0, width, height))
        geoms.append(_rect(thick, mid - bot, width - thick, mid + bot))

    elif cu == "I":
        geoms.append(_rect(0, height - thick, width, height))
        geoms.append(_rect(0, 0, width, thick))
        geoms.append(_rect(width / 2 - thick / 2, thick, width / 2 + thick / 2, height - thick))

    elif cu == "J":
        p = [
            (2 * thick, top), (right, top), (right, bot + chamfer),
            (right - chamfer, bot), (left, bot), (left, 2 * thick),
        ]
        geoms.append(_thicken(p, thick))

    elif cu == "K":
        geoms.append(_rect(0, 0, thick, height))
        geoms.append(_rect(thick, mid - bot, width - thick, mid + bot))
        geoms.append(_rect(width - thick, 0, width, mid + thick / 2))
        geoms.append(_rect(width - 2 * thick, mid + thick / 2, width - thick, height))

    elif cu == "L":
        p = [(left, height), (left, bot + chamfer), (left + chamfer, bot), (width, bot)]
        geoms.append(_thicken(p, thick))

    elif cu == "M":
        p = [(left, 0), (left, top), (right - chamfer, top), (right, top - chamfer), (right, 0)]
        geoms.append(_thicken(p, thick))
        geoms.append(_rect(width / 2 - thick / 2, thick, width / 2 + thick / 2, height - thick))

    elif cu == "N":
        p = [(left, 0), (left, top), (right - chamfer, top), (right, top - chamfer), (right, 0)]
        geoms.append(_thicken(p, thick))

    elif cu == "O":
        p = [
            (left, mid-thick), (left, top - chamfer), (left + chamfer, top),
            (right - chamfer, top), (right, top - chamfer),
            (right, bot + chamfer), (right - chamfer, bot),
            (left + chamfer, bot), (left, bot + chamfer), (left, mid-cut-thick),
        ]
        geoms.append(_thicken(p, thick))

    elif cu == "P":
        p = [
            (thick+cut, mid), (right - chamfer, mid), (right, mid + chamfer),
            (right, top - chamfer), (right - chamfer, top), (left, top), (left, 0),
        ]
        geoms.append(_thicken(p, thick))

    elif cu == "Q":
        p = [
            (width - thick - cut, 3 * thick / 2),
            (left + chamfer, 3 * thick / 2), (left, 3 * thick / 2 + chamfer),
            (left, top - chamfer), (left + chamfer, top),
            (right - chamfer, top), (right, top - chamfer), (right, thick),
        ]
        geoms.append(_thicken(p, thick))
        geoms.append(_rect(width - 2 * thick - cut, 0, width - thick - cut, thick))

    elif cu == "R":
        p = [
            (thick + cut, mid), (right, mid), (right, top - chamfer),
            (right - chamfer, top), (left, top), (left, 0),
        ]
        geoms.append(_thicken(p, thick))
        geoms.append(_rect(width - 2 * thick, 0, width - thick, mid - thick / 2))

    elif cu == "S":
        p = [
            (0, bot), (right - chamfer, bot), (right, bot + chamfer),
            (right, mid - chamfer), (right - chamfer, mid),
            (left + chamfer, mid), (left, mid + chamfer),
            (left, top - chamfer), (left + chamfer, top), (width, top),
        ]
        geoms.append(_thicken(p, thick))

    elif cu == "T":
        geoms.append(_rect(width / 2 - thick / 2, 0, width / 2 + thick / 2, height - thick))
        geoms.append(_rect(0, height - thick, width, height))

    elif cu == "U":
        p = [(right, height), (right, bot), (left + chamfer, bot), (left, bot + chamfer), (left, height)]
        geoms.append(_thicken(p, thick))

    elif cu == "V":
        p = [
            (right, height), (right, bot + width - 2 * thick),
            (left + thick, bot), (left, bot), (left, height),
        ]
        geoms.append(_thicken(p, thick))

    elif cu == "W":
        p = [(right, height), (right, bot + chamfer), (right - chamfer, bot), (left, bot), (left, height)]
        geoms.append(_thicken(p, thick))
        geoms.append(_rect(width / 2 - thick / 2, thick, width / 2 + thick / 2, height - thick))

    elif cu == "X":
        p1 = [
            (left, 0), (left, mid - width / 2 + thick / 2),
            (right, mid + width / 2 - thick / 2), (right, height),
        ]
        p2 = [
            (right, 0), (right, mid - width / 2 + thick / 2),
            (left, mid + width / 2 - thick / 2), (left, height),
        ]
        geoms.append(unary_union([_thicken(p1, thick), _thicken(p2, thick)]))

    elif cu == "Y":
        p1 = [(right, height), (right, bot + chamfer), (right - chamfer, bot), (thick, bot)]
        geoms.append(_thicken(p1, thick))
        p2 = [(left, height), (left, mid + chamfer), (left + chamfer, mid), (width - thick, mid)]
        geoms.append(_thicken(p2, thick))

    elif cu == "Z":
        p = [
            (width, bot), (left, bot),
            (left, mid - width / 2 + thick / 2),
            (right, mid + width / 2 - thick / 2),
            (right, top), (0, top),
        ]
        geoms.append(_thicken(p, thick))

    elif c == "0":
        p = [
            (left, mid-thick), (left, top - chamfer), (left + chamfer, top),
            (right - chamfer, top), (right, top - chamfer),
            (right, bot + chamfer), (right - chamfer, bot),
            (left + chamfer, bot), (left, bot + chamfer), (left, mid-cut-thick),
        ]
        geoms.append(_thicken(p, thick))
        geoms.append(_rect(width / 2 - thick / 2, mid - thick / 2, width / 2 + thick / 2, mid + thick / 2))

    elif c == "1":
        geoms.append(_rect(width / 2 - thick / 2, thick, width / 2 + thick / 2, height))
        geoms.append(_rect(0.5 * thick, height - 2 * thick, width / 2 - thick / 2, height - thick))
        geoms.append(_rect(0, 0, width, thick))

    elif c == "2":
        p = [
            (width, bot), (left, bot), (left, mid - chamfer),
            (left + chamfer, mid), (right - chamfer, mid), (right, mid + chamfer),
            (right, top - chamfer), (right - chamfer, top), (0, top),
        ]
        geoms.append(_thicken(p, thick))

    elif c == "3":
        p = [
            (0, bot), (right - chamfer, bot), (right, bot + chamfer),
            (right, top - chamfer), (right - chamfer, top), (0, top),
        ]
        geoms.append(_thicken(p, thick))
        geoms.append(_rect(thick, mid - bot, width - thick, mid + bot))

    elif c == "4":
        geoms.append(_rect(width - thick, 0, width, height))
        p2 = [(left, height), (left, mid + chamfer), (left + chamfer, mid), (width - thick, mid)]
        geoms.append(_thicken(p2, thick))

    elif c == "5":
        p = [
            (0, bot), (right - chamfer, bot), (right, bot + chamfer),
            (right, mid - chamfer), (right - chamfer, mid),
            (left, mid), (left, top), (width, top),
        ]
        geoms.append(_thicken(p, thick))

    elif c == "6":
        p = [
            (thick + cut, mid), (right - chamfer, mid), (right, mid - chamfer),
            (right, bot + chamfer), (right - chamfer, bot),
            (left + chamfer, bot), (left, bot + chamfer),
            (left, top - chamfer), (left + chamfer, top), (width - thick, top),
        ]
        geoms.append(_thicken(p, thick))

    elif c == "7":
        p = [
            (left + thick, 0),
            (left + thick, mid - width / 2 + thick / 2 + thick),
            (right, mid + width / 2 - thick / 2),
            (right, top), (0, top),
        ]
        geoms.append(_thicken(p, thick))

    elif c == "8":
        p = [
            (left, mid-thick), (left, top - chamfer), (left + chamfer, top),
            (right - chamfer, top), (right, top - chamfer), (right, mid+cut+thick),
        ]
        geoms.append(_thicken(p, thick))
        p = [
            (right, mid+thick), (right, bot + chamfer), (right - chamfer, bot),
            (left + chamfer, bot), (left, bot + chamfer), (left, mid-cut-thick),
        ]
        geoms.append(_thicken(p, thick))
        geoms.append(_rect(thick, mid - thick / 2, width - thick, mid + thick / 2))

    elif c == "9":
        p = [
            (width - thick - cut, mid), (left + chamfer, mid), (left, mid + chamfer),
            (left, top - chamfer), (left + chamfer, top),
            (right - chamfer, top), (right, top - chamfer),
            (right, bot + chamfer), (right - chamfer, bot), (thick, bot),
        ]
        geoms.append(_thicken(p, thick))

    elif c == "+":
        geoms.append(_rect(0, mid - thick / 2, width, mid + thick / 2))
        geoms.append(_rect(width / 2 - thick / 2, mid - thick / 2 - 2 * thick, width / 2 + thick / 2, mid - thick / 2))
        geoms.append(_rect(width / 2 - thick / 2, mid + thick / 2, width / 2 + thick / 2, mid + thick / 2 + 2 * thick))

    elif c == "-":
        geoms.append(_rect(0, mid - thick / 2, width, mid + thick / 2))

    elif c == "=":
        geoms.append(_rect(0, mid - 3 * thick / 2, width, mid - thick / 2))
        geoms.append(_rect(0, mid + thick / 2, width, mid + 3 * thick / 2))

    elif c == "*":
        from shapely.affinity import rotate as _rot
        cx, cy = width / 2, mid
        angle = math.degrees(math.tau / 8)
        for base in [
            _rect(0, mid - thick / 2, width, mid + thick / 2),
            _rect(width / 2 - thick / 2, mid - thick / 2 - 2 * thick, width / 2 + thick / 2, mid - thick / 2),
            _rect(width / 2 - thick / 2, mid + thick / 2, width / 2 + thick / 2, mid + thick / 2 + 2 * thick),
        ]:
            geoms.append(_rot(base, angle, origin=(cx, cy)))

    elif c == "/":
        from shapely.affinity import rotate as _rot
        cx, cy = width / 2, mid
        angle = math.degrees(math.tau / 6)
        geoms.append(_rot(_rect(-thick, mid - thick / 2, width + thick, mid + thick / 2), angle, origin=(cx, cy)))

    elif c == "^":
        p = [
            (width / 2 - 1.5 * thick, height - 2.5 * thick),
            (width / 2, height - thick),
            (width / 2 + 1.5 * thick, height - 2.5 * thick),
        ]
        geoms.append(_thicken(p, thick))

    elif c == ".":
        geoms.append(_rect(width / 2 - thick / 2, 0, width / 2 + thick / 2, thick))

    elif c == ":":
        geoms.append(_rect(width / 2 - thick / 2, thick, width / 2 + thick / 2, 2 * thick))
        geoms.append(_rect(width / 2 - thick / 2, 4 * thick, width / 2 + thick / 2, 5 * thick))

    elif c == ";":
        geoms.append(Polygon([
            (width / 2 - thick / 2, 0), (width / 2 + thick / 2, thick),
            (width / 2 + thick / 2, 2 * thick), (width / 2 - thick / 2, 2 * thick),
        ]))
        geoms.append(_rect(width / 2 - thick / 2, 4 * thick, width / 2 + thick / 2, 5 * thick))

    elif c == ",":
        geoms.append(Polygon([
            (width / 2 - thick / 2, 0), (width / 2 + thick / 2, thick),
            (width / 2 + thick / 2, 2 * thick), (width / 2 - thick / 2, 2 * thick),
        ]))

    elif c == '"':
        for dx in (-thick, thick):
            geoms.append(Polygon([
                (width / 2 - thick / 2 + dx, height - 2 * thick),
                (width / 2 + thick / 2 + dx, height - 2 * thick),
                (width / 2 + thick / 2 + dx, height),
                (width / 2 - thick / 2 + dx, height),
            ]))

    elif c == "'":
        geoms.append(Polygon([
            (width / 2 - thick / 2, height - 2 * thick),
            (width / 2 + thick / 2, height - 2 * thick),
            (width / 2 + thick / 2, height),
            (width / 2 - thick / 2, height),
        ]))

    elif c == "[":
        p = [
            (left + 4 * chamfer, bot - thick / 2), (left + chamfer, bot - thick / 2),
            (left + chamfer, top + thick / 2), (left + 4 * chamfer, top + thick / 2),
        ]
        geoms.append(_thicken(p, thick))

    elif c == "]":
        p = [
            (right - 4 * chamfer, bot - thick / 2), (right - chamfer, bot - thick / 2),
            (right - chamfer, top + thick / 2), (right - 4 * chamfer, top + thick / 2),
        ]
        geoms.append(_thicken(p, thick))

    elif c == "(":
        p = [
            (left + 4 * chamfer, bot - thick / 2), (left + 2 * chamfer, bot - thick / 2),
            (left + chamfer, bot + chamfer - thick / 2), (left + chamfer, top - chamfer + thick / 2),
            (left + 2 * chamfer, top + thick / 2), (left + 4 * chamfer, top + thick / 2),
        ]
        geoms.append(_thicken(p, thick))

    elif c == ")":
        p = [
            (right - 4 * chamfer, bot - thick / 2), (right - 2 * chamfer, bot - thick / 2),
            (right - chamfer, bot + chamfer - thick / 2), (right - chamfer, top - chamfer + thick / 2),
            (right - 2 * chamfer, top + thick / 2), (right - 4 * chamfer, top + thick / 2),
        ]
        geoms.append(_thicken(p, thick))

    elif c == "@":
        # Right-pointing triangle arrow (apex at x=width)
        geoms.append(Polygon([(0, mid - width / 2), (width, mid), (0, mid + width / 2)]))

    elif c == "#":
        # Left-pointing triangle arrow (apex at x=0)
        geoms.append(Polygon([(width, mid + width / 2), (0, mid), (width, mid - width / 2)]))

    elif c == "$":
        # Upward-pointing triangle arrow
        geoms.append(Polygon([(width / 2, mid + width / 2), (0, mid - width / 2), (width, mid - width / 2)]))

    elif c == "%":
        # Downward-pointing triangle arrow
        geoms.append(Polygon([(width / 2, mid - width / 2), (width, mid + width / 2), (0, mid + width / 2)]))

    elif c == "`":
        # Filled rounded rectangle (O shape + inner fill)
        p = [
            (left, mid), (left, top - chamfer), (left + chamfer, top),
            (right - chamfer, top), (right, top - chamfer),
            (right, bot + chamfer), (right - chamfer, bot),
            (left + chamfer, bot), (left, bot + chamfer), (left, mid),
        ]
        geoms.append(_thicken(p, thick))
        geoms.append(_rect(left + 0.5 * thick, bot + 0.5 * thick, right - 0.5 * thick, top - 0.5 * thick))

    elif c == "!":
        # Vertical stem with a dot beneath it.
        geoms.append(_rect(width / 2 - thick / 2, 2.5 * thick, width / 2 + thick / 2, height))
        geoms.append(_rect(width / 2 - thick / 2, 0, width / 2 + thick / 2, thick))

    elif c == "?":
        # Top hook curving into a centered stem, with a dot beneath.
        p = [
            (left, top - 2 * chamfer), (left, top - chamfer), (left + chamfer, top),
            (right - chamfer, top), (right, top - chamfer), (right, mid + chamfer),
            (right - chamfer, mid), (width / 2, mid), (width / 2, 2.5 * thick),
        ]
        geoms.append(_thicken(p, thick))
        geoms.append(_rect(width / 2 - thick / 2, 0, width / 2 + thick / 2, thick))

    elif c == "<":
        # Angle bracket pointing left.
        geoms.append(_thicken([(right, top), (left, mid), (right, bot)], thick))

    elif c == ">":
        # Angle bracket pointing right.
        geoms.append(_thicken([(left, top), (right, mid), (left, bot)], thick))

    elif c == "|":
        geoms.append(_rect(width / 2 - thick / 2, 0, width / 2 + thick / 2, height))

    elif c == "\\":
        # Backslash: the slash bar rotated the opposite way.
        from shapely.affinity import rotate as _rot
        cx, cy = width / 2, mid
        angle = -math.degrees(math.tau / 6)
        geoms.append(_rot(_rect(-thick, mid - thick / 2, width + thick, mid + thick / 2), angle, origin=(cx, cy)))

    elif c == "~":
        # Tilde: an S-shaped wave at mid height.
        p = [(0, mid), (width / 3, mid + 1.2 * thick), (2 * width / 3, mid - 1.2 * thick), (width, mid)]
        geoms.append(_thicken(p, thick))

    elif c == "{":
        # Left brace: top/bottom arms reach right, with a leftward pinch at mid.
        p = [
            (left + 4 * chamfer, bot - thick / 2), (left + 2 * chamfer, bot - thick / 2),
            (left + 2 * chamfer, mid - chamfer), (left + chamfer, mid),
            (left + 2 * chamfer, mid + chamfer), (left + 2 * chamfer, top + thick / 2),
            (left + 4 * chamfer, top + thick / 2),
        ]
        geoms.append(_thicken(p, thick))

    elif c == "}":
        # Right brace: mirror of "{".
        p = [
            (right - 4 * chamfer, bot - thick / 2), (right - 2 * chamfer, bot - thick / 2),
            (right - 2 * chamfer, mid - chamfer), (right - chamfer, mid),
            (right - 2 * chamfer, mid + chamfer), (right - 2 * chamfer, top + thick / 2),
            (right - 4 * chamfer, top + thick / 2),
        ]
        geoms.append(_thicken(p, thick))

    elif c == "_":
        # Underscore: a low horizontal bar on the baseline.
        geoms.append(_rect(0, 0, width, thick))

    return [g for g in geoms if not g.is_empty]


def _emit_char(glyphs, ch, x, metrics, scale=1.0, y_off=0.0):
    """Collect one glyph at (x, y_off), scaled by *scale*; return the advanced x cursor.

    Each glyph is appended as a ``(geometry, x, y, y_scale)`` tuple for `_add_geom`,
    so a whole line can be placed once its width is known. *metrics* is the
    ``(thick, height, width, space)`` tuple for the full-size font. Letters/most
    punctuation are built tall then squished vertically by *stretch*; the
    arrow/symbol set keeps its natural height.
    """
    thick, height, width, space = metrics
    stretch = .6
    if ch not in "@#$%^&*_+-=~<>":
        geoms = _char_polys(ch, thick * scale, height * scale / stretch, width * scale)
        glyphs += [(geom, x, y_off, stretch) for geom in geoms]
    else:
        geoms = _char_polys(ch, thick * scale, height * scale, width * scale)
        glyphs += [(geom, x, y_off, 1.0) for geom in geoms]
    if not geoms and not ch.isspace():
        warnings.warn(f"text: no glyph for {ch!r}; it is drawn as blank space")
    return x + width * scale + space


@gf.cell
def text(
    text: str = "abcd",
    size: float = 10.0,
    position: Coordinate = (0, 0),
    justify: str = "left",
    layer: LayerSpec = Layers.label,
) -> Component:
    """Custom geometric text component matching the pGDS font.

    Uppercase and lowercase render identically. Newlines are supported for
    multiline text.

    Subscript / superscript use a LaTeX region delimited by backticks. Outside
    backticks every character -- including ``_ ^ { }`` -- is literal, so names and
    dates like ``"2022_11_30"`` render verbatim. Inside a ``` `...` ``` region,
    ``_`` and ``^`` mark subscript / superscript: ``_x`` / ``^x`` for a single
    character and ``_{...}`` / ``^{...}`` for a multi-character span. Examples::

        "n_eff = 2.1"      -> literal (no scripting)
        "`n_{eff}` = 2.1"  -> n + eff(subscript) + " = 2.1"
        "`10^{-3}, C_2`"   -> 10 + -3(super), C + 2(sub)

    A backtick is consumed as the toggle rather than drawn, except inside a
    sub/superscript target, where it is an ordinary glyph (``` `_` ```); a region
    left open at a line break is closed implicitly. Sub/superscripts are scaled by
    ``_SCRIPT_SCALE`` (subscript on the baseline, superscript at the top).

    Parameters
    ----------
    text : str
        String to render.
    size : float
        Character height in microns.
    position : (float, float)
        (x, y) anchor position.
    justify : str
        Horizontal alignment: ``"left"``, ``"center"`` or ``"right"``.
    layer : LayerSpec
        GDS layer for the text polygons.
    """
    height = float(size)
    width = _WIDTH_RATIO * height
    thick = _THICK_RATIO * height
    space = thick
    line_height = height + 2 * thick
    metrics = (thick, height, width, space)

    justify = justify.lower()
    if justify not in ("left", "center", "right"):
        raise ValueError(f"text: justify must be 'left', 'center' or 'right', got {justify!r}")

    c = gf.Component()
    layer_idx = gf.get_layer(layer)
    dbu = c.kcl.dbu
    lines = text.split("\n")

    for line_idx, line in enumerate(lines):
        glyphs = []
        x = 0.0
        i = 0
        n = len(line)
        latex = False  # backtick region state, reset per line
        while i < n:
            ch = line[i]
            if ch == "`":            # toggle LaTeX region (consumed, never drawn)
                latex = not latex
                i += 1
                continue
            if ch == " ":
                x += width + space
                i += 1
                continue
            if latex and ch in "_^":
                # subscript on the baseline, superscript raised to the top
                y_off = 0.0 if ch == "_" else height * (1 - _SCRIPT_SCALE)
                i += 1
                if i < n and line[i] == "{":          # multi-char span
                    end = line.find("}", i + 1)
                    if end == -1:                     # unclosed: take rest of line
                        target, i = line[i + 1:], n
                    else:
                        target, i = line[i + 1:end], end + 1
                elif i < n:                           # single-char shorthand
                    target, i = line[i], i + 1
                else:                                 # trailing _/^ : nothing
                    target = ""
                for t in target:
                    if t == " ":
                        x += (width + space) * _SCRIPT_SCALE
                    else:
                        x = _emit_char(glyphs, t, x, metrics,
                                       scale=_SCRIPT_SCALE, y_off=y_off)
                continue
            # full-size glyph (literal mode, or a non-script char inside a region)
            x = _emit_char(glyphs, ch, x, metrics)
            i += 1

        line_width = x - space if x > 0 else 0.0
        dx = {"left": 0.0, "center": -line_width / 2, "right": -line_width}[justify]
        dy = -line_idx * line_height
        for geom, gx, gy, y_scale in glyphs:
            _add_geom(c, geom, gx + dx, gy + dy, layer_idx, dbu, y_scale=y_scale)

    # Snap the anchor once, so the whole string shifts by a whole number of dbu.
    c.transform(kdb.DTrans(float(position[0]), float(position[1])).to_itype(dbu))
    c.flatten()  # Merges the overlapping strokes of each glyph into one polygon.
    return c
