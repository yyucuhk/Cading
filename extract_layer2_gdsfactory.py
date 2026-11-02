"""Extract every polygon on GDS layer 2 and write it with gdsfactory.

The input GDS contains a single, flat cell made entirely from BOUNDARY records.
This small reader deliberately supports that subset, which makes the geometry
reproduction independent of gdstk/gdspy while leaving gdsfactory responsible
for creating the output GDS.

Usage:
    python extract_layer2_gdsfactory.py INPUT.gds [-o layer2_only.gds]
"""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import gdsfactory as gf


LAYER_NUMBER = 2


def _gds_real8(value: bytes) -> float:
    """Decode an eight-byte GDSII real number."""
    if value == b"\0" * 8:
        return 0.0
    sign = -1.0 if value[0] & 0x80 else 1.0
    exponent = (value[0] & 0x7F) - 64
    mantissa = int.from_bytes(value[1:], byteorder="big") / 2**56
    return sign * mantissa * 16**exponent


def read_layer_polygons(
    gds_path: str | Path, layer_number: int = LAYER_NUMBER
) -> tuple[str, list[tuple[int, list[tuple[float, float]]]]]:
    """Return ``(cell_name, [(datatype, polygon), ...])`` from a flat GDS.

    Polygon coordinates are returned in micrometres, the coordinate unit used
    by gdsfactory.
    """
    raw = Path(gds_path).read_bytes()
    offset = 0
    cell_name = "layer2_from_source"
    dbu_um = 0.001  # Standard 1 nm database unit, overridden by UNITS below.
    current: dict[str, object] | None = None
    polygons: list[tuple[int, list[tuple[float, float]]]] = []

    while offset < len(raw):
        if offset + 4 > len(raw):
            raise ValueError(f"Truncated GDS record header at byte {offset}")
        size, record_type, _data_type = struct.unpack(">HBB", raw[offset : offset + 4])
        if size < 4 or offset + size > len(raw):
            raise ValueError(f"Invalid GDS record length {size} at byte {offset}")
        payload = raw[offset + 4 : offset + size]
        offset += size

        if record_type == 0x03:  # UNITS
            # First REAL8 is the database unit expressed in user units. GDS
            # layouts conventionally use micrometres as their user unit.
            dbu_um = _gds_real8(payload[:8])
        elif record_type == 0x06:  # STRNAME
            cell_name = payload.rstrip(b"\0").decode("ascii", errors="replace")
        elif record_type == 0x08:  # BOUNDARY
            current = {}
        elif record_type == 0x0D and current is not None:  # LAYER
            current["layer"] = struct.unpack(">h", payload)[0]
        elif record_type == 0x0E and current is not None:  # DATATYPE
            current["datatype"] = struct.unpack(">h", payload)[0]
        elif record_type == 0x10 and current is not None:  # XY
            integers = struct.unpack(f">{len(payload) // 4}i", payload)
            current["points"] = [
                (integers[i] * dbu_um, integers[i + 1] * dbu_um)
                for i in range(0, len(integers), 2)
            ]
        elif record_type == 0x11 and current is not None:  # ENDEL
            if current.get("layer") == layer_number:
                points = current.get("points")
                if not isinstance(points, list) or len(points) < 4:
                    raise ValueError("Layer-2 BOUNDARY is missing valid XY data")
                # GDS boundaries repeat their first point; gdsfactory does not
                # require the closing duplicate.
                if points[0] == points[-1]:
                    points = points[:-1]
                polygons.append((int(current.get("datatype", 0)), points))
            current = None

    if not polygons:
        raise ValueError(f"No BOUNDARY polygons found on layer {layer_number}")
    return cell_name, polygons


def build_layer2_component(gds_path: str | Path) -> gf.Component:
    """Create a gdsfactory component containing only source layer 2."""
    source_name, polygons = read_layer_polygons(gds_path, LAYER_NUMBER)
    component = gf.Component(f"{source_name}_L2")
    for datatype, points in polygons:
        component.add_polygon(points, layer=(LAYER_NUMBER, datatype))
    return component


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_gds", type=Path)
    parser.add_argument("-o", "--output", type=Path, default=Path("layer2_only.gds"))
    args = parser.parse_args()

    component = build_layer2_component(args.input_gds)
    component.write_gds(args.output)
    print(f"Wrote {args.output} with {len(component.get_polygons())} Layer-2 polygons")


if __name__ == "__main__":
    main()
