"""Replace only Layer-21 Pt absorbers in 261108_optim12_Pt_abs.gds.

Each patterned absorber is a rectangular-frame lattice clipped to the
existing 30 um circular footprint.  At every (px, py) lattice point, the
outer rectangle is wx by wy and its centered hole is hx by hy.  Normally
hx = wx - 2*dx; when wx equals px, hx = wx - dx.  In all cases hy = wy - dy.
The first six slots in each vertical block are D2,D1,D3,D2,D1,D4;
the remaining two physical slots are empty.
All dimensions are micrometres.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import klayout.db as kdb


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "gds/261108_optim12_Pt_abs.gds"
OUTPUT = ROOT / "gds/261108_optim12_Pt_abs_updated.gds"
ABSORBER_LAYER = (21, 0)
FOOTPRINT_DIAMETER = 31.0
ROW_SEQUENCE = (2, 1, 3, 2, 1, 4, None, None)

DESIGNS = {
    1: dict(px=0.900, wx=0.450, py=0.450, wy=0.450, dx=0.120, dy=0.080,
            wavelength_um=3.2, absorption_percent=65),
    2: dict(px=0.350, wx=0.350, py=0.350, wy=0.350, dx=0.100, dy=0.100,
            wavelength_um=2.6, absorption_percent=70),
    3: dict(px=0.675, wx=0.450, py=0.900, wy=0.900, dx=0.125, dy=0.125,
            wavelength_um=5.0, absorption_percent=60),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def circle(cx: float, cy: float, radius: float, dbu: float) -> kdb.Region:
    points = [
        kdb.Point(
            round((cx + radius * math.cos(2 * math.pi * i / 360)) / dbu),
            round((cy + radius * math.sin(2 * math.pi * i / 360)) / dbu),
        )
        for i in range(360)
    ]
    return kdb.Region(kdb.Polygon(points))


def rectangle(cx: float, cy: float, width: float, height: float, dbu: float) -> kdb.Box:
    return kdb.Box(
        round((cx - width / 2) / dbu),
        round((cy - height / 2) / dbu),
        round((cx + width / 2) / dbu),
        round((cy + height / 2) / dbu),
    )


def patterned_absorber(cx: float, cy: float, design: dict, dbu: float) -> kdb.Region:
    radius = FOOTPRINT_DIAMETER / 2
    px, wx = design["px"], design["wx"]
    py, wy = design["py"], design["wy"]
    dx, dy = design["dx"], design["dy"]
    require(0 < wx <= px and 0 < wy <= py, "Mesh span must not exceed its period")
    hx = wx - dx if math.isclose(wx, px) else wx - 2 * dx
    hy = wy - dy
    require(hx > 0 and hy > 0, "Mesh line width leaves no central hole")

    outers = kdb.Region()
    holes = kdb.Region()
    nx = math.ceil((radius + wx) / px)
    ny = math.ceil((radius + wy) / py)
    for ix in range(-nx, nx + 1):
        x = cx + ix * px
        for iy in range(-ny, ny + 1):
            y = cy + iy * py
            outers.insert(rectangle(x, y, wx, wy, dbu))
            holes.insert(rectangle(x, y, hx, hy, dbu))
    return ((outers - holes) & circle(cx, cy, radius, dbu)).merged()


def detect_slot_centers(region: kdb.Region, dbu: float) -> tuple[list[float], list[list[float]]]:
    """Recover the existing eight columns and two eight-row blocks."""
    rows: dict[float, list[kdb.Box]] = {}
    for polygon in region.each():
        box = polygon.bbox()
        y = round(box.center().y * dbu, 3)
        rows.setdefault(y, []).append(box)
    ys = sorted(rows)
    require(len(ys) == 16, f"Expected 16 absorber rows, found {len(ys)}")

    row_x_centers = []
    for boxes in rows.values():
        intervals = sorted((box.left * dbu, box.right * dbu) for box in boxes)
        groups: list[list[float]] = []
        for left, right in intervals:
            if not groups or left - groups[-1][1] > 50:
                groups.append([left, right])
            else:
                groups[-1][1] = max(groups[-1][1], right)
        require(len(groups) == 8, f"Expected 8 absorber columns, found {len(groups)}")
        row_x_centers.append([(left + right) / 2 for left, right in groups])
    raw_xs = [sum(row[i] for row in row_x_centers) / len(row_x_centers)
              for i in range(8)]
    xs = [round(value / 100) * 100.0 for value in raw_xs]
    require(all(abs(row[i] - raw_xs[i]) < 1 for row in row_x_centers for i in range(8)),
            "Absorber column centers are inconsistent")
    require(ys[7] - ys[6] < 200 and ys[8] - ys[7] > 1000,
            "Could not identify the two eight-row blocks")
    return xs, [ys[:8], ys[8:]]


def build() -> dict:
    require(SOURCE.exists(), f"Missing source GDS: {SOURCE}")
    layout = kdb.Layout()
    layout.read(str(SOURCE))
    require(len(layout.top_cells()) == 1, "Expected one top cell")
    top = layout.top_cell()
    dbu = layout.dbu
    layer_index = layout.layer(*ABSORBER_LAYER)
    original = kdb.Region(top.shapes(layer_index)).merged()
    require(not original.is_empty(), "Layer 21 is empty")
    xs, row_blocks = detect_slot_centers(original, dbu)

    replacement = kdb.Region()
    slots = []
    for block_index, ys in enumerate(row_blocks):
        for x in xs:
            for row_index, (y, design_number) in enumerate(zip(ys, ROW_SEQUENCE, strict=True)):
                if design_number == 4:
                    geometry = circle(x, y, FOOTPRINT_DIAMETER / 2, dbu)
                elif design_number is None:
                    geometry = kdb.Region()
                else:
                    geometry = patterned_absorber(x, y, DESIGNS[design_number], dbu)
                replacement |= geometry
                slots.append(dict(block=block_index, column_x_um=x, row=row_index,
                                  center_y_um=y, design=design_number))

    top.shapes(layer_index).clear()
    top.shapes(layer_index).insert(replacement.merged())
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    options = kdb.SaveLayoutOptions()
    options.gds2_max_vertex_count = 2000
    layout.write(str(OUTPUT), options)

    check = kdb.Layout()
    check.read(str(OUTPUT))
    require(len(check.top_cells()) == 1, "Exported GDS must have one top cell")
    require((kdb.Region(check.top_cell().shapes(check.layer(101, 0))).merged()
             ^ kdb.Region(top.shapes(layout.layer(101, 0))).merged()).is_empty(),
            "Layer 101 changed during export")
    return dict(
        status="PASS",
        source_gds=str(SOURCE),
        output_gds=str(OUTPUT),
        database_unit_um=dbu,
        footprint_diameter_um=FOOTPRINT_DIAMETER,
        absorber_layer=list(ABSORBER_LAYER),
        columns_um=xs,
        row_blocks_um=row_blocks,
        row_sequence=list(ROW_SEQUENCE),
        designs=DESIGNS,
        slots=slots,
        layer101_preserved=True,
    )


if __name__ == "__main__":
    report = build()
    OUTPUT.with_suffix(".json").write_text(json.dumps(report, indent=2))
    print(f"PASS: wrote {OUTPUT}")
