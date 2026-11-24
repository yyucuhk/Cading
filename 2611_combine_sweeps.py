"""Assemble two identical 18-device sweeps in one GDS."""

from __future__ import annotations

import json
from pathlib import Path

import klayout.db as kdb


ROOT = Path(__file__).resolve().parent
GDS_DIR = ROOT / "gds"

SWEEP_GDS = GDS_DIR / "2611_design1_design2_L1400_18device_sweep.gds"
OUTPUT_GDS = GDS_DIR / "2611_duplicated_18device_sweep.gds"

SECOND_SET_X_SHIFT = 6600.0
MARKER_INSET = 150.0
DBU = 0.001
EXCLUDED_LAYERS = {(16, 0), (51, 0)}
MARKER_REMOVE_LAYERS = {(2, 0), (10, 0), (15, 0)}
MARKER_SHAPE_SOURCE_LAYER = (10, 0)
MARKER_TARGET_LAYER = (15, 0)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_top(path: Path) -> tuple[kdb.Layout, kdb.Cell]:
    require(path.exists(), f"Missing source GDS: {path}")
    layout = kdb.Layout()
    layout.read(str(path))
    tops = layout.top_cells()
    require(len(tops) == 1, f"Expected one top cell in {path.name}")
    return layout, tops[0]


def copy_flat(
    destination: kdb.Layout,
    source: kdb.Layout,
    source_top: kdb.Cell,
    name: str,
) -> kdb.Cell:
    """Copy every recursive source polygon into one reusable local cell."""
    cell = destination.create_cell(name)
    marker_regions: dict[tuple[int, int], kdb.Region] = {}
    for instance in source_top.each_inst():
        if not any(
            word in instance.cell.name.upper()
            for word in ("FRAME", "MARKER")
        ):
            continue
        for source_layer_index in source.layer_indexes():
            info = source.get_info(source_layer_index)
            layer = (info.layer, info.datatype)
            if layer not in MARKER_REMOVE_LAYERS:
                continue
            region = kdb.Region(
                instance.cell.begin_shapes_rec(source_layer_index)
            ).transformed(instance.trans)
            marker_regions.setdefault(layer, kdb.Region())
            marker_regions[layer] |= region

    outline = kdb.Region(
        source_top.begin_shapes_rec(source.layer(0, 0))
    ).merged()
    require(not outline.is_empty(), "Missing Layer-0 chip area")
    outline_center = outline.bbox().center()

    for source_layer_index in source.layer_indexes():
        region = kdb.Region(
            source_top.begin_shapes_rec(source_layer_index)
        )
        if region.is_empty():
            continue
        info = source.get_info(source_layer_index)
        layer = (info.layer, info.datatype)
        if layer in EXCLUDED_LAYERS:
            continue
        markers = marker_regions.get(layer, kdb.Region())
        region -= markers
        destination_layer = destination.layer(info.layer, info.datatype)
        cell.shapes(destination_layer).insert(region)
        if layer == MARKER_SHAPE_SOURCE_LAYER and not markers.is_empty():
            target = destination.layer(*MARKER_TARGET_LAYER)
            inset_markers = kdb.Region()
            inset_dbu = round(MARKER_INSET / DBU)
            for polygon in markers.each():
                center = polygon.bbox().center()
                dx = inset_dbu if center.x < outline_center.x else -inset_dbu
                dy = inset_dbu if center.y < outline_center.y else -inset_dbu
                inset_markers.insert(polygon.transformed(kdb.Trans(dx, dy)))
            require((inset_markers - outline).is_empty(),
                    "Inset markers must stay inside Layer 0")
            cell.shapes(target).insert(inset_markers)
    return cell


def box_um(box: kdb.Box) -> list[float]:
    return [
        box.left * DBU,
        box.bottom * DBU,
        box.right * DBU,
        box.top * DBU,
    ]


def build() -> dict:
    require(SECOND_SET_X_SHIFT > 0, "SECOND_SET_X_SHIFT must be positive")
    require(MARKER_INSET > 0, "MARKER_INSET must be positive")

    sweep_layout, sweep_top = read_top(SWEEP_GDS)
    require(abs(sweep_layout.dbu - DBU) < 1e-12,
            "Source GDS database unit must be 0.001 um")

    output = kdb.Layout()
    output.dbu = DBU
    top = output.create_cell("DUPLICATED_18DEVICE_SWEEP")
    sweep_cell = copy_flat(output, sweep_layout, sweep_top, "SWEEP_18_DEVICE")

    sweep_box = sweep_cell.bbox()
    source_center = sweep_box.center()
    half_shift_dbu = round(SECOND_SET_X_SHIFT / (2 * DBU))
    first_trans = kdb.Trans(
        -half_shift_dbu - source_center.x,
        -source_center.y,
    )
    second_trans = kdb.Trans(
        half_shift_dbu - source_center.x,
        -source_center.y,
    )
    top.insert(kdb.CellInstArray(sweep_cell.cell_index(), first_trans))
    top.insert(kdb.CellInstArray(sweep_cell.cell_index(), second_trans))

    first_box = sweep_box.transformed(first_trans)
    second_box = sweep_box.transformed(second_trans)
    actual_x_shift = (second_trans.disp.x - first_trans.disp.x) * DBU
    require(abs(actual_x_shift - SECOND_SET_X_SHIFT) < DBU,
            "Incorrect second-set x translation")
    overall_box = top.bbox()
    require(overall_box.center() == kdb.Point(0, 0),
            "Duplicated sweep must be centered at the origin")
    OUTPUT_GDS.parent.mkdir(parents=True, exist_ok=True)
    options = kdb.SaveLayoutOptions()
    options.gds2_max_vertex_count = 2000
    output.write(str(OUTPUT_GDS), options)

    report = {
        "status": "PASS",
        "output_gds": str(OUTPUT_GDS),
        "second_set_x_shift_um": actual_x_shift,
        "first_sweep_bbox_um": box_um(first_box),
        "second_sweep_bbox_um": box_um(second_box),
        "overall_bbox_um": box_um(overall_box),
        "placement": "pair centered at (0,0); second sweep is +6600 um in x",
        "marker_inset_um": MARKER_INSET,
        "marker_layer": list(MARKER_TARGET_LAYER),
        "marker_layers_removed": [
            list(layer) for layer in sorted(MARKER_REMOVE_LAYERS)
        ],
        "marker_shape_source_layer": list(MARKER_SHAPE_SOURCE_LAYER),
    }
    OUTPUT_GDS.with_suffix(".json").write_text(
        json.dumps(report, indent=2)
    )
    return report


if __name__ == "__main__":
    result = build()
    print(
        f"PASS: wrote {result['output_gds']} "
        f"(second set x shift {result['second_set_x_shift_um']:.1f} um)"
    )
