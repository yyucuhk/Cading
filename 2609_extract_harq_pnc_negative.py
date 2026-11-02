"""Extract four negative-tone PNC material geometries for COMSOL.

The source GDS contains positive-resist exposure polygons.  A physical PNC is
the narrow, modulated material gap between two neighboring exposure polygons.
This script detects the four six-device groups, extracts one representative
gap from each group, centers it at (0, 0), and writes one GDS per geometry.
"""

import sys
from pathlib import Path

import gdsfactory as gf

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from layout.core import GDS
from layout.core.layers import Layers


SOURCE_GDS = Path(
    "/Users/yuyue/Documents/Projects/Quantum transduction/Fab/2608_harq_PNC.gds"
)
SOURCE_LAYER = tuple(Layers.wg_doubleetch)  # GDS 11/0
OUTPUT_LAYER = tuple(Layers.wg_doubleetch)  # retained for COMSOL import

# The exposure polygons are about 500 um long and 4--6 um high.
MIN_EXPOSURE_LENGTH = 400.0
MIN_EXPOSURE_HEIGHT = 3.0
MAX_EXPOSURE_HEIGHT = 7.0

# Adjacent polygons in one PNC group are separated by 0.15--0.30 um.
MAX_PNC_GAP = 1.0
EXTRACTION_HALF_WIDTH = 1.0


def _recursive_polygons(component: gf.Component, layer: tuple[int, int]):
    """Return recursive integer-coordinate polygons from one GDS layer."""
    layer_index = component.kcl.layout.layer(*layer)
    iterator = component.kdb_cell.begin_shapes_rec(layer_index)
    polygons = []
    while not iterator.at_end():
        shape = iterator.shape()
        transform = iterator.trans()
        if shape.is_box():
            polygon = gf.kdb.Polygon(shape.box).transformed(transform)
        elif shape.is_polygon():
            polygon = shape.polygon.transformed(transform)
        elif shape.is_path():
            polygon = shape.path.polygon().transformed(transform)
        else:
            iterator.next()
            continue
        polygons.append(polygon)
        iterator.next()
    return polygons


def _find_pnc_groups(component: gf.Component):
    """Detect the four lower-row groups of six positive-tone polygons."""
    dbu = component.kcl.layout.dbu
    candidates = []
    for polygon in _recursive_polygons(component, SOURCE_LAYER):
        box = polygon.bbox()
        width = box.width() * dbu
        height = box.height() * dbu
        # Select the lower repeated PNC row. The companion row is equivalent.
        center_y = (box.bottom + box.top) * dbu / 2
        if (
            width >= MIN_EXPOSURE_LENGTH
            and MIN_EXPOSURE_HEIGHT <= height <= MAX_EXPOSURE_HEIGHT
            and center_y < 13620.0
        ):
            candidates.append(polygon)

    candidates.sort(key=lambda polygon: polygon.bbox().left)
    groups = []
    current = []
    for polygon in candidates:
        if current:
            previous_box = current[-1].bbox()
            gap = (polygon.bbox().left - previous_box.right) * dbu
            if gap > MAX_PNC_GAP:
                groups.append(current)
                current = []
        current.append(polygon)
    if current:
        groups.append(current)

    if len(groups) != 4 or any(len(group) != 6 for group in groups):
        counts = [len(group) for group in groups]
        raise RuntimeError(
            f"Expected four groups of six PNC exposures; detected {counts}"
        )
    return groups


def _extract_negative_pnc(
    left_polygon: gf.kdb.Polygon,
    right_polygon: gf.kdb.Polygon,
    dbu: float,
) -> tuple[gf.kdb.Region, float]:
    """Return the material gap between one adjacent exposure pair."""
    left_box = left_polygon.bbox()
    right_box = right_polygon.bbox()
    boundary_x = (left_box.right + right_box.left) // 2
    half_width_dbu = round(EXTRACTION_HALF_WIDTH / dbu)
    bottom = max(left_box.bottom, right_box.bottom)
    top = min(left_box.top, right_box.top)
    clip_box = gf.kdb.Box(
        boundary_x - half_width_dbu,
        bottom,
        boundary_x + half_width_dbu,
        top,
    )

    exposure = gf.kdb.Region(left_polygon)
    exposure.insert(right_polygon)
    material = gf.kdb.Region(clip_box) - exposure
    material.merge()

    if material.is_empty():
        raise RuntimeError("PNC complement is empty")

    material_box = material.bbox()
    center_x = (material_box.left + material_box.right) // 2
    # Center the PNC in x and place its lower/start edge exactly at y = 0.
    # This gives COMSOL imports a consistent longitudinal origin.
    material.transform(gf.kdb.Trans(-center_x, -material_box.bottom))
    nominal_gap_um = (right_box.left - left_box.right) * dbu
    return material, nominal_gap_um


def _write_region_dxf(
    region: gf.kdb.Region,
    dbu: float,
    output_path: Path,
) -> None:
    """Write a Region as closed DXF R12 polylines in micrometer units."""
    lines = [
        "0", "SECTION", "2", "HEADER",
        "9", "$INSUNITS", "70", "13",  # 13 = micrometers
        "0", "ENDSEC",
        "0", "SECTION", "2", "ENTITIES",
    ]
    polygon_count = 0
    for polygon in region.each():
        points = list(polygon.each_point_hull())
        if len(points) < 3:
            continue
        polygon_count += 1
        lines.extend(
            ["0", "POLYLINE", "8", "PNC", "66", "1", "70", "1"]
        )
        for point in points:
            lines.extend(
                [
                    "0", "VERTEX", "8", "PNC",
                    "10", f"{point.x * dbu:.9f}",
                    "20", f"{point.y * dbu:.9f}",
                    "30", "0.0",
                ]
            )
        lines.extend(["0", "SEQEND"])
    if polygon_count == 0:
        raise RuntimeError("Cannot export an empty PNC Region to DXF")
    lines.extend(["0", "ENDSEC", "0", "EOF"])
    output_path.write_text("\n".join(lines) + "\n", encoding="ascii")


def write_negative_pncs() -> list[Path]:
    """Extract and write centered negative-tone PNC GDS and DXF files."""
    if not SOURCE_GDS.is_file():
        raise FileNotFoundError(f"Source GDS not found: {SOURCE_GDS}")

    gf.gpdk.PDK.activate()
    source = gf.import_gds(SOURCE_GDS, rename_duplicated_cells=True)
    dbu = source.kcl.layout.dbu
    groups = _find_pnc_groups(source)
    output_paths = []

    for group in groups:
        # Any adjacent pair is equivalent. Use the first internal boundary.
        material, nominal_gap_um = _extract_negative_pnc(
            group[0], group[1], dbu
        )
        minimum_feature_nm = round(nominal_gap_um * 1000)
        component = gf.Component(
            f"PNC_negative_MF{minimum_feature_nm}nm"
        )
        component.kdb_cell.shapes(
            component.kcl.layout.layer(*OUTPUT_LAYER)
        ).insert(material)
        component.info["source_gds"] = str(SOURCE_GDS)
        component.info["source_layer"] = SOURCE_LAYER
        component.info["tone"] = "negative_material_geometry"
        component.info["minimum_feature_nm"] = minimum_feature_nm

        output_path = Path(GDS) / f"{component.name}.gds"
        save_options = gf.kdb.SaveLayoutOptions()
        save_options.gds2_max_vertex_count = 2000
        component.write_gds(
            output_path,
            save_options=save_options,
            with_metadata=False,
        )
        output_paths.append(output_path)

        dxf_path = Path(GDS) / f"{component.name}.dxf"
        _write_region_dxf(material, dbu, dxf_path)
        output_paths.append(dxf_path)

    return output_paths


if __name__ == "__main__":
    for path in write_negative_pncs():
        print(path)
