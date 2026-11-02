"""Check both cleaned CAD generators without assuming exactly three devices.

Run: python validate_2611_cleanup.py
Optional pre-cleanup baseline: python validate_2611_cleanup.py capture
Outputs are written only under validation_2611_cleanup/, not over production GDS.
A capture is a geometry snapshot; use it BEFORE changing code, not afterward.
Without a baseline the script checks current geometry but cannot prove preservation.
These tests are not optical simulations or fabrication design-rule checks.
"""

from __future__ import annotations

import argparse
import ast
import importlib.util
import json
from pathlib import Path

import gdsfactory as gf
import klayout.db as kdb
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "validation_2611_cleanup"
MODELS = {
    "design1": ("2611_tapered_yue_dual_resist.py", "2611_tapered_yue_dual_resist_L10_L12.gds"),
    "design2": ("2611_tapered_yue_design2.py", "2611_tapered_yue_design2.gds"),
}
REPORT = {"checks": [], "failures": []}
TOL = 0.005


def check(name, ok, **details):
    row = {"name": name, "passed": bool(ok), **details}
    REPORT["checks"].append(row)
    if not ok:
        REPORT["failures"].append(name)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "validation_report.json").write_text(json.dumps(REPORT, indent=2))
    print(json.dumps(row), flush=True)


def load(key, filename):
    spec = importlib.util.spec_from_file_location(f"cleanup_{key}", ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    gf.gpdk.PDK.activate()
    return module


def regions(layout, cell):
    result = {}
    for index in layout.layer_indices():
        reg = kdb.Region(cell.begin_shapes_rec(index)).merged()
        if not reg.is_empty():
            info = layout.get_info(index)
            result[(info.layer, info.datatype)] = reg
    return result


def read(path):
    layout = kdb.Layout()
    layout.read(str(path))
    return layout, regions(layout, layout.top_cell())


def memory(component):
    return {
        tuple(map(int, layer)): kdb.Region(
            component.kdb_cell.begin_shapes_rec(gf.get_layer(tuple(layer)))
        ).merged()
        for layer in component.layers
    }


def differences(a, b, dbu):
    return {
        f"{layer[0]}/{layer[1]}": float(
            (a.get(layer, kdb.Region()) ^ b.get(layer, kdb.Region())).area() * dbu**2
        )
        for layer in sorted(set(a) | set(b))
    }


def write(component, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    options = gf.kdb.SaveLayoutOptions()
    options.gds2_max_vertex_count = 2000
    component.write_gds(path, save_options=options, with_metadata=False)


def stable_settings(module, filename):
    tree = ast.parse((ROOT / filename).read_text())
    names = [
        n.targets[0].id
        for n in tree.body
        if isinstance(n, ast.Assign) and isinstance(n.targets[0], ast.Name)
    ]
    result = {}
    for name in names:
        if name in ("REPOSITORY_ROOT", "GDS", "EBL_SUBFIELD_SIZE"):
            continue
        value = getattr(module, name)
        try:
            result[name] = json.loads(json.dumps(value))
        except TypeError:
            pass
    return result


def device_instances(layout):
    layer = layout.find_layer(51, 0)
    if layer is None:
        return []
    return sorted(
        [
            i
            for i in layout.top_cell().each_inst()
            if not kdb.Region(i.cell.begin_shapes_rec(layer)).is_empty()
        ],
        key=lambda i: i.trans.disp.y,
    )


def run(capture=False):
    folder = OUT / ("before" if capture else "after")
    folder.mkdir(parents=True, exist_ok=True)
    baseline_available = all((OUT / "before" / f"{key}_layout.gds").exists() for key in MODELS)
    REPORT["baseline_comparison_performed"] = baseline_available and not capture
    for key, (filename, gdsname) in MODELS.items():
        m = load(key, filename)
        settings = stable_settings(m, filename)
        (folder / f"{key}_settings.json").write_text(json.dumps(settings, indent=2))
        lengths = tuple(
            dict.fromkeys([*map(float, m.MEANDER_STRAIGHT_LENGTHS), 200.0, 800.0, 1400.0])
        )
        components, snapshots = [], []
        for length in lengths:
            device = m.make_device(length)
            components.append(device)
            snap = memory(device)
            snapshots.append(snap)
            dbu = float(device.kcl.dbu)
            output = folder / f"{key}_L{length:g}.gds"
            write(device, output)
            if capture:
                continue
            label = f"{key}_L{length:g}"
            current_layout, actual = read(output)
            check(
                label + "_all_device_layers_present",
                {(2, 0), (10, 0), (11, 0), (13, 0), (14, 0), (15, 0), (16, 0), (51, 0)}
                <= set(actual),
            )
            check(label + "_no_construction_layer", (999, 0) not in actual)
            tips = m.alignment_tip_region(
                device.info["tip_centers"], device.info["connection_orientations"], dbu
            )
            check(
                label + "_alignment_ticks_only_on_layer15",
                (tips - actual[(15, 0)]).is_empty()
                and (tips & actual[(10, 0)]).is_empty()
                and len(list(tips.each())) == 4,
            )
            ridge = (actual[(15, 0)] - tips).merged()
            check(label + "_optical_ridge_connected", len(list(ridge.each())) == 1)
            check(
                label + "_ridge_protection_connected",
                len(list(actual[(16, 0)].merged().each())) == 1,
            )
            check(
                label + "_layer14_and_16_cover_ridge",
                (ridge - actual[(14, 0)]).area() * dbu**2 <= TOL
                and (ridge - actual[(16, 0)]).area() * dbu**2 <= TOL,
            )
            before = OUT / "before" / output.name
            if before.exists():
                old_layout, old = read(before)
                diff = differences(old, actual, dbu)
                check(
                    label + "_every_polygon_unchanged",
                    all(v == 0 for v in diff.values()),
                    xor_um2=diff,
                )
            for earlier, original in zip(components[:-1], snapshots[:-1], strict=True):
                diff = differences(original, memory(earlier), dbu)
                check(
                    label + f'_does_not_mutate_L{earlier.info["meander_straight_length_um"]:g}',
                    all(v == 0 for v in diff.values()),
                )
        if not capture:
            for length, original in reversed(list(zip(lengths, snapshots))):
                repeated = m.make_device(length)
                diff = differences(original, memory(repeated), float(repeated.kcl.dbu))
                check(
                    f"{key}_L{length:g}_build_order_independent", all(v == 0 for v in diff.values())
                )

        layout_path = folder / f"{key}_layout.gds"
        if capture:
            write(m.build_layout(), layout_path)
            continue
        saved = m.export_gds(layout_path)
        check(
            key + "_export_function_returns_written_path", saved == layout_path and saved.exists()
        )
        final, actual = read(saved)
        dbu = float(final.dbu)
        check(key + "_top_cell_name_preserved", final.top_cell().name == m.TOP_CELL_NAME)
        instances = device_instances(final)
        check(
            key + "_configured_device_count",
            len(instances) == len(m.MEANDER_STRAIGHT_LENGTHS),
            count=len(instances),
        )
        edges = [
            (
                kdb.Region(inst.cell.begin_shapes_rec(final.layer(11, 0)))
                .transformed(inst.trans)
                .bbox()
                .left
                * dbu
            )
            for inst in instances
        ]
        check(
            key + "_actual_left_bends_aligned",
            bool(edges) and np.allclose(edges, m.LEFT_BEND_ALIGNMENT_X, atol=dbu / 2, rtol=0),
            x_um=edges,
        )
        positions = [[i.trans.disp.x * dbu, i.trans.disp.y * dbu] for i in instances]
        expected_y = [m.BOTTOM_DEVICE_Y + i * m.DEVICE_Y_PITCH for i in range(len(instances))]
        check(
            key + "_actual_y_placement",
            np.allclose([p[1] for p in positions], expected_y, atol=dbu / 2, rtol=0),
            positions_um=positions,
        )
        before = OUT / "before" / layout_path.name
        if before.exists():
            old_layout, old = read(before)
            diff = differences(old, actual, dbu)
            check(
                key + "_complete_layout_exactly_unchanged",
                all(v == 0 for v in diff.values()),
                xor_um2=diff,
            )
            before_settings = json.loads((OUT / "before" / f"{key}_settings.json").read_text())
            altered = [
                name for name, value in before_settings.items() if settings.get(name) != value
            ]
            check(
                key + "_all_existing_parameter_values_preserved",
                not altered,
                altered_parameters=altered,
            )
        # This comparison is only a cleanup regression; future edited geometry
        # should be validated against a newly approved baseline instead.
        if baseline_available:
            stored, stored_regions = read(ROOT / "gds" / gdsname)
            diff = differences(stored_regions, actual, dbu)
            check(
                key + "_matches_approved_repository_gds",
                all(v <= TOL for v in diff.values()),
                xor_um2=diff,
            )
        saved_lengths, saved_name = m.MEANDER_STRAIGHT_LENGTHS, m.TOP_CELL_NAME
        for count in (1, 5):
            m.MEANDER_STRAIGHT_LENGTHS = tuple(200.0 + 300.0 * i for i in range(count))
            m.TOP_CELL_NAME = f"QA_{key}_count_{count}"
            test = m.build_layout()
            check(
                f"{key}_supports_{count}_devices",
                test.info["device_count"] == count and len(test.info["device_cells"]) == count,
            )
        m.MEANDER_STRAIGHT_LENGTHS = ()
        rejected = False
        try:
            m.validate_settings()
        except ValueError:
            rejected = True
        check(key + "_empty_sweep_rejected", rejected)
        m.MEANDER_STRAIGHT_LENGTHS, m.TOP_CELL_NAME = saved_lengths, saved_name
    if capture:
        print("BASELINE_CAPTURED", flush=True)
        return
    REPORT["status"] = "PASS" if not REPORT["failures"] else "FAIL"
    (OUT / "validation_report.json").write_text(json.dumps(REPORT, indent=2))
    if REPORT["failures"]:
        raise AssertionError(", ".join(REPORT["failures"]))
    print("CODE_CLEANUP_VALIDATION_PASS", len(REPORT["checks"]), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", nargs="?", choices=("verify", "capture"), default="verify")
    args = parser.parse_args()
    run(capture=args.mode == "capture")
