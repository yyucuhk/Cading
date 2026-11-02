"""Generate an interleaved Euler-meander sweep from shared layout factories.

Edit the parameter block below to configure the mask. Reusable geometry lives
in :mod:`layout.meanders`.
"""

import os
import sys
from pathlib import Path


# Allow this file to run directly from any working directory.
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import gdsfactory as gf

from layout.meanders import (
    interleaved_euler_meander_bundle,
)
from layout.core.layers import Layers
from layout.core import GDS


gf.gpdk.PDK.activate()

# Sweep settings. Each tuple becomes one compact 10-device stack.
LENGTH_GROUPS = ((5000, 11000), (7000, 9000))
N_MEANDERS = 5
WIDTH_START = 0.5
WIDTH_STEP = 0.1
CUSTOM_WIDTHS = None
N_TURNS = 1
MINIMUM_BEND_RADIUS = 120
BEND_RADIUS_STEP = 10
FACET_STRAIGHT_LENGTH = 300
FACET_BEND_RADIUS = 75
ACCESS_X_STEP = 10
INPUT_STRAIGHT_LENGTH = 300
OUTPUT_STRAIGHT_LENGTH = 300
MINIMUM_WAVEGUIDE_GAP = 2
EULER_RATIO = 1.0
NPOINTS = 720
SETS_PER_ROW = 2
SET_COLUMN_GAP = 200
SET_ROW_GAP = 200
SET_X_OFFSET = 100


def build_component() -> gf.Component:
    """Build two compact stacks, each containing two lengths × five widths."""
    if not LENGTH_GROUPS or any(not group for group in LENGTH_GROUPS):
        raise ValueError("LENGTH_GROUPS must contain non-empty length tuples")
    if SETS_PER_ROW < 1:
        raise ValueError("SETS_PER_ROW must be positive")
    if SET_COLUMN_GAP < 0 or SET_ROW_GAP < 0:
        raise ValueError("set gaps must be non-negative")

    if CUSTOM_WIDTHS is None:
        base_widths = tuple(
            round(WIDTH_START + index * WIDTH_STEP, 12)
            for index in range(N_MEANDERS)
        )
    else:
        if len(CUSTOM_WIDTHS) != N_MEANDERS:
            raise ValueError("CUSTOM_WIDTHS must contain N_MEANDERS values")
        base_widths = tuple(float(width) for width in CUSTOM_WIDTHS)

    c = gf.Component()
    x_cursor = 0.0
    y_cursor = 0.0
    row_ymax = 0.0
    maximum_length_error = 0.0
    set_grid_positions: list[list[int]] = []

    for set_index, length_group in enumerate(LENGTH_GROUPS):
        if any(length <= 0 for length in length_group):
            raise ValueError("all grouped lengths must be positive")
        if any(
            later < former
            for former, later in zip(length_group, length_group[1:])
        ):
            raise ValueError("lengths within each group must be non-decreasing")

        row_index = set_index // SETS_PER_ROW
        column_index = set_index % SETS_PER_ROW
        if column_index == 0:
            x_cursor = row_index * SET_X_OFFSET
            row_ymax = y_cursor

        # Repeat the width sweep once for every length in this set. The ten
        # resulting waveguides receive ten complementary radius levels, so
        # they remain nested without sharing or crossing a physical path.
        expanded_widths = tuple(
            width for _length in length_group for width in base_widths
        )
        first_length = float(length_group[0])
        length_offsets = tuple(
            float(length) - first_length
            for length in length_group
            for _width in base_widths
        )
        device_stack = interleaved_euler_meander_bundle(
            target_length=first_length,
            n_meanders=len(expanded_widths),
            widths=expanded_widths,
            length_step=0,
            length_offsets=length_offsets,
            n_turns=N_TURNS,
            minimum_bend_radius=MINIMUM_BEND_RADIUS,
            bend_radius_step=BEND_RADIUS_STEP,
            facet_straight_length=FACET_STRAIGHT_LENGTH,
            facet_bend_radius=FACET_BEND_RADIUS,
            access_x_step=ACCESS_X_STEP,
            input_straight_length=INPUT_STRAIGHT_LENGTH,
            output_straight_length=OUTPUT_STRAIGHT_LENGTH,
            minimum_waveguide_gap=MINIMUM_WAVEGUIDE_GAP,
            euler_ratio=EULER_RATIO,
            npoints=NPOINTS,
            allow_duplicate_widths=True,
            add_labels=False,
            layer=Layers.wg,
            label_layer=Layers.label,
        )
        reference = c.add_ref(device_stack)
        reference.move(
            (
                x_cursor - float(reference.xmin),
                y_cursor - float(reference.ymin),
            )
        )
        set_grid_positions.append([row_index, column_index])
        maximum_length_error = max(
            maximum_length_error,
            float(device_stack.info["maximum_length_error"]),
        )

        for port in reference.ports:
            c.add_port(f"S{set_index + 1}_{port.name}", port=port)

        x_cursor = float(reference.xmax) + SET_COLUMN_GAP
        row_ymax = max(row_ymax, float(reference.ymax))
        row_is_complete = column_index == SETS_PER_ROW - 1
        final_set = set_index == len(LENGTH_GROUPS) - 1
        if row_is_complete or final_set:
            y_cursor = row_ymax + SET_ROW_GAP

    c.info["device_count"] = sum(
        len(group) * len(base_widths) for group in LENGTH_GROUPS
    )
    c.info["set_count"] = len(LENGTH_GROUPS)
    c.info["devices_per_set"] = [
        len(group) * len(base_widths) for group in LENGTH_GROUPS
    ]
    c.info["length_groups"] = [list(group) for group in LENGTH_GROUPS]
    c.info["widths"] = list(base_widths)
    c.info["set_grid_positions"] = set_grid_positions
    c.info["sets_per_row"] = SETS_PER_ROW
    c.info["maximum_length_error"] = maximum_length_error
    return c


if __name__ == "__main__":
    component = build_component()
    filepath = os.path.join(GDS, "interleaved_euler_meander_sweep.gds")
    component.write_gds(filepath)
    print(f"Wrote {filepath}")
    print(f"Devices: {component.info['device_count']}")
    print(f"Maximum length error: {component.info['maximum_length_error']:.6f} um")
    component.show()
