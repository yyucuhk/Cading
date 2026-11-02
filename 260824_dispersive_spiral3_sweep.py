

import os
import sys
from pathlib import Path

# Direct execution adds only ``tapeout/yue`` to sys.path. Add the repository
# root so the shared ``layout`` package can be imported without installation.
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import gdsfactory as gf


from layout.devices.coiled import (
    coiled_same_x_facets as spiral_with_same_x_facets,
)
from layout.core.layers import Layers
from layout.core import GDS


gf.gpdk.PDK.activate()


@gf.cell
def spiral_sweep(
    lengths: tuple[float, ...] = (5000, 7000, 9000, 11000),
    widths: tuple[float, ...] = ( 0.6, 0.7, 0.8, 0.9),
    n_loops: int = 1,
    min_bend_radius: float = 120,
    radius_step: float = 5,
    end_extension_length: float = 200,
    facet_bend_radius: float = 75,
    input_straight_length: float = 1500,
    output_straight_length: float = 1500,
    x_pitch: float = 2000,
    y_pitch: float = 2500,
    second_row_x_offset: float = 1000,
    layer=Layers.wg,
    label_layer=Layers.label,
) -> gf.Component:
    """Return a labeled sweep of fixed-length spiral waveguides.

    Consecutive width pairs share one horizontal row. For each target length,
    the two widths in that pair are placed next to each other along x.
    The target lengths exclude the terminal input and output straights.
    """
    if not lengths or not widths:
        raise ValueError("lengths and widths must each contain at least one value")
    if len(widths) % 2:
        raise ValueError("widths must contain an even number of values")

    c = gf.Component()
    number_of_rows = len(widths) // 2
    for row in range(number_of_rows):
        width_pair = widths[2 * row : 2 * row + 2]
        y = (number_of_rows - 1 - row) * y_pitch
        row_x_offset = second_row_x_offset if row == 1 else 0
        for length_column, length in enumerate(lengths):
            for pair_column, width in enumerate(width_pair):
                column = 2 * length_column + pair_column
                x = column * x_pitch + row_x_offset
                xs = gf.cross_section.strip(width=width, layer=layer)
                bend90_length = float(
                    gf.components.bend_euler(
                        radius=facet_bend_radius,
                        angle=90,
                        cross_section=xs,
                    ).info["length"]
                )
                # The requested length excludes the two terminal facet straights.
                # Those straights remain in the GDS but are not counted in L5/L7/etc.
                counted_external_length = (
                    2 * end_extension_length + 2 * bend90_length
                )
                spiral_length = length - counted_external_length
                if spiral_length <= 0:
                    raise ValueError(
                        f"total length {length:g} is shorter than its "
                        f"counted external routing length {counted_external_length:g}"
                    )
                device = spiral_with_same_x_facets(
                    width=width,
                    spiral_length=spiral_length,
                    n_loops=n_loops,
                    min_bend_radius=min_bend_radius,
                    radius_step=radius_step,
                    end_extension_length=end_extension_length,
                    facet_bend_radius=facet_bend_radius,
                    input_straight_length=input_straight_length,
                    output_straight_length=output_straight_length,
                    layer=layer,
                )
                generated_counted_length = float(
                    device.info["length_excluding_input_output_straights"]
                )
                if abs(generated_counted_length - length) > 1e-3:
                    raise ValueError(
                        f"generated counted length {generated_counted_length} does not "
                        f"match target {length}"
                    )
                reference = c.add_ref(device)
                reference.move((x, y))

                label = gf.components.text(
                    text=f"L{length / 1000:g}mm_W{width:.1f}um",
                    size=35,
                    layer=label_layer,
                )
                label_reference = c.add_ref(label)
                label_reference.move((x - 350, y + 330))

    c.info["device_count"] = len(lengths) * len(widths)
    c.info["n_loops"] = n_loops
    c.info["lengths"] = list(lengths)
    c.info["widths"] = list(widths)
    c.info["row_parameter"] = "width_pair"
    c.info["column_parameter"] = "length_then_width_pair_member"
    c.info["x_pitch"] = x_pitch
    c.info["y_pitch"] = y_pitch
    c.info["second_row_x_offset"] = second_row_x_offset
    c.info["lengths_exclude_input_output_straights"] = True
    return c


if __name__ == "__main__":
    component = spiral_sweep()
    filepath = os.path.join(
        GDS, f"spiral_sweep_{component.info['device_count']}_devices.gds"
    )
    component.write_gds(filepath)
    print(f"Wrote {filepath}")
    print(f"Devices: {component.info['device_count']}")
    component.show()
