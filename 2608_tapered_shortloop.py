import os
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))
from functools import partial
import numpy as np

import tqdm
import gdsfactory as gf

from layout.core import GDS
from layout.core.layers import Layers
from layout.core.markers import calipermarker_array

from layout.pdk import load_component
from layout.coupling.tapered import (
    tapered_fiber_inout,
    tapered_fiber_coupler,
    tapered_fiber_3d,
    tapered_fiber_linear,
)


### Setup and fiducials

name = __file__.split(os.sep)[-1].split(".")[0]
print(name)

gf.gpdk.PDK.activate()

fin = gf.Component(name)

# Chip outline
y = 0
chip_shape = (4e3, 10e3)
fin.add_ref(
    gf.components.rectangle(size=(chip_shape[0], chip_shape[1]), centered=True, layer=(0, 0))
).move((0, y))

# Alignment marks in case there needs to be an MLA write or something
fin.add_ref(
    calipermarker_array(
        chip_shape=(chip_shape[0], chip_shape[1]-1000),
        arm_length=50,
        layer1=Layers.wg,
        layer2=None,    # If this layer is filled in caliper verniers will also be drawn
    )
).move((0, y))
fin.add_ref(
    calipermarker_array(
        chip_shape=(chip_shape[0]-1000, chip_shape[1]-1000),
        arm_length=80,
        layer1=Layers.wg,
        layer2=None,    # If this layer is filled in caliper verniers will also be drawn
    )
).move((0, y))

fin.add_ref(
    gf.components.rectangle(
        size=(50, 40),
        centered=True,
        layer=(10, 0),
    )
).move((0, y))

### Tapered device setup
tapers = {
    "100:NIR:NL" : partial(load_component, "taperedfiber/ln100", "tf_1300+1550nm_coopt_f494_TE.json"),
    "100:VIS:NL" : partial(load_component, "taperedfiber/ln100", "tf_780nm_f206_TE.json"),
    "250:NIR:NL" : partial(load_component, "taperedfiber/ln250", "tf_1300+1550nm_coopt_f200_TE.json"),
    "100:NIR:L" : partial(tapered_fiber_linear, width1=.98, width2=.2, length=45, taper_length=10),
    "100:VIS:L" : partial(tapered_fiber_linear, width1=.5, width2=.2, length=70),
    "250:NIR:L" : partial(tapered_fiber_linear, width1=.68, width2=.2, length=65, taper_length=40),
}
couplers = {
    key: partial(
        tapered_fiber_coupler,
        tapered_fiber_spec=partial(
            tapered_fiber_3d,
            taper_spec=spec
        ),
        label=key,
    )
    for key, spec in tapers.items()
}

N = 2
# TODO: Reshma sweep
for y in np.arange(-N, N) * 500.:
    fin.add_ref(
        tapered_fiber_inout(
            chip_length=chip_shape[0],
            tapered_fiber_spec=couplers["250:NIR:L"],
        )
    ).move((0, y))


fin.write_gds(os.path.join(GDS, f"{name}.gds"))
