


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
from layout.core import GDS


gf.gpdk.PDK.activate()


if __name__ == "__main__":
    component = spiral_with_same_x_facets()
    filepath = os.path.join(GDS, "spiral_with_same_x_facets.gds")
    component.write_gds(filepath)
    print(f"Wrote {filepath}")
    print(f"Optical length: {component.info['length']:.3f} um")
    component.show()
