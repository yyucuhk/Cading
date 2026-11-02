"""Compatibility wrapper for the shared continuous-180-degree coil factory.

New pattern code should import :func:`layout.devices.coiled.coiled` directly
and pass ``continuous_180=True``. This wrapper preserves the positional API of
the former local implementation for older Yue scripts.
"""

import sys
from pathlib import Path
from typing import Callable


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import gdsfactory as gf

from layout.core.layers import Layers
from layout.core.routing import bend_euler
from layout.devices.coiled import coiled as _shared_coiled
from layout.devices.coiled import width_profile


def coiled(
    width: "float | Callable",
    length: float,
    arm_length: float = 16e3,
    n_loops: "int | None" = None,
    arm_gap: float = 300,
    arm_spacing: float = 20.0,
    same_side: bool = False,
    turn_side: str = "right",
    bend=bend_euler,
    layer=Layers.wg,
) -> gf.Component:
    """Build the shared exact-length coil with continuous 180-degree turns."""
    return _shared_coiled(
        width=width,
        length=length,
        arm_length=arm_length,
        n_loops=n_loops,
        arm_gap=arm_gap,
        arm_spacing=arm_spacing,
        same_side=same_side,
        turn_side=turn_side,
        bend=bend,
        layer=layer,
        continuous_180=True,
    )


__all__ = ["coiled", "width_profile"]
