"""
Shared infrastructure imported by every other ``layout`` subpackage.

Importing this package installs the port helper methods (``extend``, ``bend_*``,
``flipped``, ...) onto ``gf.Port`` -- see :mod:`layout.core.ports`. ``ROOT`` and ``GDS``
are the repository root and the tape-out mask directory beneath it.
"""
from pathlib import Path

from layout.core.ports import install as _install_port_methods

_install_port_methods()

#: Repository root -- the parent of the ``layout`` package.
ROOT = Path(__file__).resolve().parents[2]

#: Where tape-out masks are written and re-imported across scripts.
GDS = ROOT / "gds"
