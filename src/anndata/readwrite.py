from __future__ import annotations

from warnings import warn

warn("Please only import from anndata, not anndata.readwrite", DeprecationWarning)

from ._io import *  # noqa: F403, E402
