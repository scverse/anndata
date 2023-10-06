from __future__ import annotations

from warnings import warn

warn("Please only import from anndata, not anndata.core", DeprecationWarning)

from ._core import *  # noqa: F403, E402
