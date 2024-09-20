from __future__ import annotations

import warnings

from anndata import io


def __getattr__(key):
    warnings.warn(
        "Importing from anndata._io is deprecated and will be removed in a future release in favor of importing from anndata.io."
        f"If you need {key} and cannot find it among our documented public imports, please open an issue.",
        FutureWarning,
    )
    return getattr(io, key)
