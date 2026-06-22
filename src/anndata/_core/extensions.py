from __future__ import annotations

import sys
import warnings

from scverse_misc import make_register_namespace_decorator

from .anndata import AnnData

__all__ = ["register_anndata_namespace"]


if sys.version_info >= (3, 12):
    register_anndata_namespace = make_register_namespace_decorator(AnnData, "adata")
else:
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            r".*argument is deprecated and will be removed in the future",
            DeprecationWarning,
        )
        register_anndata_namespace = make_register_namespace_decorator(
            AnnData, "adata", "register_anndata_namespace", "numpy"
        )
