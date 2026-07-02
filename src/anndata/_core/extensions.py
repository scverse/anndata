from __future__ import annotations

import warnings
from importlib.metadata import version

from packaging.version import Version
from scverse_misc import make_register_namespace_decorator

from .anndata import AnnData

__all__ = ["register_anndata_namespace"]


if Version(version("scverse-misc")) > Version("0.1"):
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
