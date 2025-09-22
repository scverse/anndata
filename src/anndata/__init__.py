"""Annotated multivariate observation data."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any

from ._core.anndata import AnnData
from ._core.extensions import register_anndata_namespace
from ._core.merge import concat
from ._core.raw import Raw
from ._settings import settings
from ._version import __version__
from ._warnings import (
    ExperimentalFeatureWarning,
    ImplicitModificationWarning,
    OldFormatWarning,
    WriteWarning,
)
from .io import read_h5ad, read_zarr
from .utils import module_get_attr_redirect

# Submodules need to be imported last
from . import abc, experimental, typing, io, types  # isort: skip

# We use these in tests by attribute access
from . import logging  # noqa: F401  # isort: skip

_DEPRECATED_IO = (
    "read_loom",
    "read_hdf",
    "read_excel",
    "read_umi_tools",
    "read_csv",
    "read_text",
    "read_mtx",
)
_DEPRECATED = {method: f"io.{method}" for method in _DEPRECATED_IO}


def __getattr__(attr_name: str) -> Any:
    return module_get_attr_redirect(attr_name, deprecated_mapping=_DEPRECATED)


__all__ = [
    "AnnData",
    "ExperimentalFeatureWarning",
    "ImplicitModificationWarning",
    "OldFormatWarning",
    "Raw",
    "WriteWarning",
    "__version__",
    "abc",
    "concat",
    "experimental",
    "io",
    "read_h5ad",
    "read_zarr",
    "register_anndata_namespace",
    "settings",
    "types",
    "typing",
]
