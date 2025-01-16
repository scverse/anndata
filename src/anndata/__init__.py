"""Annotated multivariate observation data."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any

try:  # See https://github.com/maresb/hatch-vcs-footgun-example
    from setuptools_scm import get_version

    __version__ = get_version(root="../..", relative_to=__file__)
except (ImportError, LookupError):
    try:
        from ._version import __version__
    except ModuleNotFoundError:
        msg = "anndata is not correctly installed. Please install it, e.g. with pip."
        raise RuntimeError(msg)

# Allowing notes to be added to exceptions. See: https://github.com/scverse/anndata/issues/868
import sys

if sys.version_info < (3, 11):
    # Backport package for exception groups
    import exceptiongroup  # noqa: F401

from ._core.anndata import AnnData
from ._core.merge import concat
from ._core.raw import Raw
from ._settings import settings
from ._warnings import (
    ExperimentalFeatureWarning,
    ImplicitModificationWarning,
    OldFormatWarning,
    WriteWarning,
)
from .io import read_h5ad, read_zarr
from .utils import module_get_attr_redirect

# Submodules need to be imported last
from . import abc, experimental, typing, io  # noqa: E402 isort: skip

# We use these in tests by attribute access
from . import logging  # noqa: F401, E402 isort: skip


def read(*args, **kwargs):
    import warnings

    warnings.warn(
        "`anndata.read` is deprecated, use `anndata.read_h5ad` instead. "
        "`ad.read` will be removed in mid 2024.",
        FutureWarning,
    )
    return read_h5ad(*args, **kwargs)


_DEPRECATED_IO = (
    "read_loom",
    "read_hdf",
    "read_excel",
    "read_umi_tools",
    "read_csv",
    "read_text",
    "read_mtx",
)
_DEPRECATED = dict((method, f"io.{method}") for method in _DEPRECATED_IO)


def __getattr__(attr_name: str) -> Any:
    return module_get_attr_redirect(attr_name, deprecated_mapping=_DEPRECATED)


__all__ = [
    # Attributes
    "__version__",
    "settings",
    # Submodules
    "abc",
    "experimental",
    "typing",
    "io",
    # Classes
    "AnnData",
    "Raw",
    # Functions
    "concat",
    "read_zarr",
    "read_h5ad",
    "read",
    # Warnings
    "OldFormatWarning",
    "WriteWarning",
    "ImplicitModificationWarning",
    "ExperimentalFeatureWarning",
]
