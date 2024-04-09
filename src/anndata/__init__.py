"""Annotated multivariate observation data."""

from __future__ import annotations

try:  # See https://github.com/maresb/hatch-vcs-footgun-example
    from setuptools_scm import get_version

    __version__ = get_version(root="../..", relative_to=__file__)
except (ImportError, LookupError):
    try:
        from ._version import __version__
    except ModuleNotFoundError:
        raise RuntimeError(
            "anndata is not correctly installed. Please install it, e.g. with pip."
        )

# Allowing notes to be added to exceptions. See: https://github.com/scverse/anndata/issues/868
import sys

if sys.version_info < (3, 11):
    # Backport package for exception groups
    import exceptiongroup  # noqa: F401

from ._core.anndata import AnnData
from ._core.merge import concat
from ._core.raw import Raw
from ._io import (
    read_csv,
    read_excel,
    read_h5ad,
    read_hdf,
    read_loom,
    read_mtx,
    read_text,
    read_umi_tools,
    read_zarr,
)
from ._warnings import (
    ExperimentalFeatureWarning,
    ImplicitModificationWarning,
    OldFormatWarning,
    WriteWarning,
)

# Experimental needs to be imported last
from . import experimental  # isort: skip

# We use these in tests by attribute access
from . import _io, logging  # noqa: F401 isort: skip


def read(*args, **kwargs):
    import warnings

    warnings.warn(
        "`anndata.read` is deprecated, use `anndata.read_h5ad` instead. "
        "`ad.read` will be removed in mid 2024.",
        FutureWarning,
    )
    return read_h5ad(*args, **kwargs)


__all__ = [
    "__version__",
    "AnnData",
    "concat",
    "Raw",
    "read_h5ad",
    "read_loom",
    "read_hdf",
    "read_excel",
    "read_umi_tools",
    "read_csv",
    "read_text",
    "read_mtx",
    "read_zarr",
    "OldFormatWarning",
    "WriteWarning",
    "ImplicitModificationWarning",
    "ExperimentalFeatureWarning",
    "experimental",
]
