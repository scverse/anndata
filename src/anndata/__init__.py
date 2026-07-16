"""Annotated multivariate observation data."""

from __future__ import annotations

from importlib.metadata import version
from typing import TYPE_CHECKING

from packaging.version import Version

if TYPE_CHECKING:
    from typing import Any

import zarr

from ._core.anndata import AnnData
from ._core.extensions import register_anndata_namespace
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
from .utils import module_get_attr_redirect, warn

# Submodules need to be imported last
from . import abc, acc, experimental, io, types, typing  # isort: skip

# We use these in tests by attribute access
from . import logging  # noqa: F401  # isort: skip

# We are going to be "guinea pigs" for this new pipeline because it should be much faster
# and we're shortchanging our users otherwise.
# So we change the pipeline if it has not been changed by the user i.e.,
# it is the old BatchedCodecPipeline.
# This pipeline fully passes ours, zarr's, and zarr's downstream CI. - Ilan

old_pipeline = zarr.config.get("codec_pipeline.path")
if "Batched" in old_pipeline and Version(version("zarr")) >= Version("3.3"):
    zarr.config.set({
        "codec_pipeline.path": "zarr.core.codec_pipeline.FusedCodecPipeline",
        "codec_pipeline.max_workers": None,
    })

__all__ = [
    "AnnData",
    "ExperimentalFeatureWarning",
    "ImplicitModificationWarning",
    "OldFormatWarning",
    "Raw",
    "WriteWarning",
    "abc",
    "acc",
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
    if attr_name == "__version__":
        from importlib.metadata import version

        msg = "`__version__` is deprecated, use `importlib.metadata.version('anndata')` instead."
        warn(msg, FutureWarning)
        return version("anndata")

    return module_get_attr_redirect(attr_name, deprecated_mapping=_DEPRECATED)
