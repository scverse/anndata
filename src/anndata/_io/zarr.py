from __future__ import annotations

import os
import warnings
from contextlib import contextmanager, nullcontext
from importlib.util import find_spec
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import zarr
from scipy import sparse
from zarr.errors import GroupNotFoundError

from .._core.anndata import AnnData
from .._settings import settings
from .._warnings import OldFormatWarning
from ..compat import _clean_uns, _from_fixed_length_strings
from ..experimental import read_dispatched, write_dispatched
from ..utils import warn
from .specs import read_elem
from .utils import _read_legacy_raw, no_write_dataset_2d, report_read_key_on_error

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any

    from zarr.core.common import AccessModeLiteral
    from zarr.storage import StoreLike

    from .._types import _GroupStorageType
    from ..typing import RWAble

from importlib.metadata import version

from packaging.version import Version


@contextmanager
def fast_zarr_context(store: StoreLike):
    # We are going to be "guinea pigs" for this new pipeline because it should be much faster
    # and we're shortchanging our users otherwise.
    # So we change the pipeline if it has not been changed by the user i.e.,
    # it is the old BatchedCodecPipeline.
    # This pipeline fully passes ours, zarr's, and zarr's downstream CI. - Ilan
    old_pipeline = zarr.config.get("codec_pipeline.path")
    use_zarr_fused = "Batched" in old_pipeline and Version(version("zarr")) >= Version(
        "3.3"
    )
    use_zarrs = find_spec("zarrs") and isinstance(
        store.store if isintance(store, zarr.Group) else store,
        str | Path | zarr.storage.LocalStore,
    )
    if use_zarrs:
        context = zarr.config.set({"codec_pipeline.path": "zarrs.ZarrsCodecPipeline"})
    elif use_zarr_fused:
        context = zarr.config.set({
            "codec_pipeline.path": "zarr.core.codec_pipeline.FusedCodecPipeline",
            "codec_pipeline.max_workers": None,
        })
    else:
        context = nullcontext()
    with (
        context,
        warnings.catch_warnings() if use_zarrs else nullcontext(),
    ):
        if use_zarrs:
            warnings.filterwarnings(
                "ignore",
                message=r".*unsupported by ZarrsCodecPipeline.*",
                category=UserWarning,
            )
        yield


@no_write_dataset_2d
def write_zarr(
    store: StoreLike,
    adata: AnnData,
    *,
    chunks: tuple[int | None, ...] | None = None,
    convert_strings_to_categoricals: bool = True,
    consolidate_metadata: bool = True,
    **ds_kwargs,
) -> None:
    """See :meth:`~anndata.AnnData.write_zarr`."""
    if convert_strings_to_categoricals:
        adata.strings_to_categoricals()
        if adata.raw is not None:
            adata.strings_to_categoricals(adata.raw.var)

    def callback(
        write_func, store, elem_name: str, elem, *, dataset_kwargs, iospec
    ) -> None:
        if (
            chunks is not None
            and not isinstance(elem, sparse.sparray | sparse.spmatrix)
            and elem_name.lstrip("/") == "X"
        ):
            dataset_kwargs = dict(dataset_kwargs, chunks=chunks)
        write_func(store, elem_name, elem, dataset_kwargs=dataset_kwargs)

    with fast_zarr_context(store):
        # TODO: Use spec writing system for this
        f = open_write_group(store)
        f.attrs.setdefault("encoding-type", "anndata")
        f.attrs.setdefault("encoding-version", "0.1.0")

        write_dispatched(f, "/", adata, callback=callback, dataset_kwargs=ds_kwargs)
        if consolidate_metadata:
            with warnings.catch_warnings():
                # Consolidated metadata will soon be a zarr convention/spec and should be safe to write.
                # There is no sense in spamming our users about this.
                # See https://github.com/zarr-developers/zarr-specs/pull/373
                warnings.filterwarnings(
                    "ignore",
                    message=r".*[Cc]onsolidated metadata.*",
                    category=UserWarning,
                )
                zarr.consolidate_metadata(f.store)


# Suffixes of paths that hold a whole store in a single file.
# `zarr.open_group` treats such a path as a directory store and finds no group in
# it, so the user has to wrap it in the matching store class themselves.
_PACKED_STORE_CLASSES: Mapping[str, str] = MappingProxyType({
    ".zip": "zarr.storage.ZipStore"
})


def _add_packed_store_note(e: BaseException, store: StoreLike) -> None:
    """Suggest the store class to use if `store` looks like a single-file store."""
    if not isinstance(store, os.PathLike | str):
        return
    path = os.fspath(store)
    if (cls_name := _PACKED_STORE_CLASSES.get(Path(path).suffix)) is not None:
        e.add_note(f"Did you mean `read_zarr({cls_name}({path!r}))`?")


def read_zarr(store: StoreLike | zarr.Group) -> AnnData:
    """\
    Read from a hierarchical Zarr array store.

    Parameters
    ----------
    store
        The filename, a :class:`~typing.MutableMapping`, or a Zarr storage class.
    """

    def callback(func, elem_name: str, elem, iospec):
        """Read with handling for backwards compat"""
        if iospec.encoding_type == "anndata" or elem_name.endswith("/"):
            attrs: dict[str, Any] = {
                k: read_dispatched(v, callback)
                for k, v in dict(elem).items()
                if not k.startswith("raw.")
            }
            return AnnData(**attrs)
        elif elem_name.startswith("/raw."):
            return None
        elif elem_name in {"/obs", "/var"}:
            return read_dataframe(elem)
        elif elem_name == "/raw":
            # Backwards compat
            return _read_legacy_raw(f, func(elem), read_dataframe, func)
        return func(elem)

    with fast_zarr_context(store):
        if isinstance(store, zarr.Group):
            f = store
        else:
            try:
                f = zarr.open_group(store, mode="r")
            except GroupNotFoundError as e:
                _add_packed_store_note(e, store)
                raise
        adata = read_dispatched(f, callback=callback)
        if not isinstance(adata, AnnData):
            msg = f"Expected an AnnData at the store root, got {type(adata).__name__}"
            raise ValueError(msg)

        # Backwards compat (should figure out which version)
        if "raw.X" in f:
            raw = AnnData(**_read_legacy_raw(f, adata.raw, read_dataframe, read_elem))
            raw.obs_names = adata.obs_names
            adata.raw = raw

        # Backwards compat for <0.7
        if isinstance(f["obs"], zarr.Array):
            _clean_uns(adata)

    return adata


@report_read_key_on_error
def read_dataset(dataset: zarr.Array):
    """Legacy method for reading datasets without encoding_type."""
    value = dataset[...]
    if not isinstance(value, np.ndarray):  # scalars have no dtype to inspect
        return value
    elif isinstance(value.dtype, str):
        pass
    elif issubclass(value.dtype.type, np.bytes_):
        value = value.astype(str).astype(object)  # bytestring -> unicode -> str
    elif len(value.dtype.descr) > 1:  # Compound dtype
        # For backwards compat, now strings are written as variable length
        value = _from_fixed_length_strings(value)
    if value.shape == () and not np.isscalar(value):
        value = value[()]
    return value


@report_read_key_on_error
def read_dataframe_legacy(dataset: zarr.Array) -> pd.DataFrame:
    """Reads old format of dataframes"""
    # NOTE: Likely that categoricals need to be removed from uns
    msg = (
        f"{dataset.name!r} was written with a very old version of AnnData. "
        "Consider rewriting it."
    )
    warn(msg, OldFormatWarning)
    df = pd.DataFrame(_from_fixed_length_strings(dataset[()]))
    df.set_index(df.columns[0], inplace=True)
    return df


@report_read_key_on_error
def read_dataframe(group: zarr.Group | zarr.Array) -> RWAble:
    """Read `obs`/`var`, which is a `DataFrame` unless it was written as a plain mapping."""
    # Fast paths
    if isinstance(group, zarr.Array):
        return read_dataframe_legacy(group)
    else:
        return read_elem(group)


def open_write_group(
    store: StoreLike, *, mode: AccessModeLiteral = "w", **kwargs
) -> zarr.Group:
    if "zarr_format" not in kwargs:
        kwargs["zarr_format"] = settings.zarr_write_format
    return zarr.open_group(store, mode=mode, **kwargs)


def is_group_consolidated(group: _GroupStorageType, *, strict: bool = True) -> bool:
    if not isinstance(group, zarr.Group):
        if strict:
            msg = f"Expected zarr.Group, got {type(group)}"
            raise TypeError(msg)
        return False
    return group.metadata.consolidated_metadata is not None
