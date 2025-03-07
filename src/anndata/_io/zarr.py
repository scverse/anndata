from __future__ import annotations

import asyncio
from pathlib import Path
from typing import TYPE_CHECKING, TypeVar
from warnings import warn

import numpy as np
import pandas as pd
import zarr
from scipy import sparse

from .._core.anndata import AnnData
from .._settings import settings
from .._warnings import OldFormatWarning
from ..compat import _clean_uns, _from_fixed_length_strings, is_zarr_v2
from ..experimental import read_dispatched_async, write_dispatched
from .specs import read_elem_async
from .specs.methods import sync_async_to_async
from .utils import _read_legacy_raw, report_read_key_on_error

if TYPE_CHECKING:
    from collections.abc import MutableMapping

    from zarr.core.common import AccessModeLiteral
    from zarr.storage import StoreLike

T = TypeVar("T")


def write_zarr(
    store: StoreLike,
    adata: AnnData,
    *,
    chunks: tuple[int, ...] | None = None,
    convert_strings_to_categoricals: bool = True,
    **ds_kwargs,
) -> None:
    if isinstance(store, Path):
        store = str(store)
    if convert_strings_to_categoricals:
        adata.strings_to_categoricals()
        if adata.raw is not None:
            adata.strings_to_categoricals(adata.raw.var)
    # TODO: Use spec writing system for this
    f = open_write_group(store)
    f.attrs.setdefault("encoding-type", "anndata")
    f.attrs.setdefault("encoding-version", "0.1.0")

    def callback(func, s, k: str, elem, dataset_kwargs, iospec):
        if (
            chunks is not None
            and not isinstance(elem, sparse.spmatrix)
            and k.lstrip("/") == "X"
        ):
            dataset_kwargs = dict(dataset_kwargs, chunks=chunks)
        func(s, k, elem, dataset_kwargs=dataset_kwargs)

    write_dispatched(f, "/", adata, callback=callback, dataset_kwargs=ds_kwargs)


def read_zarr(store: str | Path | MutableMapping | zarr.Group) -> AnnData:
    """\
    Read from a hierarchical Zarr array store.

    Parameters
    ----------
    store
        The filename, a :class:`~typing.MutableMapping`, or a Zarr storage class.
    """
    if isinstance(store, Path):
        store = str(store)

    if isinstance(store, zarr.Group):
        f = store
    else:
        f = zarr.open(store, mode="r")

    # Read with handling for backwards compat
    async def callback(func, elem_name: str, elem, iospec):
        if iospec.encoding_type == "anndata" or elem_name.endswith("/"):
            args = dict(
                await asyncio.gather(
                    *(
                        # This is covering up backwards compat in the anndata initializer
                        # In most cases we should be able to call `func(elen[k])` instead
                        sync_async_to_async(k, read_dispatched_async(elem[k], callback))
                        for k in elem.keys()
                        if not k.startswith("raw.")
                    )
                )
            )
            return AnnData(**args)
        elif elem_name.startswith("/raw."):
            return None
        elif elem_name in {"/obs", "/var"}:
            return await read_dataframe(elem)
        elif elem_name == "/raw":
            # Backwards compat
            return await _read_legacy_raw(
                f, await func(elem), read_dataframe, read_elem_async
            )
        return await func(elem)

    adata = asyncio.run(read_dispatched_async(f, callback=callback))

    # Backwards compat (should figure out which version)
    if "raw.X" in f:
        raw = AnnData(
            **asyncio.run(
                asyncio.gather(
                    _read_legacy_raw(f, adata.raw, read_dataframe, read_elem_async)
                )
            )
        )
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
    if not hasattr(value, "dtype"):
        return value
    elif isinstance(value.dtype, str):
        pass
    elif issubclass(value.dtype.type, np.bytes_):
        value = value.astype(str).astype(object)  # bytestring -> unicode -> str
    elif len(value.dtype.descr) > 1:  # Compound dtype
        # For backwards compat, now strings are written as variable length
        value = _from_fixed_length_strings(value)
    if value.shape == ():
        value = value[()]
    return value


@report_read_key_on_error
def read_dataframe_legacy(dataset: zarr.Array) -> pd.DataFrame:
    """Reads old format of dataframes"""
    # NOTE: Likely that categoricals need to be removed from uns
    warn(
        f"'{dataset.name}' was written with a very old version of AnnData. "
        "Consider rewriting it.",
        OldFormatWarning,
    )
    df = pd.DataFrame(_from_fixed_length_strings(dataset[()]))
    df.set_index(df.columns[0], inplace=True)
    return df


@report_read_key_on_error
async def read_dataframe(group: zarr.Group | zarr.Array) -> pd.DataFrame:
    # Fast paths
    if isinstance(group, zarr.Array):
        return read_dataframe_legacy(group)
    else:
        return await read_elem_async(group)


def open_write_group(
    store: StoreLike, *, mode: AccessModeLiteral = "w", **kwargs
) -> zarr.Group:
    if len({"zarr_version", "zarr_format"}.intersection(kwargs.keys())):
        msg = "Don’t specify `zarr_version` or `zarr_format` explicitly."
        raise ValueError(msg)
    kwargs["zarr_version" if is_zarr_v2() else "zarr_format"] = (
        settings.zarr_write_format
    )
    return zarr.open_group(store, mode=mode, **kwargs)
