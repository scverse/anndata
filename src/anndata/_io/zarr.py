from __future__ import annotations

import asyncio
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, TypeVar
from warnings import warn

import anyio
import numpy as np
import pandas as pd
import zarr
from scipy import sparse

from .._core.anndata import AnnData
from .._settings import settings
from .._warnings import OldFormatWarning
from ..compat import ZarrAsyncArray, _clean_uns, _from_fixed_length_strings, is_zarr_v2
from ..experimental import read_dispatched, write_dispatched
from .specs import read_elem_async
from .utils import _read_legacy_raw, contains, get, items, report_read_key_on_error

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

    async def callback(func, s, k: str, elem, dataset_kwargs, iospec):
        if (
            chunks is not None
            and not isinstance(elem, sparse.spmatrix)
            and k.lstrip("/") == "X"
        ):
            dataset_kwargs = dict(dataset_kwargs, chunks=chunks)
        await func(s, k, elem, dataset_kwargs=dataset_kwargs)

    anyio.run(
        partial(write_dispatched, dataset_kwargs=ds_kwargs), f, "/", adata, callback
    )
    if is_zarr_v2():
        zarr.convenience.consolidate_metadata(f.store)
    else:
        zarr.consolidate_metadata(f.store)


def read_zarr(store: str | Path | MutableMapping | zarr.Group) -> AnnData:
    return anyio.run(read_zarr_async, store)


async def read_zarr_async(store: str | Path | MutableMapping | zarr.Group) -> AnnData:
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
        if is_zarr_v2():
            f = zarr.open(store, mode="r")
        else:
            f = await zarr.api.asynchronous.open(store=store, mode="r")

    # Read with handling for backwards compat
    async def callback(func, elem_name: str, elem, iospec):
        if iospec.encoding_type == "anndata" or elem_name.endswith("/"):
            k_v = [(k, v) for k, v in await items(elem) if not k.startswith("raw.")]
            args = dict(
                zip(
                    (k for k, _ in k_v),
                    await asyncio.gather(
                        *(
                            # This is covering up backwards compat in the anndata initializer
                            # In most cases we should be able to call `func(elen[k])` instead
                            read_dispatched(v, callback)
                            for _, v in k_v
                        )
                    ),
                )
            )
            return AnnData(**args)
        elif elem_name.startswith("/raw."):
            return None
        elif elem_name in {"/obs", "/var"}:
            return await read_dataframe(elem)
        elif elem_name == "/raw":
            # Backwards compat
            elem_in_memory = await func(elem)
            return await _read_legacy_raw(
                f, elem_in_memory, read_dataframe, read_elem_async
            )
        return await func(elem)

    adata = await read_dispatched(f, callback=callback)
    # Backwards compat (should figure out which version)
    if await contains(f, "raw.X"):
        raw_args = await _read_legacy_raw(f, adata.raw, read_dataframe, read_elem_async)
        raw = AnnData(**raw_args)
        raw.obs_names = adata.obs_names
        adata.raw = raw
    # Backwards compat for <0.7
    if isinstance(await get(f, "obs"), zarr.AsyncArray):
        _clean_uns(adata)
    return adata


@report_read_key_on_error
async def read_dataset(dataset: zarr.Array):
    """Legacy method for reading datasets without encoding_type."""
    if is_zarr_v2():
        value = dataset[...]
    else:
        if isinstance(dataset, ZarrAsyncArray):
            value = await dataset.getitem(())
        else:
            value = dataset[()]
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
async def read_dataframe_legacy(dataset: zarr.Array) -> pd.DataFrame:
    """Reads old format of dataframes"""
    # NOTE: Likely that categoricals need to be removed from uns
    warn(
        f"'{dataset.name}' was written with a very old version of AnnData. "
        "Consider rewriting it.",
        OldFormatWarning,
    )
    if is_zarr_v2():
        data = dataset[...]
    else:
        data = await dataset._async_array.getitem(())
    df = pd.DataFrame(_from_fixed_length_strings(data))
    df.set_index(df.columns[0], inplace=True)
    return df


@report_read_key_on_error
async def read_dataframe(group: zarr.Group | zarr.Array) -> pd.DataFrame:
    # Fast paths
    if isinstance(group, zarr.Array):
        return await read_dataframe_legacy(group)
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
