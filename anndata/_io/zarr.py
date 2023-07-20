from collections.abc import MutableMapping
from pathlib import Path
from typing import TypeVar, Union
from warnings import warn

import numpy as np
from scipy import sparse
import pandas as pd
import zarr

from .._core.anndata import AnnData
from ..compat import (
    _from_fixed_length_strings,
    _clean_uns,
)
from ..experimental import read_dispatched, write_dispatched
from .utils import (
    report_read_key_on_error,
    _read_legacy_raw,
)
from .specs import read_elem, write_elem
from anndata._warnings import OldFormatWarning


T = TypeVar("T")


def write_zarr(
    store: Union[MutableMapping, str, Path],
    adata: AnnData,
    chunks=None,
    **ds_kwargs,
) -> None:
    if isinstance(store, Path):
        store = str(store)
    adata.strings_to_categoricals()
    if adata.raw is not None:
        adata.strings_to_categoricals(adata.raw.var)
    # TODO: Use spec writing system for this
    f = zarr.open(store, mode="w")
    f.attrs.setdefault("encoding-type", "anndata")
    f.attrs.setdefault("encoding-version", "0.1.0")

    def callback(func, s, k, elem, dataset_kwargs, iospec):
        if chunks is not None and not isinstance(elem, sparse.spmatrix) and k == "/X":
            func(s, k, elem, dataset_kwargs=dict(chunks=chunks, **dataset_kwargs))
        else:
            func(s, k, elem, dataset_kwargs=dataset_kwargs)

    write_dispatched(f, "/", adata, callback=callback, dataset_kwargs=ds_kwargs)


def read_zarr(store: Union[str, Path, MutableMapping, zarr.Group]) -> AnnData:
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
    def callback(func, elem_name: str, elem, iospec):
        if iospec.encoding_type == "anndata" or elem_name.endswith("/"):
            return AnnData(
                **{
                    k: read_dispatched(v, callback)
                    for k, v in elem.items()
                    if not k.startswith("raw.")
                }
            )
        elif elem_name.startswith("/raw."):
            return None
        elif elem_name in {"/obs", "/var"}:
            return read_dataframe(elem)
        elif elem_name == "/raw":
            # Backwards compat
            return _read_legacy_raw(f, func(elem), read_dataframe, func)
        return func(elem)

    adata = read_dispatched(f, callback=callback)

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
    value = dataset[...]
    if not hasattr(value, "dtype"):
        return value
    elif isinstance(value.dtype, str):
        pass
    elif issubclass(value.dtype.type, np.str_):
        value = value.astype(object)
    elif issubclass(value.dtype.type, np.string_):
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
def read_dataframe(group) -> pd.DataFrame:
    from .specs import _REGISTRY

    # Fast paths
    if isinstance(group, zarr.Array):
        return read_dataframe_legacy(group)
    else:
        return read_elem(group)
