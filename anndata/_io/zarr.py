from collections.abc import MutableMapping
from pathlib import Path
from typing import TypeVar, Union
from warnings import warn

import numpy as np
from scipy import sparse
import pandas as pd
import zarr

from .read import read_dispatched
from .write import write_dispatched

from .._core.anndata import AnnData
from ..compat import (
    _from_fixed_length_strings,
    _clean_uns,
)
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
    **dataset_kwargs,
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

    def dispatch_element(write_elem, group, key, elem):
        if key == "X" and chunks is not None and not isinstance(elem, sparse.spmatrix):
            write_elem(
                group, key, elem, dataset_kwargs=dict(chunks=chunks, **dataset_kwargs)
            )
        else:
            write_elem(group, key, elem, dataset_kwargs=dataset_kwargs)

    write_dispatched(f, adata, dispatch_element)


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

    f = zarr.open(store, mode="r")

    if "encoding-type" in f.attrs:
        return read_elem(f[""])

    # Backwards compat
    def dispatch_element(read_func, group, k, _):
        if k in {"obs", "var"}:
            return read_dataframe(group[k])
        return read_func(group[k])

    def dispatch_anndata_args(group, args):
        args["raw"] = _read_legacy_raw(
            group, args.get("raw"), read_dataframe, read_elem
        )

        if "X" in args:
            args["dtype"] = args["X"].dtype

        # Backwards compat to <0.7
        if isinstance(group["obs"], zarr.Array):
            _clean_uns(args)

        return args

    return read_dispatched(f, dispatch_element, dispatch_anndata_args)


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
