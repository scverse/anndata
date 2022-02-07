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
    if chunks is not None and not isinstance(adata.X, sparse.spmatrix):
        write_elem(
            f, "X", adata.X, dataset_kwargs=dict(chunks=chunks, **dataset_kwargs)
        )
    else:
        write_elem(f, "X", adata.X, dataset_kwargs=dataset_kwargs)
    write_elem(f, "obs", adata.obs, dataset_kwargs=dataset_kwargs)
    write_elem(f, "var", adata.var, dataset_kwargs=dataset_kwargs)
    write_elem(f, "obsm", dict(adata.obsm), dataset_kwargs=dataset_kwargs)
    write_elem(f, "varm", dict(adata.varm), dataset_kwargs=dataset_kwargs)
    write_elem(f, "obsp", dict(adata.obsp), dataset_kwargs=dataset_kwargs)
    write_elem(f, "varp", dict(adata.varp), dataset_kwargs=dataset_kwargs)
    write_elem(f, "layers", dict(adata.layers), dataset_kwargs=dataset_kwargs)
    write_elem(f, "uns", dict(adata.uns), dataset_kwargs=dataset_kwargs)
    write_elem(f, "raw", adata.raw, dataset_kwargs=dataset_kwargs)


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
    d = {}
    for k in f.keys():
        if k.startswith("raw."):
            continue
        if k in {"obs", "var"}:
            d[k] = read_dataframe(f[k])
        else:  # Base case
            d[k] = read_elem(f[k])

    d["raw"] = _read_legacy_raw(f, d.get("raw"), read_dataframe, read_elem)

    if "X" in d:
        d["dtype"] = d["X"].dtype

    _clean_uns(d)

    return AnnData(**d)


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
