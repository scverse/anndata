from collections.abc import Mapping, MutableMapping
from functools import _find_impl, singledispatch
from pathlib import Path
from types import MappingProxyType
from typing import Callable, Type, TypeVar, Union
from warnings import warn

import numpy as np
from scipy import sparse
import pandas as pd
from pandas.api.types import is_categorical_dtype
import numcodecs
import zarr

from .._core.anndata import AnnData
from .._core.raw import Raw
from ..compat import (
    _from_fixed_length_strings,
    _to_fixed_length_strings,
    _clean_uns,
)
from .utils import (
    report_read_key_on_error,
    report_write_key_on_error,
    write_attribute,
    _read_legacy_raw,
    EncodingVersions,
)
from . import WriteWarning


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
    f = zarr.open(store, mode="w")
    if chunks is not None and not isinstance(adata.X, sparse.spmatrix):
        write_attribute(f, "X", adata.X, dict(chunks=chunks, **dataset_kwargs))
    else:
        write_attribute(f, "X", adata.X, dataset_kwargs)
    write_attribute(f, "obs", adata.obs, dataset_kwargs)
    write_attribute(f, "var", adata.var, dataset_kwargs)
    write_attribute(f, "obsm", adata.obsm, dataset_kwargs)
    write_attribute(f, "varm", adata.varm, dataset_kwargs)
    write_attribute(f, "obsp", adata.obsp, dataset_kwargs)
    write_attribute(f, "varp", adata.varp, dataset_kwargs)
    write_attribute(f, "layers", adata.layers, dataset_kwargs)
    write_attribute(f, "uns", adata.uns, dataset_kwargs)
    write_attribute(f, "raw", adata.raw, dataset_kwargs)


def _write_method(cls: Type[T]) -> Callable[[zarr.Group, str, T], None]:
    return _find_impl(cls, ZARR_WRITE_REGISTRY)


@write_attribute.register(zarr.Group)
def write_attribute_zarr(f, key, value, dataset_kwargs=MappingProxyType({})):
    if key in f:
        del f[key]
    _write_method(type(value))(f, key, value, dataset_kwargs)


def write_mapping(f, key, value: Mapping, dataset_kwargs=MappingProxyType({})):
    for sub_k, sub_v in value.items():
        if not isinstance(key, str):
            warn(
                f"dict key {key} transformed to str upon writing to zarr, using "
                "string keys is recommended.",
                WriteWarning,
            )
        write_attribute(f, f"{key}/{sub_k}", sub_v, dataset_kwargs)


@report_write_key_on_error
def write_dataframe(z, key, df, dataset_kwargs=MappingProxyType({})):
    # Check arguments
    for reserved in ("__categories", "_index"):
        if reserved in df.columns:
            raise ValueError(f"{reserved!r} is a reserved name for dataframe columns.")
    group = z.create_group(key)
    group.attrs["encoding-type"] = "dataframe"
    group.attrs["encoding-version"] = EncodingVersions.dataframe.value
    group.attrs["column-order"] = list(df.columns)

    if df.index.name is not None:
        index_name = df.index.name
    else:
        index_name = "_index"
    group.attrs["_index"] = index_name

    write_series(group, index_name, df.index, dataset_kwargs)
    for colname, series in df.items():
        write_series(group, colname, series, dataset_kwargs)


@report_write_key_on_error
def write_series(group, key, series, dataset_kwargs=MappingProxyType({})):
    if series.dtype == object:
        group.create_dataset(
            key,
            shape=series.shape,
            dtype=object,
            object_codec=numcodecs.VLenUTF8(),
            **dataset_kwargs,
        )
        group[key][:] = series.values
    elif is_categorical_dtype(series):
        # This should work for categorical Index and Series
        categorical: pd.Categorical = series.values
        categories: np.ndarray = categorical.categories.values
        codes: np.ndarray = categorical.codes
        category_key = f"__categories/{key}"

        write_array(group, category_key, categories, dataset_kwargs=dataset_kwargs)
        write_array(group, key, codes, dataset_kwargs=dataset_kwargs)

        group[key].attrs["categories"] = category_key
        # Must coerce np.bool_ to bool for json writing
        group[category_key].attrs["ordered"] = bool(categorical.ordered)
    else:
        group[key] = series.values


@report_write_key_on_error
def write_not_implemented(f, key, value, dataset_kwargs=MappingProxyType({})):
    # If it’s not an array, try and make it an array. If that fails, pickle it.
    # Maybe rethink that, maybe this should just pickle,
    # and have explicit implementations for everything else
    raise NotImplementedError(
        f"Failed to write value for {key}, since a writer for type {type(value)}"
        f" has not been implemented yet."
    )


@report_write_key_on_error
def write_list(g, key, value, dataset_kwargs=MappingProxyType({})):
    write_array(g, key, np.array(value), dataset_kwargs)


@report_write_key_on_error
def write_array(g, key, value, dataset_kwargs=MappingProxyType({})):
    if value.dtype == object:
        g.create_dataset(
            key,
            shape=value.shape,
            dtype=object,
            object_codec=numcodecs.VLenUTF8(),
            **dataset_kwargs,
        )
        g[key][:] = value
    elif value.dtype.kind == "V":
        # Structured dtype
        g.create_dataset(key, data=_to_fixed_length_strings(value), **dataset_kwargs)
    else:
        g.create_dataset(key, data=value, **dataset_kwargs)


# TODO: Not working quite right
@report_write_key_on_error
def write_scalar(f, key, value, dataset_kwargs=MappingProxyType({})):
    f.create_dataset(key, data=np.array(value), **dataset_kwargs)


@report_write_key_on_error
def write_none(f, key, value, dataset_kwargs=MappingProxyType({})):
    pass


# TODO: Figure out what to do with dataset_kwargs for these
@report_write_key_on_error
def write_csr(f, key, value, dataset_kwargs=MappingProxyType({})):
    group = f.create_group(key)
    group.attrs["encoding-type"] = "csr_matrix"
    group.attrs["encoding-version"] = EncodingVersions.csr_matrix.value
    group.attrs["shape"] = value.shape
    group["data"] = value.data
    group["indices"] = value.indices
    group["indptr"] = value.indptr


@report_write_key_on_error
def write_csc(f, key, value, dataset_kwargs=MappingProxyType({})):
    group = f.create_group(key)
    group.attrs["encoding-type"] = "csc_matrix"
    group.attrs["encoding-version"] = EncodingVersions.csc_matrix.value
    group.attrs["shape"] = value.shape
    group["data"] = value.data
    group["indices"] = value.indices
    group["indptr"] = value.indptr


def write_raw(f, key, value, dataset_kwargs=MappingProxyType({})):
    group = f.create_group(key)
    group.attrs["encoding-type"] = "raw"
    group.attrs["encoding-version"] = EncodingVersions.raw.value
    group.attrs["shape"] = value.shape
    write_attribute(group, "X", value.X, dataset_kwargs)
    write_attribute(group, "var", value.var, dataset_kwargs)
    write_attribute(group, "varm", value.varm, dataset_kwargs)


ZARR_WRITE_REGISTRY = {
    type(None): write_none,
    Mapping: write_mapping,
    object: write_not_implemented,
    np.ndarray: write_array,  # Possibly merge with write_series
    list: write_list,
    pd.DataFrame: write_dataframe,
    Raw: write_raw,
    # object: write_not_implemented,
    # h5py.Dataset: write_basic,
    # type(None): write_none,
    str: write_scalar,
    float: write_scalar,
    np.floating: write_scalar,
    bool: write_scalar,
    np.bool_: write_scalar,
    int: write_scalar,
    np.integer: write_scalar,
    sparse.csr_matrix: write_csr,
    sparse.csc_matrix: write_csc,
}


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
    d = {}
    for k in f.keys():
        # Backwards compat
        if k.startswith("raw."):
            continue
        if k in {"obs", "var"}:
            d[k] = read_dataframe(f[k])
        else:  # Base case
            d[k] = read_attribute(f[k])

    d["raw"] = _read_legacy_raw(f, d.get("raw"), read_dataframe, read_attribute)

    _clean_uns(d)

    return AnnData(**d)


@singledispatch
def read_attribute(value):
    raise NotImplementedError()


@read_attribute.register(zarr.Array)
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


@read_attribute.register(zarr.Group)
@report_read_key_on_error
def read_group(group: zarr.Group):
    if "encoding-type" in group.attrs:
        enctype = group.attrs["encoding-type"]
        EncodingVersions[enctype].check(group.name, group.attrs["encoding-version"])
        if enctype == "dataframe":
            return read_dataframe(group)
        elif enctype == "csr_matrix":
            return read_csr(group)
        elif enctype == "csc_matrix":
            return read_csc(group)
        # At the moment, just treat raw as normal group
    return {k: read_attribute(group[k]) for k in group.keys()}


@report_read_key_on_error
def read_csr(group: zarr.Group) -> sparse.csr_matrix:
    return sparse.csr_matrix(
        (group["data"], group["indices"], group["indptr"]), shape=group.attrs["shape"],
    )


@report_read_key_on_error
def read_csc(group: zarr.Group) -> sparse.csc_matrix:
    return sparse.csc_matrix(
        (group["data"], group["indices"], group["indptr"]), shape=group.attrs["shape"],
    )


@report_read_key_on_error
def read_dataframe_legacy(dataset: zarr.Array) -> pd.DataFrame:
    """Reads old format of dataframes"""
    # NOTE: Likely that categoricals need to be removed from uns
    df = pd.DataFrame(_from_fixed_length_strings(dataset[()]))
    df.set_index(df.columns[0], inplace=True)
    return df


@report_read_key_on_error
def read_dataframe(group) -> pd.DataFrame:
    if isinstance(group, zarr.Array):
        return read_dataframe_legacy(group)
    columns = list(group.attrs["column-order"])
    idx_key = group.attrs["_index"]
    df = pd.DataFrame(
        {k: read_series(group[k]) for k in columns},
        index=read_series(group[idx_key]),
        columns=list(columns),
    )
    if idx_key != "_index":
        df.index.name = idx_key
    return df


@report_read_key_on_error
def read_series(dataset: zarr.Array) -> Union[np.ndarray, pd.Categorical]:
    if "categories" in dataset.attrs:
        categories = dataset.attrs["categories"]
        if isinstance(categories, str):
            categories_key = categories
            parent_name = dataset.name.rstrip(dataset.basename)
            parent = zarr.open(dataset.store)[parent_name]
            categories_dset = parent[categories_key]
            categories = categories_dset[...]
            ordered = categories_dset.attrs.get("ordered", False)
        else:
            # TODO: remove this code at some point post 0.7
            # TODO: Add tests for this
            warn(
                f"Your file {str(dataset.file.name)!r} has invalid categorical "
                "encodings due to being written from a development version of "
                "AnnData. Rewrite the file ensure you can read it in the future.",
                FutureWarning,
            )
        return pd.Categorical.from_codes(dataset[...], categories, ordered=ordered)
    else:
        return dataset[...]
