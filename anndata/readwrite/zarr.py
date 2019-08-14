from collections.abc import Mapping
from functools import _find_impl, singledispatch
from pathlib import Path
import warnings

import numpy as np
from scipy import sparse
import pandas as pd
from pandas.api.types import is_categorical_dtype, is_string_dtype

from ..core.anndata import AnnData, Raw
from ..compat import _from_fixed_length_strings, _clean_uns
from .utils import report_key_on_error
import zarr
import numcodecs

def write_zarr(store, adata, **dataset_kwargs):
    if isinstance(store, Path):
        store = str(store)
    adata.strings_to_categoricals()
    if adata.raw is not None:
        adata.strings_to_categoricals(adata.raw.var)
    f = zarr.open(store, mode='w')
    write_attribute(f, "X", adata.X, dataset_kwargs)
    write_attribute(f, "obs", adata.obs, dataset_kwargs)
    write_attribute(f, "var", adata.var, dataset_kwargs)
    write_attribute(f, "obsm", adata.obsm, dataset_kwargs)
    write_attribute(f, "varm", adata.varm, dataset_kwargs)
    write_attribute(f, "layers", adata.layers, dataset_kwargs)
    write_attribute(f, "uns", adata.uns, dataset_kwargs)
    write_attribute(f, "raw", adata.raw, dataset_kwargs)


def _write_method(cls):
    return _find_impl(cls, ZARR_WRITE_REGISTRY)


def write_attribute(f, key, value, dataset_kwargs):
    if key in f:
        del f[key]
    _write_method(type(value))(f, key, value, dataset_kwargs)

def write_mapping(f, key, value: Mapping, dataset_kwargs):
    for sub_k, sub_v in value.items():
        if not isinstance(key, str):
            warnings.warn(
                f'dict key {key} transformed to str upon writing to zarr, using '
                'string keys is recommended.'
            )
        write_attribute(f, f"{key}/{sub_k}", sub_v, dataset_kwargs)


def write_dataframe(z, k, df, dataset_kwargs):
    g = z.create_group(k)
    g.attrs["encoding-type"] = "dataframe"
    g.attrs["encoding-version"] = "0.1.0"
    g.attrs["column-order"] = list(df.columns)

    if df.index.name is not None:
        index_name = df.index.name
    else:
        index_name = "_index"
    g.attrs["_index"] = index_name
    write_series(g, index_name, df.index, dataset_kwargs)
    for colname, series in df.items():
        write_series(g, colname, series, dataset_kwargs)


def write_series(g, k, s, dataset_kwargs):
    if s.dtype == object:
        g.create_dataset(
            k,
            shape=s.shape,
            dtype=object,
            object_codec=numcodecs.VLenUTF8(),
            **dataset_kwargs
        )
        g[k][:] = s.values
    elif is_categorical_dtype(s):
        g.create_dataset(k, shape=s.shape, dtype=s.cat.codes.dtype)
        g[k][:] = s.cat.codes
        g[k].attrs["categories"] = list(s.cat.categories)
    else:
        g[k] = s.values


def write_not_implemented(f, key, value, dataset_kwargs={}):
    # If it's not an array, try and make it an array. If that fails, pickle it.
    # Maybe rethink that, maybe this should just pickle, and have explicit implementations for everything else
    raise NotImplementedError(
        f"Failed to write value for {key}, since a writer for type {type(value)}"
        f" has not been implemented yet."
    )

def write_array(f, key, value, dataset_kwargs):
    f.create_dataset(key, data=value, **dataset_kwargs)

# TODO: Not working quite right
def write_scalar(f, key, value, dataset_kwargs):
    f.create_dataset(key, data=np.array(value), **dataset_kwargs)

def write_none(f, key, value, dataset_kwargs={}):
    pass

# TODO: Figure out what to do with dataset_kwargs for these
def write_csr(f, key, value, dataset_kwargs={}):
    group = f.create_group(key)
    group.attrs["encoding-type"] = "csr_matrix"
    group.attrs["enocding-version"] = "0.1.0"
    group.attrs["shape"] = value.shape
    group["data"] = value.data
    group["indices"] = value.indices
    group["indptr"] = value.indptr

def write_csc(f, key, value, dataset_kwargs={}):
    group = f.create_group(key)
    group.attrs["encoding-type"] = "csc_matrix"
    group.attrs["enocding-version"] = "0.1.0"
    group.attrs["shape"] = value.shape
    group["data"] = value.data
    group["indices"] = value.indices
    group["indptr"] = value.indptr

def write_raw(f, key, value, dataset_kwargs={}):
    group = f.create_group(key)
    group.attrs["encoding-type"] = "raw"
    group.attrs["encoding-version"] = "0.1.0"
    group.attrs["shape"] = value.shape
    write_attribute(group, 'X', value.X, dataset_kwargs)
    write_attribute(group, 'var', value.var, dataset_kwargs)
    write_attribute(group, 'varm', value.varm, dataset_kwargs)

ZARR_WRITE_REGISTRY = {
    type(None): write_none,
    Mapping: write_mapping,
    object: write_not_implemented,
    np.ndarray: write_array,
    list: write_array,
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


def read_zarr(store):
    """Read from a hierarchical Zarr array store.

    Parameters
    ----------
    store
        The filename, a :class:`~typing.MutableMapping`, or a Zarr storage class.
    """
    if isinstance(store, Path):
        store = str(store)
    import zarr
    f = zarr.open(store, mode='r')
    d = {}
    for k in f.keys():
        # Backwards compat
        if k.startswith("raw."):
            continue
        if k in {"obs", "var"}:
            d[k] = read_dataframe(f[k])
        else:  # Base case
            d[k] = read_attribute(f[k])

    # Backwards compat
    raw = {}
    if "raw.var" in f:
        raw["var"] = read_dataframe(f["raw.var"])  # Backwards compat
    if "raw.varm" in f:
        raw["varm"] = read_attribute(f["raw.varm"])
    if "raw.X" in f:
        raw["X"] = read_attribute(f["raw.X"])
    if len(raw) > 0:
        assert "raw" not in d
        d["raw"] = raw

    _clean_uns(d)

    return AnnData(**d)


@singledispatch
def read_attribute(value):
    raise NotImplementedError()


@read_attribute.register(zarr.Array)
@report_key_on_error
def read_dataset(dataset):
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
@report_key_on_error
def read_group(group):
    if "encoding-type" in group.attrs:
        enctype = group.attrs["encoding-type"]
        if enctype == "dataframe":
            return read_dataframe(group)
        elif enctype == "csr_matrix":
            return read_csr(group)
        elif enctype == "csc_matrix":
            return read_csc(group)
    return {k: read_attribute(group[k]) for k in group.keys()}

@report_key_on_error
def read_csr(group):
    return sparse.csr_matrix(
        (group["data"], group["indices"], group["indptr"]),
        shape=group.attrs["shape"]
    )

@report_key_on_error
def read_csc(group):
    return sparse.csc_matrix(
        (group["data"], group["indices"], group["indptr"]),
        shape=group.attrs["shape"]
    )

@report_key_on_error
def read_dataframe_legacy(dataset: zarr.Array):
    """
    Reads old format of dataframes
    """
    # NOTE: Likely that categoricals need to be removed from uns
    df = pd.DataFrame(_from_fixed_length_strings(dataset[()]))
    df.set_index(df.columns[0], inplace=True)
    return df

@report_key_on_error
def read_dataframe(g):
    if isinstance(g, zarr.Array):
        return read_dataframe_legacy(g)
    df = pd.DataFrame({k: read_series(g[k]) for k in g.keys()})
    df.set_index(g.attrs["_index"], inplace=True)
    if g.attrs["_index"] == "_index":
        df.index.name = None
    if "column-order" in g.attrs:
        assert set(g.attrs["column-order"]) == set(df.columns)
        df = df[g.attrs["column-order"]]
    return df

@report_key_on_error
def read_series(d):
    if "categories" in d.attrs:
        return pd.Categorical.from_codes(
            d,
            d.attrs["categories"],
            ordered=False
        )
    else:
        return d
