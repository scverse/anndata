from pathlib import Path
from typing import Union, Tuple

import h5py as _h5py
import numpy as np
import pandas as pd
from pandas.api.types import (
    is_string_dtype,
    is_categorical_dtype
)
from scipy import sparse
from .. import h5py
from ..core.anndata import AnnData, Raw
from collections.abc import Mapping


def make_h5_cat_dtype(s):
    """This creates a hdf5 enum dtype based on a pandas categorical array."""
    # If this causes segfaults, try the alternative in the commit message
    max_code = np.max(s.codes)
    return h5py.special_dtype(enum=(
        s.codes.dtype,
        dict(zip(s.categories, np.arange(max_code + 1)))
    ))


def to_h5_dtype(a: Union[np.ndarray, pd.Series]) -> Tuple[np.ndarray, np.dtype]:
    """Given an ndarray, return tuple of hdf5 friendly array and dtype."""
    if isinstance(a, pd.Series):
        a = a.values
    dtype = a.dtype
    if is_categorical_dtype(dtype):
        return a.codes, make_h5_cat_dtype(a)
    elif is_string_dtype(dtype):
        return a, h5py.special_dtype(vlen=str)
    else:
        return a, a.dtype


def df_to_h5_recarray(df: pd.DataFrame) -> np.recarray:
    if df.index.name is None:
        names = ["index"] + list(df.columns)
    else:
        names = [df.index.name] + list(df.columns)
    arrays, formats = zip(*(to_h5_dtype(x) for x in [df.index] + [df[k] for k in df]))
    return np.rec.fromarrays(
        arrays,
        dtype={'names': names, 'formats': formats}
    )


def _to_hdf5_vlen_strings(value):
    """
    This corrects compound dtypes to work with hdf5 files.
    """
    new_dtype = []
    for dt_name, dt_type in value.dtype.descr:
        if dt_type[1] == "U":
            new_dtype.append(
                (dt_name, h5py.special_dtype(vlen=str))
            )
        else:
            new_dtype.append((dt_name, dt_type))
    return value.astype(new_dtype)


def _from_fixed_length_strings(value):
    """Convert from fixed length strings to unicode.

    For backwards compatability with older files.
    """
    new_dtype = []
    for dt in value.dtype.descr:
        dt_list = list(dt)
        dt_type = dt[1]
        if isinstance(dt_type, tuple):  # vlen strings, could probably match better
            dt_list[1] = "O"
            new_dtype.append(tuple(dt_list))
        elif issubclass(np.dtype(dt_type).type, np.string_):
            dt_list[1] = 'U{}'.format(int(dt_type[2:]))
            new_dtype.append(tuple(dt_list))
        else:
            new_dtype.append(dt)
    return value.astype(new_dtype)


def write_h5ad(filepath: Union[Path, str], adata: "AnnData", force_dense=False, dataset_kwargs={}, **kwargs):
    """Write ``.h5ad``-formatted hdf5 file.

    .. note::

        Setting compression to ``'gzip'`` can save disk space but
        will slow down writing and subsequent reading. Prior to
        v0.6.16, this was the default for parameter
        ``compression``.

    Generally, if you have sparse data that are stored as a dense
    matrix, you can dramatically improve performance and reduce
    disk space by converting to a :class:`~scipy.sparse.csr_matrix`::

        from scipy.sparse import csr_matrix
        adata.X = csr_matrix(adata.X)

    Parameters
    ----------
    adata
        AnnData object to write.
    filename
        Filename of data file. Defaults to backing file.
    compression : ``None``,  {``'gzip'``, ``'lzf'``} (default: ``None``)
        See the h5py :ref:`dataset_compression`.
    force_dense
        Write sparse data as a dense matrix.
    """
    adata.strings_to_categoricals()
    if adata.raw is not None:
        adata.strings_to_categoricals(adata.raw.var)
    dataset_kwargs = dataset_kwargs.copy()
    dataset_kwargs.update(kwargs)
    filepath = Path(filepath)
    mode = "a" if adata.isbacked else "w"
    if adata.isbacked:  # close so that we can reopen below
        adata.file.close()
    with h5py.File(filepath, mode, force_dense=force_dense) as f:
        if not (adata.isbacked and Path(adata.filename) == Path(filepath)):
            # Otherwise, X should already be up to date
            write_attribute(f, "X", adata.X, dataset_kwargs)
        write_attribute(f, "obs", adata.obs, dataset_kwargs)
        write_attribute(f, "var", adata.var, dataset_kwargs)
        write_attribute(f, "obsm", adata.obsm, dataset_kwargs)
        write_attribute(f, "varm", adata.varm, dataset_kwargs)
        write_attribute(f, "layers", adata.layers, dataset_kwargs)
        write_attribute(f, "uns", adata.uns, dataset_kwargs)
        write_attribute(f, "raw", adata.raw, dataset_kwargs)
    if adata.isbacked:
        adata.file.open(filepath, 'r+')


from functools import _find_impl, singledispatch


def _write_method(cls):
    return _find_impl(cls, H5AD_WRITE_REGISTRY)


def write_attribute(f, key, value, dataset_kwargs):
    if key in f:
        del f[key]
    _write_method(type(value))(f, key, value, dataset_kwargs)


def write_raw(f, key, value, dataset_kwargs={}):
    write_attribute(f, 'raw.X', value.X, dataset_kwargs)
    write_attribute(f, 'raw.var', value.var, dataset_kwargs)
    write_attribute(f, 'raw.varm', value.varm, dataset_kwargs)


def write_not_implemented(f, key, value, dataset_kwargs={}):
    # If it's not an array, try and make it an array. If that fails, pickle it.
    # Maybe rethink that, maybe this should just pickle, and have explicit implementations for everything else
    raise NotImplementedError(
        f"Failed to write value for {key}, since a writer for type {type(value)}"
        f" has not been implemented yet."
    )


def write_basic(f, key, value, dataset_kwargs={}):
    f.create_dataset(key, value, **dataset_kwargs)


def write_list(f, key, value, dataset_kwargs={}):
    write_array(f, key, np.array(value), dataset_kwargs)


def write_none(f, key, value, dataset_kwargs={}):
    pass

def write_scalar(f, key, value, dataset_kwargs={}):
    if "compression" in dataset_kwargs:  # Can't compress scalars, error is thrown
        dataset_kwargs = dataset_kwargs.copy()
        dataset_kwargs.pop("compression")
    write_array(f, key, np.array(value), dataset_kwargs)

def write_array(f, key, value, dataset_kwargs={}):
    # Convert unicode to fixed length strings
    if value.dtype.kind == 'U':
        value = value.astype(h5py.special_dtype(vlen=str))
    elif value.dtype.names is not None:
        value = _to_hdf5_vlen_strings(value)
    f.create_dataset(key, data=value, **dataset_kwargs)


def write_dataframe(f, key, value, dataset_kwargs={}):
    f.create_dataset(key, data=df_to_h5_recarray(value), **dataset_kwargs)
    f[key].attrs["source_type"] = "dataframe"


def write_mapping(f, key, value, dataset_kwargs={}):
    for sub_key, sub_value in value.items():
        write_attribute(f, f"{key}/{sub_key}", sub_value, dataset_kwargs)


H5AD_WRITE_REGISTRY = {
    Raw: write_raw,
    object: write_not_implemented,
    h5py.Dataset: write_basic,
    list: write_list,
    type(None): write_none,
    str: write_scalar,
    float: write_scalar,
    np.floating: write_scalar,
    bool: write_scalar,
    np.bool_: write_scalar,
    int: write_scalar,
    np.integer: write_scalar,
    np.ndarray: write_array,
    sparse.spmatrix: write_basic,
    h5py.SparseDataset: write_basic,
    pd.DataFrame: write_dataframe,
    Mapping: write_mapping
}


def read_h5ad_backed(filename, mode):
    d = {"filename": filename, "filemode": mode}
    raw = {}

    f = h5py.File(filename, mode)

    attributes = ["obsm", "varm", "uns", "layers"]
    df_attributes = ["obs", "var"]

    d.update({k: read_attribute(f[k]) for k in attributes if k in f})
    for k in df_attributes:
        if k in f:
            d[k] = read_dataframe(f[k])

    if "raw.var" in f:
        raw["var"] = read_dataframe(f["raw.var"])  # Backwards compat
    if "raw.varm" in f:
        raw["varm"] = read_attribute(f["raw.varm"])

    if len(raw) > 0:
        d["raw"] = raw
    if f.get("X", None) is not None:
        d["dtype"] = f["X"].dtype

    _clean_uns(d)

    return AnnData(**d)


# TODO: chunks
def read_h5ad(filename, backed=None, chunk_size=None):
    """Read ``.h5ad``-formatted hdf5 file.

    Parameters
    ----------
    filename
        File name of data file.
    backed : {``None``, ``'r'``, ``'r+'``}
        If ``'r'``, load :class:`~anndata.AnnData` in ``backed`` mode instead
        of fully loading it into memory (`memory` mode). If you want to modify
        backed attributes of the AnnData object, you need to choose ``'r+'``.
    chunk_size
        Used only when loading sparse dataset that is stored as dense.
        Loading iterates through chunks of the dataset of this row size
        until it reads the whole dataset.
        Higher size means higher memory consumption and higher loading speed.
    """
    if backed not in {None, False}:
        mode = backed
        if mode is True:
            mode = "r+"
        assert mode in {"r", "r+"}
        return read_h5ad_backed(filename, mode)

    attributes = ["X", "obsm", "varm", "uns", "layers"]
    df_attributes = ["obs", "var"]
    raw_attributes = ["raw.X", "raw.varm"]

    with h5py.File(filename, "r") as f:
        d = {k: read_attribute(f[k]) for k in attributes if k in f}
        raw = {k.replace("raw.", ""): read_attribute(f[k]) for k in raw_attributes if k in f}
        for k in df_attributes:
            if k in f:
                d[k] = read_dataframe(f[k])
        if "raw.var" in f:
            raw["var"] = read_dataframe(f["raw.var"])

    if len(raw) > 0:
        d["raw"] = raw
    if d.get("X", None) is not None:
        d["dtype"] = d["X"].dtype

    _clean_uns(d)  # backwards compat

    return AnnData(**d)


def _clean_uns(d):
    """Compat function for when categorical keys were stored in uns."""
    k_to_delete = []
    for k, v in d.get("uns", {}).items():
        if k.endswith('_categories'):
            k_stripped = k.replace('_categories', '')
            if isinstance(v, (str, int)):  # fix categories with a single category
                v = [v]
            for ann in ['obs', 'var']:
                if k_stripped in d[ann]:
                    d[ann][k_stripped] = pd.Categorical.from_codes(
                        codes=d[ann][k_stripped].values,
                        categories=v,
                    )
                    k_to_delete.append(k)
    for k in k_to_delete:
        del d["uns"][k]


@singledispatch
def read_attribute(value):
    raise NotImplementedError()


def read_dataframe(dataset):
    df = pd.DataFrame(_from_fixed_length_strings(dataset[()]))
    # for k, dtype in dataset.dtype.descr:
    #     if issubclass(np.dtype(dtype).type, np.string_):
    #         df[k] = df[k].astype(str)
    df.set_index(df.columns[0], inplace=True)
    dt = dataset.dtype
    for col in dt.names[1:]:
        check = _h5py.check_dtype(enum=dt[col])
        if not isinstance(check, dict):
            continue
        mapper = {v: k for k, v in check.items()}
        df[col] = pd.Categorical(df[col].map(mapper), ordered=False)
    return df


@read_attribute.register(h5py.Group)
def read_group(group: h5py.Group):
    d = dict()
    for sub_key, sub_value in group.items():
        d[sub_key] = read_attribute(sub_value)
    return d


@read_attribute.register(h5py.Dataset)
def read_dataset(dataset: h5py.Dataset):
    if dataset.attrs.get("source_type", None) == "dataframe":
        return read_dataframe(dataset)
    else:
        value = dataset[()]
        if not hasattr(value, "dtype"):
            return value
        elif isinstance(value.dtype, str):
            pass
        elif issubclass(value.dtype.type, np.string_):
            value = value.astype(str)
            if len(value) == 1:  # Backwards compat, old datasets have strings written as one element 1d arrays
                return value[0]
        elif len(value.dtype.descr) > 1:  # Compound dtype
            value = _from_fixed_length_strings(value)  # For backwards compat, now strings are written as variable length
        if value.shape == ():
            value = value[()]
        return value


@read_attribute.register(type(None))
def read_attribute_none(value):
    return None


@read_attribute.register(h5py.SparseDataset)
def read_sparse_dataset(value):
    return value.value
