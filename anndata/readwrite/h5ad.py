from pathlib import Path
from typing import Union

import h5py as _h5py
import numpy as np
import pandas as pd
from pandas.api.types import (
    is_string_dtype,
    is_categorical
)
from scipy import sparse
from .. import h5py
from ..core.anndata import AnnData, Raw
from collections.abc import Mapping


def _df_to_hdf5_recarray(df, var_len_str=True):
    """Convert a pandas dataframe into a recarray h5py can write."""
    if is_string_dtype(df.index):
        if var_len_str:
            index = df.index.values.astype(h5py.special_dtype(vlen=str))
        else:
            max_len_index = 0 if 0 in df.shape else df.index.map(len).max()
            index = df.index.values.astype('S{}'.format(max_len_index))
    else:
        index = df.index.values
    names = ['index']
    arrays = [index]
    formats = [index.dtype]
    for k in df.columns:
        names.append(k)
        if is_string_dtype(df[k]) and not is_categorical(df[k]):
            if var_len_str:
                arrays.append(df[k].values)
                formats.append(h5py.special_dtype(vlen=str))
            else:
                lengths = df[k].map(len)
                if is_categorical(lengths):
                    lengths = lengths.cat.as_ordered()
                arrays.append(df[k].values.astype('S{}'.format(lengths.max())))
                formats.append(arrays[-1].dtype)
        elif is_categorical(df[k]):
            arrays.append(df[k].cat.codes)
            formats.append(_make_h5_cat_dtype(df[k].values))
        else:
            arrays.append(df[k].values)
            formats.append(arrays[-1].dtype)
    assert len(formats) == len(arrays)
    return np.rec.fromarrays(
        arrays,
        dtype={'names': names, 'formats': formats})


def _make_h5_cat_dtype(s):
    """
    This creates a hdf5 enum dtype based on a pandas categorical array.
    """
    dt = h5py.special_dtype(enum=(
        s.codes.dtype,
        dict(zip(s.categories, np.unique(s.codes)))
    ))
    return dt


def _correct_compound_dtype(value):
    """
    This corrects compound dtypes to work with hdf5 files.
    """
    new_dtype = []
    for dt_name, dt_type in value.dtype.descr:
        if dt_type[1] == "U":
            new_dtype.append(
                (dt_name, 'S{}'.format(int(dt_type[2:]) * 4))
            )
        else:
            new_dtype.append((dt_name, dt_type))
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
    # fields = ["X", "obs", "var", "obsm", "varm", "layers", "uns"]
    adata.strings_to_categoricals()
    if adata.raw is not None:
        adata.strings_to_categoricals(adata.raw.var)
    dataset_kwargs = dataset_kwargs.copy()
    dataset_kwargs.update(kwargs)
    filepath = Path(filepath)
    mode = "a" if adata.isbacked else "w"
    if adata.isbacked:  # close so that we can reopen below
        adata.file.close()
    # we're writing to a different location than the backing file
    # - load the matrix into the memory... 
    # TODO: don't load the matrix into memory
    if adata.isbacked and filepath != adata.filename:
        X = adata.X[:]
    else:
        if adata.isbacked:
            X = adata.X
        else:
            X = adata._X  # TODO: This isn't great, what if it's a view?
            # This is so X won't be one dimensional
    with h5py.File(filepath, mode, force_dense=force_dense) as f:
        write_attribute(f, "X", X, dataset_kwargs)
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


def write_str(f, key, value, dataset_kwargs={}):
    f.create_dataset(key, np.string_(value), **dataset_kwargs)


def write_array(f, key, value, dataset_kwargs={}):
    # Convert unicode to fixed length strings
    if value.dtype.kind == 'U':
        value = value.astype(np.string_)
    elif value.dtype.names is not None:
        value = _correct_compound_dtype(value)
    # elif value.dtype.name is None:
    #     return
    f.create_dataset(key, data=value, **dataset_kwargs)


def write_dataframe(f, key, value, dataset_kwargs={}):
    f.create_dataset(key, data=_df_to_hdf5_recarray(value), **dataset_kwargs)
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
    str: write_str,
    (float, np.floating): write_basic,
    (int, np.integer): write_basic,
    (bool, np.bool_): write_basic,
    np.ndarray: write_array,
    sparse.spmatrix: write_basic,
    pd.DataFrame: write_dataframe,
    Mapping: write_mapping
}


# To allow the optional dataset_kwargs
# @dispatch(h5py.File, str, object)
# def write_attribute(f, key, value):
#     write_attribute(f, key, value, dict())

# @dispatch(h5py.File, str, Raw, dict)
# def write_attribute(f, key, value, dataset_kwargs={}):
#     write_attribute(f, 'raw.X', value.X, dataset_kwargs)
#     write_attribute(f, 'raw.var', value.var, dataset_kwargs)
#     write_attribute(f, 'raw.varm', value.varm, dataset_kwargs)

# @dispatch(h5py.File, str, object, dict)
# def write_attribute(f, key, value, dataset_kwargs={}):
#     if key in f.keys():
#         del f[key]
#     # If it's not an array, try and make it an array. If that fails, pickle it.
#     # Maybe rethink that, maybe this should just pickle, and have explicit implementations for everything else
#     raise NotImplementedError(
#         f"Failed to write value for {key}, since a writer for type {type(value)}"
#         f" has not been implemented yet."
#     )

# @dispatch(h5py.File, str, h5py.Dataset, dict)
# def write_attribute(f, key, value, dataset_kwargs={}):
#     if key in f.keys():
#         del f[key]
#     f.create_dataset(key, value, **dataset_kwargs)

# @dispatch(h5py.File, str, list, dict)
# def write_attribute(f, key, value, dataset_kwargs={}):
#     write_attribute(f, key, np.array(value), dataset_kwargs)

# @dispatch(h5py.File, str, type(None), dict)
# def write_attribute(f, key, value, dataset_kwargs={}):
#     pass

# @dispatch(h5py.File, str, str, dict)
# def write_attribute(f, key, value, dataset_kwargs={}):
#     if key in f.keys():
#         del f[key]
#     new_val = np.string_(value)
#     f.create_dataset(key, new_val, **dataset_kwargs)

# @dispatch(h5py.File, str, (float, np.floating), dict)
# def write_attribute(f, key, value, dataset_kwargs):
#     if key in f.keys():
#         del f[key]
#     f.create_dataset(key, value, **dataset_kwargs)

# @dispatch(h5py.File, str, (int, np.integer), dict)
# def write_attribute(f, key, value, dataset_kwargs={}):
#     if key in f.keys():
#         del f[key]
#     f.create_dataset(key, value, **dataset_kwargs)

# @dispatch(h5py.File, str, (bool, np.bool_), dict)
# def write_attribute(f, key, value, dataset_kwargs={}):
#     if key in f.keys():
#         del f[key]
#     f.create_dataset(key, value, **dataset_kwargs)

# @dispatch(h5py.File, str, np.ndarray, dict)
# def write_attribute(f, key, value, dataset_kwargs={}):
#     if key in f.keys() and key != "X":
#         del f[key]
#     # Convert unicode to fixed length strings
#     if value.dtype.kind == 'U':
#         value = value.astype(np.string_)
#     elif value.dtype.names is not None:
#         value = _correct_compound_dtype(value)
#     # elif value.dtype.name is None:
#     #     return
#     f.create_dataset(key, data=value, **dataset_kwargs)

# @dispatch(h5py.File, str, sparse.spmatrix, dict)
# def write_attribute(f, key, value, dataset_kwargs={}):
#     if key in f.keys() and key != "X":
#         del f[key]
#     f.create_dataset(key, data=value, **dataset_kwargs)

# @dispatch(h5py.File, str, pd.DataFrame, dict)
# def write_attribute(f, key, value, dataset_kwargs={}):
#     if key in f.keys():
#         del f[key]
#     recs = _df_to_hdf5_recarray(value)
#     f.create_dataset(key, data=recs, **dataset_kwargs)
#     f[key].attrs["source_type"] = "dataframe"

# @dispatch(h5py.File, str, Mapping, dict)
# def write_attribute(f, key, value, dataset_kwargs={}):
#     for sub_key, sub_value in value.items():
#         write_attribute(f, f"{key}/{sub_key}", sub_value, dataset_kwargs)

# TODO: chunk_size
def read_h5ad(filename, backed=None, chunk_size=None):
    d = {}
    raw = {}

    if backed is None:
        mode = "r"
    else:
        # Open in backed mode
        assert backed in ["r+", "a"]
        mode = backed
        d["filename"] = filename
        d["filemode"] = backed

    # Start reading
    f = h5py.File(filename, mode)
    if mode != "r":
        d["X"] = None
    else:
        # Handling backed objects
        d["X"] = read_attribute(f.get("X", None))
        rawX = read_attribute(f.get("raw.X", None))
        if rawX is not None:
            raw["X"] = rawX
    # Backwards compat for objects where dataframes weren't labelled
    d["obs"] = read_dataframe(f.get("obs", None))
    d["var"] = read_dataframe(f.get("var", None))
    d["obsm"] = read_attribute(f.get("obsm", None))
    d["varm"] = read_attribute(f.get("varm", None))
    d["uns"] = read_attribute(f.get("uns", None))
    d["layers"] = read_attribute(f.get("layers", None))
    # for k in ["raw.X", "raw.var", "raw.varm"]:
    for k in ["raw.var", "raw.varm"]:
        if k in f:
            raw[k.replace("raw.", "")] = read_attribute(f[k])
    # Done reading
    if mode == "r":
        f.close()

    if len(raw) > 0:
        d["raw"] = raw
    if d.get("X", None) is not None:
        d["dtype"] = d["X"].dtype

    # Backwards compat
    k_to_delete = [k for k, v in d.items() if v is None]
    for k in k_to_delete:
        del d[k]
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

    return AnnData(**d)

def read_dataframe(dataset):
    df = pd.DataFrame(dataset[()])
    df.set_index(df.columns[0], inplace=True)
    dt = dataset.dtype
    for col in dt.names[1:]:
        check = _h5py.check_dtype(enum=dt[col])
        if not isinstance(check, dict):
            continue
        mapper = {v: k for k, v in check.items()}
        df[col] = pd.Categorical(df[col].map(mapper), ordered=False)
    return df

def read_group(group: h5py.Group):
    d = dict()
    for sub_key, sub_value in group.items():
        d[sub_key] = read_attribute(sub_value)
    return d

def read_dataset(dataset: h5py.Dataset):
    if dataset.attrs.get("source_type", None) == "dataframe":
        return read_dataframe(dataset)
    else:
        value = dataset[()]
        if issubclass(value.dtype.type, np.string_):
            value = value.astype(str)
        if value.shape == ():
            value = value[()]
        return value

@singledispatch
def read_attribute(value):
    raise NotImplementedError()

@read_attribute.register(h5py.Group)
def read_attribute_group(value):
    return read_group(value)

@read_attribute.register(h5py.Dataset)
def read_attribute_dataset(value):
    return read_dataset(value)

@read_attribute.register(type(None))
def read_attribute_none(value):
    return None

@read_attribute.register(h5py.SparseDataset)
def read_attribute_sparse_dataset(value):
    return value.value

# from abc import ABC, abstractmethod

# class H5ADElement(ABC):
#     @abstractmethod
#     def match(self, element) -> bool:
#         """Check whether element is a match for these methods"""


# class H5ADReader:
#     readers
#     def __init__(self, filename):
#         self.filename = filename

#     def read(self, key):
#         if not self.file.isopen():
#             self.file = h5py.File(self.filename, "r")
#         with h5py.File(self.filename, "r") as f:
#             val = f[key]

#     def read_partial
