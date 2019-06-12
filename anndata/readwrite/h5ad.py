import h5py as _h5py

from pathlib import Path
import numpy as np
import pandas as pd
from pandas.api.types import (
    is_string_dtype,
    is_categorical
)
from scipy import sparse
from .. import h5py
# from .. import AnnData
from ..core.anndata import AnnData, Raw
from collections.abc import Mapping
from multipledispatch import dispatch


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


def write_h5ad(filepath, adata, force_dense=False, dataset_kwargs={}, **kwargs):
    fields = ["X", "obs", "var", "obsm", "varm", "layers", "uns"]
    adata.strings_to_categoricals()
    if adata.raw is not None:
        adata.strings_to_categoricals(adata.raw.var)
    dataset_kwargs = dataset_kwargs.copy()
    dataset_kwargs.update(kwargs)
    filepath = Path(filepath)
    mode = "a" if adata.isbacked else "w"
    if adata.isbacked:  # close so that we can reopen below
        adata.file.close()
    if adata.isbacked and filepath != adata.filename:
        X = adata.X[:]
    else:
        X = adata._X
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

# To allow the optional dataset_kwargs
@dispatch(h5py.File, str, object)
def write_attribute(f, key, value):
    write_attribute(f, key, value, dict())

@dispatch(h5py.File, str, Raw, dict)
def write_attribute(f, key, value, dataset_kwargs={}):
    write_attribute(f, 'raw.X', value.X, dataset_kwargs)
    write_attribute(f, 'raw.var', value.var, dataset_kwargs)
    write_attribute(f, 'raw.varm', value.varm, dataset_kwargs)

@dispatch(h5py.File, str, object, dict)
def write_attribute(f, key, value, dataset_kwargs={}):
    if key in f.keys():
        del f[key]
    # If it's not an array, try and make it an array. If that fails, pickle it.
    # Maybe rethink that, maybe this should just pickle, and have explicit implementations for everything else
    raise NotImplementedError(
        f"Failed to write value for {key}, since a writer for type {type(value)}"
        f" has not been implemented yet."
    )

@dispatch(h5py.File, str, h5py.Dataset, dict)
def write_attribute(f, key, value, dataset_kwargs={}):
    if key in f.keys():
        del f[key]
    f.create_dataset(key, value, **dataset_kwargs)

@dispatch(h5py.File, str, list, dict)
def write_attribute(f, key, value, dataset_kwargs={}):
    write_attribute(f, key, np.array(value), dataset_kwargs)

@dispatch(h5py.File, str, type(None), dict)
def write_attribute(f, key, value, dataset_kwargs={}):
    pass

@dispatch(h5py.File, str, str, dict)
def write_attribute(f, key, value, dataset_kwargs={}):
    if key in f.keys():
        del f[key]
    new_val = np.string_(value)
    f.create_dataset(key, new_val, **dataset_kwargs)

@dispatch(h5py.File, str, (float, np.floating), dict)
def write_attribute(f, key, value, dataset_kwargs):
    if key in f.keys():
        del f[key]
    f.create_dataset(key, value, **dataset_kwargs)

@dispatch(h5py.File, str, (int, np.integer), dict)
def write_attribute(f, key, value, dataset_kwargs={}):
    if key in f.keys():
        del f[key]
    f.create_dataset(key, value, **dataset_kwargs)

@dispatch(h5py.File, str, (bool, np.bool_), dict)
def write_attribute(f, key, value, dataset_kwargs={}):
    if key in f.keys():
        del f[key]
    f.create_dataset(key, value, **dataset_kwargs)

@dispatch(h5py.File, str, np.ndarray, dict)
def write_attribute(f, key, value, dataset_kwargs={}):
    if key in f.keys() and key != "X":
        del f[key]
    # Convert unicode to fixed length strings
    if value.dtype.kind == 'U':
        value = value.astype(np.string_)
    elif value.dtype.names is not None:
        value = _correct_compound_dtype(value)
    # elif value.dtype.name is None:
    #     return
    f.create_dataset(key, data=value, **dataset_kwargs)

@dispatch(h5py.File, str, sparse.spmatrix, dict)
def write_attribute(f, key, value, dataset_kwargs={}):
    if key in f.keys() and key != "X":
        del f[key]
    f.create_dataset(key, data=value, **dataset_kwargs)

@dispatch(h5py.File, str, pd.DataFrame, dict)
def write_attribute(f, key, value, dataset_kwargs={}):
    if key in f.keys():
        del f[key]
    recs = _df_to_hdf5_recarray(value)
    f.create_dataset(key, data=recs, **dataset_kwargs)
    f[key].attrs["source_type"] = "dataframe"

@dispatch(h5py.File, str, Mapping, dict)
def write_attribute(f, key, value, dataset_kwargs={}):
    for sub_key, sub_value in value.items():
        write_attribute(f, f"{key}/{sub_key}", sub_value, dataset_kwargs)

def read_h5ad(filename):
    d = {}
    with h5py.File(filename, "r") as f:
        d["X"] = read_attribute(f["X"])
        # Backwards compat for objects where dataframes weren't labelled
        d["obs"] = read_dataframe(f["obs"])
        d["var"] = read_dataframe(f["var"])
        d["obsm"] = read_attribute(f["obsm"])
        d["obsm"] = read_attribute(f["varm"])
        d["uns"] = read_attribute(f["uns"])
        d["layers"] = read_attribute(f["layers"])

def read_dataframe(dataset):
    df = pd.DataFrame(dataset[:])
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
    for sub_key, sub_value in value.items():
        d[sub_key] = read_attribute(sub_key, sub_value)
    return d

# def read_sparse_dataset(dataset: h5py.SparseDataset)

def read_dataset(dataset: h5py.Dataset):
    if dataset.attr.get("source_type", None) == "dataframe":
        return read_dataframe(dataset)

    if "source_type" in dataset.attr:
        st = dataset.attr["source_type"]
        if "source_type" == pd.DataFrame():
            return read_dataframe(dataset)
    if value.ndim == 1 and len(value) == 1 and value.dtype.names is None:
        value = value[0]

@dispatch(h5py.Group)
def read_attribute(value):
    return read_group(value)

@dispatch(h5py.Dataset)
def read_attribute(value):
    return read_dataset(value)




# class FailedMatch():
#     pass
# class Pattern(object):
#     """Checks if an object matches the criteria

#     Matches can be specified by 3 criteria

#     * type: object must be instance of this type
#     * attrs: object must have values equal to these attrs
#     * func: these functions must return true when called on the object
#     """
#     def __init__(self, *, type, attrs: dict, funcs: list):
#         self.type = type
#         self.attr = attr
#         self.funcs = funcs
    
#     def match(self, object):
#         if not isinstance(object, self.type):
#             return False
#         for attr, value in self.attrs.items():
#             obj_val = getattr(object, attr, FailedMatch)
#             if obj_val is FailedMatch:
#                 return False
#             if obj_val != value:
#                 return False
#         for func in self.funcs:
#             if func(object) == False:
#                 return False
#         return True

# HdfDataFrame = Pattern(
#     type=h5py.Datset,
#     funcs=[
#         lambda x: x.attrs["source_type"] == pd.DataFrame),
#     ]
# )

