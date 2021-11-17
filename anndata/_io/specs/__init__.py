from __future__ import annotations

from os import PathLike
from collections.abc import Mapping
from functools import singledispatch, partial, wraps
from typing import NamedTuple, Tuple, Union, Type, Callable
import typing
from types import MappingProxyType
from warnings import warn

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

import anndata as ad
from anndata import AnnData, Raw
from anndata._core.index import _normalize_indices
from anndata._core.merge import intersect_keys
from anndata._core.sparse_dataset import SparseDataset
from anndata._core import views
from anndata.compat import (
    Literal,
    OverloadedDict,
    _read_attr,
    _from_fixed_length_strings,
    _decode_structured_array,
)
from anndata._io.utils import report_write_key_on_error, check_key, H5PY_V3
from anndata._warnings import OldFormatWarning


# TODO: This probably should be replaced by a hashable Mapping due to conversion b/w "_" and "-"
class IOSpec(NamedTuple):
    encoding_type: str
    encoding_version: str


def write_spec(spec: IOSpec):
    def decorator(func: Callable):
        @wraps(func)
        def wrapper(g, k, *args, **kwargs):
            result = func(g, k, *args, **kwargs)
            g[k].attrs.setdefault("encoding-type", spec.encoding_type)
            g[k].attrs.setdefault("encoding-version", spec.encoding_version)
            return result

        return wrapper

    return decorator


class IORegistry(object):
    def __init__(self):
        self.read: typing.Mapping[IOSpec, Callable] = {}
        self.read_partial: typing.Mapping[IOSpec, Callable] = {}
        self.write: typing.Mapping[Union[Type, Tuple[Type, str]], Callable] = {}

    def register_write(
        self,
        typ: Union[type, tuple[type, str]],
        spec,
        modifiers: frozenset(str) = frozenset(),
    ):
        spec = proc_spec(spec)
        modifiers = frozenset(modifiers)

        def _register(func):
            self.write[(typ, modifiers)] = write_spec(spec)(func)
            return func

        return _register

    def get_writer(self, typ, modifiers=frozenset()):
        modifiers = frozenset(modifiers)
        return self.write[(typ, modifiers)]

    def has_writer(self, typ, modifiers):
        modifiers = frozenset(modifiers)
        return (typ, modifiers) in self.write

    def register_read(self, spec, modifiers: frozenset[str] = frozenset()):
        spec = proc_spec(spec)
        modifiers = frozenset(modifiers)

        def _register(func):
            self.read[(spec, modifiers)] = func
            return func

        return _register

    def get_reader(self, spec, modifiers=frozenset()):
        modifiers = frozenset(modifiers)
        return self.read[(spec, modifiers)]

    def has_reader(self, spec, modifiers=frozenset()):
        modifiers = frozenset(modifiers)
        return (spec, modifiers) in self.read

    def register_read_partial(self, spec, modifiers: frozenset[str] = frozenset()):
        spec = proc_spec(spec)
        modifiers = frozenset(modifiers)

        def _register(func):
            self.read_partial[(spec, modifiers)] = func
            return func

        return _register

    def get_partial_reader(self, spec, modifiers=frozenset()):
        modifiers = frozenset(modifiers)
        return self.read_partial[(spec, modifiers)]


_REGISTRY = IORegistry()


@singledispatch
def proc_spec(spec) -> IOSpec:
    raise NotImplementedError(f"proc_spec not defined for type: {type(spec)}.")


@proc_spec.register(IOSpec)
def proc_spec_spec(spec) -> IOSpec:
    return spec


@proc_spec.register(Mapping)
def proc_spec_mapping(spec) -> IOSpec:
    return IOSpec(**{k.replace("-", "_"): v for k, v in spec.items()})


def get_spec(elem: Union[h5py.Dataset, h5py.Group]) -> IOSpec:
    return proc_spec(
        {
            k: _read_attr(elem.attrs, k, "")
            for k in ["encoding-type", "encoding-version"]
        }
    )


####################
# Dispatch methods #
####################

# def is_full_slice(idx):
#     if isinstance(idx, tuple)len(idx) == 1:

#     if isinstance(idx, type(None)):
#         return True
#     elif idx is Ellipsis:
#         return True
#     elif isinstance(idx, tuple):
#         for el in idx:
#             if isinstance(el, type(None)):
#                 pass
#             elif isinstance(el, slice):
#                 if el != slice(None):
#                     return False
#             else:
#                 return False
#         return True
#     return False


def read_elem(elem, modifiers: frozenset(str) = frozenset()):
    """Read an element from an on disk store."""
    return _REGISTRY.get_reader(get_spec(elem), frozenset(modifiers))(elem)


# TODO: If all items would be read, just call normal read method
def read_elem_partial(
    elem,
    *,
    items=None,
    indices=(slice(None), slice(None)),
    modifiers: frozenset[str] = frozenset(),
):
    """Read part of an element from an on disk store."""
    return _REGISTRY.get_partial_reader(get_spec(elem), frozenset(modifiers))(
        elem, items=items, indices=indices
    )


@report_write_key_on_error
def write_elem(f: h5py.Group, k: str, elem, *args, modifiers=frozenset(), **kwargs):
    """Write an element to an on disk store."""
    if elem is None:
        return
    t = type(elem)
    if k in f:
        del f[k]
    if hasattr(elem, "dtype") and ((t, elem.dtype.kind), modifiers) in _REGISTRY.write:
        _REGISTRY.get_writer((t, elem.dtype.kind), modifiers)(
            f, k, elem, *args, **kwargs
        )
    else:
        _REGISTRY.get_writer(t, modifiers)(f, k, elem, *args, **kwargs)


################################
# Fallbacks / backwards compat #
################################

# Note: there is no need for writing in a backwards compatible format


@_REGISTRY.register_read(IOSpec("", ""))
def read_basic(elem):
    from anndata._io import h5ad

    warn(
        f"Element '{elem.name}' was written without encoding metadata.",
        OldFormatWarning,
    )

    if isinstance(elem, Mapping):
        # Backwards compat sparse arrays
        if "h5sparse_format" in elem.attrs:
            return SparseDataset(elem).to_memory()
        return {k: read_elem(v) for k, v in elem.items()}
    elif isinstance(elem, h5py.Dataset):
        return h5ad.read_dataset(elem)  # TODO: Handle legacy


@_REGISTRY.register_read_partial(IOSpec("", ""))
def read_basic_partial(elem, *, items=None, indices=(slice(None), slice(None))):
    if isinstance(elem, Mapping):
        return _read_partial(elem, items=items, indices=indices)
    elif indices != (slice(None), slice(None)):
        return elem[indices]
    else:
        return elem[()]


###########
# AnnData #
###########


def read_indices(group):
    obs_group = group["obs"]
    obs_idx_elem = obs_group[_read_attr(obs_group.attrs, "_index")]
    obs_idx = read_elem(obs_idx_elem)
    var_group = group["var"]
    var_idx_elem = var_group[_read_attr(var_group.attrs, "_index")]
    var_idx = read_elem(var_idx_elem)
    return obs_idx, var_idx


def read_partial(
    pth: PathLike,
    *,
    obs_idx=slice(None),
    var_idx=slice(None),
    X=True,
    obs=None,
    var=None,
    obsm=None,
    varm=None,
    obsp=None,
    varp=None,
    layers=None,
    uns=None,
) -> ad.AnnData:
    result = {}
    with h5py.File(pth, "r") as f:
        obs_idx, var_idx = _normalize_indices((obs_idx, var_idx), *read_indices(f))
        result["obs"] = read_elem_partial(
            f["obs"], items=obs, indices=(obs_idx, slice(None))
        )
        result["var"] = read_elem_partial(
            f["var"], items=var, indices=(var_idx, slice(None))
        )
        if X:
            result["X"] = read_elem_partial(f["X"], indices=(obs_idx, var_idx))
        else:
            result["X"] = sparse.csr_matrix((len(result["obs"]), len(result["var"])))
        if "obsm" in f:
            result["obsm"] = _read_partial(
                f["obsm"], items=obsm, indices=(obs_idx, slice(None))
            )
        if "varm" in f:
            result["varm"] = _read_partial(
                f["varm"], items=varm, indices=(var_idx, slice(None))
            )
        if "obsp" in f:
            result["obsp"] = _read_partial(
                f["obsp"], items=obsp, indices=(obs_idx, obs_idx)
            )
        if "varp" in f:
            result["varp"] = _read_partial(
                f["varp"], items=varp, indices=(var_idx, var_idx)
            )
        if "layers" in f:
            result["layers"] = _read_partial(
                f["layers"], items=layers, indices=(obs_idx, var_idx)
            )
        if "uns" in f:
            result["uns"] = _read_partial(f["uns"], items=uns)

    return ad.AnnData(**result)


def _read_partial(group, *, items=None, indices=(slice(None), slice(None))):
    if group is None:
        return None
    if items is None:
        keys = intersect_keys((group,))
    else:
        keys = intersect_keys((group, items))
    result = {}
    for k in keys:
        if isinstance(items, Mapping):
            next_items = items.get(k, None)
        else:
            next_items = None
        result[k] = read_elem_partial(group[k], items=next_items, indices=indices)
    return result


def read(pth, *, modifiers=MappingProxyType({})):
    with h5py.File(pth, "r") as f:
        results = {k: read_elem(v) for k, v in f.items()}
    return ad.AnnData(**results)


def write(adata, pth, dataset_kwargs=MappingProxyType({})):
    with h5py.File(pth, "w") as f:
        write_elem(f, "/", adata, dataset_kwargs=dataset_kwargs)


@_REGISTRY.register_write(AnnData, IOSpec("anndata", "0.1.0"))
def write_anndata(f, k, adata, dataset_kwargs=MappingProxyType({})):
    g = f.create_group(k)
    write_elem(g, "X", adata.X, dataset_kwargs=dataset_kwargs)
    write_elem(g, "obs", adata.obs, dataset_kwargs=dataset_kwargs)
    write_elem(g, "var", adata.var, dataset_kwargs=dataset_kwargs)
    write_elem(g, "obsm", dict(adata.obsm), dataset_kwargs=dataset_kwargs)
    write_elem(g, "varm", dict(adata.varm), dataset_kwargs=dataset_kwargs)
    write_elem(g, "obsp", dict(adata.obsp), dataset_kwargs=dataset_kwargs)
    write_elem(g, "varp", dict(adata.varp), dataset_kwargs=dataset_kwargs)
    write_elem(g, "layers", dict(adata.layers), dataset_kwargs=dataset_kwargs)
    write_elem(g, "uns", dict(adata.uns), dataset_kwargs=dataset_kwargs)
    write_elem(g, "raw", adata.raw, dataset_kwargs=dataset_kwargs)


@_REGISTRY.register_read(IOSpec("anndata", "0.1.0"))
@_REGISTRY.register_read(IOSpec("raw", "0.1.0"))
def read_anndata(elem):
    d = {}
    for k in [
        "X",
        "obs",
        "var",
        "obsm",
        "varm",
        "obsp",
        "varp",
        "layers",
        "uns",
        "raw",
    ]:
        if k in elem:
            d[k] = read_elem(elem[k])
    return AnnData(**d)


@_REGISTRY.register_write(Raw, IOSpec("raw", "0.1.0"))
def write_raw(f, k, raw, dataset_kwargs=MappingProxyType({})):
    g = f.create_group("raw")
    write_elem(g, "X", raw.X, dataset_kwargs=dataset_kwargs)
    write_elem(g, "var", raw.var, dataset_kwargs=dataset_kwargs)
    write_elem(g, "varm", dict(raw.varm), dataset_kwargs=dataset_kwargs)


############
# Mappings #
############


@_REGISTRY.register_read(IOSpec("dict", "0.1.0"))
def read_mapping(elem):
    return {k: read_elem(v) for k, v in elem.items()}


@_REGISTRY.register_write(OverloadedDict, IOSpec("dict", "0.1.0"))
@_REGISTRY.register_write(dict, IOSpec("dict", "0.1.0"))
def write_mapping(f, k, v, dataset_kwargs=MappingProxyType({})):
    g = f.create_group(k)
    for sub_k, sub_v in v.items():
        write_elem(g, sub_k, sub_v, dataset_kwargs=dataset_kwargs)


##############
# np.ndarray #
##############


@_REGISTRY.register_write(list, IOSpec("array", "0.2.0"))
def write_list(f, k, elem, dataset_kwargs=MappingProxyType({})):
    write_elem(f, k, np.array(elem), dataset_kwargs=dataset_kwargs)


# TODO: Is this the right behaviour for MaskedArrays?
# It's in the `AnnData.concatenate` docstring, but should we keep it?
@_REGISTRY.register_write(views.ArrayView, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(np.ndarray, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(h5py.Dataset, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(np.ma.MaskedArray, IOSpec("array", "0.2.0"))
def write_basic(f, k, elem, dataset_kwargs=MappingProxyType({})):
    """Write methods which underlying library handles nativley."""
    f.create_dataset(k, data=elem, **dataset_kwargs)


@_REGISTRY.register_read(IOSpec("array", "0.2.0"))
def read_array(elem):
    return elem[()]


@_REGISTRY.register_read_partial(IOSpec("array", "0.2.0"))
def read_array_partial(elem, *, items=None, indices=(slice(None, None))):
    return elem[indices]


# arrays of strings
@_REGISTRY.register_read(IOSpec("string-array", "0.2.0"))
def read_string_array(d):
    return read_array(d.asstr())


@_REGISTRY.register_read_partial(IOSpec("string-array", "0.2.0"))
def read_array_partial(d, items=None, indices=slice(None)):
    return read_array_partial(d.asstr(), items=items, indices=indices)


@_REGISTRY.register_write((views.ArrayView, "U"), IOSpec("string-array", "0.2.0"))
@_REGISTRY.register_write((views.ArrayView, "O"), IOSpec("string-array", "0.2.0"))
@_REGISTRY.register_write((np.ndarray, "U"), IOSpec("string-array", "0.2.0"))
@_REGISTRY.register_write((np.ndarray, "O"), IOSpec("string-array", "0.2.0"))
def write_vlen_string_array(f, k, elem, dataset_kwargs=MappingProxyType({})):
    """Write methods which underlying library handles nativley."""
    str_dtype = h5py.special_dtype(vlen=str)
    f.create_dataset(k, data=elem.astype(str_dtype), dtype=str_dtype, **dataset_kwargs)


###############
# np.recarray #
###############


def _to_hdf5_vlen_strings(value: np.ndarray) -> np.ndarray:
    """This corrects compound dtypes to work with hdf5 files."""
    new_dtype = []
    for dt_name, (dt_type, _) in value.dtype.fields.items():
        if dt_type.kind in ("U", "O"):
            new_dtype.append((dt_name, h5py.special_dtype(vlen=str)))
        else:
            new_dtype.append((dt_name, dt_type))
    return value.astype(new_dtype)


@_REGISTRY.register_read(IOSpec("rec-array", "0.2.0"))
def read_recarray(d):
    value = d[()]
    dtype = value.dtype
    value = _from_fixed_length_strings(value)
    if H5PY_V3:
        value = _decode_structured_array(value, dtype=dtype)
    return value


@_REGISTRY.register_write((np.ndarray, "V"), IOSpec("rec-array", "0.2.0"))
@_REGISTRY.register_write(np.recarray, IOSpec("rec-array", "0.2.0"))
def write_recarray(f, k, elem, dataset_kwargs=MappingProxyType({})):
    f.create_dataset(k, data=_to_hdf5_vlen_strings(elem), **dataset_kwargs)


#################
# Sparse arrays #
#################


def write_sparse_compressed(
    f, key, value, fmt: Literal["csr", "csc"], dataset_kwargs=MappingProxyType({})
):
    g = f.create_group(key)
    g.attrs["shape"] = value.shape

    # Allow resizing
    if "maxshape" not in dataset_kwargs:
        dataset_kwargs = dict(maxshape=(None,), **dataset_kwargs)

    g.create_dataset("data", data=value.data, **dataset_kwargs)
    g.create_dataset("indices", data=value.indices, **dataset_kwargs)
    g.create_dataset("indptr", data=value.indptr, **dataset_kwargs)


write_csr = partial(write_sparse_compressed, fmt="csr")
write_csc = partial(write_sparse_compressed, fmt="csc")
_REGISTRY.register_write(sparse.csr_matrix, IOSpec("csr_matrix", "0.1.0"))(write_csr)
_REGISTRY.register_write(views.SparseCSRView, IOSpec("csr_matrix", "0.1.0"))(write_csr)
_REGISTRY.register_write(sparse.csc_matrix, IOSpec("csc_matrix", "0.1.0"))(write_csc)
_REGISTRY.register_write(views.SparseCSCView, IOSpec("csc_matrix", "0.1.0"))(write_csc)


@_REGISTRY.register_write(SparseDataset, IOSpec("", "0.1.0"))
def write_sparse_dataset(f, k, elem, dataset_kwargs=MappingProxyType({})):
    write_sparse_compressed(
        f, k, elem.to_backed(), fmt=elem.format_str, dataset_kwargs=dataset_kwargs
    )
    # TODO: Cleaner way to do this
    f[k].attrs["encoding-type"] = f"{elem.format_str}_matrix"
    f[k].attrs["encoding-version"] = "0.1.0"


@_REGISTRY.register_read(IOSpec("csc_matrix", "0.1.0"))
@_REGISTRY.register_read(IOSpec("csr_matrix", "0.1.0"))
def read_sparse(elem):
    return SparseDataset(elem).to_memory()


@_REGISTRY.register_read_partial(IOSpec("csc_matrix", "0.1.0"))
@_REGISTRY.register_read_partial(IOSpec("csr_matrix", "0.1.0"))
def read_sparse_partial(elem, *, items=None, indices=(slice(None), slice(None))):
    return SparseDataset(elem)[indices]


##############
# DataFrames #
##############


@_REGISTRY.register_write(views.DataFrameView, IOSpec("dataframe", "0.2.0"))
@_REGISTRY.register_write(pd.DataFrame, IOSpec("dataframe", "0.2.0"))
def write_dataframe(f, key, df, dataset_kwargs=MappingProxyType({})):
    # Check arguments
    for reserved in ("_index",):
        if reserved in df.columns:
            raise ValueError(f"{reserved!r} is a reserved name for dataframe columns.")
    group = f.create_group(key)
    col_names = [check_key(c) for c in df.columns]
    group.attrs["column-order"] = col_names

    if df.index.name is not None:
        index_name = df.index.name
    else:
        index_name = "_index"
    group.attrs["_index"] = check_key(index_name)

    # ._values is "the best" array representation. It's the true array backing the
    # object, where `.values` is always a np.ndarray and .array is always a pandas
    # array.
    write_elem(group, index_name, df.index._values, dataset_kwargs=dataset_kwargs)
    for colname, series in df.items():
        # TODO: this should write the "true" representation of the series (i.e. the underlying array or ndarray depending)
        write_elem(group, colname, series._values, dataset_kwargs=dataset_kwargs)


@_REGISTRY.register_read(IOSpec("dataframe", "0.2.0"))
def read_dataframe(elem):
    columns = list(_read_attr(elem.attrs, "column-order"))
    idx_key = _read_attr(elem.attrs, "_index")
    df = pd.DataFrame(
        {k: read_elem(elem[k]) for k in columns},
        index=read_elem(elem[idx_key]),
        columns=list(columns),
    )
    if idx_key != "_index":
        df.index.name = idx_key
    return df


# TODO: Figure out what indices is allowed to be at each element
@_REGISTRY.register_read_partial(IOSpec("dataframe", "0.2.0"))
def read_dataframe_partial(
    elem, *, items=None, indices=(slice(None, None), slice(None, None))
):
    if items is not None:
        columns = [
            col for col in _read_attr(elem.attrs, "column-order") if col in items
        ]
    else:
        columns = list(_read_attr(elem.attrs, "column-order"))
    idx_key = _read_attr(elem.attrs, "_index")
    df = pd.DataFrame(
        {k: read_elem_partial(elem[k], indices=indices[0]) for k in columns},
        index=read_elem_partial(elem[idx_key], indices=indices[0]),
        columns=list(columns),
    )
    if idx_key != "_index":
        df.index.name = idx_key
    return df


# Backwards compat dataframe reading


@_REGISTRY.register_read({"encoding-type": "dataframe", "encoding-version": "0.1.0"})
def read_dataframe_0_1_0(elem):
    columns = _read_attr(elem.attrs, "column-order")
    idx_key = _read_attr(elem.attrs, "_index")
    df = pd.DataFrame(
        {k: read_series(elem[k]) for k in columns},
        index=read_series(elem[idx_key]),
        columns=list(columns),
    )
    if idx_key != "_index":
        df.index.name = idx_key
    return df


def read_series(dataset) -> Union[np.ndarray, pd.Categorical]:
    # For reading older dataframes
    if "categories" in dataset.attrs:
        categories_dset = dataset.parent[_read_attr(dataset.attrs, "categories")]
        categories = categories_dset[...]
        ordered = bool(_read_attr(categories_dset.attrs, "ordered", False))
        return pd.Categorical.from_codes(dataset[...], categories, ordered=ordered)
    else:
        return dataset[...]


@_REGISTRY.register_read_partial(IOSpec("dataframe", "0.1.0"))
def read_partial_dataframe_0_1_0(
    elem, *, items=None, indices=(slice(None), slice(None))
):
    if items is None:
        items = slice(None)
    else:
        items = list(items)
    return read_elem(elem)[items].iloc[indices[0]]


###############
# Categorical #
###############


@_REGISTRY.register_write(pd.Categorical, IOSpec("categorical", "0.2.0"))
def write_categorical(f, k, v, dataset_kwargs=MappingProxyType({})):
    g = f.create_group(k)
    g.attrs["ordered"] = v.ordered

    write_elem(g, "codes", v.codes, dataset_kwargs=dataset_kwargs)
    write_elem(g, "categories", v.categories._values, dataset_kwargs=dataset_kwargs)


@_REGISTRY.register_read(IOSpec("categorical", "0.2.0"))
def read_categorical(elem):
    return pd.Categorical.from_codes(
        codes=read_elem(elem["codes"]),
        categories=read_elem(elem["categories"]),
        ordered=_read_attr(elem.attrs, "ordered"),
    )


@_REGISTRY.register_read_partial(IOSpec("categorical", "0.2.0"))
def read_categorical(elem, *, items=None, indices=(slice(None),)):
    return pd.Categorical.from_codes(
        codes=read_elem_partial(elem["codes"], indices=indices),
        categories=read_elem(elem["categories"]),
        ordered=_read_attr(elem.attrs, "ordered"),
    )


###########
# Scalars #
###########


@_REGISTRY.register_read(IOSpec("numeric-scalar", "0.2.0"))
def read_scalar(elem):
    return elem[()]


def write_numeric_scalar(f, key, value, dataset_kwargs=MappingProxyType({})):
    # Can’t compress scalars, error is thrown
    dataset_kwargs = dataset_kwargs.copy()
    dataset_kwargs.pop("compression", None)
    dataset_kwargs.pop("compression_opts", None)
    f.create_dataset(key, data=np.array(value), **dataset_kwargs)


# fmt: off
for numeric_scalar_type in [
    bool, np.bool_,
    np.uint8, np.uint16, np.uint32, np.uint64,
    int, np.int8, np.int16, np.int32, np.int64,
    float, np.float16, np.float32, np.float64, np.float128,
    np.complex64, np.complex128, np.complex256,
]:
    _REGISTRY.register_write(numeric_scalar_type, IOSpec("numeric-scalar", "0.2.0"))(write_numeric_scalar)
# fmt: on


@_REGISTRY.register_read(IOSpec("string", "0.2.0"))
def read_string(elem):
    return elem.asstr()[()]


_REGISTRY.register_read(IOSpec("bytes", "0.2.0"))(read_scalar)


@_REGISTRY.register_write(np.str_, IOSpec("string", "0.2.0"))
@_REGISTRY.register_write(str, IOSpec("string", "0.2.0"))
def write_string(f, k, v, dataset_kwargs):
    dataset_kwargs = dataset_kwargs.copy()
    dataset_kwargs.pop("compression", None)
    dataset_kwargs.pop("compression_opts", None)
    f.create_dataset(
        k, data=np.array(v, dtype=h5py.string_dtype(encoding="utf-8")), **dataset_kwargs
    )


# @_REGISTRY.register_write(np.bytes_, IOSpec("bytes", "0.2.0"))
# @_REGISTRY.register_write(bytes, IOSpec("bytes", "0.2.0"))
# def write_string(f, k, v, dataset_kwargs):
#     if "compression" in dataset_kwargs:
#         dataset_kwargs = dict(dataset_kwargs)
#         dataset_kwargs.pop("compression")
#     f.create_dataset(k, data=np.array(v), **dataset_kwargs)
