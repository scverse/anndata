from __future__ import annotations

import warnings
from collections.abc import Mapping
from copy import copy
from functools import partial
from itertools import product
from types import MappingProxyType
from typing import TYPE_CHECKING
from warnings import warn

import h5py
import numpy as np
import pandas as pd
from packaging.version import Version
from scipy import sparse

import anndata as ad
from anndata import AnnData, Raw
from anndata._core import views
from anndata._core.index import _normalize_indices
from anndata._core.merge import intersect_keys
from anndata._core.sparse_dataset import _CSCDataset, _CSRDataset, sparse_dataset
from anndata._io.utils import H5PY_V3, check_key
from anndata._warnings import OldFormatWarning
from anndata.compat import (
    AwkArray,
    CupyArray,
    CupyCSCMatrix,
    CupyCSRMatrix,
    DaskArray,
    H5Array,
    H5File,
    H5Group,
    ZarrArray,
    ZarrGroup,
    _decode_structured_array,
    _from_fixed_length_strings,
    _read_attr,
    _require_group_write_dataframe,
)

from ..._settings import settings
from .registry import _REGISTRY, IOSpec, read_elem, read_elem_partial

if TYPE_CHECKING:
    from collections.abc import Callable
    from os import PathLike
    from typing import Any, Literal

    from numpy import typing as npt
    from numpy.typing import NDArray

    from anndata._types import ArrayStorageType, GroupStorageType
    from anndata.compat import SpArray
    from anndata.typing import AxisStorable, InMemoryArrayOrScalarType

    from .registry import Reader, Writer

####################
# Dask utils       #
####################

try:
    from dask.utils import SerializableLock as Lock
except ImportError:
    from threading import Lock

# to fix https://github.com/dask/distributed/issues/780
GLOBAL_LOCK = Lock()

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


def _to_cpu_mem_wrapper(write_func):
    """
    Wrapper to bring cupy types into cpu memory before writing.

    Ideally we do direct writing at some point.
    """

    def wrapper(
        f,
        k,
        cupy_val: CupyArray | CupyCSCMatrix | CupyCSRMatrix,
        *,
        _writer: Writer,
        dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
    ):
        return write_func(
            f, k, cupy_val.get(), _writer=_writer, dataset_kwargs=dataset_kwargs
        )

    return wrapper


################################
# Fallbacks / backwards compat #
################################

# Note: there is no need for writing in a backwards compatible format, maybe


@_REGISTRY.register_read(H5File, IOSpec("", ""))
@_REGISTRY.register_read(H5Group, IOSpec("", ""))
@_REGISTRY.register_read(H5Array, IOSpec("", ""))
def read_basic(
    elem: H5File | H5Group | H5Array, *, _reader: Reader
) -> dict[str, InMemoryArrayOrScalarType] | npt.NDArray | sparse.spmatrix | SpArray:
    from anndata._io import h5ad

    warn(
        f"Element '{elem.name}' was written without encoding metadata.",
        OldFormatWarning,
        stacklevel=3,
    )

    if isinstance(elem, Mapping):
        # Backwards compat sparse arrays
        if "h5sparse_format" in elem.attrs:
            return sparse_dataset(elem).to_memory()
        return {k: _reader.read_elem(v) for k, v in elem.items()}
    elif isinstance(elem, h5py.Dataset):
        return h5ad.read_dataset(elem)  # TODO: Handle legacy


@_REGISTRY.register_read(ZarrGroup, IOSpec("", ""))
@_REGISTRY.register_read(ZarrArray, IOSpec("", ""))
def read_basic_zarr(
    elem: ZarrGroup | ZarrArray, *, _reader: Reader
) -> dict[str, InMemoryArrayOrScalarType] | npt.NDArray | sparse.spmatrix | SpArray:
    from anndata._io import zarr

    warn(
        f"Element '{elem.name}' was written without encoding metadata.",
        OldFormatWarning,
        stacklevel=3,
    )

    if isinstance(elem, Mapping):
        # Backwards compat sparse arrays
        if "h5sparse_format" in elem.attrs:
            return sparse_dataset(elem).to_memory()
        return {k: _reader.read_elem(v) for k, v in elem.items()}
    elif isinstance(elem, ZarrArray):
        return zarr.read_dataset(elem)  # TODO: Handle legacy


# @_REGISTRY.register_read_partial(IOSpec("", ""))
# def read_basic_partial(elem, *, items=None, indices=(slice(None), slice(None))):
#     if isinstance(elem, Mapping):
#         return _read_partial(elem, items=items, indices=indices)
#     elif indices != (slice(None), slice(None)):
#         return elem[indices]
#     else:
#         return elem[()]


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


@_REGISTRY.register_write(ZarrGroup, AnnData, IOSpec("anndata", "0.1.0"))
@_REGISTRY.register_write(H5Group, AnnData, IOSpec("anndata", "0.1.0"))
def write_anndata(
    f: GroupStorageType,
    k: str,
    adata: AnnData,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    g = f.require_group(k)
    _writer.write_elem(g, "X", adata.X, dataset_kwargs=dataset_kwargs)
    _writer.write_elem(g, "obs", adata.obs, dataset_kwargs=dataset_kwargs)
    _writer.write_elem(g, "var", adata.var, dataset_kwargs=dataset_kwargs)
    _writer.write_elem(g, "obsm", dict(adata.obsm), dataset_kwargs=dataset_kwargs)
    _writer.write_elem(g, "varm", dict(adata.varm), dataset_kwargs=dataset_kwargs)
    _writer.write_elem(g, "obsp", dict(adata.obsp), dataset_kwargs=dataset_kwargs)
    _writer.write_elem(g, "varp", dict(adata.varp), dataset_kwargs=dataset_kwargs)
    _writer.write_elem(g, "layers", dict(adata.layers), dataset_kwargs=dataset_kwargs)
    _writer.write_elem(g, "uns", dict(adata.uns), dataset_kwargs=dataset_kwargs)
    _writer.write_elem(g, "raw", adata.raw, dataset_kwargs=dataset_kwargs)


@_REGISTRY.register_read(H5Group, IOSpec("anndata", "0.1.0"))
@_REGISTRY.register_read(H5Group, IOSpec("raw", "0.1.0"))
@_REGISTRY.register_read(H5File, IOSpec("anndata", "0.1.0"))
@_REGISTRY.register_read(H5File, IOSpec("raw", "0.1.0"))
@_REGISTRY.register_read(ZarrGroup, IOSpec("anndata", "0.1.0"))
@_REGISTRY.register_read(ZarrGroup, IOSpec("raw", "0.1.0"))
def read_anndata(elem: GroupStorageType | H5File, *, _reader: Reader) -> AnnData:
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
            d[k] = _reader.read_elem(elem[k])
    return AnnData(**d)


@_REGISTRY.register_write(H5Group, Raw, IOSpec("raw", "0.1.0"))
@_REGISTRY.register_write(ZarrGroup, Raw, IOSpec("raw", "0.1.0"))
def write_raw(
    f: GroupStorageType,
    k: str,
    raw: Raw,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    g = f.require_group(k)
    _writer.write_elem(g, "X", raw.X, dataset_kwargs=dataset_kwargs)
    _writer.write_elem(g, "var", raw.var, dataset_kwargs=dataset_kwargs)
    _writer.write_elem(g, "varm", dict(raw.varm), dataset_kwargs=dataset_kwargs)


############
# Mappings #
############


@_REGISTRY.register_read(H5Group, IOSpec("dict", "0.1.0"))
@_REGISTRY.register_read(ZarrGroup, IOSpec("dict", "0.1.0"))
def read_mapping(elem: GroupStorageType, *, _reader: Reader) -> dict[str, AxisStorable]:
    return {k: _reader.read_elem(v) for k, v in elem.items()}


@_REGISTRY.register_write(H5Group, dict, IOSpec("dict", "0.1.0"))
@_REGISTRY.register_write(ZarrGroup, dict, IOSpec("dict", "0.1.0"))
def write_mapping(
    f: GroupStorageType,
    k: str,
    v: dict[str, AxisStorable],
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    g = f.require_group(k)
    for sub_k, sub_v in v.items():
        _writer.write_elem(g, sub_k, sub_v, dataset_kwargs=dataset_kwargs)


##############
# np.ndarray #
##############


@_REGISTRY.register_write(H5Group, list, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, list, IOSpec("array", "0.2.0"))
def write_list(
    f: GroupStorageType,
    k: str,
    elem: list[AxisStorable],
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    _writer.write_elem(f, k, np.array(elem), dataset_kwargs=dataset_kwargs)


# TODO: Is this the right behavior for MaskedArrays?
# It's in the `AnnData.concatenate` docstring, but should we keep it?
@_REGISTRY.register_write(H5Group, views.ArrayView, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(H5Group, np.ndarray, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(H5Group, h5py.Dataset, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(H5Group, np.ma.MaskedArray, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, views.ArrayView, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, np.ndarray, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, h5py.Dataset, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, np.ma.MaskedArray, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, ZarrArray, IOSpec("array", "0.2.0"))
def write_basic(
    f: GroupStorageType,
    k: str,
    elem: views.ArrayView | np.ndarray | h5py.Dataset | np.ma.MaskedArray | ZarrArray,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    """Write methods which underlying library handles natively."""
    f.create_dataset(k, data=elem, **dataset_kwargs)


_REGISTRY.register_write(H5Group, CupyArray, IOSpec("array", "0.2.0"))(
    _to_cpu_mem_wrapper(write_basic)
)
_REGISTRY.register_write(ZarrGroup, CupyArray, IOSpec("array", "0.2.0"))(
    _to_cpu_mem_wrapper(write_basic)
)


@_REGISTRY.register_write(ZarrGroup, DaskArray, IOSpec("array", "0.2.0"))
def write_basic_dask_zarr(
    f: ZarrGroup,
    k: str,
    elem: DaskArray,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    import dask.array as da

    g = f.require_dataset(k, shape=elem.shape, dtype=elem.dtype, **dataset_kwargs)
    da.store(elem, g, lock=GLOBAL_LOCK)


# Adding this separately because h5py isn't serializable
# https://github.com/pydata/xarray/issues/4242
@_REGISTRY.register_write(H5Group, DaskArray, IOSpec("array", "0.2.0"))
def write_basic_dask_h5(
    f: H5Group,
    k: str,
    elem: DaskArray,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    import dask.array as da
    import dask.config as dc

    if dc.get("scheduler", None) == "dask.distributed":
        raise ValueError(
            "Cannot write dask arrays to hdf5 when using distributed scheduler"
        )

    g = f.require_dataset(k, shape=elem.shape, dtype=elem.dtype, **dataset_kwargs)
    da.store(elem, g)


@_REGISTRY.register_read(H5Array, IOSpec("array", "0.2.0"))
@_REGISTRY.register_read(ZarrArray, IOSpec("array", "0.2.0"))
@_REGISTRY.register_read(ZarrArray, IOSpec("string-array", "0.2.0"))
def read_array(elem: ArrayStorageType, *, _reader: Reader) -> npt.NDArray:
    return elem[()]


@_REGISTRY.register_read_partial(H5Array, IOSpec("array", "0.2.0"))
@_REGISTRY.register_read_partial(ZarrArray, IOSpec("string-array", "0.2.0"))
def read_array_partial(elem, *, items=None, indices=(slice(None, None))):
    return elem[indices]


@_REGISTRY.register_read_partial(ZarrArray, IOSpec("array", "0.2.0"))
def read_zarr_array_partial(elem, *, items=None, indices=(slice(None, None))):
    return elem.oindex[indices]


# arrays of strings
@_REGISTRY.register_read(H5Array, IOSpec("string-array", "0.2.0"))
def read_string_array(d: H5Array, *, _reader: Reader):
    return read_array(d.asstr(), _reader=_reader)


@_REGISTRY.register_read_partial(H5Array, IOSpec("string-array", "0.2.0"))
def read_string_array_partial(d, items=None, indices=slice(None)):
    return read_array_partial(d.asstr(), items=items, indices=indices)


@_REGISTRY.register_write(
    H5Group, (views.ArrayView, "U"), IOSpec("string-array", "0.2.0")
)
@_REGISTRY.register_write(
    H5Group, (views.ArrayView, "O"), IOSpec("string-array", "0.2.0")
)
@_REGISTRY.register_write(H5Group, (np.ndarray, "U"), IOSpec("string-array", "0.2.0"))
@_REGISTRY.register_write(H5Group, (np.ndarray, "O"), IOSpec("string-array", "0.2.0"))
def write_vlen_string_array(
    f: H5Group,
    k: str,
    elem: np.ndarray,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    """Write methods which underlying library handles nativley."""
    str_dtype = h5py.special_dtype(vlen=str)
    f.create_dataset(k, data=elem.astype(str_dtype), dtype=str_dtype, **dataset_kwargs)


@_REGISTRY.register_write(
    ZarrGroup, (views.ArrayView, "U"), IOSpec("string-array", "0.2.0")
)
@_REGISTRY.register_write(
    ZarrGroup, (views.ArrayView, "O"), IOSpec("string-array", "0.2.0")
)
@_REGISTRY.register_write(ZarrGroup, (np.ndarray, "U"), IOSpec("string-array", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, (np.ndarray, "O"), IOSpec("string-array", "0.2.0"))
def write_vlen_string_array_zarr(
    f: ZarrGroup,
    k: str,
    elem: np.ndarray,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    import numcodecs

    if Version(numcodecs.__version__) < Version("0.13"):
        msg = "Old numcodecs version detected. Please update for improved performance and stability."
        warnings.warn(msg)
        # Workaround for https://github.com/zarr-developers/numcodecs/issues/514
        if hasattr(elem, "flags") and not elem.flags.writeable:
            elem = elem.copy()

    f.create_dataset(
        k,
        shape=elem.shape,
        dtype=object,
        object_codec=numcodecs.VLenUTF8(),
        **dataset_kwargs,
    )
    f[k][:] = elem


###############
# np.recarray #
###############


def _to_hdf5_vlen_strings(value: np.ndarray) -> np.ndarray:
    """This corrects compound dtypes to work with hdf5 files."""
    new_dtype = []
    for dt_name, (dt_type, _) in value.dtype.fields.items():
        if dt_type.kind in {"U", "O"}:
            new_dtype.append((dt_name, h5py.special_dtype(vlen=str)))
        else:
            new_dtype.append((dt_name, dt_type))
    return value.astype(new_dtype)


@_REGISTRY.register_read(H5Array, IOSpec("rec-array", "0.2.0"))
@_REGISTRY.register_read(ZarrArray, IOSpec("rec-array", "0.2.0"))
def read_recarray(d: ArrayStorageType, *, _reader: Reader) -> np.recarray | npt.NDArray:
    value = d[()]
    dtype = value.dtype
    value = _from_fixed_length_strings(value)
    if H5PY_V3:
        value = _decode_structured_array(value, dtype=dtype)
    return value


@_REGISTRY.register_write(H5Group, (np.ndarray, "V"), IOSpec("rec-array", "0.2.0"))
@_REGISTRY.register_write(H5Group, np.recarray, IOSpec("rec-array", "0.2.0"))
def write_recarray(
    f: H5Group,
    k: str,
    elem: np.ndarray | np.recarray,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    f.create_dataset(k, data=_to_hdf5_vlen_strings(elem), **dataset_kwargs)


@_REGISTRY.register_write(ZarrGroup, (np.ndarray, "V"), IOSpec("rec-array", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, np.recarray, IOSpec("rec-array", "0.2.0"))
def write_recarray_zarr(
    f: ZarrGroup,
    k: str,
    elem: np.ndarray | np.recarray,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    from anndata.compat import _to_fixed_length_strings

    f.create_dataset(k, data=_to_fixed_length_strings(elem), **dataset_kwargs)


#################
# Sparse arrays #
#################


def write_sparse_compressed(
    f: GroupStorageType,
    key: str,
    value: sparse.spmatrix | SpArray,
    *,
    _writer: Writer,
    fmt: Literal["csr", "csc"],
    dataset_kwargs=MappingProxyType({}),
):
    g = f.require_group(key)
    g.attrs["shape"] = value.shape
    dataset_kwargs = dict(dataset_kwargs)
    indptr_dtype = dataset_kwargs.pop("indptr_dtype", value.indptr.dtype)

    # Allow resizing for hdf5
    if isinstance(f, H5Group) and "maxshape" not in dataset_kwargs:
        dataset_kwargs = dict(maxshape=(None,), **dataset_kwargs)

    g.create_dataset("data", data=value.data, **dataset_kwargs)
    g.create_dataset("indices", data=value.indices, **dataset_kwargs)
    g.create_dataset("indptr", data=value.indptr, dtype=indptr_dtype, **dataset_kwargs)


write_csr = partial(write_sparse_compressed, fmt="csr")
write_csc = partial(write_sparse_compressed, fmt="csc")

for store_type, (cls, spec, func) in product(
    (H5Group, ZarrGroup),
    [
        # spmatrix
        (sparse.csr_matrix, IOSpec("csr_matrix", "0.1.0"), write_csr),
        (views.SparseCSRMatrixView, IOSpec("csr_matrix", "0.1.0"), write_csr),
        (sparse.csc_matrix, IOSpec("csc_matrix", "0.1.0"), write_csc),
        (views.SparseCSCMatrixView, IOSpec("csc_matrix", "0.1.0"), write_csc),
        # sparray
        (sparse.csr_array, IOSpec("csr_matrix", "0.1.0"), write_csr),
        (views.SparseCSRArrayView, IOSpec("csr_matrix", "0.1.0"), write_csr),
        (sparse.csc_array, IOSpec("csc_matrix", "0.1.0"), write_csc),
        (views.SparseCSCArrayView, IOSpec("csc_matrix", "0.1.0"), write_csc),
        # cupy spmatrix
        (CupyCSRMatrix, IOSpec("csr_matrix", "0.1.0"), _to_cpu_mem_wrapper(write_csr)),
        (
            views.CupySparseCSRView,
            IOSpec("csr_matrix", "0.1.0"),
            _to_cpu_mem_wrapper(write_csr),
        ),
        (CupyCSCMatrix, IOSpec("csc_matrix", "0.1.0"), _to_cpu_mem_wrapper(write_csc)),
        (
            views.CupySparseCSCView,
            IOSpec("csc_matrix", "0.1.0"),
            _to_cpu_mem_wrapper(write_csc),
        ),
    ],
):
    _REGISTRY.register_write(store_type, cls, spec)(func)


@_REGISTRY.register_write(H5Group, _CSRDataset, IOSpec("", "0.1.0"))
@_REGISTRY.register_write(H5Group, _CSCDataset, IOSpec("", "0.1.0"))
@_REGISTRY.register_write(ZarrGroup, _CSRDataset, IOSpec("", "0.1.0"))
@_REGISTRY.register_write(ZarrGroup, _CSCDataset, IOSpec("", "0.1.0"))
def write_sparse_dataset(
    f: GroupStorageType,
    k: str,
    elem: _CSCDataset | _CSRDataset,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    write_sparse_compressed(
        f,
        k,
        elem._to_backed(),
        _writer=_writer,
        fmt=elem.format,
        dataset_kwargs=dataset_kwargs,
    )
    # TODO: Cleaner way to do this
    f[k].attrs["encoding-type"] = f"{elem.format}_matrix"
    f[k].attrs["encoding-version"] = "0.1.0"


@_REGISTRY.register_write(H5Group, (DaskArray, CupyArray), IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, (DaskArray, CupyArray), IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(
    H5Group, (DaskArray, CupyCSRMatrix), IOSpec("csr_matrix", "0.1.0")
)
@_REGISTRY.register_write(
    H5Group, (DaskArray, CupyCSCMatrix), IOSpec("csc_matrix", "0.1.0")
)
@_REGISTRY.register_write(
    ZarrGroup, (DaskArray, CupyCSRMatrix), IOSpec("csr_matrix", "0.1.0")
)
@_REGISTRY.register_write(
    ZarrGroup, (DaskArray, CupyCSCMatrix), IOSpec("csc_matrix", "0.1.0")
)
def write_cupy_dask_sparse(f, k, elem, _writer, dataset_kwargs=MappingProxyType({})):
    _writer.write_elem(
        f,
        k,
        elem.map_blocks(lambda x: x.get(), dtype=elem.dtype, meta=elem._meta.get()),
        dataset_kwargs=dataset_kwargs,
    )


@_REGISTRY.register_write(
    H5Group, (DaskArray, sparse.csr_matrix), IOSpec("csr_matrix", "0.1.0")
)
@_REGISTRY.register_write(
    H5Group, (DaskArray, sparse.csc_matrix), IOSpec("csc_matrix", "0.1.0")
)
@_REGISTRY.register_write(
    ZarrGroup, (DaskArray, sparse.csr_matrix), IOSpec("csr_matrix", "0.1.0")
)
@_REGISTRY.register_write(
    ZarrGroup, (DaskArray, sparse.csc_matrix), IOSpec("csc_matrix", "0.1.0")
)
def write_dask_sparse(
    f: GroupStorageType,
    k: str,
    elem: DaskArray,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    sparse_format = elem._meta.format

    def as_int64_indices(x):
        x.indptr = x.indptr.astype(np.int64, copy=False)
        x.indices = x.indices.astype(np.int64, copy=False)
        return x

    if sparse_format == "csr":
        axis = 0
    elif sparse_format == "csc":
        axis = 1
    else:
        raise NotImplementedError(
            f"Cannot write dask sparse arrays with format {sparse_format}"
        )

    def chunk_slice(start: int, stop: int) -> tuple[slice | None, slice | None]:
        result = [slice(None), slice(None)]
        result[axis] = slice(start, stop)
        return tuple(result)

    axis_chunks = elem.chunks[axis]
    chunk_start = 0
    chunk_stop = axis_chunks[0]

    _writer.write_elem(
        f,
        k,
        as_int64_indices(elem[chunk_slice(chunk_start, chunk_stop)].compute()),
        dataset_kwargs=dataset_kwargs,
    )

    disk_mtx = sparse_dataset(f[k])

    for chunk_size in axis_chunks[1:]:
        chunk_start = chunk_stop
        chunk_stop += chunk_size

        disk_mtx.append(elem[chunk_slice(chunk_start, chunk_stop)].compute())


@_REGISTRY.register_read(H5Group, IOSpec("csc_matrix", "0.1.0"))
@_REGISTRY.register_read(H5Group, IOSpec("csr_matrix", "0.1.0"))
@_REGISTRY.register_read(ZarrGroup, IOSpec("csc_matrix", "0.1.0"))
@_REGISTRY.register_read(ZarrGroup, IOSpec("csr_matrix", "0.1.0"))
def read_sparse(
    elem: GroupStorageType, *, _reader: Reader
) -> sparse.spmatrix | SpArray:
    return sparse_dataset(elem).to_memory()


@_REGISTRY.register_read_partial(H5Group, IOSpec("csc_matrix", "0.1.0"))
@_REGISTRY.register_read_partial(H5Group, IOSpec("csr_matrix", "0.1.0"))
@_REGISTRY.register_read_partial(ZarrGroup, IOSpec("csc_matrix", "0.1.0"))
@_REGISTRY.register_read_partial(ZarrGroup, IOSpec("csr_matrix", "0.1.0"))
def read_sparse_partial(elem, *, items=None, indices=(slice(None), slice(None))):
    return sparse_dataset(elem)[indices]


#################
# Awkward array #
#################


@_REGISTRY.register_write(H5Group, AwkArray, IOSpec("awkward-array", "0.1.0"))
@_REGISTRY.register_write(ZarrGroup, AwkArray, IOSpec("awkward-array", "0.1.0"))
@_REGISTRY.register_write(
    H5Group, views.AwkwardArrayView, IOSpec("awkward-array", "0.1.0")
)
@_REGISTRY.register_write(
    ZarrGroup, views.AwkwardArrayView, IOSpec("awkward-array", "0.1.0")
)
def write_awkward(
    f: GroupStorageType,
    k: str,
    v: views.AwkwardArrayView | AwkArray,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    from anndata.compat import awkward as ak

    group = f.require_group(k)
    if isinstance(v, views.AwkwardArrayView):
        # copy to remove the view attributes
        v = copy(v)
    form, length, container = ak.to_buffers(ak.to_packed(v))
    group.attrs["length"] = length
    group.attrs["form"] = form.to_json()
    for k, v in container.items():
        _writer.write_elem(group, k, v, dataset_kwargs=dataset_kwargs)


@_REGISTRY.register_read(H5Group, IOSpec("awkward-array", "0.1.0"))
@_REGISTRY.register_read(ZarrGroup, IOSpec("awkward-array", "0.1.0"))
def read_awkward(elem: GroupStorageType, *, _reader: Reader) -> AwkArray:
    from anndata.compat import awkward as ak

    form = _read_attr(elem.attrs, "form")
    length = _read_attr(elem.attrs, "length")
    container = {k: _reader.read_elem(elem[k]) for k in elem.keys()}

    return ak.from_buffers(form, int(length), container)


##############
# DataFrames #
##############


@_REGISTRY.register_write(H5Group, views.DataFrameView, IOSpec("dataframe", "0.2.0"))
@_REGISTRY.register_write(H5Group, pd.DataFrame, IOSpec("dataframe", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, views.DataFrameView, IOSpec("dataframe", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, pd.DataFrame, IOSpec("dataframe", "0.2.0"))
def write_dataframe(
    f: GroupStorageType,
    key: str,
    df: views.DataFrameView | pd.DataFrame,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    # Check arguments
    for reserved in ("_index",):
        if reserved in df.columns:
            raise ValueError(f"{reserved!r} is a reserved name for dataframe columns.")
    group = _require_group_write_dataframe(f, key, df)
    if not df.columns.is_unique:
        duplicates = list(df.columns[df.columns.duplicated()])
        raise ValueError(
            f"Found repeated column names: {duplicates}. Column names must be unique."
        )
    col_names = [check_key(c) for c in df.columns]
    group.attrs["column-order"] = col_names

    if df.index.name is not None:
        if df.index.name in col_names and not pd.Series(
            df.index, index=df.index
        ).equals(df[df.index.name]):
            raise ValueError(
                f"DataFrame.index.name ({df.index.name!r}) is also used by a column "
                "whose values are different. This is not supported. Please make sure "
                "the values are the same, or use a different name."
            )
        index_name = df.index.name
    else:
        index_name = "_index"
    group.attrs["_index"] = check_key(index_name)

    # ._values is "the best" array representation. It's the true array backing the
    # object, where `.values` is always a np.ndarray and .array is always a pandas
    # array.
    _writer.write_elem(
        group, index_name, df.index._values, dataset_kwargs=dataset_kwargs
    )
    for colname, series in df.items():
        # TODO: this should write the "true" representation of the series (i.e. the underlying array or ndarray depending)
        _writer.write_elem(
            group, colname, series._values, dataset_kwargs=dataset_kwargs
        )


@_REGISTRY.register_read(H5Group, IOSpec("dataframe", "0.2.0"))
@_REGISTRY.register_read(ZarrGroup, IOSpec("dataframe", "0.2.0"))
def read_dataframe(elem: GroupStorageType, *, _reader: Reader) -> pd.DataFrame:
    columns = list(_read_attr(elem.attrs, "column-order"))
    idx_key = _read_attr(elem.attrs, "_index")
    df = pd.DataFrame(
        {k: _reader.read_elem(elem[k]) for k in columns},
        index=_reader.read_elem(elem[idx_key]),
        columns=columns if len(columns) else None,
    )
    if idx_key != "_index":
        df.index.name = idx_key
    return df


# TODO: Figure out what indices is allowed to be at each element
@_REGISTRY.register_read_partial(H5Group, IOSpec("dataframe", "0.2.0"))
@_REGISTRY.register_read_partial(ZarrGroup, IOSpec("dataframe", "0.2.0"))
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
        columns=columns if len(columns) else None,
    )
    if idx_key != "_index":
        df.index.name = idx_key
    return df


# Backwards compat dataframe reading


@_REGISTRY.register_read(H5Group, IOSpec("dataframe", "0.1.0"))
@_REGISTRY.register_read(ZarrGroup, IOSpec("dataframe", "0.1.0"))
def read_dataframe_0_1_0(elem: GroupStorageType, *, _reader: Reader) -> pd.DataFrame:
    columns = _read_attr(elem.attrs, "column-order")
    idx_key = _read_attr(elem.attrs, "_index")
    df = pd.DataFrame(
        {k: read_series(elem[k]) for k in columns},
        index=read_series(elem[idx_key]),
        columns=columns if len(columns) else None,
    )
    if idx_key != "_index":
        df.index.name = idx_key
    return df


def read_series(dataset: h5py.Dataset) -> np.ndarray | pd.Categorical:
    # For reading older dataframes
    if "categories" in dataset.attrs:
        if isinstance(dataset, ZarrArray):
            import zarr

            parent_name = dataset.name.rstrip(dataset.basename)
            parent = zarr.open(dataset.store)[parent_name]
        else:
            parent = dataset.parent
        categories_dset = parent[_read_attr(dataset.attrs, "categories")]
        categories = read_elem(categories_dset)
        ordered = bool(_read_attr(categories_dset.attrs, "ordered", False))
        return pd.Categorical.from_codes(
            read_elem(dataset), categories, ordered=ordered
        )
    else:
        return read_elem(dataset)


@_REGISTRY.register_read_partial(H5Group, IOSpec("dataframe", "0.1.0"))
@_REGISTRY.register_read_partial(ZarrGroup, IOSpec("dataframe", "0.1.0"))
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


@_REGISTRY.register_write(H5Group, pd.Categorical, IOSpec("categorical", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, pd.Categorical, IOSpec("categorical", "0.2.0"))
def write_categorical(
    f: GroupStorageType,
    k: str,
    v: pd.Categorical,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    g = f.require_group(k)
    g.attrs["ordered"] = bool(v.ordered)

    _writer.write_elem(g, "codes", v.codes, dataset_kwargs=dataset_kwargs)
    _writer.write_elem(
        g, "categories", v.categories._values, dataset_kwargs=dataset_kwargs
    )


@_REGISTRY.register_read(H5Group, IOSpec("categorical", "0.2.0"))
@_REGISTRY.register_read(ZarrGroup, IOSpec("categorical", "0.2.0"))
def read_categorical(elem: GroupStorageType, *, _reader: Reader) -> pd.Categorical:
    return pd.Categorical.from_codes(
        codes=_reader.read_elem(elem["codes"]),
        categories=_reader.read_elem(elem["categories"]),
        ordered=bool(_read_attr(elem.attrs, "ordered")),
    )


@_REGISTRY.register_read_partial(H5Group, IOSpec("categorical", "0.2.0"))
@_REGISTRY.register_read_partial(ZarrGroup, IOSpec("categorical", "0.2.0"))
def read_partial_categorical(elem, *, items=None, indices=(slice(None),)):
    return pd.Categorical.from_codes(
        codes=read_elem_partial(elem["codes"], indices=indices),
        categories=read_elem(elem["categories"]),
        ordered=bool(_read_attr(elem.attrs, "ordered")),
    )


####################
# Pandas nullables #
####################


@_REGISTRY.register_write(
    H5Group, pd.arrays.IntegerArray, IOSpec("nullable-integer", "0.1.0")
)
@_REGISTRY.register_write(
    ZarrGroup, pd.arrays.IntegerArray, IOSpec("nullable-integer", "0.1.0")
)
@_REGISTRY.register_write(
    H5Group, pd.arrays.BooleanArray, IOSpec("nullable-boolean", "0.1.0")
)
@_REGISTRY.register_write(
    ZarrGroup, pd.arrays.BooleanArray, IOSpec("nullable-boolean", "0.1.0")
)
@_REGISTRY.register_write(
    H5Group, pd.arrays.StringArray, IOSpec("nullable-string-array", "0.1.0")
)
@_REGISTRY.register_write(
    ZarrGroup, pd.arrays.StringArray, IOSpec("nullable-string-array", "0.1.0")
)
def write_nullable(
    f: GroupStorageType,
    k: str,
    v: pd.arrays.IntegerArray | pd.arrays.BooleanArray | pd.arrays.StringArray,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    if (
        isinstance(v, pd.arrays.StringArray)
        and not settings.allow_write_nullable_strings
    ):
        msg = (
            "`anndata.settings.allow_write_nullable_strings` is False, "
            "because writing of `pd.arrays.StringArray` is new "
            "and not supported in anndata < 0.11, still use by many people. "
            "Opt-in to writing these arrays by toggling the setting to True."
        )
        raise RuntimeError(msg)
    g = f.require_group(k)
    values = (
        v.to_numpy(na_value="")
        if isinstance(v, pd.arrays.StringArray)
        else v.to_numpy(na_value=0, dtype=v.dtype.numpy_dtype)
    )
    _writer.write_elem(g, "values", values, dataset_kwargs=dataset_kwargs)
    _writer.write_elem(g, "mask", v.isna(), dataset_kwargs=dataset_kwargs)


def _read_nullable(
    elem: GroupStorageType,
    *,
    _reader: Reader,
    # BaseMaskedArray
    array_type: Callable[
        [NDArray[np.number], NDArray[np.bool_]], pd.api.extensions.ExtensionArray
    ],
) -> pd.api.extensions.ExtensionArray:
    return array_type(
        _reader.read_elem(elem["values"]),
        mask=_reader.read_elem(elem["mask"]),
    )


def _string_array(
    values: np.ndarray, mask: np.ndarray
) -> pd.api.extensions.ExtensionArray:
    """Construct a string array from values and mask."""
    arr = pd.array(values, dtype="string")
    arr[mask] = pd.NA
    return arr


_REGISTRY.register_read(H5Group, IOSpec("nullable-integer", "0.1.0"))(
    read_nullable_integer := partial(_read_nullable, array_type=pd.arrays.IntegerArray)
)
_REGISTRY.register_read(ZarrGroup, IOSpec("nullable-integer", "0.1.0"))(
    read_nullable_integer
)

_REGISTRY.register_read(H5Group, IOSpec("nullable-boolean", "0.1.0"))(
    read_nullable_boolean := partial(_read_nullable, array_type=pd.arrays.BooleanArray)
)
_REGISTRY.register_read(ZarrGroup, IOSpec("nullable-boolean", "0.1.0"))(
    read_nullable_boolean
)

_REGISTRY.register_read(H5Group, IOSpec("nullable-string-array", "0.1.0"))(
    read_nullable_string := partial(_read_nullable, array_type=_string_array)
)
_REGISTRY.register_read(ZarrGroup, IOSpec("nullable-string-array", "0.1.0"))(
    read_nullable_string
)


###########
# Scalars #
###########


@_REGISTRY.register_read(H5Array, IOSpec("numeric-scalar", "0.2.0"))
@_REGISTRY.register_read(ZarrArray, IOSpec("numeric-scalar", "0.2.0"))
def read_scalar(elem: ArrayStorageType, *, _reader: Reader) -> np.number:
    return elem[()]


def write_scalar(
    f: GroupStorageType,
    key: str,
    value,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    return f.create_dataset(key, data=np.array(value), **dataset_kwargs)


def write_hdf5_scalar(
    f: H5Group,
    key: str,
    value,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    # Can’t compress scalars, error is thrown
    dataset_kwargs = dict(dataset_kwargs)
    dataset_kwargs.pop("compression", None)
    dataset_kwargs.pop("compression_opts", None)
    f.create_dataset(key, data=np.array(value), **dataset_kwargs)


for numeric_scalar_type in [
    *(bool, np.bool_),
    *(np.uint8, np.uint16, np.uint32, np.uint64),
    *(int, np.int8, np.int16, np.int32, np.int64),
    *(float, *np.floating.__subclasses__()),
    *np.complexfloating.__subclasses__(),
]:
    _REGISTRY.register_write(
        H5Group, numeric_scalar_type, IOSpec("numeric-scalar", "0.2.0")
    )(write_hdf5_scalar)
    _REGISTRY.register_write(
        ZarrGroup, numeric_scalar_type, IOSpec("numeric-scalar", "0.2.0")
    )(write_scalar)

_REGISTRY.register_write(ZarrGroup, str, IOSpec("string", "0.2.0"))(write_scalar)
_REGISTRY.register_write(ZarrGroup, np.str_, IOSpec("string", "0.2.0"))(write_scalar)


@_REGISTRY.register_read(H5Array, IOSpec("string", "0.2.0"))
def read_hdf5_string(elem: H5Array, *, _reader: Reader) -> str:
    return elem.asstr()[()]


@_REGISTRY.register_read(ZarrArray, IOSpec("string", "0.2.0"))
def read_zarr_string(elem: ZarrArray, *, _reader: Reader) -> str:
    return str(elem[()])


_REGISTRY.register_read(H5Array, IOSpec("bytes", "0.2.0"))(read_scalar)
_REGISTRY.register_read(ZarrArray, IOSpec("bytes", "0.2.0"))(read_scalar)


@_REGISTRY.register_write(H5Group, np.str_, IOSpec("string", "0.2.0"))
@_REGISTRY.register_write(H5Group, str, IOSpec("string", "0.2.0"))
def write_string(
    f: H5Group,
    k: str,
    v: np.str_ | str,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any],
):
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
