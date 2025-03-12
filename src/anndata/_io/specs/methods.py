from __future__ import annotations

import asyncio
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
from anndata._core.sparse_dataset import _CSCDataset, _CSRDataset, sparse_dataset
from anndata._io.utils import H5PY_V3, check_key, zero_dim_array_as_scalar
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
from ...compat import is_zarr_v2
from .registry import _REGISTRY, IOSpec, read_elem

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator
    from typing import Any, Literal, TypeVar

    from numpy import typing as npt
    from numpy.typing import NDArray

    from anndata._types import ArrayStorageType, GroupStorageType
    from anndata.compat import CSArray, CSMatrix
    from anndata.typing import AxisStorable, InMemoryArrayOrScalarType

    from .registry import Reader, Writer

    T = TypeVar("T")
    C = TypeVar("C")

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


async def sync_async_to_async(
    s: T, a: asyncio.Future[C]
) -> asyncio.Future[tuple[T, C]]:
    return s, await a


################################
# Fallbacks / backwards compat #
################################

# Note: there is no need for writing in a backwards compatible format, maybe


@_REGISTRY.register_read(H5File, IOSpec("", ""))
@_REGISTRY.register_read(H5Group, IOSpec("", ""))
@_REGISTRY.register_read(H5Array, IOSpec("", ""))
async def read_basic(
    elem: H5File | H5Group | H5Array, *, _reader: Reader
) -> dict[str, InMemoryArrayOrScalarType] | npt.NDArray | CSMatrix | CSArray:
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
        return dict(
            await asyncio.gather(
                *(
                    sync_async_to_async(k, _reader.read_elem_async(v))
                    for k, v in dict(elem).items()
                )
            )
        )
    elif isinstance(elem, h5py.Dataset):
        return h5ad.read_dataset(elem)  # TODO: Handle legacy


@_REGISTRY.register_read(ZarrGroup, IOSpec("", ""))
@_REGISTRY.register_read(ZarrArray, IOSpec("", ""))
async def read_basic_zarr(
    elem: ZarrGroup | ZarrArray, *, _reader: Reader
) -> dict[str, InMemoryArrayOrScalarType] | npt.NDArray | CSMatrix | CSArray:
    from anndata._io import zarr

    warn(
        f"Element '{elem.name}' was written without encoding metadata.",
        OldFormatWarning,
        stacklevel=3,
    )
    if isinstance(elem, ZarrGroup):
        # Backwards compat sparse arrays
        if "h5sparse_format" in elem.attrs:
            return sparse_dataset(elem).to_memory()
        return dict(
            await asyncio.gather(
                *(
                    sync_async_to_async(k, _reader.read_elem_async(v))
                    for k, v in dict(elem).items()
                )
            )
        )
    elif isinstance(elem, ZarrArray):
        return zarr.read_dataset(elem)  # TODO: Handle legacy


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


@_REGISTRY.register_write(ZarrGroup, AnnData, IOSpec("anndata", "0.1.0"))
@_REGISTRY.register_write(H5Group, AnnData, IOSpec("anndata", "0.1.0"))
async def write_anndata(
    f: GroupStorageType,
    k: str,
    adata: AnnData,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    g = f.require_group(k)
    await asyncio.gather(
        _writer.write_elem_async(g, "X", adata.X, dataset_kwargs=dataset_kwargs),
        _writer.write_elem_async(g, "obs", adata.obs, dataset_kwargs=dataset_kwargs),
        _writer.write_elem_async(g, "var", adata.var, dataset_kwargs=dataset_kwargs),
        _writer.write_elem_async(
            g, "obsm", dict(adata.obsm), dataset_kwargs=dataset_kwargs
        ),
        _writer.write_elem_async(
            g, "varm", dict(adata.varm), dataset_kwargs=dataset_kwargs
        ),
        _writer.write_elem_async(
            g, "obsp", dict(adata.obsp), dataset_kwargs=dataset_kwargs
        ),
        _writer.write_elem_async(
            g, "varp", dict(adata.varp), dataset_kwargs=dataset_kwargs
        ),
        _writer.write_elem_async(
            g, "layers", dict(adata.layers), dataset_kwargs=dataset_kwargs
        ),
        _writer.write_elem_async(
            g, "uns", dict(adata.uns), dataset_kwargs=dataset_kwargs
        ),
        _writer.write_elem_async(g, "raw", adata.raw, dataset_kwargs=dataset_kwargs),
    )


@_REGISTRY.register_read(H5Group, IOSpec("anndata", "0.1.0"))
@_REGISTRY.register_read(H5Group, IOSpec("raw", "0.1.0"))
@_REGISTRY.register_read(H5File, IOSpec("anndata", "0.1.0"))
@_REGISTRY.register_read(H5File, IOSpec("raw", "0.1.0"))
@_REGISTRY.register_read(ZarrGroup, IOSpec("anndata", "0.1.0"))
@_REGISTRY.register_read(ZarrGroup, IOSpec("raw", "0.1.0"))
async def read_anndata(elem: GroupStorageType | H5File, *, _reader: Reader) -> AnnData:
    elems = [
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
    ]
    d = dict(
        await asyncio.gather(
            *(
                sync_async_to_async(k, _reader.read_elem_async(elem[k]))
                for k in elems
                if k in elem
            )
        )
    )
    return AnnData(**d)


@_REGISTRY.register_write(H5Group, Raw, IOSpec("raw", "0.1.0"))
@_REGISTRY.register_write(ZarrGroup, Raw, IOSpec("raw", "0.1.0"))
async def write_raw(
    f: GroupStorageType,
    k: str,
    raw: Raw,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    g = f.require_group(k)
    await asyncio.gather(
        _writer.write_elem_async(g, "X", raw.X, dataset_kwargs=dataset_kwargs),
        _writer.write_elem_async(g, "var", raw.var, dataset_kwargs=dataset_kwargs),
        _writer.write_elem_async(
            g, "varm", dict(raw.varm), dataset_kwargs=dataset_kwargs
        ),
    )


########
# Null #
########


@_REGISTRY.register_read(H5Array, IOSpec("null", "0.1.0"))
@_REGISTRY.register_read(ZarrArray, IOSpec("null", "0.1.0"))
async def read_null(_elem, _reader) -> None:
    return None


@_REGISTRY.register_write(H5Group, type(None), IOSpec("null", "0.1.0"))
async def write_null_h5py(f, k, _v, _writer, dataset_kwargs=MappingProxyType({})):
    f.create_dataset(k, data=h5py.Empty("f"), **dataset_kwargs)


@_REGISTRY.register_write(ZarrGroup, type(None), IOSpec("null", "0.1.0"))
async def write_null_zarr(f, k, _v, _writer, dataset_kwargs=MappingProxyType({})):
    # zarr has no first-class null dataset
    if is_zarr_v2():
        import zarr

        # zarr has no first-class null dataset
        f.create_dataset(k, data=zarr.empty(()), **dataset_kwargs)
    else:
        # TODO: why is this not actually storing the empty info with a f.empty call?
        # It fails complaining that k doesn't exist when updating the attributes.
        f.create_array(k, shape=(), dtype="bool")


############
# Mappings #
############


@_REGISTRY.register_read(H5Group, IOSpec("dict", "0.1.0"))
@_REGISTRY.register_read(ZarrGroup, IOSpec("dict", "0.1.0"))
async def read_mapping(
    elem: GroupStorageType, *, _reader: Reader
) -> dict[str, AxisStorable]:
    return dict(
        await asyncio.gather(
            *(
                sync_async_to_async(k, _reader.read_elem_async(v))
                for k, v in dict(elem).items()
            )
        )
    )


@_REGISTRY.register_write(H5Group, dict, IOSpec("dict", "0.1.0"))
@_REGISTRY.register_write(ZarrGroup, dict, IOSpec("dict", "0.1.0"))
async def write_mapping(
    f: GroupStorageType,
    k: str,
    v: dict[str, AxisStorable],
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    g = f.require_group(k)
    await asyncio.gather(
        *[
            _writer.write_elem_async(g, sub_k, sub_v, dataset_kwargs=dataset_kwargs)
            for sub_k, sub_v in v.items()
        ]
    )


##############
# np.ndarray #
##############


@_REGISTRY.register_write(H5Group, list, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, list, IOSpec("array", "0.2.0"))
async def write_list(
    f: GroupStorageType,
    k: str,
    elem: list[AxisStorable],
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    await _writer.write_elem_async(f, k, np.array(elem), dataset_kwargs=dataset_kwargs)


# TODO: Is this the right behavior for MaskedArrays?
# It's in the `AnnData.concatenate` docstring, but should we keep it?
@_REGISTRY.register_write(H5Group, views.ArrayView, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(H5Group, np.ndarray, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(H5Group, np.ma.MaskedArray, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, views.ArrayView, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, np.ndarray, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, np.ma.MaskedArray, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, ZarrArray, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, H5Array, IOSpec("array", "0.2.0"))
@zero_dim_array_as_scalar
async def write_basic(
    f: GroupStorageType,
    k: str,
    elem: views.ArrayView | np.ndarray | h5py.Dataset | np.ma.MaskedArray | ZarrArray,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    """Write methods which underlying library handles natively."""
    dataset_kwargs = dataset_kwargs.copy()
    dtype = dataset_kwargs.pop("dtype", elem.dtype)
    if isinstance(f, H5Group) or is_zarr_v2():
        f.create_dataset(k, data=elem, shape=elem.shape, dtype=dtype, **dataset_kwargs)
    else:
        f.create_array(k, shape=elem.shape, dtype=dtype, **dataset_kwargs)
        # see https://github.com/zarr-developers/zarr-python/discussions/2712
        if isinstance(elem, ZarrArray):
            await f[k]._async_array.setitem(Ellipsis, elem[...])
        else:
            f[k][...] = elem


def _iter_chunks_for_copy(
    elem: ArrayStorageType, dest: ArrayStorageType
) -> Iterator[slice | tuple[list[slice]]]:
    """
    Returns an iterator of tuples of slices for copying chunks from `elem` to `dest`.

    * If `dest` has chunks, it will return the chunks of `dest`.
    * If `dest` is not chunked, we write it in ~100MB chunks or 1000 rows, whichever is larger.
    """
    if dest.chunks and hasattr(dest, "iter_chunks"):
        return dest.iter_chunks()
    else:
        shape = elem.shape
        # Number of rows that works out to
        n_rows = max(
            ad.settings.min_rows_for_chunked_h5_copy,
            elem.chunks[0] if elem.chunks is not None else 1,
        )
        return (slice(i, min(i + n_rows, shape[0])) for i in range(0, shape[0], n_rows))


@_REGISTRY.register_write(H5Group, H5Array, IOSpec("array", "0.2.0"))
@_REGISTRY.register_write(H5Group, ZarrArray, IOSpec("array", "0.2.0"))
async def write_chunked_dense_array_to_group(
    f: GroupStorageType,
    k: str,
    elem: ArrayStorageType,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    """Write to a h5py.Dataset in chunks.

    `h5py.Group.create_dataset(..., data: h5py.Dataset)` will load all of `data` into memory
    before writing. Instead, we will write in chunks to avoid this. We don't need to do this for
    zarr since zarr handles this automatically.
    """
    dtype = dataset_kwargs.get("dtype", elem.dtype)
    kwargs = {**dataset_kwargs, "dtype": dtype}
    dest = f.create_dataset(k, shape=elem.shape, **kwargs)
    chunk_iter = _iter_chunks_for_copy(elem, dest)
    if isinstance(dest, ZarrArray) and not is_zarr_v2():
        await asyncio.gather(
            *(
                dest._async_array.setitem(chunk, elem[chunk])
                for chunk in _iter_chunks_for_copy(elem, dest)
            )
        )
    else:
        for chunk in chunk_iter:
            dest[chunk] = elem[chunk]


_REGISTRY.register_write(H5Group, CupyArray, IOSpec("array", "0.2.0"))(
    _to_cpu_mem_wrapper(write_basic)
)
_REGISTRY.register_write(ZarrGroup, CupyArray, IOSpec("array", "0.2.0"))(
    _to_cpu_mem_wrapper(write_basic)
)


@_REGISTRY.register_write(ZarrGroup, DaskArray, IOSpec("array", "0.2.0"))
async def write_basic_dask_zarr(
    f: ZarrGroup,
    k: str,
    elem: DaskArray,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    import dask.array as da

    if is_zarr_v2():
        g = f.require_dataset(k, shape=elem.shape, dtype=elem.dtype, **dataset_kwargs)
    else:
        g = f.require_array(k, shape=elem.shape, dtype=elem.dtype, **dataset_kwargs)
    da.store(elem, g, lock=GLOBAL_LOCK)


# Adding this separately because h5py isn't serializable
# https://github.com/pydata/xarray/issues/4242
@_REGISTRY.register_write(H5Group, DaskArray, IOSpec("array", "0.2.0"))
async def write_basic_dask_h5(
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
        msg = "Cannot write dask arrays to hdf5 when using distributed scheduler"
        raise ValueError(msg)

    g = f.require_dataset(k, shape=elem.shape, dtype=elem.dtype, **dataset_kwargs)
    da.store(elem, g)


@_REGISTRY.register_read(H5Array, IOSpec("array", "0.2.0"))
@_REGISTRY.register_read(ZarrArray, IOSpec("array", "0.2.0"))
@_REGISTRY.register_read(ZarrArray, IOSpec("string-array", "0.2.0"))
async def read_array(elem: ArrayStorageType, *, _reader: Reader) -> npt.NDArray:
    if not is_zarr_v2() and isinstance(elem, ZarrArray):
        return await elem._async_array.getitem(())
    return elem[()]


# arrays of strings
@_REGISTRY.register_read(H5Array, IOSpec("string-array", "0.2.0"))
async def read_string_array(d: H5Array, *, _reader: Reader):
    return await read_array(d.asstr(), _reader=_reader)


@_REGISTRY.register_write(
    H5Group, (views.ArrayView, "U"), IOSpec("string-array", "0.2.0")
)
@_REGISTRY.register_write(
    H5Group, (views.ArrayView, "O"), IOSpec("string-array", "0.2.0")
)
@_REGISTRY.register_write(H5Group, (np.ndarray, "U"), IOSpec("string-array", "0.2.0"))
@_REGISTRY.register_write(H5Group, (np.ndarray, "O"), IOSpec("string-array", "0.2.0"))
@zero_dim_array_as_scalar
async def write_vlen_string_array(
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
@zero_dim_array_as_scalar
async def write_vlen_string_array_zarr(
    f: ZarrGroup,
    k: str,
    elem: np.ndarray,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    if is_zarr_v2():
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
    else:
        from numcodecs import VLenUTF8

        dataset_kwargs = dataset_kwargs.copy()
        compressor = None
        if ad.settings.zarr_write_format == 2:
            compressor = dataset_kwargs.pop("compressor", None)

        filters, dtype = (
            ([VLenUTF8()], object)
            if ad.settings.zarr_write_format == 2
            else (None, str)
        )
        f.create_array(
            k,
            shape=elem.shape,
            dtype=dtype,
            filters=filters,
            compressor=compressor,
            **dataset_kwargs,
        )
        await f[k]._async_array.setitem(slice(None), elem)


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
async def read_recarray(
    d: ArrayStorageType, *, _reader: Reader
) -> np.recarray | npt.NDArray:
    if not is_zarr_v2() and isinstance(d, ZarrArray):
        value = await d._async_array.getitem(())
    else:
        value = d[()]
    dtype = value.dtype
    value = _from_fixed_length_strings(value)
    if H5PY_V3:
        value = _decode_structured_array(value, dtype=dtype)
    return value


@_REGISTRY.register_write(H5Group, (np.ndarray, "V"), IOSpec("rec-array", "0.2.0"))
@_REGISTRY.register_write(H5Group, np.recarray, IOSpec("rec-array", "0.2.0"))
async def write_recarray(
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
async def write_recarray_zarr(
    f: ZarrGroup,
    k: str,
    elem: np.ndarray | np.recarray,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    from anndata.compat import _to_fixed_length_strings

    elem = _to_fixed_length_strings(elem)
    if isinstance(f, H5Group) or is_zarr_v2():
        f.create_dataset(k, data=elem, shape=elem.shape, **dataset_kwargs)
    else:
        # TODO: zarr’s on-disk format v3 doesn’t support this dtype
        f.create_array(k, shape=elem.shape, dtype=elem.dtype, **dataset_kwargs)
        await f[k]._async_array.setitem(Ellipsis, elem)


#################
# Sparse arrays #
#################


async def write_sparse_compressed(
    f: GroupStorageType,
    key: str,
    value: CSMatrix | CSArray,
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
    if isinstance(f, H5Group):
        dataset_kwargs = dict(maxshape=(None,), **dataset_kwargs)
    awaitables = []
    for attr_name in ["data", "indices", "indptr"]:
        attr = getattr(value, attr_name)
        dtype = indptr_dtype if attr_name == "indptr" else attr.dtype
        awaitables.append(
            _writer.write_elem_async(
                g,
                attr_name,
                attr,
                dataset_kwargs={"dtype": dtype, **dataset_kwargs},
            )
        )
    await asyncio.gather(*awaitables)


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
async def write_sparse_dataset(
    f: GroupStorageType,
    k: str,
    elem: _CSCDataset | _CSRDataset,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    await write_sparse_compressed(
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
async def write_cupy_dask_sparse(
    f, k, elem, _writer, dataset_kwargs=MappingProxyType({})
):
    await _writer.write_elem_async(
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
async def write_dask_sparse(
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
        msg = f"Cannot write dask sparse arrays with format {sparse_format}"
        raise NotImplementedError(msg)

    def chunk_slice(start: int, stop: int) -> tuple[slice | None, slice | None]:
        result = [slice(None), slice(None)]
        result[axis] = slice(start, stop)
        return tuple(result)

    axis_chunks = elem.chunks[axis]
    chunk_start = 0
    chunk_stop = axis_chunks[0]

    await _writer.write_elem_async(
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
async def read_sparse(elem: GroupStorageType, *, _reader: Reader) -> CSMatrix | CSArray:
    return await sparse_dataset(elem).to_memory_async()


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
async def write_awkward(
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
    await asyncio.gather(
        *[
            _writer.write_elem_async(group, k, v, dataset_kwargs=dataset_kwargs)
            for k, v in container.items()
        ]
    )


@_REGISTRY.register_read(H5Group, IOSpec("awkward-array", "0.1.0"))
@_REGISTRY.register_read(ZarrGroup, IOSpec("awkward-array", "0.1.0"))
async def read_awkward(elem: GroupStorageType, *, _reader: Reader) -> AwkArray:
    from anndata.compat import awkward as ak

    form = _read_attr(elem.attrs, "form")
    length = _read_attr(elem.attrs, "length")
    container = dict(
        await asyncio.gather(
            *(
                sync_async_to_async(k, _reader.read_elem_async(elem[k]))
                for k in elem.keys()
            )
        )
    )

    return ak.from_buffers(form, int(length), container)


##############
# DataFrames #
##############


@_REGISTRY.register_write(H5Group, views.DataFrameView, IOSpec("dataframe", "0.2.0"))
@_REGISTRY.register_write(H5Group, pd.DataFrame, IOSpec("dataframe", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, views.DataFrameView, IOSpec("dataframe", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, pd.DataFrame, IOSpec("dataframe", "0.2.0"))
async def write_dataframe(
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
            msg = f"{reserved!r} is a reserved name for dataframe columns."
            raise ValueError(msg)
    group = _require_group_write_dataframe(f, key, df)
    if not df.columns.is_unique:
        duplicates = list(df.columns[df.columns.duplicated()])
        msg = f"Found repeated column names: {duplicates}. Column names must be unique."
        raise ValueError(msg)
    col_names = [check_key(c) for c in df.columns]
    group.attrs["column-order"] = col_names

    if df.index.name is not None:
        if df.index.name in col_names and not pd.Series(
            df.index, index=df.index
        ).equals(df[df.index.name]):
            msg = (
                f"DataFrame.index.name ({df.index.name!r}) is also used by a column "
                "whose values are different. This is not supported. Please make sure "
                "the values are the same, or use a different name."
            )
            raise ValueError(msg)
        index_name = df.index.name
    else:
        index_name = "_index"
    group.attrs["_index"] = check_key(index_name)

    # ._values is "the best" array representation. It's the true array backing the
    # object, where `.values` is always a np.ndarray and .array is always a pandas
    # array.
    awaitables = [
        _writer.write_elem_async(
            group, index_name, df.index._values, dataset_kwargs=dataset_kwargs
        )
    ]
    # TODO: this should write the "true" representation of the series (i.e. the underlying array or ndarray depending)
    awaitables += [
        _writer.write_elem_async(
            group, colname, series._values, dataset_kwargs=dataset_kwargs
        )
        for colname, series in df.items()
    ]
    await asyncio.gather(*awaitables)


@_REGISTRY.register_read(H5Group, IOSpec("dataframe", "0.2.0"))
@_REGISTRY.register_read(ZarrGroup, IOSpec("dataframe", "0.2.0"))
async def read_dataframe(elem: GroupStorageType, *, _reader: Reader) -> pd.DataFrame:
    columns = list(_read_attr(elem.attrs, "column-order"))
    idx_key = _read_attr(elem.attrs, "_index")
    df = pd.DataFrame(
        dict(
            await asyncio.gather(
                *(
                    sync_async_to_async(k, read_series(elem[k], _reader))
                    for k in columns
                )
            )
        ),
        index=await _reader.read_elem_async(elem[idx_key]),
        columns=columns if len(columns) else None,
    )
    if idx_key != "_index":
        df.index.name = idx_key
    return df


# Backwards compat dataframe reading


@_REGISTRY.register_read(H5Group, IOSpec("dataframe", "0.1.0"))
@_REGISTRY.register_read(ZarrGroup, IOSpec("dataframe", "0.1.0"))
async def read_dataframe_0_1_0(
    elem: GroupStorageType, *, _reader: Reader
) -> pd.DataFrame:
    columns = _read_attr(elem.attrs, "column-order")
    idx_key = _read_attr(elem.attrs, "_index")
    df = pd.DataFrame(
        dict(
            await asyncio.gather(
                *(
                    sync_async_to_async(k, read_series(elem[k], _reader))
                    for k in columns
                )
            )
        ),
        index=await read_series(elem[idx_key], _reader),
        columns=columns if len(columns) else None,
    )
    if idx_key != "_index":
        df.index.name = idx_key
    return df


async def read_series(
    dataset: h5py.Dataset, _reader: Reader
) -> np.ndarray | pd.Categorical:
    # For reading older dataframes
    if "categories" in dataset.attrs:
        if isinstance(dataset, ZarrArray):
            import zarr

            parent_name = dataset.name.rstrip(dataset.basename).strip("/")
            parent = zarr.open(dataset.store, mode="r")[parent_name]
        else:
            parent = dataset.parent
        categories_dset = parent[_read_attr(dataset.attrs, "categories")]
        categories, codes = await asyncio.gather(
            *(
                _reader.read_elem_async(categories_dset),
                _reader.read_elem_async(dataset),
            )
        )
        ordered = bool(_read_attr(categories_dset.attrs, "ordered", default=False))
        return pd.Categorical.from_codes(codes, categories, ordered=ordered)
    else:
        return await _reader.read_elem_async(dataset)


###############
# Categorical #
###############


@_REGISTRY.register_write(H5Group, pd.Categorical, IOSpec("categorical", "0.2.0"))
@_REGISTRY.register_write(ZarrGroup, pd.Categorical, IOSpec("categorical", "0.2.0"))
async def write_categorical(
    f: GroupStorageType,
    k: str,
    v: pd.Categorical,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    g = f.require_group(k)
    g.attrs["ordered"] = bool(v.ordered)

    await asyncio.gather(
        _writer.write_elem_async(g, "codes", v.codes, dataset_kwargs=dataset_kwargs),
        _writer.write_elem_async(
            g, "categories", v.categories._values, dataset_kwargs=dataset_kwargs
        ),
    )


@_REGISTRY.register_read(H5Group, IOSpec("categorical", "0.2.0"))
@_REGISTRY.register_read(ZarrGroup, IOSpec("categorical", "0.2.0"))
async def read_categorical(
    elem: GroupStorageType, *, _reader: Reader
) -> pd.Categorical:
    codes, categories = await asyncio.gather(
        *(_reader.read_elem_async(elem[k]) for k in ["codes", "categories"])
    )
    ordered = bool(_read_attr(elem.attrs, "ordered"))
    return pd.Categorical.from_codes(
        codes=codes, categories=categories, ordered=ordered
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
async def write_nullable(
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
    await asyncio.gather(
        _writer.write_elem_async(g, "values", values, dataset_kwargs=dataset_kwargs),
        _writer.write_elem_async(g, "mask", v.isna(), dataset_kwargs=dataset_kwargs),
    )


async def _read_nullable(
    elem: GroupStorageType,
    *,
    _reader: Reader,
    # BaseMaskedArray
    array_type: Callable[
        [NDArray[np.number], NDArray[np.bool_]], pd.api.extensions.ExtensionArray
    ],
) -> pd.api.extensions.ExtensionArray:
    values, mask = await asyncio.gather(
        *(_reader.read_elem_async(elem[k]) for k in ["values", "mask"])
    )
    return array_type(values, mask=mask)


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
async def read_scalar(elem: ArrayStorageType, *, _reader: Reader) -> np.number:
    # TODO: `item` ensures the return is in fact a scalar (needed after zarr v3 which now returns a 1 elem array)
    # https://github.com/zarr-developers/zarr-python/issues/2713
    if not is_zarr_v2() and isinstance(elem, ZarrArray):
        return (await elem._async_array.getitem(())).item()
    return elem[()].item()


def _remove_scalar_compression_args(dataset_kwargs: Mapping[str, Any]) -> dict:
    # Can’t compress scalars, error is thrown
    dataset_kwargs = dict(dataset_kwargs)
    for arg in (
        "compression",
        "compression_opts",
        "chunks",
        "shuffle",
        "fletcher32",
        "scaleoffset",
        "compressor",
    ):
        dataset_kwargs.pop(arg, None)
    return dataset_kwargs


async def write_scalar_zarr(
    f: ZarrGroup,
    key: str,
    value,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    # these args are ignored in v2: https://zarr.readthedocs.io/en/v2.18.4/api/hierarchy.html#zarr.hierarchy.Group.create_dataset
    # and error out in v3
    dataset_kwargs = _remove_scalar_compression_args(dataset_kwargs)
    if is_zarr_v2():
        return f.create_dataset(key, data=np.array(value), shape=(), **dataset_kwargs)
    else:
        from numcodecs import VLenUTF8

        match ad.settings.zarr_write_format, value:
            case 2, str():
                filters, dtype = [VLenUTF8()], object
            case 3, str():
                filters, dtype = None, str
            case _, _:
                filters, dtype = None, np.array(value).dtype
        a = f.create_array(
            key,
            shape=(),
            dtype=dtype,
            filters=filters,
            **dataset_kwargs,
        )
        a[...] = np.array(value)


async def write_hdf5_scalar(
    f: H5Group,
    key: str,
    value,
    *,
    _writer: Writer,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    # Can’t compress scalars, error is thrown
    dataset_kwargs = _remove_scalar_compression_args(dataset_kwargs)
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
    )(write_scalar_zarr)

_REGISTRY.register_write(ZarrGroup, str, IOSpec("string", "0.2.0"))(write_scalar_zarr)
_REGISTRY.register_write(ZarrGroup, np.str_, IOSpec("string", "0.2.0"))(
    write_scalar_zarr
)


@_REGISTRY.register_read(H5Array, IOSpec("string", "0.2.0"))
async def read_hdf5_string(elem: H5Array, *, _reader: Reader) -> str:
    return elem.asstr()[()]


@_REGISTRY.register_read(ZarrArray, IOSpec("string", "0.2.0"))
async def read_zarr_string(elem: ZarrArray, *, _reader: Reader) -> str:
    if not is_zarr_v2() and isinstance(elem, ZarrArray):
        return str(await elem._async_array.getitem(()))
    return str(elem[()])


_REGISTRY.register_read(H5Array, IOSpec("bytes", "0.2.0"))(read_scalar)
_REGISTRY.register_read(ZarrArray, IOSpec("bytes", "0.2.0"))(read_scalar)


@_REGISTRY.register_write(H5Group, np.str_, IOSpec("string", "0.2.0"))
@_REGISTRY.register_write(H5Group, str, IOSpec("string", "0.2.0"))
async def write_string(
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
# async def write_string(f, k, v, dataset_kwargs):
#     if "compression" in dataset_kwargs:
#         dataset_kwargs = dict(dataset_kwargs)
#         dataset_kwargs.pop("compression")
#     f.create_dataset(k, data=np.array(v), **dataset_kwargs)
