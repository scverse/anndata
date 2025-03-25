from __future__ import annotations

import asyncio
from copy import copy
from functools import partial
from types import MappingProxyType
from typing import TYPE_CHECKING
from warnings import warn

import numpy as np
import pandas as pd

import anndata as ad
from anndata import AnnData, Raw
from anndata._core import views
from anndata._core.aligned_mapping import AlignedView
from anndata._core.sparse_dataset import _CSCDataset, _CSRDataset, sparse_dataset
from anndata._io.specs.methods import read_series
from anndata._io.utils import check_key, zero_dim_array_as_scalar
from anndata._warnings import OldFormatWarning
from anndata.compat import (
    H5Array,
    ZarrAsyncArray,
    ZarrAsyncGroup,
    ZarrGroup,
    _from_fixed_length_strings,
    _read_attr,
)

from ..._settings import settings
from .registry import _ASYNC_REGISTRY, IOSpec, read_elem

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping
    from typing import Any, Literal, TypeVar

    import h5py
    from numpy import typing as npt
    from numpy.typing import NDArray

    from anndata._types import ArrayStorageType
    from anndata.compat import (
        AwkArray,
        CSArray,
        CSMatrix,
        CupyArray,
        CupyCSCMatrix,
        CupyCSRMatrix,
        ZarrArray,
    )
    from anndata.typing import AxisStorable, InMemoryArrayOrScalarType

    from .registry import AsyncReader, AsyncWriter

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
        _writer: AsyncWriter,
        dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
    ):
        return write_func(
            f, k, cupy_val.get(), _writer=_writer, dataset_kwargs=dataset_kwargs
        )

    return wrapper


@_ASYNC_REGISTRY.register_read(ZarrAsyncGroup, IOSpec("", ""))
@_ASYNC_REGISTRY.register_read(ZarrAsyncArray, IOSpec("", ""))
async def read_basic_zarr(
    elem: ZarrAsyncArray | ZarrAsyncGroup, *, _reader: AsyncReader
) -> dict[str, InMemoryArrayOrScalarType] | npt.NDArray | CSMatrix | CSArray:
    from anndata._io import zarr

    warn(
        f"Element '{elem.name}' was written without encoding metadata.",
        OldFormatWarning,
        stacklevel=3,
    )
    if isinstance(elem, ZarrGroup | ZarrAsyncGroup):
        # Backwards compat sparse arrays
        if "h5sparse_format" in elem.attrs:
            return await sparse_dataset(elem).to_memory_async()
        k_v = [e async for e in await elem.members()]
        return dict(
            zip(
                (k for k, _ in k_v),
                await asyncio.gather(*(_reader.read_elem(v) for _, v in k_v)),
            )
        )
    return await zarr.read_dataset(elem)  # TODO: Handle legacy


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


async def write_dict(
    f: ZarrAsyncGroup,
    key: str,
    elem,
    *,
    _writer: AsyncWriter,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
) -> None:
    if isinstance(elem, AlignedView):
        d = dict(
            zip(
                elem.keys(),
                await asyncio.gather(*(elem.getitem(k) for k in elem.keys())),
            )
        )
    else:
        d = dict(elem)
    await _writer.write_elem(f, key, d, dataset_kwargs=dataset_kwargs)


async def write_X(
    f: ZarrAsyncGroup,
    elem: AnnData | Raw,
    *,
    _writer: AsyncWriter,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    X = await elem.get_X()
    await _writer.write_elem(f, "X", X, dataset_kwargs=dataset_kwargs)


@_ASYNC_REGISTRY.register_write(ZarrAsyncGroup, AnnData, IOSpec("anndata", "0.1.0"))
async def write_anndata(
    f: ZarrAsyncGroup,
    k: str,
    adata: AnnData,
    *,
    _writer: AsyncWriter,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    g = await f.require_group(k)
    await asyncio.gather(
        write_X(g, adata, _writer=_writer, dataset_kwargs=dataset_kwargs),
        _writer.write_elem(g, "obs", adata.obs, dataset_kwargs=dataset_kwargs),
        _writer.write_elem(g, "var", adata.var, dataset_kwargs=dataset_kwargs),
        *(
            write_dict(
                g,
                elem_key,
                getattr(adata, elem_key),
                _writer=_writer,
                dataset_kwargs=dataset_kwargs,
            )
            for elem_key in ["obsm", "varm", "obsp", "varp", "layers", "uns"]
        ),
        _writer.write_elem(g, "raw", adata.raw, dataset_kwargs=dataset_kwargs),
    )


@_ASYNC_REGISTRY.register_read(ZarrAsyncGroup, IOSpec("anndata", "0.1.0"))
@_ASYNC_REGISTRY.register_read(ZarrAsyncGroup, IOSpec("raw", "0.1.0"))
async def read_anndata(elem: ZarrAsyncGroup, *, _reader: AsyncReader) -> AnnData:
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

    async def read_async(e, k):
        child = await e.get(k)
        return await _reader.read_elem(child)

    elem_keys = await elem.keys()
    elem_keys_filtered = [k async for k in elem_keys if k in elems]
    d = dict(
        zip(
            elem_keys_filtered,
            await asyncio.gather(*(read_async(elem, k) for k in elem_keys_filtered)),
        )
    )
    return AnnData(**d)


@_ASYNC_REGISTRY.register_write(ZarrAsyncGroup, Raw, IOSpec("raw", "0.1.0"))
async def write_raw(
    f: ZarrAsyncGroup,
    k: str,
    raw: Raw,
    *,
    _writer: AsyncWriter,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    g = await f.require_group(k)
    await asyncio.gather(
        write_X(g, raw, _writer=_writer, dataset_kwargs=dataset_kwargs),
        _writer.write_elem(g, "var", raw.var, dataset_kwargs=dataset_kwargs),
        write_dict(g, "varm", raw.varm, _writer=_writer, dataset_kwargs=dataset_kwargs),
    )


########
# Null #
########


@_ASYNC_REGISTRY.register_read(ZarrAsyncArray, IOSpec("null", "0.1.0"))
async def read_null(_elem, _reader) -> None:
    return None


@_ASYNC_REGISTRY.register_write(ZarrAsyncGroup, type(None), IOSpec("null", "0.1.0"))
async def write_null_zarr(f, k, _v, _writer, dataset_kwargs=MappingProxyType({})):
    await f.create_array(k, shape=(), dtype="bool")


############
# Mappings #
############


@_ASYNC_REGISTRY.register_read(ZarrAsyncGroup, IOSpec("dict", "0.1.0"))
async def read_mapping(
    elem: ZarrAsyncGroup, *, _reader: AsyncReader
) -> dict[str, AxisStorable]:
    print(elem)
    k_v = [(k, v) async for k, v in await elem.members()]
    return dict(
        zip(
            (k for k, _ in k_v),
            await asyncio.gather(*(_reader.read_elem(v) for _, v in k_v)),
        )
    )


@_ASYNC_REGISTRY.register_write(ZarrAsyncGroup, dict, IOSpec("dict", "0.1.0"))
async def write_mapping(
    f: ZarrAsyncGroup,
    k: str,
    v: dict[str, AxisStorable],
    *,
    _writer: AsyncWriter,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    g = await f.require_group(k)
    await asyncio.gather(
        *[
            _writer.write_elem(g, sub_k, sub_v, dataset_kwargs=dataset_kwargs)
            for sub_k, sub_v in v.items()
        ]
    )


##############
# np.ndarray #
##############


@_ASYNC_REGISTRY.register_write(ZarrAsyncGroup, list, IOSpec("array", "0.2.0"))
async def write_list(
    f: ZarrAsyncGroup,
    k: str,
    elem: list[AxisStorable],
    *,
    _writer: AsyncWriter,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    await _writer.write_elem(f, k, np.array(elem), dataset_kwargs=dataset_kwargs)


# TODO: Is this the right behavior for MaskedArrays?
# It's in the `AnnData.concatenate` docstring, but should we keep it?
@_ASYNC_REGISTRY.register_write(
    ZarrAsyncGroup, views.ArrayView, IOSpec("array", "0.2.0")
)
@_ASYNC_REGISTRY.register_write(ZarrAsyncGroup, np.ndarray, IOSpec("array", "0.2.0"))
@_ASYNC_REGISTRY.register_write(
    ZarrAsyncGroup, np.ma.MaskedArray, IOSpec("array", "0.2.0")
)
@_ASYNC_REGISTRY.register_write(
    ZarrAsyncGroup, ZarrAsyncArray, IOSpec("array", "0.2.0")
)
@_ASYNC_REGISTRY.register_write(ZarrAsyncGroup, H5Array, IOSpec("array", "0.2.0"))
@zero_dim_array_as_scalar
async def write_basic(
    f: ZarrAsyncGroup,
    k: str,
    elem: views.ArrayView | np.ndarray | h5py.Dataset | np.ma.MaskedArray | ZarrArray,
    *,
    _writer: AsyncWriter,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    """Write methods which underlying library handles natively."""
    dataset_kwargs = dataset_kwargs.copy()
    dtype = dataset_kwargs.pop("dtype", elem.dtype)
    arr = await f.create_array(k, shape=elem.shape, dtype=dtype, **dataset_kwargs)
    # see https://github.com/zarr-developers/zarr-python/discussions/2712
    if isinstance(elem, ZarrAsyncArray):
        elem = await elem.getitem(())
    await arr.setitem(Ellipsis, elem)


@_ASYNC_REGISTRY.register_read(ZarrAsyncArray, IOSpec("array", "0.2.0"))
@_ASYNC_REGISTRY.register_read(ZarrAsyncArray, IOSpec("string-array", "0.2.0"))
async def read_array(elem: ArrayStorageType, *, _reader: AsyncReader) -> npt.NDArray:
    return await elem.getitem(())


@_ASYNC_REGISTRY.register_write(
    ZarrAsyncGroup, (np.ndarray, "U"), IOSpec("string-array", "0.2.0")
)
@_ASYNC_REGISTRY.register_write(
    ZarrAsyncGroup, (np.ndarray, "O"), IOSpec("string-array", "0.2.0")
)
@zero_dim_array_as_scalar
async def write_vlen_string_array_zarr(
    f: ZarrAsyncGroup,
    k: str,
    elem: np.ndarray,
    *,
    _writer: AsyncWriter,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    from numcodecs import VLenUTF8

    dataset_kwargs = dataset_kwargs.copy()
    compressor = None
    if ad.settings.zarr_write_format == 2:
        compressor = dataset_kwargs.pop("compressor", None)

    filters, dtype = (
        ([VLenUTF8()], object) if ad.settings.zarr_write_format == 2 else (None, str)
    )
    arr = await f.create_array(
        k,
        shape=elem.shape,
        dtype=dtype,
        filters=filters,
        compressor=compressor,
        **dataset_kwargs,
    )
    await arr.setitem(slice(None), elem)


###############
# np.recarray #
###############


@_ASYNC_REGISTRY.register_read(ZarrAsyncArray, IOSpec("rec-array", "0.2.0"))
async def read_recarray(
    d: ZarrAsyncArray, *, _reader: AsyncReader
) -> np.recarray | npt.NDArray:
    value = await d.getitem(())
    value = _from_fixed_length_strings(value)
    return value


@_ASYNC_REGISTRY.register_write(
    ZarrAsyncGroup, (np.ndarray, "V"), IOSpec("rec-array", "0.2.0")
)
@_ASYNC_REGISTRY.register_write(
    ZarrAsyncGroup, np.recarray, IOSpec("rec-array", "0.2.0")
)
async def write_recarray_zarr(
    f: ZarrAsyncGroup,
    k: str,
    elem: np.ndarray | np.recarray,
    *,
    _writer: AsyncWriter,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    from anndata.compat import _to_fixed_length_strings

    elem = _to_fixed_length_strings(elem)
    # TODO: zarr’s on-disk format v3 doesn’t support this dtype
    arr = await f.create_array(k, shape=elem.shape, dtype=elem.dtype, **dataset_kwargs)

    await arr.setitem(Ellipsis, elem)


#################
# Sparse arrays #
#################


async def write_sparse_compressed(
    f: ZarrAsyncGroup,
    key: str,
    value: CSMatrix | CSArray,
    *,
    _writer: AsyncWriter,
    fmt: Literal["csr", "csc"],
    dataset_kwargs=MappingProxyType({}),
):
    g = await f.require_group(key)
    g.attrs["shape"] = value.shape
    dataset_kwargs = dict(dataset_kwargs)
    indptr_dtype = dataset_kwargs.pop("indptr_dtype", value.indptr.dtype)

    awaitables = []
    for attr_name in ["data", "indices", "indptr"]:
        attr = getattr(value, attr_name)
        dtype = indptr_dtype if attr_name == "indptr" else attr.dtype
        arr = await g.create_array(
            attr_name, shape=attr.shape, dtype=dtype, **dataset_kwargs
        )
        # see https://github.com/zarr-developers/zarr-python/discussions/2712
        awaitables.append(arr.setitem(Ellipsis, attr[...]))
    await asyncio.gather(*awaitables)


write_csr = partial(write_sparse_compressed, fmt="csr")
write_csc = partial(write_sparse_compressed, fmt="csc")


@_ASYNC_REGISTRY.register_write(ZarrAsyncGroup, _CSRDataset, IOSpec("", "0.1.0"))
@_ASYNC_REGISTRY.register_write(ZarrAsyncGroup, _CSCDataset, IOSpec("", "0.1.0"))
async def write_sparse_dataset(
    f: ZarrAsyncGroup,
    k: str,
    elem: _CSCDataset | _CSRDataset,
    *,
    _writer: AsyncWriter,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    await write_sparse_compressed(
        f,
        k,
        await elem.to_memory_async(),
        _writer=_writer,
        fmt=elem.format,
        dataset_kwargs=dataset_kwargs,
    )
    # TODO: Cleaner way to do this
    f[k].attrs["encoding-type"] = f"{elem.format}_matrix"
    f[k].attrs["encoding-version"] = "0.1.0"


@_ASYNC_REGISTRY.register_read(ZarrAsyncGroup, IOSpec("csc_matrix", "0.1.0"))
@_ASYNC_REGISTRY.register_read(ZarrAsyncGroup, IOSpec("csr_matrix", "0.1.0"))
async def read_sparse(
    elem: ZarrAsyncGroup, *, _reader: AsyncReader
) -> CSMatrix | CSArray:
    return await sparse_dataset(elem).to_memory_async()


#################
# Awkward array #
#################


@_ASYNC_REGISTRY.register_write(
    ZarrAsyncGroup, views.AwkwardArrayView, IOSpec("awkward-array", "0.1.0")
)
async def write_awkward(
    f: ZarrAsyncGroup,
    k: str,
    v: views.AwkwardArrayView | AwkArray,
    *,
    _writer: AsyncWriter,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    from anndata.compat import awkward as ak

    group = await f.require_group(k)
    if isinstance(v, views.AwkwardArrayView):
        # copy to remove the view attributes
        v = copy(v)
    form, length, container = ak.to_buffers(ak.to_packed(v))
    group.attrs["length"] = length
    group.attrs["form"] = form.to_json()
    await asyncio.gather(
        *[
            _writer.write_elem(group, k, v, dataset_kwargs=dataset_kwargs)
            for k, v in container.items()
        ]
    )


@_ASYNC_REGISTRY.register_read(ZarrAsyncGroup, IOSpec("awkward-array", "0.1.0"))
async def read_awkward(elem: ZarrAsyncGroup, *, _reader: AsyncReader) -> AwkArray:
    from anndata.compat import awkward as ak

    form = _read_attr(elem.attrs, "form")
    length = _read_attr(elem.attrs, "length")
    k_v = [k_v async for k_v in elem.members()]
    container = dict(
        zip(
            (k for k, _ in k_v),
            await asyncio.gather(*(_reader.read_elem(v) for _, v in k_v)),
        )
    )

    return ak.from_buffers(form, int(length), container)


##############
# DataFrames #
##############


@_ASYNC_REGISTRY.register_write(
    ZarrAsyncGroup, views.DataFrameView, IOSpec("dataframe", "0.2.0")
)
@_ASYNC_REGISTRY.register_write(
    ZarrAsyncGroup, pd.DataFrame, IOSpec("dataframe", "0.2.0")
)
async def write_dataframe(
    f: ZarrAsyncGroup,
    key: str,
    df: views.DataFrameView | pd.DataFrame,
    *,
    _writer: AsyncWriter,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    # Check arguments
    for reserved in ("_index",):
        if reserved in df.columns:
            msg = f"{reserved!r} is a reserved name for dataframe columns."
            raise ValueError(msg)
    group = await f.require_group(key)
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
        _writer.write_elem(
            group, index_name, df.index._values, dataset_kwargs=dataset_kwargs
        )
    ]
    # TODO: this should write the "true" representation of the series (i.e. the underlying array or ndarray depending)
    awaitables += [
        _writer.write_elem(
            group, colname, series._values, dataset_kwargs=dataset_kwargs
        )
        for colname, series in df.items()
    ]
    await asyncio.gather(*awaitables)


@_ASYNC_REGISTRY.register_read(ZarrAsyncGroup, IOSpec("dataframe", "0.2.0"))
async def read_dataframe(elem: ZarrAsyncGroup, *, _reader: AsyncReader) -> pd.DataFrame:
    columns = list(_read_attr(elem.attrs, "column-order"))
    idx_key = _read_attr(elem.attrs, "_index")
    column_and_index_elems = await asyncio.gather(
        *(elem.get(k) for k in columns), elem.get(idx_key)
    )
    index_elem = column_and_index_elems[-1]
    column_elems = column_and_index_elems[:-1]
    df = pd.DataFrame(
        dict(
            zip(
                columns,
                await asyncio.gather(*(read_series(e, _reader) for e in column_elems)),
            )
        ),
        index=await _reader.read_elem(index_elem),
        columns=columns if len(columns) else None,
    )
    if idx_key != "_index":
        df.index.name = idx_key
    return df


# Backwards compat dataframe reading


@_ASYNC_REGISTRY.register_read(ZarrAsyncGroup, IOSpec("dataframe", "0.1.0"))
async def read_dataframe_0_1_0(
    elem: ZarrAsyncGroup, *, _reader: AsyncReader
) -> pd.DataFrame:
    columns = _read_attr(elem.attrs, "column-order")
    idx_key = _read_attr(elem.attrs, "_index")
    column_and_index_elems = await asyncio.gather(
        *(elem.get(k) for k in columns), elem.get(idx_key)
    )
    index_elem = column_and_index_elems[-1]
    column_elems = column_and_index_elems[:-1]
    df = pd.DataFrame(
        dict(
            zip(
                columns,
                await asyncio.gather(*(read_series(e, _reader) for e in column_elems)),
            )
        ),
        index=await read_series(index_elem, _reader),
        columns=columns if len(columns) else None,
    )
    if idx_key != "_index":
        df.index.name = idx_key
    return df


###############
# Categorical #
###############


@_ASYNC_REGISTRY.register_write(
    ZarrAsyncGroup, pd.Categorical, IOSpec("categorical", "0.2.0")
)
async def write_categorical(
    f: ZarrAsyncGroup,
    k: str,
    v: pd.Categorical,
    *,
    _writer: AsyncWriter,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    g = await f.require_group(k)
    g.attrs["ordered"] = bool(v.ordered)

    await asyncio.gather(
        _writer.write_elem(g, "codes", v.codes, dataset_kwargs=dataset_kwargs),
        _writer.write_elem(
            g, "categories", v.categories._values, dataset_kwargs=dataset_kwargs
        ),
    )


@_ASYNC_REGISTRY.register_read(ZarrAsyncGroup, IOSpec("categorical", "0.2.0"))
async def read_categorical(
    elem: ZarrAsyncGroup, *, _reader: AsyncReader
) -> pd.Categorical:
    elems = await asyncio.gather(*(elem.get(k) for k in ["codes", "categories"]))
    codes, categories = await asyncio.gather(*(_reader.read_elem(e) for e in elems))
    ordered = bool(_read_attr(elem.attrs, "ordered"))
    return pd.Categorical.from_codes(
        codes=codes, categories=categories, ordered=ordered
    )


####################
# Pandas nullables #
####################


@_ASYNC_REGISTRY.register_write(
    ZarrAsyncGroup, pd.arrays.IntegerArray, IOSpec("nullable-integer", "0.1.0")
)
@_ASYNC_REGISTRY.register_write(
    ZarrAsyncGroup, pd.arrays.BooleanArray, IOSpec("nullable-boolean", "0.1.0")
)
@_ASYNC_REGISTRY.register_write(
    ZarrAsyncGroup, pd.arrays.StringArray, IOSpec("nullable-string-array", "0.1.0")
)
async def write_nullable(
    f: ZarrAsyncGroup,
    k: str,
    v: pd.arrays.IntegerArray | pd.arrays.BooleanArray | pd.arrays.StringArray,
    *,
    _writer: AsyncWriter,
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
    g = await f.require_group(k)
    values = (
        v.to_numpy(na_value="")
        if isinstance(v, pd.arrays.StringArray)
        else v.to_numpy(na_value=0, dtype=v.dtype.numpy_dtype)
    )
    await asyncio.gather(
        _writer.write_elem(g, "values", values, dataset_kwargs=dataset_kwargs),
        _writer.write_elem(g, "mask", v.isna(), dataset_kwargs=dataset_kwargs),
    )


async def _read_nullable(
    elem: ZarrAsyncGroup,
    *,
    _reader: AsyncReader,
    # BaseMaskedArray
    array_type: Callable[
        [NDArray[np.number], NDArray[np.bool_]], pd.api.extensions.ExtensionArray
    ],
) -> pd.api.extensions.ExtensionArray:
    elems = await asyncio.gather(*(elem.get(k) for k in ["values", "mask"]))
    values, mask = await asyncio.gather(*(_reader.read_elem(e) for e in elems))
    return array_type(values, mask=mask)


def _string_array(
    values: np.ndarray, mask: np.ndarray
) -> pd.api.extensions.ExtensionArray:
    """Construct a string array from values and mask."""
    arr = pd.array(values, dtype=pd.StringDtype())
    arr[mask] = pd.NA
    return arr


_ASYNC_REGISTRY.register_read(ZarrAsyncGroup, IOSpec("nullable-integer", "0.1.0"))(
    partial(_read_nullable, array_type=pd.arrays.IntegerArray)
)
_ASYNC_REGISTRY.register_read(ZarrAsyncGroup, IOSpec("nullable-boolean", "0.1.0"))(
    partial(_read_nullable, array_type=pd.arrays.BooleanArray)
)
_ASYNC_REGISTRY.register_read(ZarrAsyncGroup, IOSpec("nullable-string-array", "0.1.0"))(
    partial(_read_nullable, array_type=_string_array)
)

###########
# Scalars #
###########


@_ASYNC_REGISTRY.register_read(ZarrAsyncArray, IOSpec("numeric-scalar", "0.2.0"))
async def read_scalar(elem: ArrayStorageType, *, _reader: AsyncReader) -> np.number:
    # TODO: `item` ensures the return is in fact a scalar (needed after zarr v3 which now returns a 1 elem array)
    # https://github.com/zarr-developers/zarr-python/issues/2713
    return (await elem.getitem(())).item()


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
    _writer: AsyncWriter,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    # these args are ignored in v2: https://zarr.readthedocs.io/en/v2.18.4/api/hierarchy.html#zarr.hierarchy.Group.create_dataset
    # and error out in v3
    dataset_kwargs = _remove_scalar_compression_args(dataset_kwargs)
    from numcodecs import VLenUTF8

    match ad.settings.zarr_write_format, value:
        case 2, str():
            filters, dtype = [VLenUTF8()], object
        case 3, str():
            filters, dtype = None, str
        case _, _:
            filters, dtype = None, np.array(value).dtype
    a = await f.create_array(
        key,
        shape=(),
        dtype=dtype,
        filters=filters,
        **dataset_kwargs,
    )
    await a.setitem(Ellipsis, np.array(value))


_ASYNC_REGISTRY.register_write(ZarrAsyncGroup, str, IOSpec("string", "0.2.0"))(
    write_scalar_zarr
)
_ASYNC_REGISTRY.register_write(ZarrAsyncGroup, np.str_, IOSpec("string", "0.2.0"))(
    write_scalar_zarr
)


@_ASYNC_REGISTRY.register_read(ZarrAsyncArray, IOSpec("string", "0.2.0"))
async def read_zarr_string(elem: ZarrAsyncArray, *, _reader: AsyncReader) -> str:
    return str(await elem.getitem(()))


_ASYNC_REGISTRY.register_read(ZarrAsyncArray, IOSpec("bytes", "0.2.0"))(read_scalar)
