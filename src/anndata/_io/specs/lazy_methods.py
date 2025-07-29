from __future__ import annotations

from contextlib import contextmanager
from functools import partial, singledispatch
from pathlib import Path
from typing import TYPE_CHECKING, overload

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

import anndata as ad
from anndata._core.file_backing import filename, get_elem_name
from anndata._core.xarray import Dataset2D, requires_xarray
from anndata.abc import CSCDataset, CSRDataset
from anndata.compat import (
    NULLABLE_NUMPY_STRING_TYPE,
    DaskArray,
    H5Array,
    H5Group,
    XDataArray,
    XDataset,
    ZarrArray,
    ZarrGroup,
)

from .registry import _LAZY_REGISTRY, IOSpec

if TYPE_CHECKING:
    from collections.abc import Generator, Mapping, Sequence
    from typing import Literal, ParamSpec, TypeVar

    from anndata.experimental.backed._lazy_arrays import CategoricalArray, MaskedArray

    from ...compat import CSArray, CSMatrix, H5File
    from .registry import LazyDataStructures, LazyReader

    BlockInfo = Mapping[
        Literal[None],
        dict[str, Sequence[tuple[int, int]]],
    ]

    P = ParamSpec("P")
    R = TypeVar("R")
    D = TypeVar("D")


@overload
@contextmanager
def maybe_open_h5(
    path_or_other: Path, elem_name: str
) -> Generator[H5File, None, None]: ...
@overload
@contextmanager
def maybe_open_h5(path_or_other: D, elem_name: str) -> Generator[D, None, None]: ...
@contextmanager
def maybe_open_h5(
    path_or_other: H5File | D, elem_name: str
) -> Generator[H5File | D, None, None]:
    if not isinstance(path_or_other, Path):
        yield path_or_other
        return
    file = h5py.File(path_or_other, "r")
    try:
        yield file[elem_name]
    finally:
        file.close()


_DEFAULT_STRIDE = 1000


def compute_chunk_layout_for_axis_size(
    chunk_axis_size: int, full_axis_size: int
) -> tuple[int, ...]:
    n_strides, rest = np.divmod(full_axis_size, chunk_axis_size)
    chunk = (chunk_axis_size,) * n_strides
    if rest > 0:
        chunk += (rest,)
    return chunk


def make_dask_chunk(
    path_or_sparse_dataset: Path | D,
    elem_name: str,
    block_info: BlockInfo | None = None,
) -> CSMatrix | CSArray:
    if block_info is None:
        msg = "Block info is required"
        raise ValueError(msg)
    # We need to open the file in each task since `dask` cannot share h5py objects when using `dask.distributed`
    # https://github.com/scverse/anndata/issues/1105
    with maybe_open_h5(path_or_sparse_dataset, elem_name) as f:
        # See https://github.com/scverse/anndata/pull/2005 for why
        # should_cache_indptr is False.
        # The prupose of caching the indptr was when the dataset is reused
        # which is in general the case but is not here.  Hence
        # caching it on every access to the dataset here is quite costly.
        mtx = (
            ad.io.sparse_dataset(f, should_cache_indptr=False)
            if isinstance(f, H5Group)
            else f
        )
        idx = tuple(
            slice(start, stop) for start, stop in block_info[None]["array-location"]
        )
        chunk = mtx[idx]
    return chunk


@singledispatch
def get_chunksize(obj) -> tuple[int, ...]:
    if hasattr(obj, "chunks"):
        return obj.chunks
    msg = "object of type {type(obj)} has no recognized chunks"
    raise ValueError(msg)


@_LAZY_REGISTRY.register_read(H5Group, IOSpec("csc_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(H5Group, IOSpec("csr_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("csc_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("csr_matrix", "0.1.0"))
def read_sparse_as_dask(
    elem: H5Group | ZarrGroup,
    *,
    _reader: LazyReader,
    chunks: tuple[int, ...] | None = None,  # only tuple[int, int] is supported here
) -> DaskArray:
    import dask.array as da

    path_or_sparse_dataset = (
        Path(filename(elem))
        if isinstance(elem, H5Group)
        else ad.io.sparse_dataset(elem, should_cache_indptr=False)
    )
    elem_name = get_elem_name(elem)
    shape: tuple[int, int] = tuple(elem.attrs["shape"])
    if isinstance(path_or_sparse_dataset, CSRDataset | CSCDataset):
        dtype = path_or_sparse_dataset.dtype
    else:
        dtype = elem["data"].dtype
    is_csc: bool = elem.attrs["encoding-type"] == "csc_matrix"

    stride: int = _DEFAULT_STRIDE
    major_dim, minor_dim = (1, 0) if is_csc else (0, 1)
    if chunks is not None:
        if len(chunks) != 2:
            msg = "`chunks` must be a tuple of two integers"
            raise ValueError(msg)
        if chunks[minor_dim] not in {shape[minor_dim], -1, None}:
            msg = (
                "Only the major axis can be chunked. "
                f"Try setting chunks to {((-1, _DEFAULT_STRIDE) if is_csc else (_DEFAULT_STRIDE, -1))}"
            )
            raise ValueError(msg)
        stride = (
            chunks[major_dim]
            if chunks[major_dim] not in {None, -1}
            else shape[major_dim]
        )

    shape_minor, shape_major = shape if is_csc else shape[::-1]
    chunks_major = compute_chunk_layout_for_axis_size(stride, shape_major)
    chunks_minor = (shape_minor,)
    chunk_layout = (
        (chunks_minor, chunks_major) if is_csc else (chunks_major, chunks_minor)
    )
    memory_format = sparse.csc_matrix if is_csc else sparse.csr_matrix
    make_chunk = partial(make_dask_chunk, path_or_sparse_dataset, elem_name)
    da_mtx = da.map_blocks(
        make_chunk,
        dtype=dtype,
        chunks=chunk_layout,
        meta=memory_format((0, 0), dtype=dtype),
    )
    return da_mtx


def resolve_chunks(
    elem: H5Array | ZarrArray,
    chunks_arg: tuple[int, ...] | None,
    shape: tuple[int, ...],
) -> tuple[int, ...]:
    shape = tuple(elem.shape)
    if chunks_arg is not None:
        # None and -1 on a given axis indicate that one should use the shape
        # in `dask`'s semantics.
        return tuple(
            c if c not in {None, -1} else s
            for c, s in zip(chunks_arg, shape, strict=True)
        )
    elif elem.chunks is None:  # h5 unchunked
        return tuple(min(_DEFAULT_STRIDE, s) for s in shape)
    return elem.chunks


@_LAZY_REGISTRY.register_read(H5Array, IOSpec("string-array", "0.2.0"))
def read_h5_string_array(
    elem: H5Array,
    *,
    _reader: LazyReader,
    chunks: tuple[int] | None = None,
) -> DaskArray:
    import dask.array as da

    from anndata._io.h5ad import read_dataset

    chunks = resolve_chunks(elem, chunks, tuple(elem.shape))
    return da.from_array(read_dataset(elem), chunks=chunks)


@_LAZY_REGISTRY.register_read(H5Array, IOSpec("array", "0.2.0"))
def read_h5_array(
    elem: H5Array, *, _reader: LazyReader, chunks: tuple[int, ...] | None = None
) -> DaskArray:
    import dask.array as da

    path = Path(elem.file.filename)
    elem_name: str = elem.name
    shape = tuple(elem.shape)
    dtype = elem.dtype
    chunks = resolve_chunks(elem, chunks, shape)

    chunk_layout = tuple(
        compute_chunk_layout_for_axis_size(chunks[i], shape[i])
        for i in range(len(shape))
    )

    make_chunk = partial(make_dask_chunk, path, elem_name)
    return da.map_blocks(
        make_chunk, dtype=dtype, chunks=chunk_layout, meta=np.array([])
    )


@_LAZY_REGISTRY.register_read(ZarrArray, IOSpec("string-array", "0.2.0"))
@_LAZY_REGISTRY.register_read(ZarrArray, IOSpec("array", "0.2.0"))
def read_zarr_array(
    elem: ZarrArray, *, _reader: LazyReader, chunks: tuple[int, ...] | None = None
) -> DaskArray:
    import dask.array as da

    return da.from_zarr(elem, chunks=chunks)


def _gen_xarray_dict_iterator_from_elems(
    elem_dict: dict[str, LazyDataStructures],
    dim_name: str,
    index: np.NDArray,
) -> Generator[tuple[str, XDataArray], None, None]:
    from anndata.experimental.backed._lazy_arrays import CategoricalArray, MaskedArray

    from ...compat import XDataArray
    from ...compat import xarray as xr

    for k, v in elem_dict.items():
        if isinstance(v, DaskArray) and k != dim_name:
            data_array = XDataArray(v, coords=[index], dims=[dim_name], name=k)
        elif isinstance(v, CategoricalArray | MaskedArray) and k != dim_name:
            variable = xr.Variable(
                data=xr.core.indexing.LazilyIndexedArray(v), dims=[dim_name]
            )
            data_array = XDataArray(
                variable,
                coords=[index],
                dims=[dim_name],
                name=k,
                attrs={
                    "base_path_or_zarr_group": v.base_path_or_zarr_group,
                    "elem_name": v.elem_name,
                    "is_nullable_string": isinstance(v, MaskedArray)
                    and v.dtype == NULLABLE_NUMPY_STRING_TYPE,
                },
            )
        elif k == dim_name:
            data_array = XDataArray(
                index, coords=[index], dims=[dim_name], name=dim_name
            )
        else:
            msg = f"Could not read {k}: {v} from into xarray Dataset2D"
            raise ValueError(msg)
        yield k, data_array


DUMMY_RANGE_INDEX_KEY = "_anndata_dummy_range_index"


@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("dataframe", "0.2.0"))
@_LAZY_REGISTRY.register_read(H5Group, IOSpec("dataframe", "0.2.0"))
@requires_xarray
def read_dataframe(
    elem: H5Group | ZarrGroup,
    *,
    _reader: LazyReader,
    use_range_index: bool = False,
    chunks: tuple[int] | None = None,
) -> Dataset2D:
    elem_dict = {
        k: _reader.read_elem(elem[k], chunks=chunks)
        for k in [*elem.attrs["column-order"], elem.attrs["_index"]]
    }
    # If we use a range index, the coord axis needs to have the special dim name
    # which is used below as well.
    if not use_range_index:
        dim_name = elem.attrs["_index"]
        # no sense in reading this in multiple times
        index = elem_dict[dim_name].compute()
    else:
        dim_name = DUMMY_RANGE_INDEX_KEY
        index = pd.RangeIndex(len(elem_dict[elem.attrs["_index"]])).astype("str")
    elem_xarray_dict = dict(
        _gen_xarray_dict_iterator_from_elems(elem_dict, dim_name, index)
    )
    if use_range_index:
        elem_xarray_dict[DUMMY_RANGE_INDEX_KEY] = XDataArray(
            index,
            coords=[index],
            dims=[DUMMY_RANGE_INDEX_KEY],
            name=DUMMY_RANGE_INDEX_KEY,
        )
    ds = Dataset2D(XDataset(elem_xarray_dict))
    ds.is_backed = True
    # We ensure the indexing_key attr always points to the true index
    # so that the roundtrip works even for the `use_range_index` `True` case
    ds.true_index_dim = elem.attrs["_index"]
    return ds


@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("categorical", "0.2.0"))
@_LAZY_REGISTRY.register_read(H5Group, IOSpec("categorical", "0.2.0"))
@requires_xarray
def read_categorical(
    elem: H5Group | ZarrGroup,
    *,
    _reader: LazyReader,
) -> CategoricalArray:
    from anndata.experimental.backed._lazy_arrays import CategoricalArray

    base_path_or_zarr_group = (
        Path(filename(elem)) if isinstance(elem, H5Group) else elem
    )
    elem_name = get_elem_name(elem)
    return CategoricalArray(
        codes=elem["codes"],
        categories=elem["categories"],
        ordered=elem.attrs["ordered"],
        base_path_or_zarr_group=base_path_or_zarr_group,
        elem_name=elem_name,
    )


@requires_xarray
def read_nullable(
    elem: H5Group | ZarrGroup,
    *,
    encoding_type: Literal[
        "nullable-integer", "nullable-boolean", "nullable-string-array"
    ],
    _reader: LazyReader,
) -> MaskedArray:
    from anndata.experimental.backed._lazy_arrays import MaskedArray

    base_path_or_zarr_group = (
        Path(filename(elem)) if isinstance(elem, H5Group) else elem
    )
    elem_name = get_elem_name(elem)
    return MaskedArray(
        values=elem["values"],
        mask=elem.get("mask", None),
        dtype_str=encoding_type,
        base_path_or_zarr_group=base_path_or_zarr_group,
        elem_name=elem_name,
    )


for dtype in ["integer", "boolean", "string-array"]:
    for group_type in [ZarrGroup, H5Group]:
        _LAZY_REGISTRY.register_read(group_type, IOSpec(f"nullable-{dtype}", "0.1.0"))(
            partial(read_nullable, encoding_type=f"nullable-{dtype}")
        )
