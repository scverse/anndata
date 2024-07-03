from __future__ import annotations

from contextlib import contextmanager
from functools import singledispatch
from pathlib import Path, PurePosixPath
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Literal, overload

import h5py
import numpy as np
from scipy import sparse

import anndata as ad
from anndata.compat import H5Array, H5Group, ZarrArray, ZarrGroup

from .registry import _LAZY_REGISTRY, IOSpec

if TYPE_CHECKING:
    from collections.abc import Mapping


@overload
def make_block_indexer(
    *,
    is_csc: Literal[True],
    stride: int,
    shape: tuple[int, int],
    block_id: tuple[int, int],
) -> tuple[slice, slice]: ...
@overload
def make_block_indexer(
    *,
    is_csc: Literal[False],
    stride: int,
    shape: tuple[int, int],
    block_id: tuple[int, int],
) -> tuple[slice]: ...


def make_block_indexer(
    *, is_csc: bool, stride: int, shape: tuple[int, int], block_id: tuple[int, int]
) -> tuple[slice, slice] | tuple[slice]:
    index1d = slice(
        block_id[is_csc] * stride,
        min((block_id[is_csc] * stride) + stride, shape[0]),
    )
    if is_csc:
        return (slice(None), index1d)
    return (index1d,)


@contextmanager
def maybe_open_h5(path_or_group: Path | ZarrGroup, elem_name: str):
    if not isinstance(path_or_group, Path):
        yield path_or_group
        return
    file = h5py.File(path_or_group, "r")
    try:
        yield file[elem_name]
    finally:
        file.close()


_DEFAULT_STRIDE = 1000


def compute_chunk_layout_for_axis_shape(
    chunk_axis_shape: int, full_axis_shape: int
) -> tuple[int, ...]:
    n_strides, rest = np.divmod(full_axis_shape, chunk_axis_shape)
    chunk = (chunk_axis_shape,) * n_strides
    if rest > 0:
        chunk += (rest,)
    return chunk


@singledispatch
def get_elem_name(x):
    raise NotImplementedError(f"Not implemented for {type(x)}")


@get_elem_name.register(H5Group)
def _(x):
    return x.name


@get_elem_name.register(ZarrGroup)
def _(x):
    return PurePosixPath(x.path).name


@_LAZY_REGISTRY.register_read(H5Group, IOSpec("csc_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(H5Group, IOSpec("csr_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("csc_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("csr_matrix", "0.1.0"))
def read_sparse_as_dask(
    elem: H5Group | ZarrGroup,
    _reader,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    import dask.array as da

    path_or_group = Path(elem.file.filename) if isinstance(elem, H5Group) else elem
    elem_name = get_elem_name(elem)
    shape: tuple[int, int] = elem.attrs["shape"]
    dtype = elem["data"].dtype
    is_csc: bool = elem.attrs["encoding-type"] == "csc_matrix"

    chunks = dataset_kwargs.get("chunks", None)
    stride: int = _DEFAULT_STRIDE
    if chunks is not None:
        if len(chunks) != 2:
            raise ValueError("`chunks` must be a tuple of two integers")
        if chunks[int(not is_csc)] != shape[int(not is_csc)]:
            raise ValueError("Only the major axis can be chunked")
        stride = chunks[int(is_csc)]

    def make_dask_chunk(block_id: tuple[int, int]):
        # We need to open the file in each task since `dask` cannot share h5py objects when using `dask.distributed`
        # https://github.com/scverse/anndata/issues/1105
        with maybe_open_h5(path_or_group, elem_name) as f:
            mtx = ad.experimental.sparse_dataset(f)
            index = make_block_indexer(
                is_csc=is_csc, stride=stride, shape=shape, block_id=block_id
            )
            chunk = mtx[index]
        return chunk

    shape_minor, shape_major = shape if is_csc else shape[::-1]
    chunks_major = compute_chunk_layout_for_axis_shape(stride, shape_major)
    chunks_minor = (shape_minor,)
    chunk_layout = (
        (chunks_minor, chunks_major) if is_csc else (chunks_major, chunks_minor)
    )
    memory_format = sparse.csc_matrix if is_csc else sparse.csr_matrix
    da_mtx = da.map_blocks(
        make_dask_chunk,
        dtype=dtype,
        chunks=chunk_layout,
        meta=memory_format((0, 0), dtype=np.float32),
    )
    return da_mtx


@_LAZY_REGISTRY.register_read(H5Array, IOSpec("array", "0.2.0"))
def read_h5_array(
    elem,
    _reader,
    chunks: tuple[int] | None = None,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    import dask.array as da

    path = Path(elem.file.filename)
    elem_name = elem.name
    shape = elem.shape
    dtype = elem.dtype
    chunks: tuple[int, ...] = dataset_kwargs.get(
        "chunks", (_DEFAULT_STRIDE,) * len(shape)
    )

    def make_dask_chunk(block_id: tuple[int, int]):
        with maybe_open_h5(path, elem_name) as f:
            idx = ()
            for i in range(len(shape)):
                start = block_id[i] * chunks[i]
                stop = min(((block_id[i] * chunks[i]) + chunks[i]), shape[i])
                idx += (slice(start, stop),)
            return f[idx]

    chunk_layout = tuple(
        compute_chunk_layout_for_axis_shape(chunks[i], shape[i])
        for i in range(len(shape))
    )

    return da.map_blocks(
        make_dask_chunk,
        dtype=dtype,
        chunks=chunk_layout,
    )


@_LAZY_REGISTRY.register_read(ZarrArray, IOSpec("array", "0.2.0"))
def read_zarr_array(
    elem, _reader, dataset_kwargs: Mapping[str, Any] = MappingProxyType({})
):
    import dask.array as da

    return da.from_zarr(elem)
