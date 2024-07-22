from __future__ import annotations

from contextlib import contextmanager
from functools import wraps
from pathlib import Path
from typing import TYPE_CHECKING, Literal, ParamSpec, TypeVar, Union

import h5py
import numpy as np
from scipy import sparse

import anndata as ad
from anndata._core.file_backing import filename, get_elem_name
from anndata.compat import H5Array, H5Group, ZarrArray, ZarrGroup

from .registry import _LAZY_REGISTRY, IOSpec

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Concatenate

    from anndata.compat import DaskArray

    from .registry import DaskReader


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


P = ParamSpec("P")
R = TypeVar("R")
BlockInfo = dict[
    Literal[None],
    dict[str, Union[tuple[int, ...], list[tuple[int, ...]]]],
]


def require_block_info(
    f: Callable[Concatenate[BlockInfo, P], R],
) -> Callable[Concatenate[BlockInfo | None, P], R]:
    @wraps(f)
    def wrapper(
        block_info: BlockInfo | None = None, *args: P.args, **kwargs: P.kwargs
    ) -> R:
        if block_info is None:
            msg = "Block info is required"
            raise ValueError(msg)
        return f(block_info, *args, **kwargs)

    return wrapper


@_LAZY_REGISTRY.register_read(H5Group, IOSpec("csc_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(H5Group, IOSpec("csr_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("csc_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("csr_matrix", "0.1.0"))
def read_sparse_as_dask(
    elem: H5Group | ZarrGroup,
    *,
    _reader: DaskReader,
    chunks: tuple[int, ...] | None = None,
) -> DaskArray:
    import dask.array as da

    path_or_group = Path(filename(elem)) if isinstance(elem, H5Group) else elem
    elem_name = get_elem_name(elem)
    shape: tuple[int, int] = tuple(elem.attrs["shape"])
    dtype = elem["data"].dtype
    is_csc: bool = elem.attrs["encoding-type"] == "csc_matrix"

    stride: int = _DEFAULT_STRIDE
    if chunks is not None:
        if len(chunks) != 2:
            raise ValueError("`chunks` must be a tuple of two integers")
        if chunks[int(not is_csc)] != shape[int(not is_csc)]:
            raise ValueError(
                "Only the major axis can be chunked. "
                f"Try setting chunks to {((-1, _DEFAULT_STRIDE) if is_csc else (_DEFAULT_STRIDE, -1))}"
            )
        stride = chunks[int(is_csc)]

    @require_block_info
    def make_dask_chunk(block_info: BlockInfo):
        # We need to open the file in each task since `dask` cannot share h5py objects when using `dask.distributed`
        # https://github.com/scverse/anndata/issues/1105
        with maybe_open_h5(path_or_group, elem_name) as f:
            mtx = ad.experimental.sparse_dataset(f)
            array_location = block_info[None]["array-location"]
            index = (
                slice(array_location[0][0], array_location[0][1]),
                slice(array_location[1][0], array_location[1][1]),
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
        meta=memory_format((0, 0), dtype=dtype),
    )
    return da_mtx


@_LAZY_REGISTRY.register_read(H5Array, IOSpec("array", "0.2.0"))
def read_h5_array(
    elem: H5Array, *, _reader: DaskReader, chunks: tuple[int, ...] | None = None
) -> DaskArray:
    import dask.array as da

    path = Path(elem.file.filename)
    elem_name = elem.name
    shape = tuple(elem.shape)
    dtype = elem.dtype
    chunks: tuple[int, ...] = (
        chunks if chunks is not None else (_DEFAULT_STRIDE,) * len(shape)
    )

    @require_block_info
    def make_dask_chunk(block_info: BlockInfo):
        with maybe_open_h5(path, elem_name) as f:
            idx = ()
            for i in range(len(shape)):
                array_location = block_info[None]["array-location"][i]
                idx += (slice(array_location[0], array_location[1]),)
            return f[idx]

    chunk_layout = tuple(
        compute_chunk_layout_for_axis_shape(chunks[i], shape[i])
        for i in range(len(shape))
    )

    return da.map_blocks(make_dask_chunk, dtype=dtype, chunks=chunk_layout)


@_LAZY_REGISTRY.register_read(ZarrArray, IOSpec("array", "0.2.0"))
def read_zarr_array(
    elem: ZarrArray, *, _reader: DaskReader, chunks: tuple[int, ...] | None = None
) -> DaskArray:
    chunks: tuple[int, ...] = chunks if chunks is not None else elem.chunks
    import dask.array as da

    return da.from_zarr(elem, chunks=chunks)
