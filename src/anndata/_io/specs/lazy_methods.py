from __future__ import annotations

from contextlib import contextmanager
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, overload

import h5py
import numpy as np
from scipy import sparse

import anndata as ad
from anndata.abc import CSCDataset, CSRDataset

from ..._core.file_backing import filename, get_elem_name
from ...compat import H5Array, H5Group, ZarrArray, ZarrGroup
from .registry import _LAZY_REGISTRY, IOSpec

if TYPE_CHECKING:
    from collections.abc import Generator, Mapping, Sequence
    from typing import Literal, ParamSpec, TypeVar

    from ...compat import DaskArray, H5File, SpArray
    from .registry import DaskReader

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


def compute_chunk_layout_for_axis_shape(
    chunk_axis_shape: int, full_axis_shape: int
) -> tuple[int, ...]:
    n_strides, rest = np.divmod(full_axis_shape, chunk_axis_shape)
    chunk = (chunk_axis_shape,) * n_strides
    if rest > 0:
        chunk += (rest,)
    return chunk


def make_dask_chunk(
    path_or_sparse_dataset: Path | D,
    elem_name: str,
    block_info: BlockInfo | None = None,
) -> sparse.csr_matrix | sparse.csc_matrix | SpArray:
    if block_info is None:
        msg = "Block info is required"
        raise ValueError(msg)
    # We need to open the file in each task since `dask` cannot share h5py objects when using `dask.distributed`
    # https://github.com/scverse/anndata/issues/1105
    with maybe_open_h5(path_or_sparse_dataset, elem_name) as f:
        mtx = ad.io.sparse_dataset(f) if isinstance(f, H5Group) else f
        idx = tuple(
            slice(start, stop) for start, stop in block_info[None]["array-location"]
        )
        chunk = mtx[idx]
    return chunk


@_LAZY_REGISTRY.register_read(H5Group, IOSpec("csc_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(H5Group, IOSpec("csr_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("csc_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("csr_matrix", "0.1.0"))
def read_sparse_as_dask(
    elem: H5Group | ZarrGroup,
    *,
    _reader: DaskReader,
    chunks: tuple[int, ...] | None = None,  # only tuple[int, int] is supported here
) -> DaskArray:
    import dask.array as da

    path_or_sparse_dataset = (
        Path(filename(elem))
        if isinstance(elem, H5Group)
        else ad.io.sparse_dataset(elem)
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
            raise ValueError("`chunks` must be a tuple of two integers")
        if chunks[minor_dim] not in {shape[minor_dim], -1, None}:
            raise ValueError(
                "Only the major axis can be chunked. "
                f"Try setting chunks to {((-1, _DEFAULT_STRIDE) if is_csc else (_DEFAULT_STRIDE, -1))}"
            )
        stride = (
            chunks[major_dim]
            if chunks[major_dim] not in {None, -1}
            else shape[major_dim]
        )

    shape_minor, shape_major = shape if is_csc else shape[::-1]
    chunks_major = compute_chunk_layout_for_axis_shape(stride, shape_major)
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


@_LAZY_REGISTRY.register_read(H5Array, IOSpec("array", "0.2.0"))
def read_h5_array(
    elem: H5Array, *, _reader: DaskReader, chunks: tuple[int, ...] | None = None
) -> DaskArray:
    import dask.array as da

    path = Path(elem.file.filename)
    elem_name: str = elem.name
    shape = tuple(elem.shape)
    dtype = elem.dtype
    chunks: tuple[int, ...] = (
        tuple(c if c not in {None, -1} else shape[i] for i, c in enumerate(chunks))
        if chunks is not None
        else (_DEFAULT_STRIDE,) * len(shape)
    )

    chunk_layout = tuple(
        compute_chunk_layout_for_axis_shape(chunks[i], shape[i])
        for i in range(len(shape))
    )

    make_chunk = partial(make_dask_chunk, path, elem_name)
    return da.map_blocks(make_chunk, dtype=dtype, chunks=chunk_layout)


@_LAZY_REGISTRY.register_read(ZarrArray, IOSpec("array", "0.2.0"))
def read_zarr_array(
    elem: ZarrArray, *, _reader: DaskReader, chunks: tuple[int, ...] | None = None
) -> DaskArray:
    chunks: tuple[int, ...] = chunks if chunks is not None else elem.chunks
    import dask.array as da

    return da.from_zarr(elem, chunks=chunks)
