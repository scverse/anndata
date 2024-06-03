from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse

import anndata as ad
from anndata.compat import H5Array, H5Group, ZarrArray, ZarrGroup

from .registry import _LAZY_REGISTRY, IOSpec


def make_index(is_csc, stride, shape, block_id):
    index = (
        slice(
            block_id[is_csc] * stride,
            min((block_id[is_csc] * stride) + stride, shape[0]),
        ),
    )
    if is_csc:
        return (slice(None, None, None),) + index
    return index


@contextmanager
def maybe_open_h5(filename_or_elem: str | ZarrGroup, elem_name: str):
    if isinstance(filename_or_elem, str):
        file = h5py.File(filename_or_elem, "r")
        try:
            yield file[elem_name]
        finally:
            file.close()
    else:
        try:
            yield filename_or_elem
        finally:
            pass


@_LAZY_REGISTRY.register_read(H5Group, IOSpec("csc_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(H5Group, IOSpec("csr_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("csc_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("csr_matrix", "0.1.0"))
def read_sparse_as_dask(elem, _reader, stride: int = 100):
    import dask.array as da

    filename_or_elem = elem.file.filename if isinstance(elem, H5Group) else elem
    elem_name = elem.name if isinstance(elem, H5Group) else Path(elem.path).name
    shape = elem.attrs["shape"]
    dtype = elem["data"].dtype
    is_csc = elem.attrs["encoding-type"] == "csc_matrix"

    def make_dask_chunk(block_id=None):
        # We need to open the file in each task since `dask` cannot share h5py objects when using `dask.distributed`
        # https://github.com/scverse/anndata/issues/1105
        with maybe_open_h5(filename_or_elem, elem_name) as f:
            mtx = ad.experimental.sparse_dataset(f)
            index = make_index(is_csc, stride, shape, block_id)
            chunk = mtx[index]
        return chunk

    chunks = [None, None]
    major_index = int(is_csc)
    minor_index = (is_csc + 1) % 2
    chunks[minor_index] = (shape[minor_index],)
    chunks[major_index] = (stride,) * (shape[major_index] // stride) + (
        shape[major_index] % stride,
    )
    memory_format = [sparse.csr_matrix, sparse.csc_matrix][major_index]
    da_mtx = da.map_blocks(
        make_dask_chunk,
        dtype=dtype,
        chunks=chunks,
        meta=memory_format((0, 0), dtype=np.float32),
    )
    return da_mtx


@_LAZY_REGISTRY.register_read(H5Array, IOSpec("array", "0.2.0"))
def read_h5_array(elem, _reader, chunk_size: int = 1000):
    import dask.array as da

    if not hasattr(elem, "chunks") or elem.chunks is None:
        return da.from_array(elem, chunks=(chunk_size,) * len(elem.shape))
    return da.from_array(elem)


@_LAZY_REGISTRY.register_read(ZarrArray, IOSpec("array", "0.2.0"))
def read_zarr_array(elem, _reader):
    import dask.array as da

    return da.from_zarr(elem)
