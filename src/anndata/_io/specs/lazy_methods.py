from __future__ import annotations

import dask.array as da
import h5py
import numpy as np
from scipy import sparse

import anndata as ad
from anndata.compat import H5Array, H5Group, ZarrArray, ZarrGroup

from .registry import _LAZY_REGISTRY, IOSpec

# TODO: settings
stride = 100
h5_chunks = 1000


def make_dask_array(is_csc, shape, make_dask_chunk, dtype):
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


@_LAZY_REGISTRY.register_read(H5Group, IOSpec("csc_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(H5Group, IOSpec("csr_matrix", "0.1.0"))
def read_sparse_as_dask_h5(elem, _reader):
    filename = elem.file.filename
    elem_name = elem.name
    shape = elem.attrs["shape"]
    dtype = elem["data"].dtype
    is_csc = elem.attrs["encoding-type"] == "csc_matrix"

    def make_dask_chunk(block_id=None):
        # We need to open the file in each task since `dask` cannot share h5py objects when using `dask.distributed`
        # https://github.com/scverse/anndata/issues/1105
        with h5py.File(filename, "r") as f:
            mtx = ad.experimental.sparse_dataset(f[elem_name])
            index = make_index(is_csc, stride, shape, block_id)
            chunk = mtx[*index]
        return chunk

    return make_dask_array(is_csc, shape, make_dask_chunk, dtype)


@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("csc_matrix", "0.1.0"))
@_LAZY_REGISTRY.register_read(ZarrGroup, IOSpec("csr_matrix", "0.1.0"))
def read_sparse_as_dask_zarr(elem, _reader):
    shape = elem.attrs["shape"]
    dtype = elem["data"].dtype
    is_csc = elem.attrs["encoding-type"] == "csc_matrix"

    def make_dask_chunk(block_id=None):
        mtx = ad.experimental.sparse_dataset(elem)
        index = make_index(is_csc, stride, shape, block_id)
        return mtx[*index]

    return make_dask_array(is_csc, shape, make_dask_chunk, dtype)


@_LAZY_REGISTRY.register_read(H5Array, IOSpec("array", "0.2.0"))
def read_h5_array(elem, _reader):
    if not hasattr(elem, "chunks") or elem.chunks is None:
        return da.from_array(elem, chunks=(h5_chunks,) * len(elem.shape))
    return da.from_array(elem)


@_LAZY_REGISTRY.register_read(ZarrArray, IOSpec("array", "0.2.0"))
def read_zarr_array(elem, _reader):
    return da.from_zarr(elem)
