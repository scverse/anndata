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
    index1d = slice(
        block_id[is_csc] * stride,
        min((block_id[is_csc] * stride) + stride, shape[0]),
    )
    if is_csc:
        return (slice(None, None, None), index1d)
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

    n_strides, rest = np.divmod(shape[major_index], stride)
    chunks_major = (stride,) * n_strides + (rest,)
    chunks_minor = (shape[minor_index],)
    chunks = (chunks_minor, chunks_major) if is_csc else (chunks_major, chunks_minor)
    memory_format = sparse.csc_matrix if is_csc else sparse.csr_matrix
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
