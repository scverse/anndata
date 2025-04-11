"""\
This file contains tests for conversion made during io.
"""

from __future__ import annotations

import h5py
import numpy as np
import pytest
from scipy import sparse

import anndata as ad
from anndata.compat import CSMatrix
from anndata.tests.helpers import GEN_ADATA_NO_XARRAY_ARGS, assert_equal, gen_adata


@pytest.fixture(
    params=[sparse.csr_matrix, sparse.csc_matrix, np.array],
    ids=["scipy-csr", "scipy-csc", "np-array"],
)
def mtx_format(request):
    return request.param


@pytest.fixture(
    params=[sparse.csr_matrix, sparse.csc_matrix],
    ids=["scipy-csr", "scipy-csc"],
)
def spmtx_format(request):
    return request.param


@pytest.fixture(params=[("raw/X",), ("X",), ("X", "raw/X")])
def to_convert(request):
    return request.param


def test_sparse_to_dense_disk(tmp_path, mtx_format, to_convert):
    mem_pth = tmp_path / "orig.h5ad"
    dense_from_mem_pth = tmp_path / "dense_mem.h5ad"
    dense_from_disk_pth = tmp_path / "dense_disk.h5ad"
    mem = gen_adata((50, 50), mtx_format, **GEN_ADATA_NO_XARRAY_ARGS)
    mem.raw = mem.copy()

    mem.write_h5ad(mem_pth)
    disk = ad.read_h5ad(mem_pth, backed="r")

    mem.write_h5ad(dense_from_mem_pth, as_dense=to_convert)
    disk.write_h5ad(dense_from_disk_pth, as_dense=to_convert)

    with h5py.File(dense_from_mem_pth, "r") as f:
        for k in to_convert:
            assert isinstance(f[k], h5py.Dataset)
    with h5py.File(dense_from_disk_pth, "r") as f:
        for k in to_convert:
            assert isinstance(f[k], h5py.Dataset)

    for backed in [None, "r"]:
        from_mem = ad.read_h5ad(dense_from_mem_pth, backed=backed)
        from_disk = ad.read_h5ad(dense_from_disk_pth, backed=backed)
        assert_equal(mem, from_mem)
        assert_equal(mem, from_disk)
        assert_equal(disk, from_mem)
        assert_equal(disk, from_disk)


def test_sparse_to_dense_inplace(tmp_path, spmtx_format):
    pth = tmp_path / "adata.h5ad"
    orig = gen_adata((50, 50), spmtx_format, **GEN_ADATA_NO_XARRAY_ARGS)
    orig.raw = orig.copy()
    orig.write(pth)
    backed = ad.read_h5ad(pth, backed="r+")
    backed.write(as_dense=("X", "raw/X"))
    new = ad.read_h5ad(pth)

    assert_equal(orig, new)
    assert_equal(backed, new)

    assert isinstance(new.X, np.ndarray)
    assert isinstance(new.raw.X, np.ndarray)
    assert isinstance(orig.X, spmtx_format)
    assert isinstance(orig.raw.X, spmtx_format)
    assert isinstance(backed.X, h5py.Dataset)
    assert isinstance(backed.raw.X, h5py.Dataset)


def test_sparse_to_dense_errors(tmp_path):
    adata = ad.AnnData(X=sparse.random(50, 50, format="csr"))
    adata.layers["like_X"] = adata.X.copy()
    with pytest.raises(ValueError, match=r"Cannot specify writing"):
        adata.write_h5ad(tmp_path / "failure.h5ad", as_dense=("raw/X",))
    with pytest.raises(NotImplementedError):
        adata.write_h5ad(tmp_path / "failure.h5ad", as_dense=("raw", "X"))
    with pytest.raises(NotImplementedError):
        adata.write_h5ad(tmp_path / "failure.h5ad", as_dense=("layers/like_X",))


def test_dense_to_sparse_memory(tmp_path, spmtx_format, to_convert):
    dense_path = tmp_path / "dense.h5ad"
    orig = gen_adata((50, 50), np.array, **GEN_ADATA_NO_XARRAY_ARGS)
    orig.raw = orig.copy()
    orig.write_h5ad(dense_path)
    assert not isinstance(orig.X, CSMatrix)
    assert not isinstance(orig.raw.X, CSMatrix)

    curr = ad.read_h5ad(dense_path, as_sparse=to_convert, as_sparse_fmt=spmtx_format)

    if "X" in to_convert:
        assert isinstance(curr.X, spmtx_format)
    if "raw/X" in to_convert:
        assert isinstance(curr.raw.X, spmtx_format)

    assert_equal(orig, curr)


def test_dense_to_sparse_errors(tmp_path):
    dense_pth = tmp_path / "dense.h5ad"
    adata = ad.AnnData(X=np.ones((50, 50)))
    adata.layers["like_X"] = adata.X.copy()
    adata.write(dense_pth)

    with pytest.raises(NotImplementedError):
        ad.read_h5ad(dense_pth, as_sparse=("X",), as_sparse_fmt=sparse.coo_matrix)
    with pytest.raises(NotImplementedError):
        ad.read_h5ad(dense_pth, as_sparse=("layers/like_X",))
