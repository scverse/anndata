"""
This file contains tests for conversion made during io.
"""
import h5py
import numpy as np
import pytest
from scipy import sparse

import anndata as ad
from anndata.tests.helpers import gen_adata, assert_equal


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


def test_sparse_to_dense_disk(tmp_path, mtx_format):
    mem_pth = tmp_path / "orig.h5ad"
    dense_from_mem_pth = tmp_path / "dense_mem.h5ad"
    dense_from_disk_pth = tmp_path / "dense_disk.h5ad"
    mem = gen_adata((50, 50), mtx_format)

    mem.write_h5ad(mem_pth)
    disk = ad.read_h5ad(mem_pth, backed="r")

    mem.write_h5ad(dense_from_mem_pth, force_dense=True)
    disk.write_h5ad(dense_from_disk_pth, force_dense=True)

    with h5py.File(dense_from_mem_pth, "r") as f:
        assert isinstance(f["X"], h5py.Dataset)
    with h5py.File(dense_from_disk_pth, "r") as f:
        assert isinstance(f["X"], h5py.Dataset)

    for backed in [None, "r"]:
        from_mem = ad.read_h5ad(dense_from_mem_pth, backed=backed)
        from_disk = ad.read_h5ad(dense_from_disk_pth, backed=backed)
        assert_equal(mem, from_mem)
        assert_equal(mem, from_disk)
        assert_equal(disk, from_mem)
        assert_equal(disk, from_disk)


def test_dense_to_sparse_memory(tmp_path, spmtx_format):
    dense_path = tmp_path / "dense.h5ad"
    orig = gen_adata((50, 50), np.array)
    orig.write_h5ad(dense_path)

    curr = ad.read_h5ad(dense_path, dense_X_as_sparse=spmtx_format)

    assert_equal(orig, curr)
