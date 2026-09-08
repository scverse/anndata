from __future__ import annotations

import numpy as np
import pytest
import zarr
from scipy import sparse
from zarr.storage import MemoryStore

import anndata as ad
from anndata import AnnData, Raw
from anndata._core.sparse_dataset import sparse_dataset
from anndata.compat import CupyCSRMatrix
from anndata.tests.helpers import assert_equal


@pytest.mark.gpu
def test_gpu():
    """
    For testing that the gpu mark works
    """
    import cupy  # This test shouldn't run if cupy isn't installed

    cupy.ones(1)


@pytest.mark.gpu
def test_adata_raw_gpu():
    import cupy as cp
    from cupyx.scipy import sparse as cupy_sparse

    adata = AnnData(
        X=cupy_sparse.random(500, 50, density=0.01, format="csr", dtype=cp.float32)
    )
    adata.raw = adata.copy()
    assert isinstance(adata.raw.X, sparse.csr_matrix)


@pytest.mark.gpu
def test_raw_gpu():
    import cupy as cp
    from cupyx.scipy import sparse as cupy_sparse

    adata = AnnData(
        X=cupy_sparse.random(500, 50, density=0.01, format="csr", dtype=cp.float32)
    )
    araw = Raw(adata)
    assert isinstance(araw.X, sparse.csr_matrix)


@pytest.mark.gpu
@pytest.mark.parametrize(
    "index",
    [
        pytest.param(
            slice(0, 25, 2),
            marks=pytest.mark.xfail(
                reason="see https://github.com/zarr-developers/zarr-python/issues/3640"
            ),
        ),
        pytest.param(
            np.array([0, 10, 20]),
            marks=pytest.mark.xfail(
                reason="see https://github.com/zarr-developers/zarr-python/issues/3640"
            ),
        ),
        slice(None),
        slice(3, 10),
    ],
)
def test_get_with_zarr_gpu(index: slice | np.ndarray):
    adata = AnnData(X=sparse.random(50, 100, format="csr"))
    store = MemoryStore()
    # compressor None because there are no GPU compressors right now
    ad.io.write_zarr(store, adata, compressor=None)
    g = zarr.open_group(store, mode="r")
    adata = AnnData(X=sparse_dataset(g["X"]))
    expected = AnnData(X=ad.io.read_elem(g["X"]))
    assert isinstance(adata.X[index], sparse.csr_matrix)
    with zarr.config.enable_gpu():
        assert isinstance(adata.X[index], CupyCSRMatrix)
        assert_equal(expected.X[index], adata.X[index].get())
    assert isinstance(adata.X[index], sparse.csr_matrix)
