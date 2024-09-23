from __future__ import annotations

import pytest
from scipy import sparse

from anndata import AnnData, Raw


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
