from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
import zarr
from scipy import sparse

from anndata import AnnData, Raw
from anndata._core.sparse_dataset import sparse_dataset
from anndata.compat import CupyCSRMatrix

if TYPE_CHECKING:
    from pathlib import Path


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
def test_get_with_zarr_gpu(tmp_path: Path):
    adata = AnnData(X=sparse.random(50, 100, format="csr"))
    zarr_path = tmp_path / "gpu_adata.zarr"
    adata.write_zarr(zarr_path)
    g = zarr.open_group(zarr_path, mode="r")
    adata = AnnData(X=sparse_dataset(g["X"]))
    with zarr.config.enable_gpu():
        assert isinstance(adata.X[...], CupyCSRMatrix)
