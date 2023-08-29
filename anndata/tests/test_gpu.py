import pytest
from anndata import AnnData, Raw
from scipy.sparse import isspmatrix_csr

@pytest.mark.gpu
def test_gpu():
    """
    For testing that the gpu mark works
    """
    import cupy  # This test shouldn't run if cupy isn't installed

    cupy.ones(1)

@pytest.mark.gpu
def test_adata_raw_gpu():
    from cupyx.scipy.sparse import random
    import cupy as cp
    adata = AnnData(X = random(500,50,density=0.01, format='csr', dtype=cp.float32))
    adata.raw = adata
    assert isspmatrix_csr(adata.raw.X)

@pytest.mark.gpu
def test_raw_gpu():
    from cupyx.scipy.sparse import random
    import cupy as cp
    adata = AnnData(X = random(500,50,density=0.01, format='csr', dtype=cp.float32))
    araw = Raw(adata)
    assert isspmatrix_csr(araw.X)
