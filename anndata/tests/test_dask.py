"""
For tests using dask
"""
import anndata as ad
import pandas as pd
from anndata._core.anndata import AnnData
import pytest
from anndata.tests.helpers import assert_equal, asarray

pytest.importorskip("dask.array")


def test_dask_X_view():
    import dask.array as da

    M, N = 50, 30
    adata = ad.AnnData(
        obs=pd.DataFrame(index=[f"cell{i:02}" for i in range(M)]),
        var=pd.DataFrame(index=[f"gene{i:02}" for i in range(N)]),
    )
    adata.X = da.ones((M, N))
    view = adata[:30]
    view.copy()


@pytest.fixture(params=["h5ad", "zarr"])
def diskfmt(request):
    return request.param


def test_dask_write(tmp_path, diskfmt):
    import dask.array as da

    pth = tmp_path / f"test_write.{diskfmt}"
    write = lambda x, y: getattr(x, f"write_{diskfmt}")(y)
    read = lambda x: getattr(ad, f"read_{diskfmt}")(x)
    M, N = 50, 30
    adata = ad.AnnData(
        obs=pd.DataFrame(index=[f"cell{i:02}" for i in range(M)]),
        var=pd.DataFrame(index=[f"gene{i:02}" for i in range(N)]),
    )
    adata.X = da.random.randint((M, N))
    write(adata, pth)
    result = read(pth)
    assert_equal(adata, result)


@pytest.mark.parametrize(
    "shape,chunks",
    [
        [(200, 100), (100, 100)],
        [(200, 100), (100, 100)],
        [(20, 10), (1, 1)],
        [(20, 10), (1, 1)],
    ],
)
def test_assign_X(shape, chunks):
    """Check if assignment works"""
    import dask.array as da
    import numpy as np

    darr = da.from_array(np.ones(shape), chunks=chunks)
    adata = AnnData(darr)
    adata_copy = adata.copy()
    adata.X = -1 * da.from_array(np.ones(shape), chunks=chunks)
    assert_equal(asarray(adata.X), -1 * np.ones(shape))
    assert_equal(asarray(adata_copy.X), np.ones(shape))


@pytest.mark.parametrize(
    "shape,chunks,idx",
    [
        [(200, 100), (100, 100), (2, 3)],
        [(200, 100), (100, 100), (50, 73)],
        [(20, 10), (1, 1), (0, 0)],
        [(20, 10), (1, 1), (4, 3)],
    ],
)
def test_idx_2d(shape, chunks, idx):
    """Check if setting the given indices work"""
    import dask.array as da
    import numpy as np

    darr = da.from_array(np.ones(shape), chunks=chunks)
    adata = AnnData(darr)
    adata_copy = adata.copy()
    adata.X = -1 * da.from_array(np.ones(shape), chunks=chunks)

    assert_equal(adata.X[idx[0], idx[1]], -1)
    assert_equal(adata_copy.X[idx[0], idx[1]], 1)
