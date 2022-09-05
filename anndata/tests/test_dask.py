"""
For tests using dask
"""
import anndata as ad
import pandas as pd
from anndata._core.anndata import AnnData
import numpy.testing as npt
import pytest

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

# TODO: Appropriate path to write in testing.
# TODO: Check written object


def test_dask_write_h5ad():
    import dask.array as da

    M, N = 50, 30
    adata = ad.AnnData(
        obs=pd.DataFrame(index=[f"cell{i:02}" for i in range(M)]),
        var=pd.DataFrame(index=[f"gene{i:02}" for i in range(N)]),
    )
    adata.X = da.ones((M, N))
    view = adata[:30]
    view.copy()
    view.write_h5ad("test.h5ad")

# TODO: Appropriate path to write in testing.
# TODO: Check written object


def test_dask_write():
    import dask.array as da

    M, N = 50, 30
    adata = ad.AnnData(
        obs=pd.DataFrame(index=[f"cell{i:02}" for i in range(M)]),
        var=pd.DataFrame(index=[f"gene{i:02}" for i in range(N)]),
    )
    adata.X = da.ones((M, N))
    view = adata[:30]
    view.copy()
    view.write("test.h5ad")


@pytest.mark.parametrize(
    "shape,chunks",
    [
        [(200,100),(100,100)],
        [(200,100),(100,100)],
        [(20,10),(1,1)],
        [(20,10),(1,1)],
    ]
)
def test_assign(shape,chunks,idx):
    """Check if setting the given indices work"""
    import dask.array as da
    import numpy as np

    darr = da.from_array(np.ones(shape),chunks=chunks)
    adata = AnnData(darr)
    adata_copy = adata.copy()
    adata.X = -1 * da.from_array(np.ones(shape),chunks=chunks)
    npt.assert_array_equal(adata.X,-1 * np.ones(shape))
    npt.assert_allclose(adata_copy.X,np.ones(shape))

    npt.assert_equal(adata.X[idx[0],idx[1]] , -1)
    npt.assert_equal(adata_copy.X[idx[0],idx[1]] , 1)


@pytest.mark.parametrize(
    "shape,chunks,idx",
    [
        [(200,100),(100,100),(2,3)],
        [(200,100),(100,100),(50,73)],
        [(20,10),(1,1),(0,0)],
        [(20,10),(1,1),(4,3)],
    ]

)
def test_idx_2d(shape,chunks,idx):
    """Check if setting the given indices work"""
    import dask.array as da
    import numpy as np

    darr = da.from_array(np.ones(shape),chunks=chunks)
    adata = AnnData(darr)
    adata_copy = adata.copy()
    adata.X = -1 * da.from_array(np.ones(shape),chunks=chunks)

    npt.assert_equal(adata.X[idx[0],idx[1]] , -1)
    npt.assert_equal(adata_copy.X[idx[0],idx[1]] , 1)
