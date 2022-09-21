"""
For tests using dask
"""
import anndata as ad
import pandas as pd
from anndata._core.anndata import AnnData
import pytest
from anndata.tests.helpers import assert_equal, asarray

pytest.importorskip("dask.array")


@pytest.fixture(
    params=[
        [(2000, 1000), (100, 100)],
        [(200, 100), (100, 100)],
        [(200, 100), (100, 100)],
        [(20, 10), (1, 1)],
        [(20, 10), (1, 1)],
    ]
)
def sizes(request):
    return request.param


@pytest.fixture(params=["h5ad", "zarr"])
def diskfmt(request):
    return request.param


@pytest.fixture
def adata(sizes):
    import numpy as np
    import dask.array as da

    (M, N), chunks = sizes
    X = da.random.random((M, N), chunks=chunks)
    obs = pd.DataFrame(
        dict(batch=np.array(["a", "b"])[np.random.randint(0, 2, M)]),
        index=[f"cell{i:03d}" for i in range(M)],
    )
    var = pd.DataFrame(index=[f"gene{i:03d}" for i in range(N)])

    return AnnData(X, obs=obs, var=var)


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


def test_dask_write(adata, tmp_path, diskfmt):
    import dask.array as da

    pth = tmp_path / f"test_write.{diskfmt}"
    write = lambda x, y: getattr(x, f"write_{diskfmt}")(y)
    read = lambda x: getattr(ad, f"read_{diskfmt}")(x)

    M, N = adata.X.shape
    adata.obsm["a"] = da.ones((M, 10))
    adata.varm["a"] = da.ones((N, 10))

    orig = adata
    write(orig, pth)
    curr = read(pth)
    assert da.equal(curr.varm["a"], orig.varm["a"]).all()
    assert da.equal(curr.obsm["a"], orig.obsm["a"]).all()
    # assert_equal(asarray(orig), asarray(curr)) # find way to compare


def test_assign_X(adata):
    """Check if assignment works"""
    import dask.array as da
    import numpy as np

    adata.X = da.ones(adata.X.shape)

    # This won't work since the setter converts the data as ndarray
    adata_copy = adata.copy()

    adata.X = -1 * da.ones(adata.X.shape)
    assert_equal(asarray(adata.X), -1 * np.ones(adata.X.shape))
    assert_equal(asarray(adata_copy.X), np.ones(adata.X.shape))


# TODO: Views


def test_assert_equal_dask_arrays():

    import dask.array as da

    a = da.from_array([[1, 2, 3], [4, 5, 6]])
    b = da.from_array([[1, 2, 3], [4, 5, 6]])

    assert_equal(a, b)

    c = da.ones(10, dtype="int32")
    d = da.ones(10, dtype="int64")
    assert_equal(c, d)
