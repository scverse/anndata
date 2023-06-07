"""
For tests using dask
"""
import anndata as ad
import pandas as pd
from anndata._core.anndata import AnnData
import pytest
from anndata.tests.helpers import (
    as_dense_dask_array,
    GEN_ADATA_DASK_ARGS,
    gen_adata,
    assert_equal,
)
from anndata.compat import DaskArray

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
    import dask.array as da
    import numpy as np

    (M, N), chunks = sizes
    X = da.random.random((M, N), chunks=chunks)
    obs = pd.DataFrame(
        {"batch": np.random.choice(["a", "b"], M)},
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
    import numpy as np

    pth = tmp_path / f"test_write.{diskfmt}"
    write = lambda x, y: getattr(x, f"write_{diskfmt}")(y)
    read = lambda x: getattr(ad, f"read_{diskfmt}")(x)

    M, N = adata.X.shape
    adata.obsm["a"] = da.random.random((M, 10))
    adata.obsm["b"] = da.random.random((M, 10))
    adata.varm["a"] = da.random.random((N, 10))

    orig = adata
    write(orig, pth)
    curr = read(pth)

    with pytest.raises(Exception):
        assert_equal(curr.obsm["a"], curr.obsm["b"])

    assert_equal(curr.varm["a"], orig.varm["a"])
    assert_equal(curr.obsm["a"], orig.obsm["a"])

    assert isinstance(curr.X, np.ndarray)
    assert isinstance(curr.obsm["a"], np.ndarray)
    assert isinstance(curr.varm["a"], np.ndarray)
    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["a"], DaskArray)
    assert isinstance(orig.varm["a"], DaskArray)


def test_dask_to_memory_check_array_types(adata, tmp_path, diskfmt):
    import dask.array as da
    import numpy as np

    pth = tmp_path / f"test_write.{diskfmt}"
    write = lambda x, y: getattr(x, f"write_{diskfmt}")(y)
    read = lambda x: getattr(ad, f"read_{diskfmt}")(x)

    M, N = adata.X.shape
    adata.obsm["a"] = da.random.random((M, 10))
    adata.obsm["b"] = da.random.random((M, 10))
    adata.varm["a"] = da.random.random((N, 10))

    orig = adata
    write(orig, pth)
    curr = read(pth)

    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["a"], DaskArray)
    assert isinstance(orig.varm["a"], DaskArray)

    mem = orig.to_memory()

    with pytest.raises(Exception):
        assert_equal(curr.obsm["a"], curr.obsm["b"])

    assert_equal(curr.varm["a"], orig.varm["a"])
    assert_equal(curr.obsm["a"], orig.obsm["a"])
    assert_equal(mem.obsm["a"], orig.obsm["a"])
    assert_equal(mem.varm["a"], orig.varm["a"])

    assert isinstance(curr.X, np.ndarray)
    assert isinstance(curr.obsm["a"], np.ndarray)
    assert isinstance(curr.varm["a"], np.ndarray)
    assert isinstance(mem.X, np.ndarray)
    assert isinstance(mem.obsm["a"], np.ndarray)
    assert isinstance(mem.varm["a"], np.ndarray)
    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["a"], DaskArray)
    assert isinstance(orig.varm["a"], DaskArray)


def test_dask_to_memory_copy_check_array_types(adata, tmp_path, diskfmt):
    import dask.array as da
    import numpy as np

    pth = tmp_path / f"test_write.{diskfmt}"
    write = lambda x, y: getattr(x, f"write_{diskfmt}")(y)
    read = lambda x: getattr(ad, f"read_{diskfmt}")(x)

    M, N = adata.X.shape
    adata.obsm["a"] = da.random.random((M, 10))
    adata.obsm["b"] = da.random.random((M, 10))
    adata.varm["a"] = da.random.random((N, 10))

    orig = adata
    write(orig, pth)
    curr = read(pth)

    mem = orig.to_memory(copy=True)

    with pytest.raises(Exception):
        assert_equal(curr.obsm["a"], curr.obsm["b"])

    assert_equal(curr.varm["a"], orig.varm["a"])
    assert_equal(curr.obsm["a"], orig.obsm["a"])
    assert_equal(mem.obsm["a"], orig.obsm["a"])
    assert_equal(mem.varm["a"], orig.varm["a"])

    assert isinstance(curr.X, np.ndarray)
    assert isinstance(curr.obsm["a"], np.ndarray)
    assert isinstance(curr.varm["a"], np.ndarray)
    assert isinstance(mem.X, np.ndarray)
    assert isinstance(mem.obsm["a"], np.ndarray)
    assert isinstance(mem.varm["a"], np.ndarray)
    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["a"], DaskArray)
    assert isinstance(orig.varm["a"], DaskArray)


def test_dask_copy_check_array_types(adata):
    import dask.array as da

    M, N = adata.X.shape
    adata.obsm["a"] = da.random.random((M, 10))
    adata.obsm["b"] = da.random.random((M, 10))
    adata.varm["a"] = da.random.random((N, 10))

    orig = adata
    curr = adata.copy()

    with pytest.raises(Exception):
        assert_equal(curr.obsm["a"], curr.obsm["b"])

    assert_equal(curr.varm["a"], orig.varm["a"])
    assert_equal(curr.obsm["a"], orig.obsm["a"])

    assert isinstance(curr.X, DaskArray)
    assert isinstance(curr.obsm["a"], DaskArray)
    assert isinstance(curr.varm["a"], DaskArray)
    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["a"], DaskArray)
    assert isinstance(orig.varm["a"], DaskArray)


def test_assign_X(adata):
    """Check if assignment works"""
    import dask.array as da
    import numpy as np
    from anndata.compat import DaskArray

    adata.X = da.ones(adata.X.shape)
    prev_type = type(adata.X)
    adata_copy = adata.copy()

    adata.X = -1 * da.ones(adata.X.shape)
    assert prev_type is DaskArray and type(adata_copy.X) is DaskArray
    assert_equal(adata.X, -1 * np.ones(adata.X.shape))
    assert_equal(adata_copy.X, np.ones(adata.X.shape))


# Test if dask arrays turn into numpy arrays after to_memory is called
def test_dask_to_memory_unbacked():
    import numpy as np

    orig = gen_adata((15, 10), X_type=as_dense_dask_array, **GEN_ADATA_DASK_ARGS)
    orig.uns = {"da": {"da": as_dense_dask_array(np.ones(12))}}

    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["da"], DaskArray)
    assert isinstance(orig.layers["da"], DaskArray)
    assert isinstance(orig.varm["da"], DaskArray)
    assert isinstance(orig.uns["da"]["da"], DaskArray)

    curr = orig.to_memory()

    assert_equal(orig, curr)
    assert isinstance(curr.X, np.ndarray)
    assert isinstance(curr.obsm["da"], np.ndarray)
    assert isinstance(curr.varm["da"], np.ndarray)
    assert isinstance(curr.layers["da"], np.ndarray)
    assert isinstance(curr.uns["da"]["da"], np.ndarray)
    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["da"], DaskArray)
    assert isinstance(orig.layers["da"], DaskArray)
    assert isinstance(orig.varm["da"], DaskArray)
    assert isinstance(orig.uns["da"]["da"], DaskArray)


# Test if dask arrays turn into numpy arrays after to_memory is called
def test_dask_to_memory_copy_unbacked():
    import numpy as np

    orig = gen_adata((15, 10), X_type=as_dense_dask_array, **GEN_ADATA_DASK_ARGS)
    orig.uns = {"da": {"da": as_dense_dask_array(np.ones(12))}}

    curr = orig.to_memory(copy=True)

    assert_equal(orig, curr)
    assert isinstance(curr.X, np.ndarray)
    assert isinstance(curr.obsm["da"], np.ndarray)
    assert isinstance(curr.varm["da"], np.ndarray)
    assert isinstance(curr.layers["da"], np.ndarray)
    assert isinstance(curr.uns["da"]["da"], np.ndarray)
    assert isinstance(orig.X, DaskArray)
    assert isinstance(orig.obsm["da"], DaskArray)
    assert isinstance(orig.layers["da"], DaskArray)
    assert isinstance(orig.varm["da"], DaskArray)
    assert isinstance(orig.uns["da"]["da"], DaskArray)


def test_to_memory_raw():
    import numpy as np
    import dask.array as da

    orig = gen_adata((20, 10), **GEN_ADATA_DASK_ARGS)
    orig.X = da.ones((20, 10))

    with_raw = orig[:, ::2].copy()
    with_raw.raw = orig.copy()

    assert isinstance(with_raw.raw.X, DaskArray)
    assert isinstance(with_raw.raw.varm["da"], DaskArray)

    curr = with_raw.to_memory()

    assert isinstance(with_raw.raw.X, DaskArray)
    assert isinstance(with_raw.raw.varm["da"], DaskArray)
    assert isinstance(curr.raw.X, np.ndarray)
    assert isinstance(curr.raw.varm["da"], np.ndarray)


def test_to_memory_copy_raw():
    import numpy as np
    import dask.array as da

    orig = gen_adata((20, 10), **GEN_ADATA_DASK_ARGS)
    orig.X = da.ones((20, 10))

    with_raw = orig[:, ::2].copy()
    with_raw.raw = orig.copy()

    assert isinstance(with_raw.raw.X, DaskArray)
    assert isinstance(with_raw.raw.varm["da"], DaskArray)

    curr = with_raw.to_memory(copy=True)

    assert isinstance(with_raw.raw.X, DaskArray)
    assert isinstance(with_raw.raw.varm["da"], DaskArray)
    assert isinstance(curr.raw.X, np.ndarray)
    assert isinstance(curr.raw.varm["da"], np.ndarray)
