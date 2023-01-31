import pytest

import anndata as ad

pytest.importorskip("pytest_memray")

# ------------------------------------------------------------------------------
# Some test data
# ------------------------------------------------------------------------------


@pytest.fixture(params=["layers", "obsm", "varm"])
def mapping_name(request):
    return request.param


@pytest.fixture(params=["obs", "var"])
def attr_name(request):
    return request.param


# ------------------------------------------------------------------------------
# The test functions
# ------------------------------------------------------------------------------


# Does some stuff so that dask can cache the
# subclasscheck before the run.
@pytest.fixture
def alloc_cache():
    import dask.array as da

    N = 2**6
    size = ((N, N), (N, N))

    adata = ad.AnnData(
        da.random.random(*size),
        **{
            "layers": dict(m=da.random.random(*size)),
            "obsm": dict(m=da.random.random(*size)),
            "obs": dict(m=da.random.random((N))),
            "var": dict(m=da.random.random((N))),
            "varm": dict(m=da.random.random(*size)),
        },
    )
    subset = adata[:10, :][:, :10]
    for mn in ["varm", "obsm", "layers"]:
        m = getattr(subset, mn)["m"]
        m[0, 0] = 100
    _ = adata.to_memory(copy=False)


# Based on previous results expected to allocate ~74 kb
# This tests our assumption which is used on test_modify_view_mapping_component_memory
@pytest.mark.usefixtures("alloc_cache")
@pytest.mark.limit_memory("80 KB")
def test_size_of_view(mapping_name):
    import dask.array as da

    N = 2**6
    size = ((N, N), (N, N))

    adata = ad.AnnData(
        da.random.random(*size),
        **{mapping_name: dict(m=da.random.random(*size))},
    )
    _ = adata.to_memory(copy=False)


# Normally should expect something around 80 kbs
# Pandas does some indexing stuff that requires 256kb sometimes
# since the array we allocated would be 4mb this should be ok
@pytest.mark.usefixtures("alloc_cache")
@pytest.mark.limit_memory("340 KB")
def test_modify_view_mapping_component_memory(mapping_name):
    import dask.array as da

    N = 2**6
    M = 2**9

    size = ((M, M), (M, M))

    adata = ad.AnnData(
        da.random.random(*size),
        **{mapping_name: dict(m=da.random.random(*size))},
    )
    subset = adata[:N, :][:, :N]
    assert subset.is_view
    m = getattr(subset, mapping_name)["m"]
    m[0, 0] = 100


# Normally should expect something around 80 kbs
# Pandas does some indexing stuff that requires 256kb sometimes
# since the array we allocated would be 4mb this should be ok
@pytest.mark.usefixtures("alloc_cache")
@pytest.mark.limit_memory("340 KB")
def test_modify_view_X_memory(mapping_name):
    import dask.array as da

    N = 2**6
    M = 2**9

    size = ((M, M), (M, M))

    adata = ad.AnnData(
        da.random.random(*size),
        **{mapping_name: dict(m=da.random.random(*size))},
    )
    subset = adata[:N, :][:, :N]
    assert subset.is_view
    m = subset.X
    m[0, 0] = 100


# Normally should expect something around 80 kbs
# Pandas does some indexing stuff that requires 256kb sometimes
# since the array we allocated would be 4mb this should be ok


@pytest.mark.usefixtures("alloc_cache")
@pytest.mark.limit_memory("340 KB")
def test_modify_view_mapping_obs_var_memory(attr_name):
    import dask.array as da

    N = 2**6
    M = 2**9

    size = ((M, M), (M, M))

    adata = ad.AnnData(
        da.random.random(*size),
        **{attr_name: dict(m=da.random.random(M))},
    )
    subset = adata[:N, :][:, :N]
    assert subset.is_view
    m = getattr(subset, attr_name)['m']
    m[0] = 100
