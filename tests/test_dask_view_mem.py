from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

import anndata as ad

if TYPE_CHECKING:
    import pandas as pd

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


@pytest.fixture(params=[True, False])
def give_chunks(request):
    return request.param


# ------------------------------------------------------------------------------
# The test functions
# ------------------------------------------------------------------------------


# Does some stuff so that dask can cache the
# subclasscheck before the run.
@pytest.fixture
def _alloc_cache():
    import dask.array as da

    N = 2**6
    size = ((N, N), (N, N))

    adata = ad.AnnData(
        da.random.random(*size),
        layers=dict(m=da.random.random(*size)),
        obsm=dict(m=da.random.random(*size)),
        obs=dict(m=da.random.random(N)),
        var=dict(m=da.random.random(N)),
        varm=dict(m=da.random.random(*size)),
    )
    subset = adata[:10, :][:, :10]
    for mn in ["varm", "obsm", "layers"]:
        m = getattr(subset, mn)["m"]
        m[0, 0] = 100
    _ = adata.to_memory(copy=False)


# Theoretically this is expected to allocate:
# N*N*4 bytes per matrix (we have 2).
# N*4 bytes per index (we have 1).
# N*N*(2**3) + N*(2**2) bytes
# N*N*(2**3) + N*(2**2) bytes
# 2**19 + 2**10
# if we put a 2 factor on 2**19
# the results seems more accurate with the experimental results
# For example from dask.random we allocate 1mb
@pytest.mark.usefixtures("_alloc_cache")
@pytest.mark.limit_memory("1.5 MB")
def test_size_of_view(mapping_name, give_chunks):
    import dask.array as da

    N = 2**8
    size = ((N, N), (N, N)) if give_chunks else ((N, N), "auto")

    adata = ad.AnnData(
        da.random.random(*size),
        **{mapping_name: dict(m=da.random.random(*size))},
    )
    _ = adata.to_memory(copy=False)


# Normally should expect something around 90 kbs
# Pandas does some indexing stuff that requires more sometimes
# since the array we allocated would be 4mb for both arrays + 2mb
# Thus, if we allocated it all it should at least have 6mb
# experimentally we should at least have 10mb
# for index this should be ok
@pytest.mark.usefixtures("_alloc_cache")
@pytest.mark.limit_memory("1.5 MB")
def test_modify_view_mapping_component_memory(mapping_name, give_chunks):
    import dask.array as da

    N = 2**8
    M = 2**9

    size = ((M, M), (M, M)) if give_chunks else ((M, M), "auto")

    adata = ad.AnnData(
        da.random.random(*size),
        **{mapping_name: dict(m=da.random.random(*size))},
    )
    subset = adata[:N, :N]
    assert subset.is_view
    m = getattr(subset, mapping_name)["m"]
    m[0, 0] = 100


# Normally should expect something around 90 kbs
# Pandas does some indexing stuff that requires more sometimes
# since the array we allocated would be 4mb for both arrays + 2mb
# Thus, if we allocated it all it should at least have 6mb
# experimentally we should at least have 10mb
# for index this should be ok
@pytest.mark.usefixtures("_alloc_cache")
@pytest.mark.limit_memory("1.5 MB")
def test_modify_view_X_memory(mapping_name, give_chunks):
    import dask.array as da

    N = 2**8
    M = 2**9

    size = ((M, M), (M, M)) if give_chunks else ((M, M), "auto")

    adata = ad.AnnData(
        da.random.random(*size),
        **{mapping_name: dict(m=da.random.random(*size))},
    )
    subset = adata[:N, :N]
    assert subset.is_view
    m = subset.X
    with pytest.warns(
        ad.ImplicitModificationWarning,
        match=r"Trying to modify attribute `.X` of view, initializing view as actual.",
    ):
        m[0, 0] = 100


# Normally should expect something around 90 kbs
# Pandas does some indexing stuff that requires more sometimes
# since the array we allocated would be 4mb for both arrays + 2mb
# Thus, if we allocated it all it should at least have 6mb
# experimentally we should at least have 10mb
# for index this should be ok
@pytest.mark.usefixtures("_alloc_cache")
@pytest.mark.limit_memory("1.5 MB")
def test_modify_view_mapping_obs_var_memory(attr_name, give_chunks):
    import dask.array as da

    N = 2**8
    M = 2**9

    size = ((M, M), (M, M)) if give_chunks else ((M, M), "auto")

    adata = ad.AnnData(
        da.random.random(*size),
        **{attr_name: dict(m=da.random.random(M))},
    )
    subset = adata[:N, :N]
    assert subset.is_view
    m: pd.Series = getattr(subset, attr_name)["m"]
    m.iloc[0] = 100
