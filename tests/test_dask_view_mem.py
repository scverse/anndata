from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

import anndata as ad

if TYPE_CHECKING:
    from typing import Literal

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
def _alloc_cache() -> None:
    import dask.array as da

    n = 2**6
    size = ((n, n), (n, n))

    adata = ad.AnnData(
        da.random.random(*size),
        layers=dict(m=da.random.random(*size)),
        obsm=dict(m=da.random.random(*size)),
        obs=dict(m=da.random.random(n)),
        var=dict(m=da.random.random(n)),
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
# As of 2025.09.* dask, this needs a bit more than the previous 1.5mb.
# TODO: Why?
@pytest.mark.usefixtures("_alloc_cache")
@pytest.mark.limit_memory("2.2 MB")
def test_size_of_view(
    *, mapping_name: Literal["layers", "obsm", "varm"], give_chunks: bool
) -> None:
    import dask.array as da

    n = 2**8
    size = ((n, n), (n, n)) if give_chunks else ((n, n), "auto")

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
def test_modify_view_mapping_component_memory(
    *, mapping_name: Literal["layers", "obsm", "varm"], give_chunks: bool
) -> None:
    import dask.array as da

    m, n = 2**9, 2**8

    size = ((m, m), (m, m)) if give_chunks else ((m, m), "auto")

    adata = ad.AnnData(
        da.random.random(*size),
        **{mapping_name: dict(m=da.random.random(*size))},
    )
    subset = adata[:n, :n]
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
def test_modify_view_x_memory(
    *, mapping_name: Literal["layers", "obsm", "varm"], give_chunks: bool
) -> None:
    import dask.array as da

    m, n = 2**9, 2**8

    size = ((m, m), (m, m)) if give_chunks else ((m, m), "auto")

    adata = ad.AnnData(
        da.random.random(*size),
        **{mapping_name: dict(m=da.random.random(*size))},
    )
    subset = adata[:n, :n]
    assert subset.is_view
    x = subset.X
    with pytest.warns(
        ad.ImplicitModificationWarning,
        match=r"Trying to modify attribute `.X` of view, initializing view as actual.",
    ):
        x[0, 0] = 100


# Normally should expect something around 90 kbs
# Pandas does some indexing stuff that requires more sometimes
# since the array we allocated would be 4mb for both arrays + 2mb
# Thus, if we allocated it all it should at least have 6mb
# experimentally we should at least have 10mb
# for index this should be ok
@pytest.mark.usefixtures("_alloc_cache")
@pytest.mark.limit_memory("1.5 MB")
def test_modify_view_mapping_obs_var_memory(
    *, attr_name: Literal["obs", "var"], give_chunks: bool
) -> None:
    import dask.array as da

    m, n = 2**9, 2**8

    size = ((m, m), (m, m)) if give_chunks else ((m, m), "auto")

    adata = ad.AnnData(
        da.random.random(*size),
        **{attr_name: dict(m=da.random.random(m))},
    )
    subset = adata[:n, :n]
    assert subset.is_view
    m: pd.Series = getattr(subset, attr_name)["m"]
    m.iloc[0] = 100
