from __future__ import annotations

from importlib.util import find_spec
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from anndata._core.anndata import AnnData
from anndata.experimental import read_backed
from anndata.tests.helpers import (
    AccessTrackingStore,
    as_dense_dask_array,
    assert_equal,
    gen_adata,
    gen_typed_df,
)


@pytest.fixture(
    params=[sparse.csr_matrix, sparse.csc_matrix, np.array, as_dense_dask_array],
    ids=["scipy-csr", "scipy-csc", "np-array", "dask_array"],
)
def mtx_format(request):
    return request.param


@pytest.fixture(params=[sparse.csr_matrix, sparse.csc_matrix])
def sparse_format(request):
    return request.param


@pytest.fixture(params=["zarr", "h5ad"])
def dskfmt(request):
    return request.param


needs_xarray = pytest.mark.skipif(
    not find_spec("xarray"), reason="Xarray is not installed"
)


@needs_xarray
def test_access_count_obs_var(tmp_path, mtx_format):
    base_pth = Path(tmp_path)
    orig_pth = base_pth / "orig.zarr"
    M = 1000000  # forces zarr to chunk `obs` columns multiple ways - that way 1 access to `int64` below is actually only one access
    N = 5
    obs_names = pd.Index(f"cell{i}" for i in range(M))
    var_names = pd.Index(f"gene{i}" for i in range(N))
    obs = gen_typed_df(M, obs_names)
    var = gen_typed_df(N, var_names)
    orig = AnnData(
        obs=obs,
        var=var,
        X=mtx_format(np.random.binomial(100, 0.005, (M, N)).astype(np.float32)),
    )
    orig.write_zarr(orig_pth)
    store = AccessTrackingStore(orig_pth)
    remote = read_backed(store)
    store.initialize_key_trackers(
        ["obs/cat/codes", "obs/cat/categories", "obs/int64", "var/int64", "X"]
    )
    # a series of methods that should __not__ read in any data
    remote.X  # the initial (non-subset) access to `X` should not read in data
    remote.shape
    remote.var
    remote.obs
    remote.obs["int64"]
    remote.obs["int64"]
    remote.obs["cat"]
    assert store.get_access_count("obs/int64") == 0, store.get_subkeys_accessed(
        "obs/int64"
    )
    assert (
        store.get_access_count("obs/cat/categories") == 0
    ), store.get_subkeys_accessed("obs/cat/categories")
    # This should only cause categories to be read in once
    remote.obs["cat"].dtype
    remote.obs["cat"].dtype
    remote.obs["cat"].dtype
    assert (
        store.get_access_count("obs/cat/categories") == 1
    ), store.get_subkeys_accessed("obs/cat/categories")
    subset = remote[
        (remote.obs["cat"] == "a").data, :
    ]  # `.data` for xarray, but should we handle internally?
    subset.obs["int64"]
    sub_subset = subset[0:10, :]
    sub_subset.obs["int64"]
    assert store.get_access_count("X") == 0, store.get_subkeys_accessed("X")
    assert store.get_access_count("obs/int64") == 0, store.get_subkeys_accessed(
        "obs/int64"
    )
    assert store.get_access_count("var/int64") == 0, store.get_subkeys_accessed(
        "var/int64"
    )
    # all codes read in for subset
    assert store.get_access_count("obs/cat/codes") == 4, store.get_subkeys_accessed(
        "obs/cat/codes"
    )
    remote[0:10, :].obs["int64"][0:10].compute()
    assert store.get_access_count("obs/int64") == 1, store.get_subkeys_accessed(
        "obs/int64"
    )
    # one for 0, .zmetadata handles .zarray
    assert store.get_access_count("var/int64") == 0, store.get_subkeys_accessed(
        "var/int64"
    )  # never accessed


@needs_xarray
def test_to_memory(tmp_path, mtx_format, dskfmt):
    adata = gen_adata((1000, 1000), mtx_format)
    base_pth = Path(tmp_path)
    orig_pth = base_pth / f"orig.{dskfmt}"
    write = lambda x: getattr(x, f"write_{dskfmt}")(orig_pth)
    write(adata)
    remote = read_backed(orig_pth)
    remote_to_memory = remote.to_memory()
    assert_equal(remote_to_memory, adata)


@needs_xarray
def test_view_to_memory(tmp_path, mtx_format, dskfmt):
    adata = gen_adata((1000, 1000), mtx_format)
    base_pth = Path(tmp_path)
    orig_pth = base_pth / f"orig.{dskfmt}"
    write = lambda x: getattr(x, f"write_{dskfmt}")(orig_pth)
    write(adata)
    remote = read_backed(orig_pth)
    subset_obs = adata.obs["obs_cat"] == "a"
    assert_equal(adata[subset_obs, :], remote[subset_obs, :].to_memory())

    subset_var = adata.var["var_cat"] == "a"
    assert_equal(adata[:, subset_var], remote[:, subset_var].to_memory())


@needs_xarray
def test_view_of_view_to_memory(tmp_path, mtx_format, dskfmt):
    adata = gen_adata((1000, 1000), mtx_format)
    base_pth = Path(tmp_path)
    orig_pth = base_pth / f"orig.{dskfmt}"
    write = lambda x: getattr(x, f"write_{dskfmt}")(orig_pth)
    write(adata)
    remote = read_backed(orig_pth)
    subset_obs = (adata.obs["obs_cat"] == "a") | (adata.obs["obs_cat"] == "b")
    subsetted_adata = adata[subset_obs, :]
    subset_subset_obs = subsetted_adata.obs["obs_cat"] == "b"
    subsetted_subsetted_adata = subsetted_adata[subset_subset_obs, :]
    assert_equal(
        subsetted_subsetted_adata,
        remote[subset_obs, :][subset_subset_obs, :].to_memory(),
    )

    subset_var = (adata.var["var_cat"] == "a") | (adata.var["var_cat"] == "b")
    subsetted_adata = adata[:, subset_var]
    subset_subset_var = subsetted_adata.var["var_cat"] == "b"
    subsetted_subsetted_adata = subsetted_adata[:, subset_subset_var]
    assert_equal(
        subsetted_subsetted_adata,
        remote[:, subset_var][:, subset_subset_var].to_memory(),
    )


@needs_xarray
def test_unconsolidated(tmp_path, mtx_format):
    adata = gen_adata((1000, 1000), mtx_format)
    base_pth = Path(tmp_path)
    orig_pth = base_pth / "orig.zarr"
    write = lambda x: getattr(x, "write_zarr")(orig_pth)
    write(adata)
    (Path(orig_pth) / ".zmetadata").unlink()
    store = AccessTrackingStore(orig_pth)
    store.initialize_key_trackers(["obs/.zgroup", ".zgroup"])
    with pytest.warns(UserWarning, match=r"Did not read zarr as consolidated"):
        remote = read_backed(store)
    remote_to_memory = remote.to_memory()
    assert_equal(remote_to_memory, adata)
    assert store.get_access_count("obs/.zgroup") == 1, store.get_subkeys_accessed(
        "obs/.zgroup"
    )
