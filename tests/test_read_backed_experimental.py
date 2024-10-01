from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import anndata as ad
from anndata import AnnData
from anndata.experimental import read_lazy
from anndata.tests.helpers import (
    AccessTrackingStore,
    as_dense_dask_array,
    assert_equal,
    gen_adata,
    gen_typed_df,
)

if TYPE_CHECKING:
    from typing import Literal


@pytest.fixture(
    params=[sparse.csr_matrix, sparse.csc_matrix, np.array, as_dense_dask_array],
    ids=["scipy-csr", "scipy-csc", "np-array", "dask_array"],
    scope="session",
)
def mtx_format(request):
    return request.param


@pytest.fixture(params=["zarr", "h5ad"], scope="session")
def dskfmt(request):
    return request.param


@pytest.fixture(params=[True, False], scope="session")
def load_annotation_index(request):
    return request.param


@pytest.fixture(scope="session")
def adata_remote_orig(
    tmp_path_factory, dskfmt: str, mtx_format, load_annotation_index: bool
) -> tuple[AnnData, AnnData]:
    """Create remote fixtures, one without a range index and the other with"""
    orig_path = tmp_path_factory.mktemp(f"orig.{dskfmt}")
    orig = gen_adata((1000, 1000), mtx_format)
    orig.write_zarr(orig_path)
    return read_lazy(orig_path, load_annotation_index=load_annotation_index), orig


@pytest.fixture
def adata_remote_with_store_tall_skinny(
    tmp_path_factory, mtx_format
) -> tuple[AnnData, AccessTrackingStore]:
    orig_path = tmp_path_factory.mktemp("orig.zarr")
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
    orig.write_zarr(orig_path)
    store = AccessTrackingStore(orig_path)
    remote = read_lazy(store)
    return remote, store


pytestmark = pytest.mark.skipif(
    not find_spec("xarray"), reason="Xarray is not installed"
)


def test_access_count_obs_var(adata_remote_with_store_tall_skinny):
    remote, store = adata_remote_with_store_tall_skinny
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
    store.assert_access_count("obs/int64", 0)
    store.assert_access_count("obs/cat/categories", 0)
    subset = remote[remote.obs["cat"] == "a", :]
    subset.obs["int64"]
    sub_subset = subset[0:10, :]
    sub_subset.obs["int64"]
    sub_subset.var["int64"]
    store.assert_access_count("X", 0)
    store.assert_access_count("obs/int64", 0)
    store.assert_access_count("var/int64", 0)
    # all codes read in for subset (from 4 chunks)
    store.assert_access_count("obs/cat/codes", 4)
    remote[0:10, :].obs["int64"][0:10].compute()
    store.assert_access_count("obs/int64", 1)
    # .zmetadata handles .zarray so simple access does not cause any read
    store.assert_access_count("var/int64", 0)


def test_access_count_index(adata_remote_with_store_tall_skinny):
    _, store = adata_remote_with_store_tall_skinny
    store.reset_key_trackers()
    store.initialize_key_trackers(["obs/_index"])
    read_lazy(store, load_annotation_index=False)
    store.assert_access_count("obs/_index", 0)
    read_lazy(store)
    # 16 is number of chunks
    store.assert_access_count("obs/_index", 16)


def test_access_count_dtype(adata_remote_with_store_tall_skinny):
    remote, store = adata_remote_with_store_tall_skinny
    store.initialize_key_trackers(["obs/cat/categories"])
    store.assert_access_count("obs/cat/categories", 0)
    # This should only cause categories to be read in once
    remote.obs["cat"].dtype
    remote.obs["cat"].dtype
    remote.obs["cat"].dtype
    store.assert_access_count("obs/cat/categories", 1)


def test_to_memory(adata_remote_orig):
    remote, orig = adata_remote_orig
    remote_to_memory = remote.to_memory()
    assert_equal(remote_to_memory, orig)


def test_view_to_memory(adata_remote_orig):
    remote, orig = adata_remote_orig
    subset_obs = orig.obs["obs_cat"] == "a"
    assert_equal(orig[subset_obs, :], remote[subset_obs, :].to_memory())

    subset_var = orig.var["var_cat"] == "a"
    assert_equal(orig[:, subset_var], remote[:, subset_var].to_memory())


def test_view_of_view_to_memory(adata_remote_orig):
    remote, orig = adata_remote_orig
    subset_obs = (orig.obs["obs_cat"] == "a") | (orig.obs["obs_cat"] == "b")
    subsetted_adata = orig[subset_obs, :]
    subset_subset_obs = subsetted_adata.obs["obs_cat"] == "b"
    subsetted_subsetted_adata = subsetted_adata[subset_subset_obs, :]
    assert_equal(
        subsetted_subsetted_adata,
        remote[subset_obs, :][subset_subset_obs, :].to_memory(),
    )

    subset_var = (orig.var["var_cat"] == "a") | (orig.var["var_cat"] == "b")
    subsetted_adata = orig[:, subset_var]
    subset_subset_var = subsetted_adata.var["var_cat"] == "b"
    subsetted_subsetted_adata = subsetted_adata[:, subset_subset_var]
    assert_equal(
        subsetted_subsetted_adata,
        remote[:, subset_var][:, subset_subset_var].to_memory(),
    )


def test_unconsolidated(tmp_path, mtx_format):
    adata = gen_adata((1000, 1000), mtx_format)
    orig_pth = tmp_path / "orig.zarr"
    adata.write_zarr(orig_pth)
    (orig_pth / ".zmetadata").unlink()
    store = AccessTrackingStore(orig_pth)
    store.initialize_key_trackers(["obs/.zgroup", ".zgroup"])
    with pytest.warns(UserWarning, match=r"Did not read zarr as consolidated"):
        remote = read_lazy(store)
    remote_to_memory = remote.to_memory()
    assert_equal(remote_to_memory, adata)
    store.assert_access_count("obs/.zgroup", 1)


@pytest.mark.parametrize("join", ["outer", "inner"])
@pytest.mark.parametrize("are_vars_different", [True, False])
def test_concat(
    tmp_path, join: Literal["outer", "inner"], are_vars_different: bool, mtx_format
):
    from anndata.experimental.backed._compat import Dataset

    lazy_adatas = []
    adatas = []
    stores: list[AccessTrackingStore] = []
    var_indices = []
    M = 1000
    N = 50
    n_datasets = 2
    for dataset_index in range(n_datasets):
        orig_path = tmp_path / f"orig_{dataset_index}.zarr"
        orig_path.mkdir()
        obs_names = pd.Index(f"cell_{dataset_index}_{i}" for i in range(M))
        var_names = pd.Index(
            f"gene_{i}{f'_{dataset_index}_ds' if are_vars_different and (i % 2) else ''}"
            for i in range(N)
        )
        var_indices.append(var_names)
        obs = gen_typed_df(M, obs_names)
        var = gen_typed_df(N, var_names)
        orig = AnnData(
            obs=obs,
            var=var,
            X=mtx_format(np.random.binomial(100, 0.005, (M, N)).astype(np.float32)),
        )
        orig.write_zarr(orig_path)
        store = AccessTrackingStore(orig_path)
        store.initialize_key_trackers(["obs/int64", "var/int64", "X"])
        lazy_adatas += [read_lazy(store)]
        adatas += [orig]
        stores += [store]
    concated_remote = ad.concat(lazy_adatas, join=join)
    assert isinstance(concated_remote.obs, Dataset)
    for i in range(n_datasets):
        stores[i].assert_access_count("obs/int64", 0)
        stores[i].assert_access_count("X", 0)
        stores[i].assert_access_count("var/int64", 0)
    concatenated_memory = ad.concat(adatas, join=join)
    # account for differences

    # name is lost normally, should fix
    obs_memory = concatenated_memory.obs
    obs_memory.index.name = "obs_names"

    # remote has object dtype, need to convert back for integers booleans etc.
    def correct_extension_dtype_differences(remote: pd.DataFrame, memory: pd.DataFrame):
        for col in memory.columns:
            dtype = memory[col].dtype
            if pd.api.types.is_extension_array_dtype(dtype):
                remote[col] = remote[col].astype(dtype)
        return remote, memory

    assert_equal(
        *correct_extension_dtype_differences(
            concated_remote.obs.to_pandas(), obs_memory
        )
    )
    # check non-different variables, taken from first annotation.  all others are null so incomparable
    pd_index = pd.Index(filter(lambda x: not x.endswith("ds"), var_indices[0]))
    var_df = adatas[0][:, pd_index].var.copy()
    var_df.index.name = "var_names"
    remote_df_corrected, _ = correct_extension_dtype_differences(
        concated_remote[:, pd_index].var.to_pandas(), var_df
    )
    #  TODO:xr.merge always upcasts to float due to NA and you can't downcast?
    for col in remote_df_corrected.columns:
        dtype = remote_df_corrected[col].dtype
        if dtype in [np.float64, np.float32]:
            var_df[col] = var_df[col].astype(dtype)
    assert_equal(remote_df_corrected, var_df)

    assert_equal(concated_remote.X, concatenated_memory.X)
