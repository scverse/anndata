from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest

from anndata import AnnData
from anndata.compat import DaskArray
from anndata.experimental import read_lazy
from anndata.tests.helpers import (
    GEN_ADATA_NO_XARRAY_ARGS,
    AccessTrackingStore,
    assert_equal,
    gen_adata,
)

from .conftest import ANNDATA_ELEMS

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path

    from anndata._types import AnnDataElem

pytestmark = pytest.mark.skipif(not find_spec("xarray"), reason="xarray not installed")


@pytest.mark.parametrize(
    ("elem_key", "sub_key"),
    [
        ("raw", "X"),
        ("obs", "cat"),
        ("obs", "int64"),
        *((elem_name, None) for elem_name in ANNDATA_ELEMS),
    ],
)
def test_access_count_elem_access(
    remote_store_tall_skinny: AccessTrackingStore,
    adata_remote_tall_skinny: AnnData,
    elem_key: AnnDataElem,
    sub_key: str,
    simple_subset_func: Callable[[AnnData], AnnData],
):
    full_path = f"{elem_key}/{sub_key}" if sub_key is not None else elem_key
    remote_store_tall_skinny.initialize_key_trackers({full_path, "X"})
    # a series of methods that should __not__ read in any data
    elem = getattr(simple_subset_func(adata_remote_tall_skinny), elem_key)
    if sub_key is not None:
        getattr(elem, sub_key)
    remote_store_tall_skinny.assert_access_count(full_path, 0)
    remote_store_tall_skinny.assert_access_count("X", 0)


def test_access_count_subset(
    remote_store_tall_skinny: AccessTrackingStore,
    adata_remote_tall_skinny: AnnData,
):
    non_obs_elem_names = filter(lambda e: e != "obs", ANNDATA_ELEMS)
    remote_store_tall_skinny.initialize_key_trackers(
        ["obs/cat/codes", *non_obs_elem_names]
    )
    adata_remote_tall_skinny[adata_remote_tall_skinny.obs["cat"] == "a", :]
    # all codes read in for subset (from 1 chunk)
    remote_store_tall_skinny.assert_access_count("obs/cat/codes", 1)
    for elem_name in non_obs_elem_names:
        remote_store_tall_skinny.assert_access_count(elem_name, 0)


def test_access_count_subset_column_compute(
    remote_store_tall_skinny: AccessTrackingStore,
    adata_remote_tall_skinny: AnnData,
):
    remote_store_tall_skinny.initialize_key_trackers(["obs/int64"])
    adata_remote_tall_skinny[adata_remote_tall_skinny.shape[0] // 2, :].obs[
        "int64"
    ].compute()
    # two chunks needed for 0:10 subset
    remote_store_tall_skinny.assert_access_count("obs/int64", 1)


def test_access_count_index(
    remote_store_tall_skinny: AccessTrackingStore,
):
    remote_store_tall_skinny.initialize_key_trackers(["obs/_index"])
    read_lazy(remote_store_tall_skinny, load_annotation_index=False)
    remote_store_tall_skinny.assert_access_count("obs/_index", 0)
    read_lazy(remote_store_tall_skinny)
    # 4 is number of chunks
    remote_store_tall_skinny.assert_access_count("obs/_index", 4)


def test_access_count_dtype(
    remote_store_tall_skinny: AccessTrackingStore,
    adata_remote_tall_skinny: AnnData,
):
    remote_store_tall_skinny.initialize_key_trackers(["obs/cat/categories"])
    remote_store_tall_skinny.assert_access_count("obs/cat/categories", 0)
    # This should only cause categories to be read in once
    adata_remote_tall_skinny.obs["cat"].dtype  # noqa: B018
    adata_remote_tall_skinny.obs["cat"].dtype  # noqa: B018
    adata_remote_tall_skinny.obs["cat"].dtype  # noqa: B018
    remote_store_tall_skinny.assert_access_count("obs/cat/categories", 1)


def test_uns_uses_dask(adata_remote: AnnData):
    assert isinstance(adata_remote.uns["nested"]["nested_further"]["array"], DaskArray)


def test_to_memory(adata_remote: AnnData, adata_orig: AnnData):
    remote_to_memory = adata_remote.to_memory()
    assert_equal(remote_to_memory, adata_orig)


def test_access_counts_obsm_df(tmp_path: Path):
    adata = AnnData(
        X=np.array(np.random.rand(100, 20)),
    )
    adata.obsm["df"] = pd.DataFrame(
        {"col1": np.random.rand(100), "col2": np.random.rand(100)},
        index=adata.obs_names,
    )
    adata.write_zarr(tmp_path)
    store = AccessTrackingStore(tmp_path)
    store.initialize_key_trackers(["obsm/df"])
    read_lazy(store, load_annotation_index=False)
    store.assert_access_count("obsm/df", 0)


def test_view_to_memory(adata_remote: AnnData, adata_orig: AnnData):
    obs_cats = adata_orig.obs["obs_cat"].cat.categories
    subset_obs = adata_orig.obs["obs_cat"] == obs_cats[0]
    assert_equal(adata_orig[subset_obs, :], adata_remote[subset_obs, :].to_memory())

    var_cats = adata_orig.var["var_cat"].cat.categories
    subset_var = adata_orig.var["var_cat"] == var_cats[0]
    assert_equal(adata_orig[:, subset_var], adata_remote[:, subset_var].to_memory())


def test_view_of_view_to_memory(adata_remote: AnnData, adata_orig: AnnData):
    cats_obs = adata_orig.obs["obs_cat"].cat.categories
    subset_obs = (adata_orig.obs["obs_cat"] == cats_obs[0]) | (
        adata_orig.obs["obs_cat"] == cats_obs[1]
    )
    subsetted_adata = adata_orig[subset_obs, :]
    subset_subset_obs = subsetted_adata.obs["obs_cat"] == cats_obs[1]
    subsetted_subsetted_adata = subsetted_adata[subset_subset_obs, :]
    assert_equal(
        subsetted_subsetted_adata,
        adata_remote[subset_obs, :][subset_subset_obs, :].to_memory(),
    )

    cats_var = adata_orig.var["var_cat"].cat.categories
    subset_var = (adata_orig.var["var_cat"] == cats_var[0]) | (
        adata_orig.var["var_cat"] == cats_var[1]
    )
    subsetted_adata = adata_orig[:, subset_var]
    subset_subset_var = subsetted_adata.var["var_cat"] == cats_var[1]
    subsetted_subsetted_adata = subsetted_adata[:, subset_subset_var]
    assert_equal(
        subsetted_subsetted_adata,
        adata_remote[:, subset_var][:, subset_subset_var].to_memory(),
    )


def test_unconsolidated(tmp_path: Path, mtx_format):
    adata = gen_adata((1000, 1000), mtx_format, **GEN_ADATA_NO_XARRAY_ARGS)
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
