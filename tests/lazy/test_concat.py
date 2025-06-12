from __future__ import annotations

from functools import reduce
from importlib.util import find_spec
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest

import anndata as ad
from anndata._core.file_backing import to_memory
from anndata.experimental import read_lazy
from anndata.tests.helpers import GEN_ADATA_NO_XARRAY_ARGS, assert_equal, gen_adata

from .conftest import ANNDATA_ELEMS, get_key_trackers_for_columns_on_axis

pytestmark = pytest.mark.skipif(not find_spec("xarray"), reason="xarray not installed")

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path
    from typing import Literal

    from numpy.typing import NDArray

    from anndata import AnnData
    from anndata._types import AnnDataElem, Join_T
    from anndata.tests.helpers import AccessTrackingStore


def unify_extension_dtypes(
    remote: pd.DataFrame, memory: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    For concatenated lazy datasets, we send the extension arrays through dask
    But this means we lose the pandas dtype, so this function corrects that.

    Parameters
    ----------
    remote
        The dataset that comes from the concatenated lazy operation
    memory
        The in-memory, "correct" version

    Returns
    -------
    The two dataframes unified
    """
    for col in memory.columns:
        dtype = memory[col].dtype
        if pd.api.types.is_extension_array_dtype(dtype):
            remote[col] = remote[col].astype(dtype)
    return remote, memory


@pytest.mark.parametrize("join", ["outer", "inner"])
@pytest.mark.parametrize(
    ("elem_key", "sub_key"),
    [
        ("obs", "cat"),
        ("obs", "int64"),
        *((elem_name, None) for elem_name in ANNDATA_ELEMS),
    ],
)
def test_concat_access_count(
    adatas_for_concat: list[AnnData],
    stores_for_concat: list[AccessTrackingStore],
    lazy_adatas_for_concat: list[AnnData],
    join: Join_T,
    elem_key: AnnDataElem,
    sub_key: str,
    simple_subset_func: Callable[[AnnData], AnnData],
):
    # track all elems except categories from categoricals because they must be read in for concatenation
    # due to the dtype check on the elements (which causes `categories` to be read in)
    non_categorical_columns = (
        f"{elem}/{col}" if "cat" not in col else f"{elem}/{col}/codes"
        for elem in ["obs", "var"]
        for col in adatas_for_concat[0].obs.columns
    )
    category_columns = (
        f"{elem}/{col}/categories"
        for elem in ["obs", "var"]
        for col in adatas_for_concat[0].obs.columns
        if "cat" in col
    )
    non_obs_var_keys = filter(lambda e: e not in {"obs", "var"}, ANNDATA_ELEMS)
    zero_access_count_keys = [*non_categorical_columns, *non_obs_var_keys]
    keys_to_track = [*zero_access_count_keys, *category_columns]
    for store in stores_for_concat:
        store.initialize_key_trackers(keys_to_track)
    concated_remote = ad.concat(lazy_adatas_for_concat, join=join)
    # a series of methods that should __not__ read in any data
    elem = getattr(simple_subset_func(concated_remote), elem_key)
    if sub_key is not None:
        getattr(elem, sub_key)
    for store in stores_for_concat:
        for elem in zero_access_count_keys:
            store.assert_access_count(elem, 0)
        for elem in category_columns:
            # once for .zarray, once for the actual data
            store.assert_access_count(elem, 2)


def test_concat_to_memory_obs_access_count(
    adatas_for_concat: list[AnnData],
    stores_for_concat: list[AccessTrackingStore],
    lazy_adatas_for_concat: list[AnnData],
    join: Join_T,
    simple_subset_func: Callable[[AnnData], AnnData],
):
    """This test ensures that only the necessary chunks are accessed in `to_memory` call after a subsetting operation"""
    concated_remote = simple_subset_func(ad.concat(lazy_adatas_for_concat, join=join))
    concated_remote_subset = simple_subset_func(concated_remote)
    n_datasets = len(adatas_for_concat)
    obs_keys_to_track = get_key_trackers_for_columns_on_axis(
        adatas_for_concat[0], "obs"
    )
    for store in stores_for_concat:
        store.initialize_key_trackers(obs_keys_to_track)
    concated_remote_subset.to_memory()
    # check access count for the stores - only the first should be accessed when reading into memory
    for col in obs_keys_to_track:
        stores_for_concat[0].assert_access_count(col, 1)
        for i in range(1, n_datasets):
            # if the shapes are the same, data was read in to bring the object into memory; otherwise, not
            stores_for_concat[i].assert_access_count(
                col, concated_remote_subset.shape[0] == concated_remote.shape[0]
            )


def test_concat_to_memory_obs(
    adatas_for_concat: list[AnnData],
    lazy_adatas_for_concat: list[AnnData],
    join: Join_T,
    simple_subset_func: Callable[[AnnData], AnnData],
):
    concatenated_memory = simple_subset_func(ad.concat(adatas_for_concat, join=join))
    concated_remote = simple_subset_func(ad.concat(lazy_adatas_for_concat, join=join))
    assert_equal(
        *unify_extension_dtypes(to_memory(concated_remote.obs), concatenated_memory.obs)
    )


def test_concat_to_memory_obs_dtypes(
    lazy_adatas_for_concat: list[AnnData],
    join: Join_T,
):
    concated_remote = ad.concat(lazy_adatas_for_concat, join=join)
    # check preservation of non-categorical dtypes on the concat axis
    assert concated_remote.obs["int64"].dtype == "int64"
    assert concated_remote.obs["uint8"].dtype == "uint8"
    assert concated_remote.obs["nullable-int"].dtype == "int32"
    assert concated_remote.obs["float64"].dtype == "float64"
    assert concated_remote.obs["bool"].dtype == "bool"
    assert concated_remote.obs["nullable-bool"].dtype == "bool"


def test_concat_to_memory_var(
    var_indices_for_concat: list[pd.Index],
    adatas_for_concat: list[AnnData],
    stores_for_concat: list[AccessTrackingStore],
    lazy_adatas_for_concat: list[AnnData],
    join: Join_T,
    simple_subset_func: Callable[[AnnData], AnnData],
    *,
    are_vars_different: bool,
):
    """\
    The goal of this test to ensure that the various `join` operations work as expected
    under various scenarios.

    We test two things here: first, we take all the overlapping indices for var.
    Then if the underlying vars are different and this is an outer join
    (i.e., there are non-overlapping indices), we take the unique indices from one of the dataframes.
    We then check if the var dataframe subsetted from lazily-concatenated object and put into memory
    is the same as the underlying anndata object that created it, up to some corrections.

    We also test for key access counts to ensure that data was not taken from the var df of other
    on-disk anndata objects that might be different i.e., in the case of an outer join.
    """
    concated_remote = simple_subset_func(ad.concat(lazy_adatas_for_concat, join=join))
    var_keys_to_track = get_key_trackers_for_columns_on_axis(
        adatas_for_concat[0], "var"
    )
    for store in stores_for_concat:
        store.initialize_key_trackers(var_keys_to_track)
    # check non-different variables, taken from first annotation.
    pd_index_overlapping = reduce(pd.Index.intersection, var_indices_for_concat)
    var_df_overlapping = adatas_for_concat[0][:, pd_index_overlapping].var.copy()
    test_cases = [(pd_index_overlapping, var_df_overlapping, 0)]
    if are_vars_different and join == "outer":
        # check a set of unique variables from the first object since we only take from there if different
        pd_index_only_ds_0 = pd.Index(
            filter(lambda x: "0_ds" in x, var_indices_for_concat[1])
        )
        var_df_only_ds_0 = adatas_for_concat[0][:, pd_index_only_ds_0].var.copy()
        test_cases.append((pd_index_only_ds_0, var_df_only_ds_0, 0))
    for pd_index, var_df, store_idx in test_cases:
        remote_df = to_memory(concated_remote[:, pd_index].var)
        remote_df_corrected, _ = unify_extension_dtypes(remote_df, var_df)
        # NOTE: xr.merge always upcasts to float due to NA and you can't downcast?
        for col in remote_df_corrected.columns:
            dtype = remote_df_corrected[col].dtype
            if dtype in [np.float64, np.float32]:
                var_df[col] = var_df[col].astype(dtype)
        assert_equal(remote_df_corrected, var_df)
        for key in var_keys_to_track:
            stores_for_concat[store_idx].assert_access_count(key, 1)
            for store in stores_for_concat:
                if store != stores_for_concat[store_idx]:
                    store.assert_access_count(key, 0)
        stores_for_concat[store_idx].reset_key_trackers()


@pytest.mark.xdist_group("dask")
def test_concat_data_with_cluster_to_memory(
    adata_remote: AnnData, join: Join_T, local_cluster_addr: str
) -> None:
    import dask.distributed as dd

    with dd.Client(local_cluster_addr):
        ad.concat([adata_remote, adata_remote], join=join).to_memory()


@pytest.mark.parametrize(
    "index",
    [
        pytest.param(
            slice(50, 150),
            id="slice",
        ),
        pytest.param(
            np.arange(95, 105),
            id="consecutive integer array",
        ),
        pytest.param(
            np.random.randint(80, 110, 5),
            id="random integer array",
        ),
        pytest.param(
            np.random.choice([True, False], 200),
            id="boolean array",
        ),
        pytest.param(slice(None), id="full slice"),
        pytest.param("a", id="categorical_subset"),
        pytest.param(None, id="No index"),
    ],
)
def test_concat_data_subsetting(
    adata_remote: AnnData,
    adata_orig: AnnData,
    join: Join_T,
    index: slice | NDArray | Literal["a"] | None,
):
    from anndata._core.xarray import Dataset2D

    remote_concatenated = ad.concat([adata_remote, adata_remote], join=join)
    if index is not None:
        if np.isscalar(index) and index == "a":
            index = remote_concatenated.obs["obs_cat"] == "a"
        remote_concatenated = remote_concatenated[index]
    orig_concatenated = ad.concat([adata_orig, adata_orig], join=join)
    if index is not None:
        orig_concatenated = orig_concatenated[index]
    in_memory_remote_concatenated = remote_concatenated.to_memory()
    corrected_remote_obs, corrected_memory_obs = unify_extension_dtypes(
        in_memory_remote_concatenated.obs, orig_concatenated.obs
    )
    assert isinstance(remote_concatenated.obs, Dataset2D)
    assert_equal(corrected_remote_obs, corrected_memory_obs)
    assert_equal(in_memory_remote_concatenated.X, orig_concatenated.X)
    assert (
        in_memory_remote_concatenated.var_names.tolist()
        == orig_concatenated.var_names.tolist()
    )


@pytest.mark.parametrize(
    ("attr", "key"),
    (
        pytest.param(param[0], param[1], id="-".join(map(str, param)))
        for param in [("obs", None), ("var", None), ("obsm", "df"), ("varm", "df")]
    ),
)
def test_concat_df_ds_mixed_types(
    adata_remote: AnnData,
    adata_orig: AnnData,
    join: Join_T,
    attr: str,
    key: str | None,
    *,
    load_annotation_index: bool,
):
    def with_elem_in_memory(adata: AnnData, attr: str, key: str | None) -> AnnData:
        parent_elem = getattr(adata, attr)
        if key is not None:
            getattr(adata, attr)[key] = to_memory(parent_elem[key])
            return adata
        setattr(adata, attr, to_memory(parent_elem))
        return adata

    if not load_annotation_index:
        pytest.skip(
            "Testing for mixed types is independent of the axis since the indices always have to match."
        )
    remote = with_elem_in_memory(adata_remote, attr, key)
    in_memory_concatenated = ad.concat([adata_orig, adata_orig], join=join)
    mixed_concatenated = ad.concat([remote, adata_orig], join=join)
    assert_equal(mixed_concatenated, in_memory_concatenated)


def test_concat_bad_mixed_types(tmp_path: Path):
    orig = gen_adata((100, 200), np.array, **GEN_ADATA_NO_XARRAY_ARGS)
    orig.write_zarr(tmp_path)
    remote = read_lazy(tmp_path)
    orig.obsm["df"] = orig.obsm["array"]
    with pytest.raises(ValueError, match=r"Cannot concatenate a Dataset2D*"):
        ad.concat([remote, orig], join="outer")
