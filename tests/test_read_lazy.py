from __future__ import annotations

import typing
from contextlib import nullcontext
from importlib.util import find_spec
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import anndata as ad
from anndata import AnnData
from anndata._types import ANNDATA_ELEMS
from anndata.compat import DaskArray
from anndata.experimental import read_lazy
from anndata.tests.helpers import (
    AccessTrackingStore,
    as_dense_dask_array,
    assert_equal,
    gen_adata,
    gen_typed_df,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Generator
    from pathlib import Path
    from typing import Literal

    from numpy.typing import NDArray
pytestmark = pytest.mark.skipif(
    not find_spec("xarray"), reason="Xarray is not installed"
)

ANNDATA_ELEMS_LIST = typing.get_args(ANNDATA_ELEMS)


@pytest.fixture(
    params=[sparse.csr_matrix, sparse.csc_matrix, np.array, as_dense_dask_array],
    ids=["scipy-csr", "scipy-csc", "np-array", "dask_array"],
    scope="session",
)
def mtx_format(request):
    return request.param


@pytest.fixture(
    params=[True, False], ids=["vars_different", "vars_same"], scope="session"
)
def are_vars_different(request):
    return request.param


@pytest.fixture(params=["zarr", "h5ad"], scope="session")
def diskfmt(request):
    return request.param


@pytest.fixture(params=[True, False], scope="session")
def load_annotation_index(request):
    return request.param


@pytest.fixture(params=["outer", "inner"], scope="session")
def join(request):
    return request.param


@pytest.fixture(
    params=[
        pytest.param(lambda x: x, id="full"),
        pytest.param(lambda x: x[0:10, :], id="subset"),
    ],
    scope="session",
)
def simple_subset_func(request):
    return request.param


# TODO: why does `read_lazy().to_memory()` cause `Dataset2D.to_memory()` to lose index name in
# multi-threaded tests when only opened once i.e., without this Callable?
@pytest.fixture(scope="session")
def adata_remote_orig_with_path(
    tmp_path_factory,
    diskfmt: str,
    mtx_format,
    worker_id: str = "serial",
) -> tuple[AnnData, AnnData]:
    """Create remote fixtures, one without a range index and the other with"""
    file_name = f"orig_{worker_id}.{diskfmt}"
    if diskfmt == "h5ad":
        orig_path = tmp_path_factory.mktemp("h5ad_file_dir") / file_name
    else:
        orig_path = tmp_path_factory.mktemp(file_name)
    orig = gen_adata((1000, 1100), mtx_format)
    orig.raw = orig.copy()
    getattr(orig, f"write_{diskfmt}")(orig_path)
    return orig_path, orig


@pytest.fixture
def adata_remote_orig(
    adata_remote_orig_with_path: tuple[Path, AnnData], load_annotation_index: bool
) -> tuple[AnnData, AnnData]:
    orig_path, orig = adata_remote_orig_with_path
    return read_lazy(orig_path, load_annotation_index=load_annotation_index), orig


@pytest.fixture(scope="session")
def adata_remote_with_store_tall_skinny_path(
    tmp_path_factory,
    mtx_format,
    worker_id: str = "serial",
) -> Path:
    orig_path = tmp_path_factory.mktemp(f"orig_{worker_id}.zarr")
    M = 100_000  # forces zarr to chunk `obs` columns multiple ways - that way 1 access to `int64` below is actually only one access
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
    orig.raw = orig.copy()
    orig.write_zarr(orig_path)
    return orig_path


@pytest.fixture(scope="session")
def adatas_paths_var_indices_for_concatenation(
    tmp_path_factory, are_vars_different: bool, worker_id: str = "serial"
) -> tuple[list[AnnData], list[Path], list[pd.Index]]:
    adatas = []
    var_indices = []
    paths = []
    M = 1000
    N = 50
    n_datasets = 3
    for dataset_index in range(n_datasets):
        orig_path = tmp_path_factory.mktemp(f"orig_{worker_id}_{dataset_index}.zarr")
        paths.append(orig_path)
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
            X=np.random.binomial(100, 0.005, (M, N)).astype(np.float32),
        )
        orig.write_zarr(orig_path)
        adatas.append(orig)
    return adatas, paths, var_indices


@pytest.fixture
def concatenation_objects(
    adatas_paths_var_indices_for_concatenation,
) -> tuple[list[AnnData], list[pd.Index], list[AccessTrackingStore], list[AnnData]]:
    adatas, paths, var_indices = adatas_paths_var_indices_for_concatenation
    stores = [AccessTrackingStore(path) for path in paths]
    lazys = [read_lazy(store) for store in stores]
    return adatas, var_indices, stores, lazys


@pytest.fixture
def adata_remote_with_store_tall_skinny(
    adata_remote_with_store_tall_skinny_path: Path,
) -> tuple[AnnData, AccessTrackingStore]:
    store = AccessTrackingStore(adata_remote_with_store_tall_skinny_path)
    remote = read_lazy(store)
    return remote, store


def get_key_trackers_for_columns_on_axis(
    adata: AnnData, axis: Literal["obs", "var"]
) -> Generator[str, None, None]:
    """Generate keys for tracking, using `codes` from categorical columns instead of the column name

    Parameters
    ----------
    adata
        Object to get keys from
    axis
        Axis to get keys from

    Yields
    ------
        Keys for tracking
    """
    for col in getattr(adata, axis).columns:
        yield f"{axis}/{col}" if "cat" not in col else f"{axis}/{col}/codes"


@pytest.mark.parametrize(
    ("elem_key", "sub_key"),
    [
        ("raw", "X"),
        ("obs", "cat"),
        ("obs", "int64"),
        *((elem_name, None) for elem_name in ANNDATA_ELEMS_LIST),
    ],
)
def test_access_count_elem_access(
    adata_remote_with_store_tall_skinny: tuple[AnnData, AccessTrackingStore],
    elem_key: ANNDATA_ELEMS,
    sub_key: str,
    simple_subset_func: Callable[[AnnData], AnnData],
):
    remote, store = adata_remote_with_store_tall_skinny
    full_path = f"{elem_key}/{sub_key}" if sub_key is not None else elem_key
    store.initialize_key_trackers({full_path, "X"})
    # a series of methods that should __not__ read in any data
    elem = getattr(simple_subset_func(remote), elem_key)
    if sub_key is not None:
        getattr(elem, sub_key)
    store.assert_access_count(full_path, 0)
    store.assert_access_count("X", 0)


def test_access_count_subset(
    adata_remote_with_store_tall_skinny: tuple[AnnData, AccessTrackingStore],
):
    remote, store = adata_remote_with_store_tall_skinny
    non_obs_elem_names = filter(lambda e: e != "obs", ANNDATA_ELEMS_LIST)
    store.initialize_key_trackers(["obs/cat/codes", *non_obs_elem_names])
    remote[remote.obs["cat"] == "a", :]
    # all codes read in for subset (from 1 chunk)
    store.assert_access_count("obs/cat/codes", 1)
    for elem_name in non_obs_elem_names:
        store.assert_access_count(elem_name, 0)


def test_access_count_subset_column_compute(
    adata_remote_with_store_tall_skinny: tuple[AnnData, AccessTrackingStore],
):
    remote, store = adata_remote_with_store_tall_skinny
    store.initialize_key_trackers(["obs/int64"])
    remote[remote.shape[0] // 2, :].obs["int64"].compute()
    # two chunks needed for 0:10 subset
    store.assert_access_count("obs/int64", 1)


def test_access_count_index(
    adata_remote_with_store_tall_skinny: tuple[AnnData, AccessTrackingStore],
):
    _, store = adata_remote_with_store_tall_skinny
    store.initialize_key_trackers(["obs/_index"])
    read_lazy(store, load_annotation_index=False)
    store.assert_access_count("obs/_index", 0)
    read_lazy(store)
    # 4 is number of chunks
    store.assert_access_count("obs/_index", 4)


def test_access_count_dtype(
    adata_remote_with_store_tall_skinny: tuple[AnnData, AccessTrackingStore],
):
    remote, store = adata_remote_with_store_tall_skinny
    store.initialize_key_trackers(["obs/cat/categories"])
    store.assert_access_count("obs/cat/categories", 0)
    # This should only cause categories to be read in once
    remote.obs["cat"].dtype
    remote.obs["cat"].dtype
    remote.obs["cat"].dtype
    store.assert_access_count("obs/cat/categories", 1)


def test_uns_uses_dask(adata_remote_orig: tuple[AnnData, AnnData]):
    remote, _ = adata_remote_orig
    assert isinstance(remote.uns["nested"]["nested_further"]["array"], DaskArray)


def test_to_memory(adata_remote_orig: tuple[AnnData, AnnData]):
    remote, orig = adata_remote_orig
    remote_to_memory = remote.to_memory()
    assert_equal(remote_to_memory, orig)


def test_view_to_memory(adata_remote_orig: tuple[AnnData, AnnData]):
    remote, orig = adata_remote_orig
    subset_obs = orig.obs["obs_cat"] == "a"
    assert_equal(orig[subset_obs, :], remote[subset_obs, :].to_memory())

    subset_var = orig.var["var_cat"] == "a"
    assert_equal(orig[:, subset_var], remote[:, subset_var].to_memory())


def test_view_of_view_to_memory(adata_remote_orig: tuple[AnnData, AnnData]):
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


def test_unconsolidated(tmp_path: Path, mtx_format):
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


# remote has object dtype, need to convert back for integers booleans etc.
def correct_extension_dtype_differences(remote: pd.DataFrame, memory: pd.DataFrame):
    for col in memory.columns:
        dtype = memory[col].dtype
        if pd.api.types.is_extension_array_dtype(dtype):
            remote[col] = remote[col].astype(dtype)
    return remote, memory


ANNDATA_ELEMS_LIST = typing.get_args(ANNDATA_ELEMS)


@pytest.mark.parametrize("join", ["outer", "inner"])
@pytest.mark.parametrize(
    ("elem_key", "sub_key"),
    [
        ("obs", "cat"),
        ("obs", "int64"),
        *((elem_name, None) for elem_name in ANNDATA_ELEMS_LIST),
    ],
)
def test_concat_access_count(
    concatenation_objects: tuple[
        list[AnnData], list[pd.Index], list[AccessTrackingStore], list[AnnData]
    ],
    join: Literal["outer", "inner"],
    elem_key: ANNDATA_ELEMS,
    sub_key: str,
    simple_subset_func: Callable[[AnnData], AnnData],
):
    adatas, _, stores, lazy_adatas = concatenation_objects
    # track all elems except codes because they must be read in for concatenation
    non_categorical_columns = (
        f"{elem}/{col}" if "cat" not in col else f"{elem}/{col}/codes"
        for elem in ["obs", "var"]
        for col in adatas[0].obs.columns
    )
    non_obs_var_keys = filter(lambda e: e not in {"obs", "var"}, ANNDATA_ELEMS_LIST)
    keys_to_track = [*non_categorical_columns, *non_obs_var_keys]
    for store in stores:
        store.initialize_key_trackers(keys_to_track)
    concated_remote = ad.concat(lazy_adatas, join=join)
    # a series of methods that should __not__ read in any data
    elem = getattr(simple_subset_func(concated_remote), elem_key)
    if sub_key is not None:
        getattr(elem, sub_key)
    for store in stores:
        for elem in keys_to_track:
            store.assert_access_count(elem, 0)


def test_concat_to_memory_obs_access_count(
    concatenation_objects: tuple[
        list[AnnData], list[pd.Index], list[AccessTrackingStore], list[AnnData]
    ],
    join: Literal["outer", "inner"],
    simple_subset_func: Callable[[AnnData], AnnData],
):
    adatas, _, stores, lazy_adatas = concatenation_objects
    concated_remote = simple_subset_func(ad.concat(lazy_adatas, join=join))
    concated_remote_subset = simple_subset_func(concated_remote)
    n_datasets = len(adatas)
    obs_keys_to_track = get_key_trackers_for_columns_on_axis(adatas[0], "obs")
    for store in stores:
        store.initialize_key_trackers(obs_keys_to_track)
    concated_remote_subset.to_memory()
    # check access count for the stores - only the first should be accessed when reading into memory
    for col in obs_keys_to_track:
        stores[0].assert_access_count(col, 1)
        for i in range(1, n_datasets):
            # if the shapes are the same, data was read in to bring the object into memory; otherwise, not
            stores[i].assert_access_count(
                col, concated_remote_subset.shape[0] == concated_remote.shape[0]
            )


def test_concat_to_memory_obs(
    concatenation_objects: tuple[
        list[AnnData], list[pd.Index], list[AccessTrackingStore], list[AnnData]
    ],
    join: Literal["outer", "inner"],
    simple_subset_func: Callable[[AnnData], AnnData],
):
    adatas, _, _, lazy_adatas = concatenation_objects
    concatenated_memory = simple_subset_func(ad.concat(adatas, join=join))
    concated_remote = simple_subset_func(ad.concat(lazy_adatas, join=join))
    # TODO: name is lost normally, should fix
    obs_memory = concatenated_memory.obs
    obs_memory.index.name = "obs_names"
    assert_equal(
        *correct_extension_dtype_differences(
            concated_remote.obs.to_pandas(), concatenated_memory.obs
        )
    )


def test_concat_to_memory_obs_dtypes(
    concatenation_objects: tuple[
        list[AnnData], list[pd.Index], list[AccessTrackingStore], list[AnnData]
    ],
    join: Literal["outer", "inner"],
):
    _, _, _, lazy_adatas = concatenation_objects
    concated_remote = ad.concat(lazy_adatas, join=join)
    # check preservation of non-categorical dtypes on the concat axis
    assert concated_remote.obs["int64"].dtype == "int64"
    assert concated_remote.obs["uint8"].dtype == "uint8"
    assert concated_remote.obs["nullable-int"].dtype == "int32"
    assert concated_remote.obs["float64"].dtype == "float64"
    assert concated_remote.obs["bool"].dtype == "bool"
    assert concated_remote.obs["nullable-bool"].dtype == "bool"


def test_concat_to_memory_var(
    concatenation_objects: tuple[
        list[AnnData], list[pd.Index], list[AccessTrackingStore], list[AnnData]
    ],
    join: Literal["outer", "inner"],
    are_vars_different: bool,
    simple_subset_func: Callable[[AnnData], AnnData],
):
    adatas, var_indices, stores, lazy_adatas = concatenation_objects
    concated_remote = simple_subset_func(ad.concat(lazy_adatas, join=join))
    var_keys_to_track = get_key_trackers_for_columns_on_axis(adatas[0], "var")
    for store in stores:
        store.initialize_key_trackers(var_keys_to_track)
    # check non-different variables, taken from first annotation.
    pd_index_overlapping = pd.Index(
        filter(lambda x: not x.endswith("ds"), var_indices[0])
    )
    var_df_overlapping = adatas[0][:, pd_index_overlapping].var.copy()
    test_cases = [(pd_index_overlapping, var_df_overlapping, 0)]
    if are_vars_different and join == "outer":
        # check a set of unique variables from the first object since we only take from there if different
        pd_index_only_ds_0 = pd.Index(filter(lambda x: "0_ds" in x, var_indices[1]))
        var_df_only_ds_0 = adatas[0][:, pd_index_only_ds_0].var.copy()
        test_cases.append((pd_index_only_ds_0, var_df_only_ds_0, 0))
    for pd_index, var_df, store_idx in test_cases:
        var_df.index.name = "var_names"
        remote_df = concated_remote[:, pd_index].var.to_pandas()
        remote_df_corrected, _ = correct_extension_dtype_differences(remote_df, var_df)
        #  TODO:xr.merge always upcasts to float due to NA and you can't downcast?
        for col in remote_df_corrected.columns:
            dtype = remote_df_corrected[col].dtype
            if dtype in [np.float64, np.float32]:
                var_df[col] = var_df[col].astype(dtype)
        assert_equal(remote_df_corrected, var_df)
        for key in var_keys_to_track:
            stores[store_idx].assert_access_count(key, 1)
            for store in stores:
                if store != stores[store_idx]:
                    store.assert_access_count(key, 0)
        stores[store_idx].reset_key_trackers()


@pytest.mark.parametrize(
    "index",
    [
        pytest.param(
            slice(500, 1500),
            id="slice",
        ),
        pytest.param(
            np.arange(950, 1050),
            id="consecutive integer array",
        ),
        pytest.param(
            np.random.choice(np.arange(800, 1100), 500),
            id="random integer array",
        ),
        pytest.param(
            np.random.choice([True, False], 2000),
            id="boolean array",
        ),
        pytest.param(slice(None), id="full slice"),
        pytest.param("a", id="categorical_subset"),
        pytest.param(None, id="No index"),
    ],
)
def test_concat_full_and_subsets(
    adata_remote_orig: tuple[AnnData, AnnData],
    join,
    index: slice | NDArray | Literal["a"] | None,
    load_annotation_index: bool,
):
    from anndata.experimental.backed._compat import Dataset2D

    remote, orig = adata_remote_orig

    maybe_warning_context = (
        pytest.warns(UserWarning, match=r"Concatenating with a pandas numeric")
        if not load_annotation_index
        else nullcontext()
    )
    with maybe_warning_context:
        remote_concatenated = ad.concat([remote, remote], join=join)
    if index is not None:
        if np.isscalar(index) and index == "a":
            index = remote_concatenated.obs["obs_cat"] == "a"
        remote_concatenated = remote_concatenated[index]
    orig_concatenated = ad.concat([orig, orig], join=join)
    if index is not None:
        orig_concatenated = orig_concatenated[index]
    in_memory_remote_concatenated = remote_concatenated.to_memory()
    corrected_remote_obs, corrected_memory_obs = correct_extension_dtype_differences(
        in_memory_remote_concatenated.obs, orig_concatenated.obs
    )
    assert isinstance(remote_concatenated.obs, Dataset2D)
    assert_equal(corrected_remote_obs, corrected_memory_obs)
    assert_equal(in_memory_remote_concatenated.X, orig_concatenated.X)
    assert (
        in_memory_remote_concatenated.var_names.tolist()
        == orig_concatenated.var_names.tolist()
    )
