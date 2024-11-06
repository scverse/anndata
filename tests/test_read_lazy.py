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
from anndata._core.file_backing import to_memory
from anndata._types import AnnDataElem
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

    from anndata._types import Join_T

pytestmark = pytest.mark.skipif(
    not find_spec("xarray"), reason="Xarray is not installed"
)

ANNDATA_ELEMS = typing.get_args(AnnDataElem)


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


@pytest.fixture(
    params=[True, False],
    scope="session",
    ids=["load-annotation-index", "dont-load-annotation-index"],
)
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


@pytest.fixture(scope="session")
def adata_remote_orig_with_path(
    tmp_path_factory,
    diskfmt: str,
    mtx_format,
    worker_id: str = "serial",
) -> tuple[Path, AnnData]:
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
def adata_remote(
    adata_remote_orig_with_path: tuple[Path, AnnData], load_annotation_index: bool
) -> AnnData:
    orig_path, _ = adata_remote_orig_with_path
    return read_lazy(orig_path, load_annotation_index=load_annotation_index)


@pytest.fixture
def adata_orig(adata_remote_orig_with_path: tuple[Path, AnnData]) -> AnnData:
    _, orig = adata_remote_orig_with_path
    return orig


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
def var_indices_for_concat(
    adatas_paths_var_indices_for_concatenation,
) -> list[pd.Index]:
    _, _, var_indices = adatas_paths_var_indices_for_concatenation
    return var_indices


@pytest.fixture
def adatas_for_concat(
    adatas_paths_var_indices_for_concatenation,
) -> list[AnnData]:
    adatas, _, _ = adatas_paths_var_indices_for_concatenation
    return adatas


@pytest.fixture
def stores_for_concat(
    adatas_paths_var_indices_for_concatenation,
) -> list[AccessTrackingStore]:
    _, paths, _ = adatas_paths_var_indices_for_concatenation
    return [AccessTrackingStore(path) for path in paths]


@pytest.fixture
def lazy_adatas_for_concat(
    stores_for_concat,
) -> list[AnnData]:
    return [read_lazy(store) for store in stores_for_concat]


@pytest.fixture
def adata_remote_with_store_tall_skinny(
    adata_remote_with_store_tall_skinny_path: Path,
) -> tuple[AnnData, AccessTrackingStore]:
    store = AccessTrackingStore(adata_remote_with_store_tall_skinny_path)
    remote = read_lazy(store)
    return remote, store


@pytest.fixture
def remote_store_tall_skinny(
    adata_remote_with_store_tall_skinny_path: Path,
) -> AccessTrackingStore:
    return AccessTrackingStore(adata_remote_with_store_tall_skinny_path)


@pytest.fixture
def adata_remote_tall_skinny(
    remote_store_tall_skinny: AccessTrackingStore,
) -> AnnData:
    remote = read_lazy(remote_store_tall_skinny)
    return remote


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
    adata_remote_tall_skinny.obs["cat"].dtype
    adata_remote_tall_skinny.obs["cat"].dtype
    adata_remote_tall_skinny.obs["cat"].dtype
    remote_store_tall_skinny.assert_access_count("obs/cat/categories", 1)


def test_uns_uses_dask(adata_remote: AnnData):
    assert isinstance(adata_remote.uns["nested"]["nested_further"]["array"], DaskArray)


def test_to_memory(adata_remote: AnnData, adata_orig: AnnData):
    remote_to_memory = adata_remote.to_memory()
    assert_equal(remote_to_memory, adata_orig)


def test_view_to_memory(adata_remote: AnnData, adata_orig: AnnData):
    subset_obs = adata_orig.obs["obs_cat"] == "a"
    assert_equal(adata_orig[subset_obs, :], adata_remote[subset_obs, :].to_memory())

    subset_var = adata_orig.var["var_cat"] == "a"
    assert_equal(adata_orig[:, subset_var], adata_remote[:, subset_var].to_memory())


def test_view_of_view_to_memory(adata_remote: AnnData, adata_orig: AnnData):
    subset_obs = (adata_orig.obs["obs_cat"] == "a") | (adata_orig.obs["obs_cat"] == "b")
    subsetted_adata = adata_orig[subset_obs, :]
    subset_subset_obs = subsetted_adata.obs["obs_cat"] == "b"
    subsetted_subsetted_adata = subsetted_adata[subset_subset_obs, :]
    assert_equal(
        subsetted_subsetted_adata,
        adata_remote[subset_obs, :][subset_subset_obs, :].to_memory(),
    )

    subset_var = (adata_orig.var["var_cat"] == "a") | (adata_orig.var["var_cat"] == "b")
    subsetted_adata = adata_orig[:, subset_var]
    subset_subset_var = subsetted_adata.var["var_cat"] == "b"
    subsetted_subsetted_adata = subsetted_adata[:, subset_subset_var]
    assert_equal(
        subsetted_subsetted_adata,
        adata_remote[:, subset_var][:, subset_subset_var].to_memory(),
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


ANNDATA_ELEMS = typing.get_args(AnnDataElem)


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
    # track all elems except codes because they must be read in for concatenation
    non_categorical_columns = (
        f"{elem}/{col}" if "cat" not in col else f"{elem}/{col}/codes"
        for elem in ["obs", "var"]
        for col in adatas_for_concat[0].obs.columns
    )
    non_obs_var_keys = filter(lambda e: e not in {"obs", "var"}, ANNDATA_ELEMS)
    keys_to_track = [*non_categorical_columns, *non_obs_var_keys]
    for store in stores_for_concat:
        store.initialize_key_trackers(keys_to_track)
    concated_remote = ad.concat(lazy_adatas_for_concat, join=join)
    # a series of methods that should __not__ read in any data
    elem = getattr(simple_subset_func(concated_remote), elem_key)
    if sub_key is not None:
        getattr(elem, sub_key)
    for store in stores_for_concat:
        for elem in keys_to_track:
            store.assert_access_count(elem, 0)


def test_concat_to_memory_obs_access_count(
    adatas_for_concat: list[AnnData],
    stores_for_concat: list[AccessTrackingStore],
    lazy_adatas_for_concat: list[AnnData],
    join: Join_T,
    simple_subset_func: Callable[[AnnData], AnnData],
):
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
    # TODO: name is lost normally, should fix
    obs_memory = concatenated_memory.obs
    obs_memory.index.name = "obs_names"
    assert_equal(
        *unify_extension_dtypes(
            concated_remote.obs.to_pandas(), concatenated_memory.obs
        )
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
    are_vars_different: bool,
    simple_subset_func: Callable[[AnnData], AnnData],
):
    concated_remote = simple_subset_func(ad.concat(lazy_adatas_for_concat, join=join))
    var_keys_to_track = get_key_trackers_for_columns_on_axis(
        adatas_for_concat[0], "var"
    )
    for store in stores_for_concat:
        store.initialize_key_trackers(var_keys_to_track)
    # check non-different variables, taken from first annotation.
    pd_index_overlapping = pd.Index(
        filter(lambda x: not x.endswith("ds"), var_indices_for_concat[0])
    )
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
        var_df.index.name = "var_names"
        remote_df = concated_remote[:, pd_index].var.to_pandas()
        remote_df_corrected, _ = unify_extension_dtypes(remote_df, var_df)
        #  TODO:xr.merge always upcasts to float due to NA and you can't downcast?
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


def test_concat_data_with_cluster_to_memory(
    adata_remote: AnnData,
    join: Join_T,
):
    import dask.distributed as dd

    with (
        dd.LocalCluster(n_workers=1, threads_per_worker=1) as cluster,
        dd.Client(cluster),
    ):
        ad.concat([adata_remote, adata_remote], join=join).to_memory()


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
            np.random.randint(800, 1100, 500),
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
def test_concat_data(
    adata_remote: AnnData,
    adata_orig: AnnData,
    join: Join_T,
    index: slice | NDArray | Literal["a"] | None,
    load_annotation_index: bool,
):
    from anndata.experimental.backed._compat import Dataset2D

    with (
        pytest.warns(UserWarning, match=r"Concatenating with a pandas numeric")
        if not load_annotation_index
        else nullcontext()
    ):
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
    load_annotation_index: bool,
    join: Join_T,
    attr: str,
    key: str | None,
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


def test_concat_bad_mixed_types(tmp_path: str):
    orig = gen_adata((100, 200), np.array)
    orig.write_zarr(tmp_path)
    remote = read_lazy(tmp_path)
    orig.obsm["df"] = orig.obsm["array"]
    with pytest.raises(ValueError, match=r"Cannot concatenate a Dataset2D*"):
        ad.concat([remote, orig], join="outer")
