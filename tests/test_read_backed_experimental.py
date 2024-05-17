from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import zarr
from scipy import sparse

from anndata._core.anndata import AnnData
from anndata.experimental import read_backed
from anndata.experimental.backed._lazy_arrays import (
    CategoricalArray,
    MaskedArray,
)
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


@pytest.fixture()
def categorical_lazy_arr(tmp_path_factory):
    base_path = tmp_path_factory.getbasetemp()
    z = zarr.open_group(base_path, mode="w")
    z["codes"] = np.array([0, 1, 0, 1, 1, 2, 2, 1, 2, 0, 1, 1, 1, 2, 1, 2])
    z["categories"] = np.array(["foo", "bar", "jazz"])
    z.attrs["ordered"] = False
    z = zarr.open(base_path)
    return CategoricalArray(
        codes=z["codes"],
        categories=z["categories"],
        ordered=z.attrs["ordered"],
        drop_unused_categories=True,
    )


@pytest.fixture()
def nullable_boolean_lazy_arr(tmp_path_factory):
    base_path = tmp_path_factory.getbasetemp()
    z = zarr.open_group(base_path, mode="w")
    z["values"] = np.array(
        [
            True,
            False,
            True,
            False,
            False,
            True,
            False,
            False,
            True,
            True,
            False,
            False,
            False,
            True,
            False,
            True,
        ]
    )
    z["mask"] = np.array(
        [
            True,
            True,
            True,
            True,
            True,
            False,
            False,
            True,
            False,
            True,
            True,
            True,
            True,
            False,
            True,
            False,
        ]
    )
    z = zarr.open(base_path)
    return MaskedArray(values=z["values"], mask=z["mask"], dtype_str="nullable-boolean")


@pytest.fixture()
def nullable_boolean_lazy_arr_no_mask(tmp_path_factory):
    base_path = tmp_path_factory.getbasetemp()
    z = zarr.open_group(base_path, mode="w")
    z["values"] = np.array(
        [
            True,
            False,
            True,
            False,
            False,
            True,
            False,
            False,
            True,
            True,
            False,
            False,
            False,
            True,
            False,
            True,
        ]
    )
    z = zarr.open(base_path)
    return MaskedArray(values=z["values"], mask=None, dtype_str="nullable-boolean")


@pytest.fixture()
def nullable_integer_lazy_arr(tmp_path_factory):
    base_path = tmp_path_factory.getbasetemp()
    z = zarr.open_group(base_path, mode="w")
    z["values"] = np.array([0, 1, 0, 1, 1, 2, 2, 1, 2, 0, 1, 1, 1, 2, 1, 2])
    z["mask"] = np.array(
        [
            True,
            True,
            True,
            True,
            True,
            False,
            False,
            True,
            False,
            True,
            True,
            True,
            True,
            False,
            True,
            False,
        ]
    )
    z = zarr.open(base_path)
    return MaskedArray(values=z["values"], mask=z["mask"], dtype_str="nullable-integer")


@pytest.fixture()
def nullable_integer_lazy_arr_no_mask(tmp_path_factory):
    base_path = tmp_path_factory.getbasetemp()
    z = zarr.open_group(base_path, mode="w")
    z["values"] = np.array([0, 1, 0, 1, 1, 2, 2, 1, 2, 0, 1, 1, 1, 2, 1, 2])
    z = zarr.open(base_path)
    return MaskedArray(values=z["values"], mask=None, dtype_str="nullable-integer")


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
    store.set_key_trackers(["obs/int64", "var/int64", "obs/cat/codes", "X"])
    remote = read_backed(store)
    # a series of methods that should __not__ read in any data
    remote.X  # the initial (non-subset) access to `X` should not read in data
    remote.shape
    remote.var
    remote.obs
    remote.obs["int64"]
    remote.var["int64"]
    assert store.get_access_count("obs/cat/codes") == 0, store.get_subkeys_accessed(
        "obs/cat/codes"
    )
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


def test_to_memory(tmp_path, mtx_format, dskfmt):
    adata = gen_adata((1000, 1000), mtx_format)
    base_pth = Path(tmp_path)
    orig_pth = base_pth / f"orig.{dskfmt}"
    write = lambda x: getattr(x, f"write_{dskfmt}")(orig_pth)
    write(adata)
    remote = read_backed(orig_pth)
    remote_to_memory = remote.to_memory()
    assert_equal(remote_to_memory, adata)


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


def test_lazy_categorical_array_properties(categorical_lazy_arr):
    assert len(categorical_lazy_arr[0:3]) == 3
    assert type(categorical_lazy_arr[0:3]) == pd.Categorical
    assert len(categorical_lazy_arr[...]) == len(categorical_lazy_arr)
    assert type(categorical_lazy_arr[...]) == pd.Categorical


def test_lazy_categorical_array_equality(categorical_lazy_arr):
    assert (categorical_lazy_arr[0] == "foo").all()
    assert (categorical_lazy_arr[3:5] == "bar").all()
    assert (categorical_lazy_arr == "foo").any()


def test_lazy_categorical_array_subset_subset(categorical_lazy_arr):
    subset_susbet = categorical_lazy_arr[0:10][5:10]
    assert len(subset_susbet) == 5
    assert type(subset_susbet) == pd.Categorical
    assert (
        subset_susbet[...]
        == pd.Categorical.from_codes(
            codes=[2, 2, 1, 2, 0],
            categories=["foo", "bar", "jazz"],
            ordered=False,
        ).remove_unused_categories()
    ).all()


def test_nullable_boolean_array_properties(nullable_boolean_lazy_arr):
    assert len(nullable_boolean_lazy_arr[0:3]) == 3
    assert type(nullable_boolean_lazy_arr[0:3]) == pd.arrays.BooleanArray
    assert len(nullable_boolean_lazy_arr[...]) == len(nullable_boolean_lazy_arr)
    assert type(nullable_boolean_lazy_arr[...]) == pd.arrays.BooleanArray


def test_nullable_boolean_array_equality(nullable_boolean_lazy_arr):
    assert (nullable_boolean_lazy_arr[0] == pd.NA).all()
    assert (nullable_boolean_lazy_arr[3:5] == pd.NA).all()
    assert (nullable_boolean_lazy_arr[5:7] == np.array([True, False])).all()


def test_nullable_boolean_array_subset_subset(nullable_boolean_lazy_arr):
    subset_susbet = nullable_boolean_lazy_arr[0:10][5:10]
    assert len(subset_susbet) == 5
    assert type(subset_susbet) == pd.arrays.BooleanArray
    assert (
        subset_susbet[...]
        == pd.arrays.BooleanArray(
            values=np.array([True, False, False, True, True]),
            mask=np.array([False, False, True, False, True]),
        )
    ).all()


def test_nullable_boolean_array_no_mask_equality(nullable_boolean_lazy_arr_no_mask):
    assert nullable_boolean_lazy_arr_no_mask[0] == True  # noqa
    assert (nullable_boolean_lazy_arr_no_mask[3:5] == False).all()  # noqa
    assert (nullable_boolean_lazy_arr_no_mask[5:7] == np.array([True, False])).all()


def test_nullable_boolean_array_no_mask_subset_subset(
    nullable_boolean_lazy_arr_no_mask,
):
    subset_susbet = nullable_boolean_lazy_arr_no_mask[0:10][5:10]
    assert len(subset_susbet) == 5
    assert type(subset_susbet) == pd.arrays.BooleanArray
    assert (
        subset_susbet[...]
        == pd.array(
            np.array([True, False, False, True, True]),
        )
    ).all()


def test_nullable_integer_array_properties(nullable_integer_lazy_arr):
    assert len(nullable_integer_lazy_arr[0:3]) == 3
    assert type(nullable_integer_lazy_arr[0:3]) == pd.arrays.IntegerArray
    assert len(nullable_integer_lazy_arr[...]) == len(nullable_integer_lazy_arr)
    assert type(nullable_integer_lazy_arr[...]) == pd.arrays.IntegerArray


def test_nullable_integer_array_equality(nullable_integer_lazy_arr):
    assert (nullable_integer_lazy_arr[0] == pd.NA).all()
    assert (nullable_integer_lazy_arr[3:5] == pd.NA).all()
    assert (nullable_integer_lazy_arr[5:7] == np.array([2, 2])).all()


def test_nullable_integer_array_subset_subset(nullable_integer_lazy_arr):
    subset_susbet = nullable_integer_lazy_arr[0:10][5:10]
    assert len(subset_susbet) == 5
    assert type(subset_susbet) == pd.arrays.IntegerArray
    assert (
        subset_susbet[...]
        == pd.arrays.IntegerArray(
            values=np.array([2, 2, 1, 2, 0]),
            mask=np.array([False, False, True, False, True]),
        )
    ).all()


def test_nullable_integer_array_no_mask_equality(nullable_integer_lazy_arr_no_mask):
    assert (nullable_integer_lazy_arr_no_mask[0] == pd.NA).all()
    assert (nullable_integer_lazy_arr_no_mask[3:5] == 1).all()
    assert (nullable_integer_lazy_arr_no_mask[5:7] == np.array([2, 2])).all()


def test_nullable_integer_array_no_mask_subset_subset(
    nullable_integer_lazy_arr_no_mask,
):
    subset_susbet = nullable_integer_lazy_arr_no_mask[0:10][5:10]
    assert len(subset_susbet) == 5
    assert type(subset_susbet) == pd.arrays.IntegerArray
    assert (
        subset_susbet[...]
        == pd.array(
            np.array([2, 2, 1, 2, 0]),
        )
    ).all()
