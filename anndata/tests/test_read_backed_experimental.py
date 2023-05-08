from pathlib import Path

import pytest
import numpy as np
import pandas as pd
from scipy import sparse
import zarr

from anndata.tests.helpers import (
    as_dense_dask_array,
    gen_adata,
    subset_func,
)
from anndata.experimental.read_backed import read_backed, LazyCategoricalArray, LazyMaskedArray
from anndata.utils import asarray


@pytest.fixture(
    params=[sparse.csr_matrix, sparse.csc_matrix, np.array, as_dense_dask_array],
    ids=["scipy-csr", "scipy-csc", "np-array", "dask_array"],
)
def mtx_format(request):
    return request.param


@pytest.fixture(params=[sparse.csr_matrix, sparse.csc_matrix])
def sparse_format(request):
    return request.param


@pytest.fixture()
def categorical_zarr_group(tmp_path_factory):
    base_path = tmp_path_factory.getbasetemp()
    z = zarr.open_group(base_path, mode="w")
    z["codes"] =np.array([0, 1, 0, 1, 1, 2, 2, 1, 2, 0, 1, 1, 1, 2, 1, 2])
    z["categories"] = np.array(["foo", "bar", "jazz"])
    z.attrs["ordered"] = False
    z = zarr.open(base_path)
    return z

@pytest.fixture()
def nullable_boolean_zarr_group(tmp_path_factory):
    base_path = tmp_path_factory.getbasetemp()
    z = zarr.open_group(base_path, mode="w")
    z["values"] =np.array([True, False, True, False, False, True, False, False, True, True, False, False, False, True, False, True])
    z["mask"] = np.array([True, True, True, True, True, False, False, True, False, True, True, True, True, False, True, False])
    z = zarr.open(base_path)
    return z

@pytest.fixture()
def nullable_boolean_zarr_group_no_mask(tmp_path_factory):
    base_path = tmp_path_factory.getbasetemp()
    z = zarr.open_group(base_path, mode="w")
    z["values"] = np.array([True, False, True, False, False, True, False, False, True, True, False, False, False, True, False, True])
    z = zarr.open(base_path)
    return z

@pytest.fixture()
def nullable_integer_zarr_group(tmp_path_factory):
    base_path = tmp_path_factory.getbasetemp()
    z = zarr.open_group(base_path, mode="w")
    z["values"] = np.array([0, 1, 0, 1, 1, 2, 2, 1, 2, 0, 1, 1, 1, 2, 1, 2])
    z["mask"] = np.array([True, True, True, True, True, False, False, True, False, True, True, True, True, False, True, False])
    z = zarr.open(base_path)
    return z

@pytest.fixture()
def nullable_integer_zarr_group_no_mask(tmp_path_factory):
    base_path = tmp_path_factory.getbasetemp()
    z = zarr.open_group(base_path, mode="w")
    z["values"] = np.array([0, 1, 0, 1, 1, 2, 2, 1, 2, 0, 1, 1, 1, 2, 1, 2])
    z = zarr.open(base_path)
    return z


def test_read_write_X(tmp_path, mtx_format):
    base_pth = Path(tmp_path)
    orig_pth = base_pth / "orig.zarr"
    # remote_pth = base_pth / "backed.zarr"

    orig = gen_adata((1000, 1000), mtx_format)
    orig.write_zarr(orig_pth)

    remote = read_backed(orig_pth)
    # remote.write_zarr(remote_pth) # need to implement writing!

    assert np.all(asarray(orig.X) == asarray(remote.X))
    assert (orig.obs == remote.obs.to_df()[orig.obs.columns]).all().all()
    assert (orig.var == remote.var.to_df()[orig.var.columns]).all().all()
    assert (orig.obsm["array"] == remote.obsm["array"].compute()).all()


def test_read_write_full(tmp_path, mtx_format):
    adata = gen_adata((1000, 1000), mtx_format)
    base_pth = Path(tmp_path)
    orig_pth = base_pth / "orig.zarr"
    adata.write_zarr(orig_pth)
    remote = read_backed(orig_pth)
    assert np.all(asarray(adata.X) == asarray(remote.X))
    assert (adata.obs == remote.obs.to_df()[adata.obs.columns]).all().all()
    assert (adata.var == remote.var.to_df()[adata.var.columns]).all().all()
    assert (adata.obsm["array"] == remote.obsm["array"].compute()).all()


def test_read_write_view(tmp_path, mtx_format):
    adata = gen_adata((1000, 1000), mtx_format)
    base_pth = Path(tmp_path)
    orig_pth = base_pth / "orig.zarr"
    adata.write_zarr(orig_pth)
    remote = read_backed(orig_pth)
    subset = adata.obs["obs_cat"] == "a"
    assert np.all(asarray(adata[subset, :].X) == asarray(remote[subset, :].X))
    assert (
        (adata[subset, :].obs == remote[subset, :].obs.to_df()[adata.obs.columns])
        .all()
        .all()
    )
    assert (
        (adata[subset, :].var == remote[subset, :].var.to_df()[adata.var.columns])
        .all()
        .all()
    )
    assert (
        adata[subset, :].obsm["array"] == remote[subset, :].obsm["array"].compute()
    ).all()


def test_read_write_view_of_view(tmp_path, mtx_format):
    adata = gen_adata((1000, 1000), mtx_format)
    base_pth = Path(tmp_path)
    orig_pth = base_pth / "orig.zarr"
    adata.write_zarr(orig_pth)
    remote = read_backed(orig_pth)
    subset = (adata.obs["obs_cat"] == "a") | (adata.obs["obs_cat"] == "b")
    subsetted_adata = adata[subset, :]
    subset_subset = subsetted_adata.obs["obs_cat"] == "b"
    subsetted_subsetted_adata = subsetted_adata[subset_subset, :]
    assert np.all(
        asarray(subsetted_subsetted_adata.X)
        == asarray(remote[subset, :][subset_subset, :].X)
    )
    assert (
        (
            subsetted_subsetted_adata.obs
            == remote[subset, :][subset_subset, :].obs.to_df()[adata.obs.columns]
        )
        .all()
        .all()
    )
    assert (
        (
            subsetted_subsetted_adata.var
            == remote[subset, :][subset_subset, :].var.to_df()[adata.var.columns]
        )
        .all()
        .all()
    )
    assert (
        subsetted_subsetted_adata.obsm["array"]
        == remote[subset, :][subset_subset, :].obsm["array"].compute()
    ).all()


def test_lazy_categorical_array_properties(categorical_zarr_group):
    arr = LazyCategoricalArray(categorical_zarr_group)
    assert len(arr[0:3]) == 3
    assert type(arr[0:3]) == pd.Categorical
    assert len(arr[()]) == len(arr)
    assert type(arr[()]) == pd.Categorical


def test_lazy_categorical_array_equality(categorical_zarr_group):
    arr = LazyCategoricalArray(categorical_zarr_group)
    assert (arr[0] == "foo").all()
    assert (arr[3:5] == "bar").all()
    assert (arr == "foo").any()


def test_lazy_categorical_array_subset_subset(categorical_zarr_group):
    arr = LazyCategoricalArray(categorical_zarr_group)
    subset_susbet = arr[0:10][5:10]
    assert len(subset_susbet) == 5
    assert type(subset_susbet) == pd.Categorical
    assert (
        subset_susbet[()]
        == pd.Categorical.from_codes(
            codes=[2, 2, 1, 2, 0],
            categories=["foo", "bar", "jazz"],
            ordered=False,
        ).remove_unused_categories()
    ).all()


def test_nullable_boolean_array_properties(nullable_boolean_zarr_group):
    arr = LazyMaskedArray(nullable_boolean_zarr_group, 'nullable-boolean')
    assert len(arr[0:3]) == 3
    assert type(arr[0:3]) == pd.arrays.BooleanArray
    assert len(arr[()]) == len(arr)
    assert type(arr[()]) == pd.arrays.BooleanArray


def test_nullable_boolean_array_equality(nullable_boolean_zarr_group):
    arr = LazyMaskedArray(nullable_boolean_zarr_group, 'nullable-boolean')
    assert (arr[0] == pd.NA).all()
    assert (arr[3:5] == pd.NA).all()
    assert (arr[5:7] == np.array([True, False])).all()


def test_nullable_boolean_array_subset_subset(nullable_boolean_zarr_group):
    arr = LazyMaskedArray(nullable_boolean_zarr_group, 'nullable-boolean')
    subset_susbet = arr[0:10][5:10]
    assert len(subset_susbet) == 5
    assert type(subset_susbet) == pd.arrays.BooleanArray
    assert (
        subset_susbet[()]
        == pd.arrays.BooleanArray(
            values=np.array([True, False, False, True, True]),
            mask=np.array([False, False, True, False, True]),
        )
    ).all()

def test_nullable_boolean_array_no_mask_equality(nullable_boolean_zarr_group_no_mask):
    arr = LazyMaskedArray(nullable_boolean_zarr_group_no_mask, 'nullable-boolean')
    assert (arr[0] == True).all()
    assert (arr[3:5] == False).all()
    assert (arr[5:7] == np.array([True, False])).all()


def test_nullable_boolean_array_no_mask_subset_subset(nullable_boolean_zarr_group_no_mask):
    arr = LazyMaskedArray(nullable_boolean_zarr_group_no_mask, 'nullable-boolean')
    subset_susbet = arr[0:10][5:10]
    assert len(subset_susbet) == 5
    assert type(subset_susbet) == pd.arrays.BooleanArray
    assert (
        subset_susbet[()]
        == pd.array(np.array([True, False, False, True, True]),)
    ).all()

def test_nullable_integer_array_properties(nullable_integer_zarr_group):
    arr = LazyMaskedArray(nullable_integer_zarr_group, 'nullable-integer')
    assert len(arr[0:3]) == 3
    assert type(arr[0:3]) == pd.arrays.IntegerArray
    assert len(arr[()]) == len(arr)
    assert type(arr[()]) == pd.arrays.IntegerArray


def test_nullable_integer_array_equality(nullable_integer_zarr_group):
    arr = LazyMaskedArray(nullable_integer_zarr_group, 'nullable-integer')
    assert (arr[0] == pd.NA).all()
    assert (arr[3:5] == pd.NA).all()
    assert (arr[5:7] == np.array([2, 2])).all()


def test_nullable_integer_array_subset_subset(nullable_integer_zarr_group):
    arr = LazyMaskedArray(nullable_integer_zarr_group, 'nullable-integer')
    subset_susbet = arr[0:10][5:10]
    assert len(subset_susbet) == 5
    assert type(subset_susbet) == pd.arrays.IntegerArray
    assert (
        subset_susbet[()]
        == pd.arrays.IntegerArray(
            values=np.array([2, 2, 1, 2, 0]),
            mask=np.array([False, False, True, False, True]),
        )
    ).all()

def test_nullable_integer_array_no_mask_equality(nullable_integer_zarr_group_no_mask):
    arr = LazyMaskedArray(nullable_integer_zarr_group_no_mask, 'nullable-integer')
    assert (arr[0] == pd.NA).all()
    assert (arr[3:5] == 1).all()
    assert (arr[5:7] == np.array([2, 2])).all()


def test_nullable_integer_array_no_mask_subset_subset(nullable_integer_zarr_group_no_mask):
    arr = LazyMaskedArray(nullable_integer_zarr_group_no_mask, 'nullable-integer')
    subset_susbet = arr[0:10][5:10]
    assert len(subset_susbet) == 5
    assert type(subset_susbet) == pd.arrays.IntegerArray
    assert (
        subset_susbet[()]
        == pd.array(np.array([2, 2, 1, 2, 0]),)
    ).all()



    
