from pathlib import Path

import pytest
import numpy as np
import pandas as pd
from scipy import sparse
import zarr
from anndata._core.anndata import AnnData

from anndata.tests.helpers import (
    as_dense_dask_array,
    gen_adata,
    gen_typed_df,
)
from anndata.experimental.read_backed import read_backed, LazyCategoricalArray, LazyMaskedArray
from anndata.utils import asarray

from zarr import DirectoryStore
import traceback
class AccessTrackingStore(DirectoryStore):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._access_count = {}

    def __getitem__(self, key):
        print(key)
        for tracked in self._access_count:
            if tracked in key:
                self._access_count[tracked] += 1
        return super().__getitem__(key)

    def get_access_count(self, key):
        return self._access_count[key]

    def set_key_trackers(self, keys_to_track):
        for k in keys_to_track:
            self._access_count[k] = 0


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


def test_access_count_obs_var(tmp_path, mtx_format):
    base_pth = Path(tmp_path)
    orig_pth = base_pth / "orig.zarr"
    M = 1000000 # forces zarr to chunk `obs` columns multiple ways - that way 1 access to `int64` below is actually only one access
    N = 5
    obs_names = pd.Index(f"cell{i}" for i in range(M))
    var_names = pd.Index(f"gene{i}" for i in range(N))
    obs = gen_typed_df(M, obs_names)
    var = gen_typed_df(N, var_names)
    orig = AnnData(obs=obs, var=var, X=mtx_format(np.random.binomial(100, 0.005, (M, N)).astype(np.float32)))
    orig.write_zarr(orig_pth)
    store = AccessTrackingStore(orig_pth)
    store.set_key_trackers(["obs/int64", "var/int64", "obs/cat/codes"])
    remote = read_backed(store)
    # a series of methods that should __not__ read in any data
    remote.X
    remote.shape
    remote.var
    remote.obs
    remote.obs['int64']
    remote.var['int64']
    # only the `cat` should be read in
    subset = remote[remote.obs['cat'] == 'a', :]
    subset.obs['int64']
    assert store.get_access_count("obs/int64") == 0
    assert store.get_access_count("var/int64") == 0
    assert store.get_access_count("obs/cat/codes") == 4 # entire thing needs to be read in for subset
    remote[0:10, :].obs['int64'][0:10].compute()
    assert store.get_access_count("obs/int64") == 1 # one for 0, .zmetadata handles .zarray
    assert store.get_access_count("var/int64") == 0 # never accessed

def test_access_count_obsp_varp_layers(tmp_path):
    base_pth = Path(tmp_path)
    orig_pth = base_pth / "orig.zarr"
    M = 1000
    N = 1000
    orig = gen_adata((M, N), mtx_format)
    orig.write_zarr(orig_pth)
    store = AccessTrackingStore(orig_pth)
    store.set_key_trackers(["obsp", "varp", "layers"])
    remote = read_backed(store)
    # these operations should not read in any data
    subset = remote[0:10, 500:600]
    subset.obsp['array']
    subset.varp['array']
    subset.layers['array']
    assert store.get_access_count("obsp") == 0
    assert store.get_access_count("varp") == 0
    assert store.get_access_count("layers") == 0


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



    
