from pathlib import Path

import pytest
import numpy as np
import pandas as pd
from scipy import sparse
import zarr

import anndata as ad
from anndata.tests.helpers import (
    as_dense_dask_array,
    subset_func,
)
from anndata.experimental.read_remote import read_remote, LazyCategoricalArray
from anndata.utils import asarray

subset_func2 = subset_func
# -------------------------------------------------------------------------------
# Some test data
# -------------------------------------------------------------------------------


@pytest.fixture
def adata():
    X_list = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
    ]  # data matrix of shape n_obs x n_vars
    X = np.array(X_list)
    obs_dict = dict(  # annotation of observations / rows
        row_names=["name1", "name2", "name3"],  # row annotation
        oanno1=["cat1", "cat2", "cat2"],  # categorical annotation
        oanno2=["o1", "o2", "o3"],  # string annotation
        oanno3=[2.1, 2.2, 2.3],  # float annotation
    )
    var_dict = dict(vanno1=[3.1, 3.2, 3.3])  # annotation of variables / columns
    uns_dict = dict(  # unstructured annotation
        oanno1_colors=["#000000", "#FFFFFF"], uns2=["some annotation"]
    )
    return ad.AnnData(
        X,
        obs=obs_dict,
        var=var_dict,
        uns=uns_dict,
        obsm=dict(o1=np.zeros((X.shape[0], 10))),
        varm=dict(v1=np.ones((X.shape[1], 20))),
        layers=dict(float=X.astype(float), sparse=sparse.csr_matrix(X)),
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


@pytest.fixture()
def categorical_zarr_group(tmp_path_factory):
    base_path = tmp_path_factory.getbasetemp()
    z = zarr.open_group(base_path, mode="w")
    z["codes"] = [0, 1, 0, 1, 1, 2, 2]
    z["categories"] = ["foo", "bar", "jazz"]
    z.attrs["ordered"] = False
    z = zarr.open(base_path)
    return z


def test_read_write_X(tmp_path, mtx_format):
    base_pth = Path(tmp_path)
    orig_pth = base_pth / "orig.zarr"
    # remote_pth = base_pth / "backed.zarr"

    orig = ad.AnnData(mtx_format(asarray(sparse.random(10, 10, format="csr"))))
    orig.write_zarr(orig_pth)

    remote = read_remote(orig_pth)
    # remote.write_zarr(remote_pth) # need to implement writing!

    assert np.all(asarray(orig.X) == asarray(remote.X))
    assert (orig.obs == remote.obs.to_df()).all().all()
    assert (orig.var == remote.var.to_df()).all().all()


def test_read_write_full(adata, tmp_path):
    base_pth = Path(tmp_path)
    orig_pth = base_pth / "orig.zarr"
    adata.write_zarr(orig_pth)
    remote = read_remote(orig_pth)
    assert np.all(asarray(adata.X) == asarray(remote.X))
    assert (adata.obs == remote.obs.to_df()).all().all()
    assert (adata.var == remote.var.to_df()).all().all()


def test_read_write_view(adata, tmp_path):
    base_pth = Path(tmp_path)
    orig_pth = base_pth / "orig.zarr"
    adata.write_zarr(orig_pth)
    remote = read_remote(orig_pth)
    subset = adata.obs["oanno1"] == "cat1"
    assert np.all(asarray(adata[subset, :].X) == asarray(remote[subset, :].X))
    assert (adata[subset, :].obs == remote[subset, :].obs.to_df()).all().all()
    assert (adata[subset, :].var == remote[subset, :].var.to_df()).all().all()


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
