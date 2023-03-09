from pathlib import Path
import re

import joblib
import pytest
import numpy as np
from scipy import sparse

import anndata as ad
from anndata.tests.helpers import (
    as_dense_dask_array,
    GEN_ADATA_DASK_ARGS,
    gen_adata,
    assert_equal,
    subset_func,
)
from anndata.experimental.read_remote import read_remote
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


def test_read_write_X(tmp_path, mtx_format):
    base_pth = Path(tmp_path)
    orig_pth = base_pth / "orig.zarr"
    # remote_pth = base_pth / "backed.zarr"

    orig = ad.AnnData(mtx_format(asarray(sparse.random(10, 10, format="csr"))))
    orig.write_zarr(orig_pth)

    remote = read_remote(orig_pth)
    # remote.write_zarr(remote_pth) # need to implement writing!

    assert np.all(asarray(orig.X) == asarray(remote.X))
