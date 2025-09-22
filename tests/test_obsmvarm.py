from __future__ import annotations

from functools import partial

import joblib
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from anndata import AnnData
from anndata.compat import CupyArray
from anndata.tests.helpers import as_cupy, get_multiindex_columns_df

M, N = (100, 100)


@pytest.fixture(
    params=[
        pytest.param(
            partial(as_cupy, typ=CupyArray), id="cupy_array", marks=pytest.mark.gpu
        ),
        pytest.param(np.array, id="numpy_array"),
    ],
    ids=["cupy", "numpy"],
)
def array_type(request):
    return request.param


@pytest.fixture
def adata():
    X = np.zeros((M, N))
    obs = pd.DataFrame(
        dict(batch=np.array(["a", "b"])[np.random.randint(0, 2, M)]),
        index=[f"cell{i:03d}" for i in range(N)],
    )
    var = pd.DataFrame(index=[f"gene{i:03d}" for i in range(N)])
    return AnnData(X, obs=obs, var=var)


def test_assignment_dict(adata: AnnData):
    d_obsm = dict(
        a=pd.DataFrame(
            dict(a1=np.ones(M), a2=[f"a{i}" for i in range(M)]),
            index=adata.obs_names,
        ),
        b=np.zeros((M, 2)),
    )
    d_varm = dict(
        a=pd.DataFrame(
            dict(a1=np.ones(N), a2=[f"a{i}" for i in range(N)]),
            index=adata.var_names,
        ),
        b=np.zeros((N, 2)),
    )
    adata.obsm = d_obsm
    for k, v in d_obsm.items():
        assert np.all(adata.obsm[k] == v)
    adata.varm = d_varm
    for k, v in d_varm.items():
        assert np.all(adata.varm[k] == v)


def test_setting_ndarray(adata: AnnData):
    adata.obsm["a"] = np.ones((M, 10))
    adata.varm["a"] = np.ones((N, 10))
    assert np.all(adata.obsm["a"] == np.ones((M, 10)))
    assert np.all(adata.varm["a"] == np.ones((N, 10)))

    h = joblib.hash(adata)
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.obsm["b"] = np.ones((int(M / 2), 10))
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.obsm["b"] = np.ones((int(M * 2), 10))
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.varm["b"] = np.ones((int(N / 2), 10))
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.varm["b"] = np.ones((int(N * 2), 10))
    assert h == joblib.hash(adata)


def test_setting_dataframe(adata: AnnData):
    obsm_df = pd.DataFrame(dict(b_1=np.ones(M), b_2=["a"] * M), index=adata.obs_names)
    varm_df = pd.DataFrame(dict(b_1=np.ones(N), b_2=["a"] * N), index=adata.var_names)

    adata.obsm["b"] = obsm_df
    assert np.all(adata.obsm["b"] == obsm_df)
    adata.varm["b"] = varm_df
    assert np.all(adata.varm["b"] == varm_df)

    bad_obsm_df = obsm_df.copy()
    bad_obsm_df.reset_index(inplace=True)
    with pytest.raises(ValueError, match=r"index does not match.*obs names"):
        adata.obsm["c"] = bad_obsm_df

    bad_varm_df = varm_df.copy()
    bad_varm_df.reset_index(inplace=True)
    with pytest.raises(ValueError, match=r"index does not match.*var names"):
        adata.varm["c"] = bad_varm_df


def test_setting_sparse(adata: AnnData):
    obsm_sparse = sparse.random(M, 100, format="csr")
    adata.obsm["a"] = obsm_sparse
    assert not np.any((adata.obsm["a"] != obsm_sparse).data)

    varm_sparse = sparse.random(N, 100, format="csr")
    adata.varm["a"] = varm_sparse
    assert not np.any((adata.varm["a"] != varm_sparse).data)

    h = joblib.hash(adata)

    bad_obsm_sparse = sparse.random(M * 2, M, format="csr")
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.obsm["b"] = bad_obsm_sparse

    bad_varm_sparse = sparse.random(N * 2, N, format="csr")
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.varm["b"] = bad_varm_sparse

    assert h == joblib.hash(adata)


def test_setting_daskarray(adata: AnnData):
    import dask.array as da

    adata.obsm["a"] = da.ones((M, 10))
    adata.varm["a"] = da.ones((N, 10))
    assert da.all(adata.obsm["a"] == da.ones((M, 10)))
    assert da.all(adata.varm["a"] == da.ones((N, 10)))
    assert isinstance(adata.obsm["a"], da.Array)
    assert isinstance(adata.varm["a"], da.Array)

    h = joblib.hash(adata)
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.obsm["b"] = da.ones((int(M / 2), 10))
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.obsm["b"] = da.ones((int(M * 2), 10))
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.varm["b"] = da.ones((int(N / 2), 10))
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.varm["b"] = da.ones((int(N * 2), 10))
    assert h == joblib.hash(adata)


def test_shape_error(adata: AnnData):
    with pytest.raises(
        ValueError,
        match=(
            r"Value passed for key 'b' is of incorrect shape\. "
            r"Values of obsm must match dimensions \('obs',\) of parent\. "
            r"Value had shape \(101,\) while it should have had \(100,\)\."
        ),
    ):
        adata.obsm["b"] = np.zeros((adata.shape[0] + 1, adata.shape[0]))


def test_error_set_multiindex_df(adata: AnnData):
    df = get_multiindex_columns_df((adata.shape[0], 20))
    with pytest.raises(ValueError, match=r"MultiIndex columns are not supported"):
        adata.obsm["df"] = df


def test_1d_declaration(array_type):
    adata = AnnData(np.ones((5, 20)), obsm={"1d-array": array_type(np.ones(5))})
    assert adata.obsm["1d-array"].shape == (5, 1)


def test_1d_set(adata, array_type):
    adata.varm["1d-array"] = array_type(np.ones(adata.shape[1]))
    assert adata.varm["1d-array"].shape == (adata.shape[1], 1)
