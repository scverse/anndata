import joblib
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import anndata

M, N = (100, 100)


@pytest.fixture
def adata():
    X = np.zeros((M, N))
    obs = pd.DataFrame(
        dict(batch=np.array(["a", "b"])[np.random.randint(0, 2, M)]),
        index=[f"cell{i:03d}" for i in range(N)],
    )
    var = pd.DataFrame(index=[f"gene{i:03d}" for i in range(N)])
    return anndata.AnnData(X, obs=obs, var=var)


def test_assigmnent_dict(adata):
    d_obsm = dict(
        a=pd.DataFrame(
            dict(a1=np.ones(M), a2=[f"a{i}" for i in range(M)]), index=adata.obs_names,
        ),
        b=np.zeros((M, 2)),
    )
    d_varm = dict(
        a=pd.DataFrame(
            dict(a1=np.ones(N), a2=[f"a{i}" for i in range(N)]), index=adata.var_names,
        ),
        b=np.zeros((N, 2)),
    )
    adata.obsm = d_obsm
    for k, v in d_obsm.items():
        assert np.all(adata.obsm[k] == v)
    adata.varm = d_varm
    for k, v in d_varm.items():
        assert np.all(adata.varm[k] == v)


def test_setting_ndarray(adata):
    adata.obsm["a"] = np.ones((M, 10))
    adata.varm["a"] = np.ones((N, 10))
    assert np.all(adata.obsm["a"] == np.ones((M, 10)))
    assert np.all(adata.varm["a"] == np.ones((N, 10)))

    h = joblib.hash(adata)
    with pytest.raises(ValueError):
        adata.obsm["b"] = np.ones((int(M / 2), 10))
    with pytest.raises(ValueError):
        adata.obsm["b"] = np.ones((int(M * 2), 10))
    with pytest.raises(ValueError):
        adata.varm["b"] = np.ones((int(N / 2), 10))
    with pytest.raises(ValueError):
        adata.varm["b"] = np.ones((int(N * 2), 10))
    assert h == joblib.hash(adata)


def test_setting_dataframe(adata):
    obsm_df = pd.DataFrame(dict(b_1=np.ones(M), b_2=["a"] * M), index=adata.obs_names)
    varm_df = pd.DataFrame(dict(b_1=np.ones(N), b_2=["a"] * N), index=adata.var_names)

    adata.obsm["b"] = obsm_df
    assert np.all(adata.obsm["b"] == obsm_df)
    adata.varm["b"] = varm_df
    assert np.all(adata.varm["b"] == varm_df)

    bad_obsm_df = obsm_df.copy()
    bad_obsm_df.reset_index(inplace=True)
    with pytest.raises(ValueError):
        adata.obsm["c"] = bad_obsm_df

    bad_varm_df = varm_df.copy()
    bad_varm_df.reset_index(inplace=True)
    with pytest.raises(ValueError):
        adata.varm["c"] = bad_varm_df


def test_setting_sparse(adata):
    obsm_sparse = sparse.random(M, 100)
    adata.obsm["a"] = obsm_sparse
    assert not np.any((adata.obsm["a"] != obsm_sparse).data)

    varm_sparse = sparse.random(N, 100)
    adata.varm["a"] = varm_sparse
    assert not np.any((adata.varm["a"] != varm_sparse).data)

    h = joblib.hash(adata)

    bad_obsm_sparse = sparse.random(M * 2, M)
    with pytest.raises(ValueError):
        adata.obsm["b"] = bad_obsm_sparse

    bad_varm_sparse = sparse.random(N * 2, N)
    with pytest.raises(ValueError):
        adata.varm["b"] = bad_varm_sparse

    assert h == joblib.hash(adata)
