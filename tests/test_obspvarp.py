# TODO: These tests should share code with test_layers, and test_obsmvarm
from __future__ import annotations

import warnings

import joblib
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from anndata import AnnData
from anndata.tests.helpers import gen_typed_df_t2_size
from anndata.utils import asarray

M, N = (200, 100)


@pytest.fixture
def adata():
    X = np.zeros((M, N))
    obs = pd.DataFrame(
        dict(batch=np.array(["a", "b"])[np.random.randint(0, 2, M)]),
        index=[f"cell{i:03d}" for i in range(M)],
    )
    var = pd.DataFrame(index=[f"gene{i:03d}" for i in range(N)])
    return AnnData(X, obs=obs, var=var)


def test_assigmnent_dict(adata: AnnData):
    d_obsp = dict(
        a=pd.DataFrame(np.ones((M, M)), columns=adata.obs_names, index=adata.obs_names),
        b=np.zeros((M, M)),
        c=sparse.random(M, M, format="csr"),
    )
    d_varp = dict(
        a=pd.DataFrame(np.ones((N, N)), columns=adata.var_names, index=adata.var_names),
        b=np.zeros((N, N)),
        c=sparse.random(N, N, format="csr"),
    )
    adata.obsp = d_obsp
    for k, v in d_obsp.items():
        assert np.all(asarray(adata.obsp[k]) == asarray(v))
    adata.varp = d_varp
    for k, v in d_varp.items():
        assert np.all(asarray(adata.varp[k]) == asarray(v))


def test_setting_ndarray(adata: AnnData):
    adata.obsp["a"] = np.ones((M, M))
    adata.varp["a"] = np.ones((N, N))
    assert np.all(adata.obsp["a"] == np.ones((M, M)))
    assert np.all(adata.varp["a"] == np.ones((N, N)))

    h = joblib.hash(adata)
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.obsp["b"] = np.ones((int(M / 2), M))
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.obsp["b"] = np.ones((M, int(M * 2)))
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.varp["b"] = np.ones((int(N / 2), 10))
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.varp["b"] = np.ones((N, int(N * 2)))
    assert h == joblib.hash(adata)


def test_setting_sparse(adata: AnnData):
    obsp_sparse = sparse.random(M, M, format="csr")
    adata.obsp["a"] = obsp_sparse
    assert not np.any((adata.obsp["a"] != obsp_sparse).data)

    varp_sparse = sparse.random(N, N, format="csr")
    adata.varp["a"] = varp_sparse
    assert not np.any((adata.varp["a"] != varp_sparse).data)

    h = joblib.hash(adata)

    bad_obsp_sparse = sparse.random(M * 2, M, format="csr")
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.obsp["b"] = bad_obsp_sparse

    bad_varp_sparse = sparse.random(N * 2, N, format="csr")
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.varp["b"] = bad_varp_sparse

    assert h == joblib.hash(adata)


@pytest.mark.parametrize(("field", "dim"), [("obsp", M), ("varp", N)])
@pytest.mark.parametrize(
    ("df", "homogenous", "dtype"),
    [
        (lambda dim: gen_typed_df_t2_size(dim, dim), True, np.object_),
        (lambda dim: pd.DataFrame(np.random.randn(dim, dim)), False, np.floating),
    ],
    ids=["heterogeneous", "homogeneous"],
)
def test_setting_dataframe(adata: AnnData, field, dim, homogenous, df, dtype):
    if homogenous:
        with pytest.warns(UserWarning, match=rf"{field.title()} 'df'.*dtype object"):
            getattr(adata, field)["df"] = df(dim)
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            getattr(adata, field)["df"] = df(dim)
    assert isinstance(getattr(adata, field)["df"], np.ndarray)
    assert np.issubdtype(getattr(adata, field)["df"].dtype, dtype)


def test_setting_daskarray(adata: AnnData):
    import dask.array as da

    adata.obsp["a"] = da.ones((M, M))
    adata.varp["a"] = da.ones((N, N))
    assert da.all(adata.obsp["a"] == da.ones((M, M)))
    assert da.all(adata.varp["a"] == da.ones((N, N)))
    assert isinstance(adata.obsp["a"], da.Array)
    assert isinstance(adata.varp["a"], da.Array)

    h = joblib.hash(adata)
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.obsp["b"] = da.ones((int(M / 2), M))
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.obsp["b"] = da.ones((M, int(M * 2)))
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.varp["b"] = da.ones((int(N / 2), 10))
    with pytest.raises(ValueError, match=r"incorrect shape"):
        adata.varp["b"] = da.ones((N, int(N * 2)))
    assert h == joblib.hash(adata)


def test_shape_error(adata: AnnData):
    with pytest.raises(
        ValueError,
        match=(
            r"Value passed for key 'a' is of incorrect shape\. "
            r"Values of obsp must match dimensions \('obs', 'obs'\) of parent\. "
            r"Value had shape \(201, 200\) while it should have had \(200, 200\)\."
        ),
    ):
        adata.obsp["a"] = np.zeros((adata.shape[0] + 1, adata.shape[0]))
