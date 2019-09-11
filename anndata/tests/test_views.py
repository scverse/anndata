from operator import ne, mul

import joblib
import numpy as np
from scipy import sparse
import pandas as pd
import pytest

import anndata as ad

from anndata.tests.helpers import (
    gen_adata,
    subset_func,
    slice_subset,
    single_subset,
    asarray,
    assert_equal,
)

# -------------------------------------------------------------------------------
# Some test data
# -------------------------------------------------------------------------------

# data matrix of shape n_obs x n_vars
X_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
# annotation of observations / rows
obs_dict = dict(
    row_names=['name1', 'name2', 'name3'],  # row annotation
    oanno1=['cat1', 'cat2', 'cat2'],  # categorical annotation
    oanno2=['o1', 'o2', 'o3'],  # string annotation
    oanno3=[2.1, 2.2, 2.3],  # float annotation
)
# annotation of variables / columns
var_dict = dict(vanno1=[3.1, 3.2, 3.3])
# unstructured annotation
uns_dict = dict(oanno1_colors=['#000000', '#FFFFFF'], uns2=['some annotation'])


subset_func2 = subset_func


@pytest.fixture
def adata():
    adata = ad.AnnData(np.zeros((100, 100)))
    adata.obsm['o'] = np.zeros((100, 50))
    adata.varm['o'] = np.zeros((100, 50))
    return adata


@pytest.fixture(params=[asarray, sparse.csr_matrix, sparse.csc_matrix])
def adata_parameterized(request):
    return gen_adata(shape=(200, 300), X_type=request.param)


@pytest.fixture(
    params=[np.array, sparse.csr_matrix, sparse.csc_matrix],
    ids=["np_array", "scipy_csr", "scipy_csc"],
)
def matrix_type(request):
    return request.param


@pytest.fixture(params=["layers", "obsm", "varm"])
def mapping_name(request):
    return request.param


# -------------------------------------------------------------------------------
# The test functions
# -------------------------------------------------------------------------------


def test_views():
    X = np.array(X_list)
    adata = ad.AnnData(
        X, obs=obs_dict, var=var_dict, uns=uns_dict, dtype='int32'
    )

    assert adata[:, 0].isview
    assert adata[:, 0].X.tolist() == np.reshape([1, 4, 7], (3, 1)).tolist()

    adata[:2, 0].X = [0, 0]

    assert adata[:, 0].X.tolist() == np.reshape([0, 0, 7], (3, 1)).tolist()

    adata_subset = adata[:2, [0, 1]]

    assert adata_subset.isview
    # now transition to actual object
    adata_subset.obs['foo'] = range(2)
    assert not adata_subset.isview

    assert adata_subset.obs['foo'].tolist() == list(range(2))


def test_modify_view_component(matrix_type, mapping_name):
    adata = ad.AnnData(
        np.zeros((10, 10)),
        **{mapping_name: {"m": matrix_type(asarray(sparse.random(10, 10)))}},
    )
    init_hash = joblib.hash(adata)

    subset = adata[:5, :][:, :5]
    assert subset.isview
    m = getattr(subset, mapping_name)["m"]
    m[0, 0] = 100
    assert not subset.isview
    assert getattr(subset, mapping_name)["m"][0, 0] == 100

    assert init_hash == joblib.hash(adata)


# These tests could probably be condensed into a fixture based test for obsm and varm
def test_set_obsm_key(adata):
    init_hash = joblib.hash(adata)

    orig_obsm_val = adata.obsm['o'].copy()
    subset_obsm = adata[:50]
    assert subset_obsm.isview
    subset_obsm.obsm['o'] = np.ones((50, 20))
    assert not subset_obsm.isview
    assert np.all(adata.obsm["o"] == orig_obsm_val)

    assert init_hash == joblib.hash(adata)


def test_set_varm_key(adata):
    init_hash = joblib.hash(adata)

    orig_varm_val = adata.varm['o'].copy()
    subset_varm = adata[:, :50]
    assert subset_varm.isview
    subset_varm.varm['o'] = np.ones((50, 20))
    assert not subset_varm.isview
    assert np.all(adata.varm["o"] == orig_varm_val)

    assert init_hash == joblib.hash(adata)


def test_set_obs(adata, subset_func):
    init_hash = joblib.hash(adata)

    subset = adata[subset_func(adata.obs_names), :]

    new_obs = pd.DataFrame(
        {"a": np.ones(subset.n_obs), "b": np.ones(subset.n_obs)},
        index=subset.obs_names,
    )

    assert subset.isview
    subset.obs = new_obs
    assert not subset.isview
    assert np.all(subset.obs == new_obs)

    assert joblib.hash(adata) == init_hash


def test_set_var(adata, subset_func):
    init_hash = joblib.hash(adata)

    subset = adata[:, subset_func(adata.var_names)]

    new_var = pd.DataFrame(
        dict(a=np.ones(subset.n_vars), b=np.ones(subset.n_vars)),
        index=subset.var_names,
    )

    assert subset.isview
    subset.var = new_var
    assert not subset.isview
    assert np.all(subset.var == new_var)

    assert joblib.hash(adata) == init_hash


def test_set_obsm(adata):
    init_hash = joblib.hash(adata)

    dim0_size = np.random.randint(2, adata.shape[0] - 1)
    dim1_size = np.random.randint(1, 99)
    orig_obsm_val = adata.obsm["o"].copy()
    subset_idx = np.random.choice(adata.obs_names, dim0_size, replace=False)

    subset = adata[subset_idx, :]
    assert subset.isview
    subset.obsm = {"o": np.ones((dim0_size, dim1_size))}
    assert not subset.isview
    assert np.all(orig_obsm_val == adata.obsm["o"])  # Checking for mutation
    assert np.all(subset.obsm["o"] == np.ones((dim0_size, dim1_size)))

    subset = adata[subset_idx, :]
    subset_hash = joblib.hash(subset)
    with pytest.raises(ValueError):
        subset.obsm = {"o": np.ones((dim0_size + 1, dim1_size))}
    with pytest.raises(ValueError):
        subset.varm = {"o": np.ones((dim0_size - 1, dim1_size))}
    assert subset_hash == joblib.hash(subset)

    assert init_hash == joblib.hash(
        adata
    )  # Only modification have been made to a view


def test_set_varm(adata):
    init_hash = joblib.hash(adata)

    dim0_size = np.random.randint(2, adata.shape[1] - 1)
    dim1_size = np.random.randint(1, 99)
    orig_varm_val = adata.varm["o"].copy()
    subset_idx = np.random.choice(adata.var_names, dim0_size, replace=False)

    subset = adata[:, subset_idx]
    assert subset.isview
    subset.varm = {"o": np.ones((dim0_size, dim1_size))}
    assert not subset.isview
    assert np.all(orig_varm_val == adata.varm["o"])  # Checking for mutation
    assert np.all(subset.varm["o"] == np.ones((dim0_size, dim1_size)))

    subset = adata[:, subset_idx]
    subset_hash = joblib.hash(subset)
    with pytest.raises(ValueError):
        subset.varm = {"o": np.ones((dim0_size + 1, dim1_size))}
    with pytest.raises(ValueError):
        subset.varm = {"o": np.ones((dim0_size - 1, dim1_size))}
    assert subset_hash == joblib.hash(
        subset
    )  # subset should not be changed by failed setting

    assert init_hash == joblib.hash(adata)


# TODO: Determine if this is the intended behavior, or just the behaviour we've had for a while
def test_not_set_subset_X(matrix_type, subset_func):
    adata = ad.AnnData(matrix_type(asarray(sparse.random(20, 20))))
    init_hash = joblib.hash(adata)
    orig_X_val = adata.X.copy()
    while True:
        subset_idx = slice_subset(adata.obs_names)
        if len(adata[subset_idx, :]) > 2:
            break
    subset = adata[subset_idx, :]

    subset = adata[:, subset_idx]

    internal_idx = subset_func(np.arange(subset.X.shape[1]))
    assert subset.isview
    subset.X[:, internal_idx] = 1
    assert not subset.isview
    assert not np.any(asarray(adata.X != orig_X_val))

    assert init_hash == joblib.hash(adata)


def test_set_scalar_subset_X(matrix_type, subset_func):
    adata = ad.AnnData(matrix_type(np.zeros((10, 10))))
    orig_X_val = adata.X.copy()
    subset_idx = slice_subset(adata.obs_names)

    adata_subset = adata[subset_idx, :]

    adata_subset.X = 1

    assert adata_subset.isview
    assert np.all(asarray(adata[subset_idx, :].X) == 1)

    assert asarray((orig_X_val != adata.X)).sum() == mul(*adata_subset.shape)


# TODO: Use different kind of subsetting for adata and view
def test_set_subset_obsm(adata, subset_func):
    init_hash = joblib.hash(adata)
    orig_obsm_val = adata.obsm["o"].copy()

    while True:
        subset_idx = slice_subset(adata.obs_names)
        if len(adata[subset_idx, :]) > 2:
            break
    subset = adata[subset_idx, :]

    internal_idx = subset_func(np.arange(subset.obsm["o"].shape[0]))
    assert subset.isview
    subset.obsm["o"][internal_idx] = 1
    assert not subset.isview
    assert np.all(adata.obsm["o"] == orig_obsm_val)

    assert init_hash == joblib.hash(adata)


def test_set_subset_varm(adata, subset_func):
    init_hash = joblib.hash(adata)
    orig_varm_val = adata.varm["o"].copy()

    while True:
        subset_idx = slice_subset(adata.var_names)
        if (adata[:, subset_idx]).shape[1] > 2:
            break
    subset = adata[:, subset_idx]

    internal_idx = subset_func(np.arange(subset.varm["o"].shape[0]))
    assert subset.isview
    subset.varm["o"][internal_idx] = 1
    assert not subset.isview
    assert np.all(adata.varm["o"] == orig_varm_val)

    assert init_hash == joblib.hash(adata)


@pytest.mark.parametrize('attr', ["obsm", "varm", "obsp", "varp", "layers"])
def test_view_failed_delitem(attr):
    adata = gen_adata((10, 10))
    view = adata[5:7, :][:, :5]
    adata_hash = joblib.hash(adata)
    view_hash = joblib.hash(view)

    with pytest.raises(KeyError):
        getattr(view, attr).__delitem__("not a key")

    assert view.isview
    assert adata_hash == joblib.hash(adata)
    assert view_hash == joblib.hash(view)


@pytest.mark.parametrize('attr', ["obsm", "varm", "obsp", "varp", "layers"])
def test_view_delitem(attr):
    adata = gen_adata((10, 10))
    getattr(adata, attr)["to_delete"] = np.ones((10, 10))
    assert (
        type(getattr(adata, attr)["to_delete"]) is np.ndarray
    )  # Shouldn't be a subclass, should be an ndarray
    view = adata[5:7, :][:, :5]
    adata_hash = joblib.hash(adata)
    view_hash = joblib.hash(view)

    getattr(view, attr).__delitem__("to_delete")

    assert not view.isview
    assert "to_delete" not in getattr(view, attr)
    assert "to_delete" in getattr(adata, attr)
    assert adata_hash == joblib.hash(adata)
    assert view_hash != joblib.hash(view)


def test_layers_view():
    X = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    L = np.array([[10, 11, 12], [13, 14, 15], [16, 17, 18]])
    real_adata = ad.AnnData(X)
    real_adata.layers["L"] = L
    view_adata = real_adata[1:, 1:]
    real_hash = joblib.hash(real_adata)
    view_hash = joblib.hash(view_adata)

    assert view_adata.isview

    with pytest.raises(ValueError):
        view_adata.layers["L2"] = L + 2

    assert view_adata.isview  # Failing to set layer item makes adata not view
    assert real_hash == joblib.hash(real_adata)
    assert view_hash == joblib.hash(view_adata)

    view_adata.layers["L2"] = L[1:, 1:] + 2

    assert not view_adata.isview
    assert real_hash == joblib.hash(real_adata)
    assert view_hash != joblib.hash(view_adata)


# TODO: This can be flaky. Make that stop
def test_view_of_view(matrix_type, subset_func, subset_func2):
    adata = gen_adata((30, 15), X_type=matrix_type)
    adata.raw = adata
    if subset_func is single_subset:
        pytest.xfail("Other subset generating functions have trouble with this")
    var_s1 = subset_func(adata.var_names, min_size=4)
    var_view1 = adata[:, var_s1]
    var_s2 = subset_func2(var_view1.var_names)
    var_view2 = var_view1[:, var_s2]
    assert var_view2._adata_ref is adata
    obs_s1 = subset_func(adata.obs_names, min_size=4)
    obs_view1 = adata[obs_s1, :]
    obs_s2 = subset_func2(obs_view1.obs_names)
    assert adata[obs_s1, :][:, var_s1][obs_s2, :]._adata_ref is adata

    view_of_actual_copy = (
        adata[:, var_s1].copy()[obs_s1, :].copy()[:, var_s2].copy()
    )
    view_of_view_copy = adata[:, var_s1][obs_s1, :][:, var_s2].copy()

    assert_equal(view_of_actual_copy, view_of_view_copy, exact=True)


def test_view_of_view_modification():
    adata = ad.AnnData(np.zeros((10, 10)))
    adata[0, :][:, 5:].X = np.ones(5)
    assert np.all(adata.X[0, 5:] == np.ones(5))
    adata[[1, 2], :][:, [1, 2]].X = np.ones((2, 2))
    assert np.all(adata.X[1:3, 1:3] == np.ones((2, 2)))

    adata.X = sparse.csr_matrix(adata.X)
    adata[0, :][:, 5:].X = np.ones(5) * 2
    assert np.all(asarray(adata.X)[0, 5:] == np.ones(5) * 2)
    adata[[1, 2], :][:, [1, 2]].X = np.ones((2, 2)) * 2
    assert np.all(asarray(adata.X)[1:3, 1:3] == np.ones((2, 2)) * 2)


def test_double_index(subset_func, subset_func2):
    adata = gen_adata((10, 10))
    obs_subset = subset_func(adata.obs_names)
    var_subset = subset_func2(adata.var_names)
    v1 = adata[obs_subset, var_subset]
    v2 = adata[obs_subset, :][:, var_subset]

    assert np.all(asarray(v1.X) == asarray(v2.X))
    assert np.all(v1.obs == v2.obs)
    assert np.all(v1.var == v2.var)
