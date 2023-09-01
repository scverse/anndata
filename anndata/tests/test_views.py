from copy import deepcopy
from operator import mul

import joblib
import numpy as np
from scipy import sparse
import pandas as pd
import pytest

import anndata as ad
from anndata._core.index import _normalize_index
from anndata._core.views import ArrayView, SparseCSRView, SparseCSCView
from anndata.compat import CupyCSCMatrix
from anndata.utils import asarray
from anndata.tests.helpers import (
    gen_adata,
    subset_func,
    slice_subset,
    single_subset,
    assert_equal,
    GEN_ADATA_DASK_ARGS,
    BASE_MATRIX_PARAMS,
    DASK_MATRIX_PARAMS,
    CUPY_MATRIX_PARAMS,
)
from dask.base import tokenize, normalize_token


# ------------------------------------------------------------------------------
# Some test data
# ------------------------------------------------------------------------------

# data matrix of shape n_obs x n_vars
X_list = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
# annotation of observations / rows
obs_dict = dict(
    row_names=["name1", "name2", "name3"],  # row annotation
    oanno1=["cat1", "cat2", "cat2"],  # categorical annotation
    oanno2=["o1", "o2", "o3"],  # string annotation
    oanno3=[2.1, 2.2, 2.3],  # float annotation
)
# annotation of variables / columns
var_dict = dict(vanno1=[3.1, 3.2, 3.3])
# unstructured annotation
uns_dict = dict(oanno1_colors=["#000000", "#FFFFFF"], uns2=["some annotation"])

subset_func2 = subset_func


class NDArraySubclass(np.ndarray):
    def view(self, dtype=None, typ=None):
        return self


@pytest.fixture
def adata():
    adata = ad.AnnData(np.zeros((100, 100)))
    adata.obsm["o"] = np.zeros((100, 50))
    adata.varm["o"] = np.zeros((100, 50))
    return adata


@pytest.fixture(
    params=BASE_MATRIX_PARAMS + DASK_MATRIX_PARAMS + CUPY_MATRIX_PARAMS,
)
def matrix_type(request):
    return request.param


@pytest.fixture(params=BASE_MATRIX_PARAMS + DASK_MATRIX_PARAMS)
def matrix_type_no_gpu(request):
    return request.param


@pytest.fixture(params=BASE_MATRIX_PARAMS)
def matrix_type_base(request):
    return request.param


@pytest.fixture(params=["layers", "obsm", "varm"])
def mapping_name(request):
    return request.param


# ------------------------------------------------------------------------------
# The test functions
# ------------------------------------------------------------------------------


def test_views():
    X = np.array(X_list, dtype="int32")
    adata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)

    assert adata[:, 0].is_view
    assert adata[:, 0].X.tolist() == np.reshape([1, 4, 7], (3, 1)).tolist()

    adata[:2, 0].X = [0, 0]

    assert adata[:, 0].X.tolist() == np.reshape([0, 0, 7], (3, 1)).tolist()

    adata_subset = adata[:2, [0, 1]]

    assert adata_subset.is_view
    # now transition to actual object
    adata_subset.obs["foo"] = range(2)
    assert not adata_subset.is_view

    assert adata_subset.obs["foo"].tolist() == list(range(2))


def test_view_subset_shapes():
    adata = gen_adata((20, 10), **GEN_ADATA_DASK_ARGS)

    view = adata[:, ::2]
    assert view.var.shape == (5, 8)
    assert {k: v.shape[0] for k, v in view.varm.items()} == {k: 5 for k in view.varm}


def test_modify_view_component(matrix_type, mapping_name):
    adata = ad.AnnData(
        np.zeros((10, 10)),
        **{mapping_name: dict(m=matrix_type(asarray(sparse.random(10, 10))))},
    )
    init_hash = joblib.hash(adata)

    subset = adata[:5, :][:, :5]
    assert subset.is_view
    m = getattr(subset, mapping_name)["m"]
    m[0, 0] = 100
    assert not subset.is_view
    assert getattr(subset, mapping_name)["m"][0, 0] == 100

    assert init_hash == joblib.hash(adata)


@pytest.mark.parametrize("attr", ["obsm", "varm"])
def test_set_obsm_key(adata, attr):
    init_hash = joblib.hash(adata)

    orig_val = getattr(adata, attr)["o"].copy()
    subset = adata[:50] if attr == "obsm" else adata[:, :50]

    assert subset.is_view

    with pytest.warns(ad.ImplicitModificationWarning, match=rf".*\.{attr}\['o'\].*"):
        getattr(subset, attr)["o"] = new_val = np.ones((50, 20))

    assert not subset.is_view
    assert np.all(getattr(adata, attr)["o"] == orig_val)
    assert np.any(getattr(subset, attr)["o"] == new_val)

    assert init_hash == joblib.hash(adata)


def test_set_obs(adata, subset_func):
    init_hash = joblib.hash(adata)

    subset = adata[subset_func(adata.obs_names), :]

    new_obs = pd.DataFrame(
        dict(a=np.ones(subset.n_obs), b=np.ones(subset.n_obs)),
        index=subset.obs_names,
    )

    assert subset.is_view
    subset.obs = new_obs
    assert not subset.is_view
    assert np.all(subset.obs == new_obs)

    assert joblib.hash(adata) == init_hash


def test_set_var(adata, subset_func):
    init_hash = joblib.hash(adata)

    subset = adata[:, subset_func(adata.var_names)]

    new_var = pd.DataFrame(
        dict(a=np.ones(subset.n_vars), b=np.ones(subset.n_vars)),
        index=subset.var_names,
    )

    assert subset.is_view
    subset.var = new_var
    assert not subset.is_view
    assert np.all(subset.var == new_var)

    assert joblib.hash(adata) == init_hash


def test_drop_obs_column():
    adata = ad.AnnData(np.array(X_list, dtype="int32"), obs=obs_dict)

    subset = adata[:2]
    assert subset.is_view
    # returns a copy of obs
    assert subset.obs.drop(columns=["oanno1"]).columns.tolist() == ["oanno2", "oanno3"]
    assert subset.is_view
    # would modify obs, so it should actualize subset and not modify adata
    subset.obs.drop(columns=["oanno1"], inplace=True)
    assert not subset.is_view
    assert subset.obs.columns.tolist() == ["oanno2", "oanno3"]

    assert adata.obs.columns.tolist() == ["oanno1", "oanno2", "oanno3"]


def test_set_obsm(adata):
    init_hash = joblib.hash(adata)

    dim0_size = np.random.randint(2, adata.shape[0] - 1)
    dim1_size = np.random.randint(1, 99)
    orig_obsm_val = adata.obsm["o"].copy()
    subset_idx = np.random.choice(adata.obs_names, dim0_size, replace=False)

    subset = adata[subset_idx, :]
    assert subset.is_view
    subset.obsm = dict(o=np.ones((dim0_size, dim1_size)))
    assert not subset.is_view
    assert np.all(orig_obsm_val == adata.obsm["o"])  # Checking for mutation
    assert np.all(subset.obsm["o"] == np.ones((dim0_size, dim1_size)))

    subset = adata[subset_idx, :]
    subset_hash = joblib.hash(subset)
    with pytest.raises(ValueError):
        subset.obsm = dict(o=np.ones((dim0_size + 1, dim1_size)))
    with pytest.raises(ValueError):
        subset.varm = dict(o=np.ones((dim0_size - 1, dim1_size)))
    assert subset_hash == joblib.hash(subset)

    # Only modification have been made to a view
    assert init_hash == joblib.hash(adata)


def test_set_varm(adata):
    init_hash = joblib.hash(adata)

    dim0_size = np.random.randint(2, adata.shape[1] - 1)
    dim1_size = np.random.randint(1, 99)
    orig_varm_val = adata.varm["o"].copy()
    subset_idx = np.random.choice(adata.var_names, dim0_size, replace=False)

    subset = adata[:, subset_idx]
    assert subset.is_view
    subset.varm = dict(o=np.ones((dim0_size, dim1_size)))
    assert not subset.is_view
    assert np.all(orig_varm_val == adata.varm["o"])  # Checking for mutation
    assert np.all(subset.varm["o"] == np.ones((dim0_size, dim1_size)))

    subset = adata[:, subset_idx]
    subset_hash = joblib.hash(subset)
    with pytest.raises(ValueError):
        subset.varm = dict(o=np.ones((dim0_size + 1, dim1_size)))
    with pytest.raises(ValueError):
        subset.varm = dict(o=np.ones((dim0_size - 1, dim1_size)))
    # subset should not be changed by failed setting
    assert subset_hash == joblib.hash(subset)
    assert init_hash == joblib.hash(adata)


# TODO: Determine if this is the intended behavior,
#       or just the behaviour we’ve had for a while
def test_not_set_subset_X(matrix_type_base, subset_func):
    adata = ad.AnnData(matrix_type_base(asarray(sparse.random(20, 20))))
    init_hash = joblib.hash(adata)
    orig_X_val = adata.X.copy()
    while True:
        subset_idx = slice_subset(adata.obs_names)
        if len(adata[subset_idx, :]) > 2:
            break
    subset = adata[subset_idx, :]

    subset = adata[:, subset_idx]

    internal_idx = _normalize_index(
        subset_func(np.arange(subset.X.shape[1])), subset.var_names
    )
    assert subset.is_view
    subset.X[:, internal_idx] = 1
    assert not subset.is_view
    assert not np.any(asarray(adata.X != orig_X_val))

    assert init_hash == joblib.hash(adata)


@normalize_token.register(ad.AnnData)
def tokenize_anndata(adata: ad.AnnData):
    res = []
    if adata.X is not None:
        res.append(tokenize(adata.X))
    res.extend([tokenize(adata.obs), tokenize(adata.var)])
    for attr in ["obsm", "varm", "obsp", "varp", "layers"]:
        elem = getattr(adata, attr)
        res.append(tokenize(list(elem.items())))
    res.append(joblib.hash(adata.uns))
    if adata.raw is not None:
        res.append(tokenize(adata.raw.to_adata()))
    return tuple(res)


# TODO: Determine if this is the intended behavior,
#       or just the behaviour we’ve had for a while
def test_not_set_subset_X_dask(matrix_type_no_gpu, subset_func):
    adata = ad.AnnData(matrix_type_no_gpu(asarray(sparse.random(20, 20))))
    init_hash = tokenize(adata)
    orig_X_val = adata.X.copy()
    while True:
        subset_idx = slice_subset(adata.obs_names)
        if len(adata[subset_idx, :]) > 2:
            break
    subset = adata[subset_idx, :]

    subset = adata[:, subset_idx]

    internal_idx = _normalize_index(
        subset_func(np.arange(subset.X.shape[1])), subset.var_names
    )
    assert subset.is_view
    subset.X[:, internal_idx] = 1
    assert not subset.is_view
    assert not np.any(asarray(adata.X != orig_X_val))

    assert init_hash == tokenize(adata)


def test_set_scalar_subset_X(matrix_type, subset_func):
    adata = ad.AnnData(matrix_type(np.zeros((10, 10))))
    orig_X_val = adata.X.copy()
    subset_idx = slice_subset(adata.obs_names)

    adata_subset = adata[subset_idx, :]

    adata_subset.X = 1

    assert adata_subset.is_view
    assert np.all(asarray(adata[subset_idx, :].X) == 1)
    if isinstance(adata.X, CupyCSCMatrix):
        # Comparison broken for CSC matrices
        # https://github.com/cupy/cupy/issues/7757
        assert asarray((orig_X_val.tocsr() != adata.X.tocsr())).sum() == mul(
            *adata_subset.shape
        )
    else:
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

    internal_idx = _normalize_index(
        subset_func(np.arange(subset.obsm["o"].shape[0])), subset.obs_names
    )

    assert subset.is_view
    subset.obsm["o"][internal_idx] = 1
    assert not subset.is_view
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

    internal_idx = _normalize_index(
        subset_func(np.arange(subset.varm["o"].shape[0])), subset.var_names
    )

    assert subset.is_view
    subset.varm["o"][internal_idx] = 1
    assert not subset.is_view
    assert np.all(adata.varm["o"] == orig_varm_val)

    assert init_hash == joblib.hash(adata)


@pytest.mark.parametrize("attr", ["obsm", "varm", "obsp", "varp", "layers"])
def test_view_failed_delitem(attr):
    adata = gen_adata((10, 10), **GEN_ADATA_DASK_ARGS)
    view = adata[5:7, :][:, :5]
    adata_hash = joblib.hash(adata)
    view_hash = joblib.hash(view)

    with pytest.raises(KeyError):
        getattr(view, attr).__delitem__("not a key")

    assert view.is_view
    assert adata_hash == joblib.hash(adata)
    assert view_hash == joblib.hash(view)


@pytest.mark.parametrize("attr", ["obsm", "varm", "obsp", "varp", "layers"])
def test_view_delitem(attr):
    adata = gen_adata((10, 10), **GEN_ADATA_DASK_ARGS)
    getattr(adata, attr)["to_delete"] = np.ones((10, 10))
    # Shouldn’t be a subclass, should be an ndarray
    assert type(getattr(adata, attr)["to_delete"]) is np.ndarray
    view = adata[5:7, :][:, :5]
    adata_hash = joblib.hash(adata)
    view_hash = joblib.hash(view)

    with pytest.warns(
        ad.ImplicitModificationWarning, match=rf".*\.{attr}\['to_delete'\].*"
    ):
        getattr(view, attr).__delitem__("to_delete")

    assert not view.is_view
    assert "to_delete" not in getattr(view, attr)
    assert "to_delete" in getattr(adata, attr)
    assert adata_hash == joblib.hash(adata)
    assert view_hash != joblib.hash(view)


@pytest.mark.parametrize(
    "attr", ["X", "obs", "var", "obsm", "varm", "obsp", "varp", "layers", "uns"]
)
def test_view_delattr(attr, subset_func):
    base = gen_adata((10, 10), **GEN_ADATA_DASK_ARGS)
    orig_hash = tokenize(base)
    subset = base[subset_func(base.obs_names), subset_func(base.var_names)]
    empty = ad.AnnData(obs=subset.obs[[]], var=subset.var[[]])

    delattr(subset, attr)

    assert not subset.is_view
    # Should now have same value as default
    assert_equal(getattr(subset, attr), getattr(empty, attr))
    assert orig_hash == tokenize(base)  # Original should not be modified


@pytest.mark.parametrize(
    "attr", ["obs", "var", "obsm", "varm", "obsp", "varp", "layers", "uns"]
)
def test_view_setattr_machinery(attr, subset_func, subset_func2):
    # Tests that setting attributes on a view doesn't mess anything up too bad
    adata = gen_adata((10, 10), **GEN_ADATA_DASK_ARGS)
    view = adata[subset_func(adata.obs_names), subset_func2(adata.var_names)]

    actual = view.copy()
    setattr(view, attr, getattr(actual, attr))
    assert_equal(actual, view, exact=True)


def test_layers_view():
    X = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    L = np.array([[10, 11, 12], [13, 14, 15], [16, 17, 18]])
    real_adata = ad.AnnData(X)
    real_adata.layers["L"] = L
    view_adata = real_adata[1:, 1:]
    real_hash = joblib.hash(real_adata)
    view_hash = joblib.hash(view_adata)

    assert view_adata.is_view

    with pytest.raises(ValueError):
        view_adata.layers["L2"] = L + 2

    assert view_adata.is_view  # Failing to set layer item makes adata not view
    assert real_hash == joblib.hash(real_adata)
    assert view_hash == joblib.hash(view_adata)

    view_adata.layers["L2"] = L[1:, 1:] + 2

    assert not view_adata.is_view
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

    view_of_actual_copy = adata[:, var_s1].copy()[obs_s1, :].copy()[:, var_s2].copy()

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
    adata = gen_adata((10, 10), **GEN_ADATA_DASK_ARGS)
    obs_subset = subset_func(adata.obs_names)
    var_subset = subset_func2(adata.var_names)
    v1 = adata[obs_subset, var_subset]
    v2 = adata[obs_subset, :][:, var_subset]

    assert np.all(asarray(v1.X) == asarray(v2.X))
    assert np.all(v1.obs == v2.obs)
    assert np.all(v1.var == v2.var)


def test_view_retains_ndarray_subclass():
    adata = ad.AnnData(np.zeros((10, 10)))
    adata.obsm["foo"] = np.zeros((10, 5)).view(NDArraySubclass)

    view = adata[:5, :]

    assert isinstance(view.obsm["foo"], NDArraySubclass)
    assert view.obsm["foo"].shape == (5, 5)


def test_modify_uns_in_copy():
    # https://github.com/scverse/anndata/issues/571
    adata = ad.AnnData(np.ones((5, 5)), uns={"parent": {"key": "value"}})
    adata_copy = adata[:3].copy()
    adata_copy.uns["parent"]["key"] = "new_value"
    assert adata.uns["parent"]["key"] != adata_copy.uns["parent"]["key"]


@pytest.mark.parametrize("index", [-101, 100, (slice(None), -101), (slice(None), 100)])
def test_invalid_scalar_index(adata, index):
    # https://github.com/scverse/anndata/issues/619
    with pytest.raises(IndexError, match=r".*index.* out of range\."):
        _ = adata[index]


@pytest.mark.parametrize("obs", [False, True])
@pytest.mark.parametrize("index", [-100, -50, -1])
def test_negative_scalar_index(adata, index: int, obs: bool):
    pos_index = index + (adata.n_obs if obs else adata.n_vars)

    if obs:
        adata_pos_subset = adata[pos_index]
        adata_neg_subset = adata[index]
    else:
        adata_pos_subset = adata[:, pos_index]
        adata_neg_subset = adata[:, index]

    np.testing.assert_array_equal(
        adata_pos_subset.obs_names, adata_neg_subset.obs_names
    )
    np.testing.assert_array_equal(
        adata_pos_subset.var_names, adata_neg_subset.var_names
    )


def test_viewness_propagation_nan():
    """Regression test for https://github.com/scverse/anndata/issues/239"""
    adata = ad.AnnData(np.random.random((10, 10)))
    adata = adata[:, [0, 2, 4]]
    v = adata.X.var(axis=0)
    assert not isinstance(v, ArrayView), type(v).mro()
    # this used to break
    v[np.isnan(v)] = 0


def test_viewness_propagation_allclose(adata):
    """Regression test for https://github.com/scverse/anndata/issues/191"""
    adata.varm["o"][4:10] = np.tile(np.nan, (10 - 4, adata.varm["o"].shape[1]))
    a = adata[:50].copy()
    b = adata[:50]

    # .copy() turns view to ndarray, so this was fine:
    assert np.allclose(a.varm["o"], b.varm["o"].copy(), equal_nan=True)
    # Next line triggered the mutation:
    assert np.allclose(a.varm["o"], b.varm["o"], equal_nan=True)
    # Showing that the mutation didn’t happen:
    assert np.allclose(a.varm["o"], b.varm["o"].copy(), equal_nan=True)


@pytest.mark.parametrize("spmat", [sparse.csr_matrix, sparse.csc_matrix])
def test_deepcopy_subset(adata, spmat: type):
    adata.obsp["arr"] = np.zeros((adata.n_obs, adata.n_obs))
    adata.obsp["spmat"] = spmat((adata.n_obs, adata.n_obs))

    adata = deepcopy(adata[:10].copy())

    assert isinstance(adata.obsp["arr"], np.ndarray)
    assert not isinstance(adata.obsp["arr"], ArrayView)
    np.testing.assert_array_equal(adata.obsp["arr"].shape, (10, 10))

    assert isinstance(adata.obsp["spmat"], spmat)
    assert not isinstance(
        adata.obsp["spmat"],
        SparseCSRView if spmat is sparse.csr_matrix else SparseCSCView,
    )
    np.testing.assert_array_equal(adata.obsp["spmat"].shape, (10, 10))


# https://github.com/scverse/anndata/issues/680
@pytest.mark.parametrize("array_type", [asarray, sparse.csr_matrix, sparse.csc_matrix])
@pytest.mark.parametrize("attr", ["X", "layers", "obsm", "varm", "obsp", "varp"])
def test_view_mixin_copies_data(adata, array_type: type, attr):
    N = 100
    adata = ad.AnnData(
        obs=pd.DataFrame(index=np.arange(N)), var=pd.DataFrame(index=np.arange(N))
    )

    X = array_type(sparse.eye(N, N).multiply(np.arange(1, N + 1)))
    if attr == "X":
        adata.X = X
    else:
        getattr(adata, attr)["arr"] = X

    view = adata[:50]

    if attr == "X":
        arr_view = view.X
    else:
        arr_view = getattr(view, attr)["arr"]

    arr_view_copy = arr_view.copy()

    if sparse.issparse(X):
        assert not np.shares_memory(arr_view.indices, arr_view_copy.indices)
        assert not np.shares_memory(arr_view.indptr, arr_view_copy.indptr)
        assert not np.shares_memory(arr_view.data, arr_view_copy.data)

        arr_view_copy.data[0] = -5
        assert not np.array_equal(arr_view_copy.data, arr_view.data)
    else:
        assert not np.shares_memory(arr_view, arr_view_copy)

        arr_view_copy[0, 0] = -5
        assert not np.array_equal(arr_view_copy, arr_view)


def test_copy_X_dtype():
    adata = ad.AnnData(sparse.eye(50, dtype=np.float64, format="csr"))
    adata_c = adata[::2].copy()
    assert adata_c.X.dtype == adata.X.dtype
