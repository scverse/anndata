"""Tests for the attribute .X"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import anndata as ad
from anndata import AnnData
from anndata._warnings import ImplicitModificationWarning
from anndata.tests.helpers import GEN_ADATA_NO_XARRAY_ARGS, assert_equal, gen_adata
from anndata.utils import asarray

UNLABELLED_ARRAY_TYPES = [
    pytest.param(sparse.csr_matrix, id="csr"),
    pytest.param(sparse.csc_matrix, id="csc"),
    pytest.param(sparse.csr_array, id="csr_array"),
    pytest.param(sparse.csc_array, id="csc_array"),
    pytest.param(asarray, id="ndarray"),
]
SINGULAR_SHAPES = [
    pytest.param(shape, id=str(shape)) for shape in [(1, 10), (10, 1), (1, 1)]
]


@pytest.mark.parametrize("shape", SINGULAR_SHAPES)
@pytest.mark.parametrize("orig_array_type", UNLABELLED_ARRAY_TYPES)
@pytest.mark.parametrize("new_array_type", UNLABELLED_ARRAY_TYPES)
def test_setter_singular_dim(shape, orig_array_type, new_array_type):
    # https://github.com/scverse/anndata/issues/500
    adata = gen_adata(shape, X_type=orig_array_type)
    to_assign = new_array_type(np.ones(shape))
    adata.X = to_assign
    np.testing.assert_equal(asarray(adata.X), 1)
    assert isinstance(adata.X, type(to_assign))


def test_repeat_indices_view():
    adata = gen_adata((10, 10), X_type=np.asarray)
    subset = adata[[0, 0, 1, 1], :]
    mat = np.array([np.ones(adata.shape[1]) * i for i in range(4)])
    with pytest.warns(
        FutureWarning,
        match=r"You are attempting to set `X` to a matrix on a view which has non-unique indices",
    ):
        subset.X = mat


@pytest.mark.parametrize("orig_array_type", UNLABELLED_ARRAY_TYPES)
@pytest.mark.parametrize("new_array_type", UNLABELLED_ARRAY_TYPES)
def test_setter_view(orig_array_type, new_array_type):
    adata = gen_adata((10, 10), X_type=orig_array_type)
    orig_X = adata.X
    to_assign = new_array_type(np.ones((9, 9)))
    if isinstance(orig_X, np.ndarray) and sparse.issparse(to_assign):
        # https://github.com/scverse/anndata/issues/500
        pytest.xfail("Cannot set a dense array with a sparse array")
    view = adata[:9, :9]
    view.X = to_assign
    np.testing.assert_equal(asarray(view.X), np.ones((9, 9)))
    assert isinstance(view.X, type(orig_X))


###############################
# Tests for `adata.X is None` #
###############################


def test_set_x_is_none():
    # test setter and getter
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]), dict(o1=[1, 2], o2=[3, 4]))
    adata.X = None
    assert adata.X is None


def test_del_set_equiv_X():
    """Tests that `del adata.X` is equivalent to `adata.X = None`"""
    # test setter and deleter
    orig = gen_adata((10, 10))
    copy = orig.copy()

    del orig.X
    copy.X = None

    assert orig.X is None
    assert_equal(orig, copy)

    # Check that deleting again is still fine
    del orig.X
    assert orig.X is None


@pytest.mark.parametrize(
    ("obs", "var", "shape_expected"),
    [
        pytest.param(dict(obs_names=["1", "2"]), None, (2, 0), id="obs"),
        pytest.param(None, dict(var_names=["a", "b"]), (0, 2), id="var"),
        pytest.param(
            dict(obs_names=["1", "2", "3"]),
            dict(var_names=["a", "b"]),
            (3, 2),
            id="both",
        ),
    ],
)
def test_init_x_as_none_shape_from_obs_var(obs, var, shape_expected):
    adata = AnnData(None, obs, var)
    assert adata.X is None
    assert adata.shape == shape_expected


def test_init_x_as_none_explicit_shape():
    shape = (3, 5)
    adata = AnnData(None, uns=dict(test=np.array((3, 3))), shape=shape)
    assert adata.X is None
    assert adata.shape == shape


@pytest.mark.parametrize("shape", [*SINGULAR_SHAPES, pytest.param((5, 3), id="(5, 3)")])
def test_transpose_with_X_as_none(shape):
    adata = gen_adata(shape, X_type=lambda x: None)
    adataT = adata.transpose()
    assert_equal(adataT.shape, shape[::-1])
    assert_equal(adataT.obsp.keys(), adata.varp.keys())
    assert_equal(adataT.T, adata)


def test_copy():
    adata = AnnData(
        None,
        obs=pd.DataFrame(index=[f"cell{i:03}" for i in range(100)]),
        var=pd.DataFrame(index=[f"gene{i:03}" for i in range(200)]),
    )
    assert_equal(adata.copy(), adata)


def test_copy_view():
    adata = AnnData(
        None,
        obs=pd.DataFrame(index=[f"cell{i:03}" for i in range(100)]),
        var=pd.DataFrame(index=[f"gene{i:03}" for i in range(200)]),
    )
    v = adata[::-2, ::-2]
    assert_equal(v.copy(), v)


############
# IO tests #
############


def test_io_missing_X(tmp_path, diskfmt):
    file_pth = tmp_path / f"x_none_adata.{diskfmt}"
    write = lambda obj, pth: getattr(obj, f"write_{diskfmt}")(pth)
    read = lambda pth: getattr(ad, f"read_{diskfmt}")(pth)

    adata = gen_adata((20, 30), **GEN_ADATA_NO_XARRAY_ARGS)
    del adata.X

    write(adata, file_pth)
    from_disk = read(file_pth)

    assert_equal(from_disk, adata)


def test_set_dense_x_view_from_sparse():
    x = np.zeros((100, 30))
    x1 = np.ones((100, 30))
    orig = ad.AnnData(x)
    view = orig[:30]
    with (
        pytest.warns(
            UserWarning,
            match=r"Trying to set a dense array with a sparse array on a view",
        ),
        pytest.warns(
            ImplicitModificationWarning, match=r"Modifying `X` on a view results"
        ),
    ):
        view.X = sparse.csr_matrix(x1[:30])
    assert_equal(view.X, x1[:30])
    assert_equal(orig.X[:30], x1[:30])  # change propagates through
    assert_equal(orig.X[30:], x[30:])  # change propagates through


def test_fail_on_non_csr_csc_matrix():
    X = sparse.eye(100, format="coo")
    with pytest.raises(
        ValueError,
        match=r"Only CSR and CSC.*",
    ):
        ad.AnnData(X=X)
