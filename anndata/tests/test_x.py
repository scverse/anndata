"""Tests for the attribute .X"""
import numpy as np
from scipy import sparse

from anndata import AnnData
from anndata.utils import asarray

import pytest

from anndata.tests.helpers import gen_adata, assert_equal

UNLABELLED_ARRAY_TYPES = [
    pytest.param(sparse.csr_matrix, id="csr"),
    pytest.param(sparse.csc_matrix, id="csc"),
    pytest.param(asarray, id="ndarray"),
]
SINGULAR_SHAPES = [
    pytest.param(shape, id=str(shape)) for shape in [(1, 10), (10, 1), (1, 1)]
]


@pytest.mark.parametrize("shape", SINGULAR_SHAPES)
@pytest.mark.parametrize("orig_array_type", UNLABELLED_ARRAY_TYPES)
@pytest.mark.parametrize("new_array_type", UNLABELLED_ARRAY_TYPES)
def test_setter_singular_dim(shape, orig_array_type, new_array_type):
    # https://github.com/theislab/anndata/issues/500
    adata = gen_adata(shape, X_type=orig_array_type)
    adata.X = new_array_type(np.ones(shape))
    np.testing.assert_equal(asarray(adata.X), 1)


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


def test_init_X_as_none():
    # test initialiser
    shape = (3, 5)
    adata = AnnData(None, uns=dict(test=np.array((3, 3))), shape=shape)
    assert adata.X is None
    assert adata.shape == shape


@pytest.mark.parametrize("shape", SINGULAR_SHAPES + [pytest.param((5, 3), id="(5, 3)")])
def test_transpose_with_X_as_none(shape):
    adata = gen_adata(shape, X_type=lambda x: None)
    adataT = adata.transpose()
    assert_equal(adataT.shape, shape[::-1])
    assert_equal(adataT.obsp.keys(), adata.varp.keys())
    assert_equal(adataT.T, adata)
