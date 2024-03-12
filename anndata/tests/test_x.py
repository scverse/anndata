"""Tests for the attribute .X"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import anndata as ad
from anndata import AnnData
from anndata.tests.helpers import assert_equal, gen_adata
from anndata.utils import asarray

UNLABELLED_ARRAY_TYPES = [
    pytest.param(sparse.csr_matrix, id="csr"),
    pytest.param(sparse.csc_matrix, id="csc"),
    pytest.param(asarray, id="ndarray"),
]
SINGULAR_SHAPES = [
    pytest.param(shape, id=str(shape)) for shape in [(1, 10), (10, 1), (1, 1)]
]


@pytest.fixture(params=["h5ad", "zarr"])
def diskfmt(request):
    return request.param


@pytest.mark.parametrize("shape", SINGULAR_SHAPES)
@pytest.mark.parametrize("orig_array_type", UNLABELLED_ARRAY_TYPES)
@pytest.mark.parametrize("new_array_type", UNLABELLED_ARRAY_TYPES)
def test_setter_singular_dim(shape, orig_array_type, new_array_type):
    # https://github.com/scverse/anndata/issues/500
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


@pytest.mark.parametrize("shape", SINGULAR_SHAPES + [pytest.param((5, 3), id="(5, 3)")])
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

    adata = gen_adata((20, 30))
    del adata.X

    write(adata, file_pth)
    from_disk = read(file_pth)

    assert_equal(from_disk, adata)
