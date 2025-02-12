from __future__ import annotations

import warnings
from importlib.util import find_spec

import numpy as np
import pandas as pd
import pytest
from numba.core.errors import NumbaDeprecationWarning

from anndata import AnnData, ImplicitModificationWarning, read_h5ad
from anndata.io import read_loom
from anndata.tests.helpers import gen_typed_df_t2_size

X_ = np.arange(12).reshape((3, 4))
L = np.arange(12).reshape((3, 4)) + 12


@pytest.fixture(params=[X_, None])
def X(request):
    return request.param


def test_creation(X: np.ndarray | None):
    adata = AnnData(X=X, layers=dict(L=L.copy()))

    assert list(adata.layers.keys()) == ["L"]
    assert "L" in adata.layers
    assert "X" not in adata.layers
    assert "some_other_thing" not in adata.layers
    assert (adata.layers["L"] == L).all()
    assert adata.shape == L.shape


def test_views():
    adata = AnnData(X=X_, layers=dict(L=L.copy()))
    adata_view = adata[1:, 1:]

    assert adata_view.layers.is_view
    assert adata_view.layers.parent_mapping == adata.layers

    assert adata_view.layers.keys() == adata.layers.keys()
    assert (adata_view.layers["L"] == adata.layers["L"][1:, 1:]).all()

    adata.layers["S"] = X_

    assert adata_view.layers.keys() == adata.layers.keys()
    assert (adata_view.layers["S"] == adata.layers["S"][1:, 1:]).all()

    with pytest.warns(ImplicitModificationWarning):
        adata_view.layers["T"] = X_[1:, 1:]

    assert not adata_view.layers.is_view
    assert not adata_view.is_view


@pytest.mark.parametrize(
    ("df", "homogenous", "dtype"),
    [
        (lambda: gen_typed_df_t2_size(*X_.shape), True, np.object_),
        (lambda: pd.DataFrame(X_**2), False, np.int_),
    ],
)
def test_set_dataframe(homogenous, df, dtype):
    adata = AnnData(X_)
    if homogenous:
        with pytest.warns(UserWarning, match=r"Layer 'df'.*dtype object"):
            adata.layers["df"] = df()
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            adata.layers["df"] = df()
    assert isinstance(adata.layers["df"], np.ndarray)
    assert np.issubdtype(adata.layers["df"].dtype, dtype)


def test_readwrite(X: np.ndarray | None, backing_h5ad):
    adata = AnnData(X=X, layers=dict(L=L.copy()))
    adata.write(backing_h5ad)
    adata_read = read_h5ad(backing_h5ad)

    assert adata.layers.keys() == adata_read.layers.keys()
    assert (adata.layers["L"] == adata_read.layers["L"]).all()


@pytest.mark.skipif(find_spec("loompy") is None, reason="loompy not installed")
def test_readwrite_loom(tmp_path):
    loom_path = tmp_path / "test.loom"
    adata = AnnData(X=X_, layers=dict(L=L.copy()))

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=NumbaDeprecationWarning)
        # loompy uses “is” for ints
        warnings.filterwarnings("ignore", category=SyntaxWarning)
        warnings.filterwarnings(
            "ignore",
            message=r"datetime.datetime.utcnow\(\) is deprecated",
            category=DeprecationWarning,
        )
        adata.write_loom(loom_path)
    adata_read = read_loom(loom_path, X_name="")

    assert adata.layers.keys() == adata_read.layers.keys()
    assert (adata.layers["L"] == adata_read.layers["L"]).all()


def test_backed():
    # backed mode for layers isn’t implemented, layers stay in memory
    pass


def test_copy():
    adata = AnnData(X=X_, layers=dict(L=L.copy()))
    bdata = adata.copy()
    # check that we don’t create too many references
    assert bdata._layers is bdata.layers._data
    # check that we have a copy
    adata.layers["L"] += 10
    assert np.all(adata.layers["L"] != bdata.layers["L"])  # 201


def test_shape_error():
    adata = AnnData(X=X_)
    with pytest.raises(
        ValueError,
        match=(
            r"Value passed for key 'L' is of incorrect shape\. "
            r"Values of layers must match dimensions \('obs', 'var'\) of parent\. "
            r"Value had shape \(4, 4\) while it should have had \(3, 4\)\."
        ),
    ):
        adata.layers["L"] = np.zeros((X_.shape[0] + 1, X_.shape[1]))
