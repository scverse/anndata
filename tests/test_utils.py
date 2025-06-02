from __future__ import annotations

from itertools import repeat

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import anndata as ad
from anndata.tests.helpers import gen_typed_df
from anndata.utils import adapt_vars_like, make_index_unique


def test_make_index_unique():
    index = pd.Index(["val", "val", "val-1", "val-1"])
    with pytest.warns(UserWarning):
        result = make_index_unique(index)
    expected = pd.Index(["val", "val-2", "val-1", "val-1-1"])
    assert list(expected) == list(result)
    assert result.is_unique


def test_adata_unique_indices():
    m, n = (10, 20)
    obs_index = pd.Index(repeat("a", m), name="obs")
    var_index = pd.Index(repeat("b", n), name="var")

    adata = ad.AnnData(
        X=sparse.random(m, n, format="csr"),
        obs=gen_typed_df(m, index=obs_index),
        var=gen_typed_df(n, index=var_index),
        obsm={"df": gen_typed_df(m, index=obs_index)},
        varm={"df": gen_typed_df(n, index=var_index)},
    )

    pd.testing.assert_index_equal(adata.obsm["df"].index, adata.obs_names)
    pd.testing.assert_index_equal(adata.varm["df"].index, adata.var_names)

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    assert adata.obs_names.name == "obs"
    assert adata.var_names.name == "var"

    assert len(pd.unique(adata.obs_names)) == m
    assert len(pd.unique(adata.var_names)) == n

    pd.testing.assert_index_equal(adata.obsm["df"].index, adata.obs_names)
    pd.testing.assert_index_equal(adata.varm["df"].index, adata.var_names)

    v = adata[:5, :5]

    assert v.obs_names.name == "obs"
    assert v.var_names.name == "var"

    pd.testing.assert_index_equal(v.obsm["df"].index, v.obs_names)
    pd.testing.assert_index_equal(v.varm["df"].index, v.var_names)


@pytest.mark.parametrize(
    ("source", "target", "expected_X"),
    [
        pytest.param(
            ad.AnnData(X=np.ones((1, 3)), var=pd.DataFrame(index=["a", "b", "c"])),
            ad.AnnData(
                X=np.array([[1, 2, 3]]), var=pd.DataFrame(index=["a", "b", "c"])
            ),
            np.array([[1, 2, 3]]),
            id="exact_match",
        ),
        pytest.param(
            ad.AnnData(X=np.ones((1, 3)), var=pd.DataFrame(index=["a", "b", "c"])),
            ad.AnnData(
                X=np.array([[3, 2, 1]]), var=pd.DataFrame(index=["c", "b", "a"])
            ),
            np.array([[1, 2, 3]]),
            id="different_order",
        ),
    ],
)
def test_adapt_vars(source, target, expected_X):
    output = adapt_vars_like(source, target)
    np.testing.assert_array_equal(output.X, expected_X)
    assert list(output.var_names) == list(source.var_names)


@pytest.mark.parametrize(
    ("source", "target", "fill_value", "expected_X"),
    [
        pytest.param(
            ad.AnnData(X=np.ones((1, 2)), var=pd.DataFrame(index=["g1", "g2"])),
            ad.AnnData(X=np.array([[7, 8]]), var=pd.DataFrame(index=["g3", "g4"])),
            0.5,
            np.array([[0.5, 0.5]]),
            id="no_shared_genes",
        ),
        pytest.param(
            ad.AnnData(X=np.ones((1, 3)), var=pd.DataFrame(index=["g1", "g2", "g3"])),
            ad.AnnData(X=np.array([[1, 3]]), var=pd.DataFrame(index=["g1", "g3"])),
            -1,
            np.array([[1, -1, 3]]),
            id="missing_genes",
        ),
    ],
)
def test_adapt_vars_with_fill_value(source, target, fill_value, expected_X):
    output = adapt_vars_like(source, target, fill_value=fill_value)
    np.testing.assert_array_equal(output.X, expected_X)
    assert list(output.var_names) == list(source.var_names)


def test_adapt_vars_target_X_none():
    source = ad.AnnData(
        X=np.ones((2, 2)),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    target = ad.AnnData(
        X=None,
        var=pd.DataFrame(index=["g2", "g3"]),
        obs=pd.DataFrame(index=["cell1", "cell2"]),
    )
    output = adapt_vars_like(source, target, fill_value=-1)
    assert output.X is None
    assert list(output.var_names) == list(source.var_names)
