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


def test_adapt_vars_exact_match():
    # Test that adapt_vars_like works when the source and target have the same var names
    # and the same number of variables
    source = ad.AnnData(X=np.ones((1, 3)), var=pd.DataFrame(index=["a", "b", "c"]))
    target = ad.AnnData(
        X=np.array([[1, 2, 3]]), var=pd.DataFrame(index=["a", "b", "c"])
    )
    output = adapt_vars_like(source, target)
    np.testing.assert_array_equal(output.X, target.X)
    assert (output.var.index == target.var.index).all()


def test_adapt_vars_different_order():
    # Test that adapt_vars_like works when the source and target have the same var names
    # but in a different order
    source = ad.AnnData(X=np.ones((1, 3)), var=pd.DataFrame(index=["a", "b", "c"]))
    target = ad.AnnData(
        X=np.array([[3, 2, 1]]), var=pd.DataFrame(index=["c", "b", "a"])
    )
    output = adapt_vars_like(source, target)
    np.testing.assert_array_equal(output.X, [[1, 2, 3]])


def test_adapt_vars_none_X_raises():
    source = ad.AnnData(X=np.ones((1, 2)), var=pd.DataFrame(index=["g1", "g2"]))
    target = ad.AnnData(X=None, var=pd.DataFrame(index=["g1", "g2"]))
    with pytest.raises(ValueError, match="target.X is None"):
        adapt_vars_like(source, target)


def test_adapt_vars_no_shared_genes():
    source = ad.AnnData(X=np.ones((1, 2)), var=pd.DataFrame(index=["g1", "g2"]))
    target = ad.AnnData(X=np.array([[7, 8]]), var=pd.DataFrame(index=["g3", "g4"]))
    output = adapt_vars_like(source, target, fill_value=0.5)
    np.testing.assert_array_equal(output.X, [[0.5, 0.5]])


def test_adapt_vars_missing_genes():
    source = ad.AnnData(X=np.ones((1, 3)), var=pd.DataFrame(index=["g1", "g2", "g3"]))
    target = ad.AnnData(X=np.array([[1, 3]]), var=pd.DataFrame(index=["g1", "g3"]))
    output = adapt_vars_like(source, target, fill_value=-1)
    np.testing.assert_array_equal(output.X, [[1, -1, 3]])
