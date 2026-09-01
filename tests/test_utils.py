from __future__ import annotations

from itertools import repeat

import pandas as pd
import pytest
from scipy import sparse

import anndata as ad
from anndata.tests.helpers import gen_typed_df
from anndata.utils import make_index_unique


def test_make_index_unique() -> None:
    index = pd.Index(["val", "val", pd.NA, "val-1", "val-1", pd.NA], dtype="string")
    with pytest.warns(
        UserWarning, match=r"Suffix used.*index values difficult to interpret"
    ):
        result = make_index_unique(index)
    expected = pd.Index(
        ["val", "val-2", pd.NA, "val-1", "val-1-1", pd.NA], dtype="string"
    )
    pd.testing.assert_index_equal(result, expected)
    assert result[result.notna()].is_unique


def assert_df_index_equal(df: object, expected: pd.Index) -> None:
    assert isinstance(df, pd.DataFrame)
    pd.testing.assert_index_equal(df.index, expected)


def test_adata_unique_indices() -> None:
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

    assert_df_index_equal(adata.obsm["df"], adata.obs_names)
    assert_df_index_equal(adata.varm["df"], adata.var_names)

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    assert adata.obs_names.name == "obs"
    assert adata.var_names.name == "var"

    assert len(pd.unique(adata.obs_names)) == m
    assert len(pd.unique(adata.var_names)) == n

    assert_df_index_equal(adata.obsm["df"], adata.obs_names)
    assert_df_index_equal(adata.varm["df"], adata.var_names)

    v = adata[:5, :5]

    assert v.obs_names.name == "obs"
    assert v.var_names.name == "var"

    assert_df_index_equal(v.obsm["df"], v.obs_names)
    assert_df_index_equal(v.varm["df"], v.var_names)
