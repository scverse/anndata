import pandas as pd
from scipy import sparse
from itertools import repeat
import pytest

import anndata as ad
from anndata.utils import make_index_unique, concatenate_uns
from anndata.tests.helpers import gen_typed_df


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

def test_concatenate_uns():
    visium_1_uns = {
        "spatial": {
            "library1": {
                "images": {"hires": 0, "lowres": 1},
                "scalefactors": {"hires_scale": 0, "lowres_scale": 1},
            }
        },
        "pca": 0,
    }
    visium_2_uns = {
        "spatial": {
            "library2": {"images": {"hires": 2}, "scalefactors": {"hires_scale": 2}}
        },
        "pca": 2,
    }
    chromium_1_uns = {"pca": 3}

    all_adata_uns = [visium_1_uns, visium_2_uns, chromium_1_uns]

    uns = concatenate_uns(all_adata_uns, merge_uns="if_clean")
    assert "pca" not in uns.keys()

    uns = concatenate_uns(all_adata_uns, merge_uns="always")
    assert uns["pca"] == 3
    assert len(uns["spatial"].keys()) == 2

    uns = concatenate_uns(all_adata_uns, merge_uns="never")
    assert isinstance(uns, dict)
    assert len(uns) == 0