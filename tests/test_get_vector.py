from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import anndata as ad


def test_amgibuous_keys():
    """Tests that an error is raised if obs_vector or var_vector is ambiguous."""
    var_keys = ["The", "quick", "brown", "fox", "jumps", "over", "the", "lazy", "dog"]
    obs_keys = [
        "Lorem",
        "ipsum",
        "dolor",
        "sit",
        "amet",
        "consectetur",
        "adipiscing",
        "elit",
    ]
    adata = ad.AnnData(
        X=sparse.random(len(obs_keys), len(var_keys), format="csr"),
        layers={"layer": sparse.random(len(obs_keys), len(var_keys), format="csr")},
        obs=pd.DataFrame(
            np.random.randn(len(obs_keys), len(obs_keys) + len(var_keys)),
            index=obs_keys,
            columns=obs_keys + var_keys,
        ),
        var=pd.DataFrame(
            np.random.randn(len(var_keys), len(obs_keys) + len(var_keys)),
            index=var_keys,
            columns=var_keys + obs_keys,
        ),
    )

    adata.raw = adata.copy()

    for k in var_keys:
        # These are mostly to check that the test is working
        assert k in adata.var_names
        assert k in adata.obs.columns
        # Now the actual checks:
        with pytest.raises(ValueError, match=r".*var_names.*obs\.columns.*"):
            adata.obs_vector(k)
        with pytest.raises(ValueError, match=r".*var_names.*obs\.columns.*"):
            adata.obs_vector(k, layer="layer")

        # Should uniquely select column from in adata.var
        assert list(adata.var[k]) == list(adata.var_vector(k))
        assert list(adata.var[k]) == list(adata.var_vector(k, layer="layer"))

        assert list(adata.raw.var[k]) == list(adata.raw.var_vector(k))

    for k in obs_keys:
        assert k in adata.obs_names
        assert k in adata.var.columns
        with pytest.raises(ValueError, match=r".*obs_names.*var\.columns"):
            adata.var_vector(k)
        with pytest.raises(ValueError, match=r".*obs_names.*var\.columns"):
            adata.var_vector(k, layer="layer")

        assert list(adata.obs[k]) == list(adata.obs_vector(k))
        assert list(adata.obs[k]) == list(adata.obs_vector(k, layer="layer"))

        with pytest.raises(ValueError, match=r".*obs_names.*var\.columns*"):
            adata.raw.var_vector(k)
