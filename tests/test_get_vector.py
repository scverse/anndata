from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import anndata as ad

pytestmark = [
    pytest.mark.filterwarnings("ignore:Use anndata.acc.A instead of:FutureWarning"),
]

OBS_KEYS = [
    *["Lorem", "ipsum", "dolor", "sit", "amet"],
    *["consectetur", "adipiscing", "elit"],
]
VAR_KEYS = ["The", "quick", "brown", "fox", "jumps", "over", "the", "lazy", "dog"]


@pytest.fixture(scope="session")
def adata() -> ad.AnnData:
    adata = ad.AnnData(
        X=sparse.random(len(OBS_KEYS), len(VAR_KEYS), format="csr"),
        layers=dict(layer=sparse.random(len(OBS_KEYS), len(VAR_KEYS), format="csr")),
        obs=pd.DataFrame(
            np.random.randn(len(OBS_KEYS), len(OBS_KEYS) + len(VAR_KEYS)),
            index=OBS_KEYS,
            columns=[*OBS_KEYS, *VAR_KEYS],
        ),
        var=pd.DataFrame(
            np.random.randn(len(VAR_KEYS), len(OBS_KEYS) + len(VAR_KEYS)),
            index=VAR_KEYS,
            columns=[*VAR_KEYS, *OBS_KEYS],
        ),
    )
    adata.raw = adata.copy()
    return adata


@pytest.mark.parametrize("key", OBS_KEYS)
def test_amgibuous_keys_obs(adata: ad.AnnData, key: str) -> None:
    """Tests that an error is raised if obs_vector is ambiguous."""
    assert key in adata.obs_names
    assert key in adata.var.columns

    with pytest.raises(ValueError, match=r".*obs_names.*var\.columns"):
        adata.var_vector(key)
    with pytest.raises(ValueError, match=r".*obs_names.*var\.columns"):
        adata.var_vector(key, layer="layer")

    assert list(adata.obs[key]) == list(adata.obs_vector(key))
    assert list(adata.obs[key]) == list(adata.obs_vector(key, layer="layer"))

    with pytest.raises(ValueError, match=r".*obs_names.*var\.columns*"):
        adata.raw.var_vector(key)


@pytest.mark.parametrize("key", VAR_KEYS)
def test_amgibuous_keys_var(adata: ad.AnnData, key: str) -> None:
    """Tests that an error is raised if var_vector is ambiguous."""
    assert key in adata.var_names
    assert key in adata.obs.columns

    with pytest.raises(ValueError, match=r".*var_names.*obs\.columns.*"):
        adata.obs_vector(key)
    with pytest.raises(ValueError, match=r".*var_names.*obs\.columns.*"):
        adata.obs_vector(key, layer="layer")

    assert list(adata.var[key]) == list(adata.var_vector(key))
    assert list(adata.var[key]) == list(adata.var_vector(key, layer="layer"))

    assert list(adata.raw.var[key]) == list(adata.raw.var_vector(key))
