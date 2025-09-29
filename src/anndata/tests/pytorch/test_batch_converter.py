from __future__ import annotations

import numpy as np
import pytest
import torch

import anndata as ad
from anndata.experimental.pytorch import AnnLoader, batch_dict_converter


def _make_dummy_adata(n_obs: int = 10, n_vars: int = 4):
    X = np.random.rand(n_obs, n_vars).astype(np.float32)
    obs = {"group": np.arange(n_obs) % 2}
    return ad.AnnData(X=X, obs=obs)


@pytest.mark.parametrize("num_workers", [0, 2])
def test_batch_converter_returns_dict(num_workers):
    adata = _make_dummy_adata(16, 3)
    loader = AnnLoader(
        adata,
        batch_size=4,
        batch_converter=batch_dict_converter,
        num_workers=num_workers,
    )
    batch = next(iter(loader))

    assert isinstance(batch, dict)
    assert "x" in batch and isinstance(batch["x"], torch.Tensor)
    assert batch["x"].shape == (4, 3)
    assert "group" in batch
