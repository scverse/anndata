from __future__ import annotations

import numpy as np
import pytest

import anndata as ad

pytest.importorskip("torch")

import torch

from anndata.experimental.pytorch import AnnLoader, batch_dict_converter


def _make_dummy_adata(n_obs: int = 10, n_vars: int = 4):
    X = np.random.rand(n_obs, n_vars).astype(np.float32)
    obs = {"group": np.arange(n_obs) % 2, "category": ["A", "B"] * (n_obs // 2)}
    return ad.AnnData(X=X, obs=obs)


def test_batch_converter_default_behavior():
    """Test that without batch_converter, we get AnnCollectionView."""
    adata = _make_dummy_adata(16, 3)
    loader = AnnLoader(adata, batch_size=4)
    batch = next(iter(loader))

    # Should be AnnCollectionView without converter
    assert hasattr(batch, "X")
    assert hasattr(batch, "obs")


def test_batch_converter_returns_dict():
    """Test that batch_dict_converter returns proper dict format."""
    adata = _make_dummy_adata(16, 3)
    loader = AnnLoader(
        adata,
        batch_size=4,
        batch_converter=batch_dict_converter,
        num_workers=0,  # Single-threaded for now
    )
    batch = next(iter(loader))

    assert isinstance(batch, dict)
    assert "x" in batch
    assert isinstance(batch["x"], torch.Tensor)
    assert batch["x"].shape == (4, 3)
    assert "group" in batch
    assert isinstance(batch["group"], torch.Tensor)
    assert "category" in batch


def test_custom_batch_converter():
    """Test that custom batch converters work."""
    adata = _make_dummy_adata(8, 2)

    def custom_converter(batch):
        result = batch_dict_converter(batch)
        result["custom_field"] = torch.tensor([42])
        result["x_sum"] = result["x"].sum()
        return result

    loader = AnnLoader(
        adata,
        batch_size=4,
        batch_converter=custom_converter,
        num_workers=0,
    )
    batch = next(iter(loader))

    assert isinstance(batch, dict)
    assert "x" in batch
    assert "group" in batch
    assert "custom_field" in batch
    assert "x_sum" in batch
    assert batch["custom_field"].item() == 42
    assert isinstance(batch["x_sum"], torch.Tensor)


def test_batch_converter_multiprocessing_works():
    """Test that batch converter now works with num_workers > 0."""
    adata = _make_dummy_adata(16, 3)
    loader = AnnLoader(
        adata,
        batch_size=4,
        batch_converter=batch_dict_converter,
        num_workers=2,
    )

    # This should now work with our multiprocessing fix
    batch = next(iter(loader))
    assert isinstance(batch, dict)
    assert "x" in batch
    assert batch["x"].shape == (4, 3)
