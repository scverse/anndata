from __future__ import annotations

import numpy as np
import pytest

import anndata as ad

pytest.importorskip("torch")

from anndata.experimental.pytorch import AnnLoader


@pytest.fixture
def adata():
    return ad.AnnData(X=np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]))


def test_annloader_default_device(adata):
    """AnnLoader with default device='cpu' produces CPU tensors."""
    loader = AnnLoader(adata, batch_size=2)
    batch = next(iter(loader))
    assert batch.X.device.type == "cpu"


def test_annloader_explicit_cpu_device(adata):
    """AnnLoader with explicit device='cpu' produces CPU tensors."""
    loader = AnnLoader(adata, batch_size=2, device="cpu")
    batch = next(iter(loader))
    assert batch.X.device.type == "cpu"


def test_annloader_use_cuda_deprecation_warning(adata):
    """Passing use_cuda emits a FutureWarning."""
    with pytest.warns(FutureWarning, match="use_cuda.*deprecated"):
        # use_cuda=False should still emit warning (parameter was explicitly passed)
        AnnLoader(adata, batch_size=2, use_cuda=False)


def test_annloader_use_cuda_and_device_conflict(adata):
    """Passing both use_cuda and device raises ValueError."""
    with pytest.raises(ValueError, match="Cannot specify both"):
        AnnLoader(adata, batch_size=2, use_cuda=True, device="cuda")
