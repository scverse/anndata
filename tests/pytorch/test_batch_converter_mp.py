from __future__ import annotations

import numpy as np
import pytest

import anndata as ad
from anndata.experimental.pytorch import AnnLoader, batch_dict_converter


@pytest.mark.filterwarnings("ignore:Series.__getitem__")
def test_worker_safe_batch_converter():
    """AnnLoader should work with num_workers > 0 when batch_converter is supplied."""
    adata = ad.AnnData(X=np.random.rand(32, 4).astype(np.float32))

    loader = AnnLoader(
        adata,
        batch_size=8,
        batch_converter=batch_dict_converter,
        num_workers=2,
    )

    batch = next(iter(loader))
    assert isinstance(batch, dict)
    assert batch["x"].shape == (8, 4)


def test_default_collate_fails_with_anncollection():
    """Sanity-check that vanilla DataLoader still fails, documenting why fix is needed."""
    from torch.utils.data import DataLoader

    from anndata.experimental.multi_files import AnnCollection

    adata1 = ad.AnnData(X=np.random.rand(4, 3).astype(np.float32))
    adata2 = ad.AnnData(X=np.random.rand(4, 3).astype(np.float32))
    coll = AnnCollection([adata1, adata2])
    dataset = coll  # AnnCollection implements __getitem__/__len__

    failing_loader = DataLoader(dataset, batch_size=4, num_workers=2)
    with pytest.raises(TypeError):
        next(iter(failing_loader))
