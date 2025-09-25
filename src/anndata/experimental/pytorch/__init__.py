from __future__ import annotations

from ._ann_dataset import (
    TORCH_AVAILABLE,
    AnnDataset,
)
from ._annloader import AnnLoader

# Import transforms if torch is available
if TORCH_AVAILABLE:
    from . import transforms
    from .transforms import (
        Compose,
        Transform,
        get_obs_key,
    )

    __all__ = [
        "TORCH_AVAILABLE",
        "AnnDataset",
        "AnnLoader",
        "Compose",
        "Transform",
        "get_obs_key",
        "transforms",
    ]
else:
    __all__ = [
        "TORCH_AVAILABLE",
        "AnnDataset",
        "AnnLoader",
    ]
