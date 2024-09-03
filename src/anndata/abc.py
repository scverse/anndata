from __future__ import annotations

from ._core.sparse_dataset import AbstractCSCDataset as CSCDataset
from ._core.sparse_dataset import AbstractCSRDataset as CSRDataset

__all__ = [
    "CSRDataset",
    "CSCDataset",
]
