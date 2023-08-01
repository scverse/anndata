from __future__ import annotations

from .multi_files import AnnCollection
from .pytorch import AnnLoader

from anndata._io.specs import read_elem, write_elem, IOSpec
from anndata._core.sparse_dataset import sparse_dataset, CSRDataset, CSCDataset
from ._dispatch_io import read_dispatched, write_dispatched
from .merge import concat_on_disk

__all__ = [
    "AnnCollection",
    "AnnLoader",
    "read_elem",
    "write_elem",
    "read_dispatched",
    "write_dispatched",
    "IOSpec",
    "concat_on_disk",
    "sparse_dataset",
    "CSRDataset",
    "CSCDataset",
]
