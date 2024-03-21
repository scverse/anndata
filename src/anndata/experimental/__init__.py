from __future__ import annotations

from anndata._core.sparse_dataset import CSCDataset, CSRDataset, sparse_dataset
from anndata._io.specs import IOSpec, read_elem, write_elem

from ._dispatch_io import read_dispatched, write_dispatched
from .merge import concat_on_disk
from .multi_files import AnnCollection
from .pytorch import AnnLoader

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
