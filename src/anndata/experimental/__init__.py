from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

from anndata._core.sparse_dataset import CSCDataset, CSRDataset, sparse_dataset
from anndata._io.specs import IOSpec, read_elem, read_elem_as_dask, write_elem

from .._types import Read, ReadCallback, StorageType, Write, WriteCallback
from ._dispatch_io import read_dispatched, write_dispatched
from .merge import concat_on_disk
from .multi_files import AnnCollection
from .pytorch import AnnLoader

if TYPE_CHECKING:
    from typing import Any

deprecated = ["CSRDataset", "CSCDataset", "read_elem", "write_elem"]


def __getattr__(key: str) -> Any:
    if key in deprecated:
        msg = f"Importing {key} from `anndata.experimental` is deprecated. Import from `anndata` directly."
        warnings.warn(msg, FutureWarning)
    return globals()[key]


__all__ = [
    "AnnCollection",
    "AnnLoader",
    "read_elem",
    "write_elem",
    "read_elem_as_dask",
    "read_dispatched",
    "write_dispatched",
    "IOSpec",
    "concat_on_disk",
    "sparse_dataset",
    "CSRDataset",
    "CSCDataset",
    "Read",
    "Write",
    "ReadCallback",
    "WriteCallback",
    "StorageType",
]
