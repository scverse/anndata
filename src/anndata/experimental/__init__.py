from __future__ import annotations

from anndata._core.sparse_dataset import CSCDataset, CSRDataset, sparse_dataset
from anndata._io.specs import IOSpec, read_elem, write_elem

from .._core.storage import StorageType
from .._types import InMemoryElem as _InMemoryElem
from .._types import Read, ReadCallback, Write, WriteCallback
from .._types import RWAble as _RWAble
from .._types import RWAbleDict as _RWAbleDict
from .._types import RWAbleList as _RWAbleList
from ._dispatch_io import read_dispatched, write_dispatched
from .merge import concat_on_disk
from .multi_files import AnnCollection
from .pytorch import AnnLoader

# Sphinx can’t find data docstrings when objects are re-exported
InMemoryElem = _InMemoryElem
"""An in-memory element that can be read and written."""
RWAble = _RWAble
"""A serializable object."""
RWAbleDict = _RWAbleDict
"""A dict containing serializable objects."""
RWAbleList = _RWAbleList
"""A list containing serializable objects."""

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
    "InMemoryElem",
    "Read",
    "RWAbleDict",
    "RWAbleList",
    "RWAble",
    "Write",
    "ReadCallback",
    "WriteCallback",
    "StorageType",
]
