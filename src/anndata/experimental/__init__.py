from __future__ import annotations

from anndata._core.sparse_dataset import sparse_dataset
from anndata._io.specs import IOSpec, read_elem_as_dask

from .._core.storage import StorageType
from .._types import InMemoryElem as _InMemoryElem
from .._types import Read, ReadCallback, Write, WriteCallback
from ._dispatch_io import read_dispatched, write_dispatched
from .merge import concat_on_disk
from .multi_files import AnnCollection
from .pytorch import AnnLoader

# Sphinx can’t find data docstrings when objects are re-exported
InMemoryElem = _InMemoryElem
"""An in-memory element that can be read and written, including an :class:`anndata.AnnData` objects."""

__all__ = [
    "AnnCollection",
    "AnnLoader",
    "read_elem_as_dask",
    "read_dispatched",
    "write_dispatched",
    "IOSpec",
    "concat_on_disk",
    "sparse_dataset",
    "InMemoryElem",
    "Read",
    "Write",
    "ReadCallback",
    "WriteCallback",
    "StorageType",
]
