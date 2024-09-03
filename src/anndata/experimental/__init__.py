from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import anndata

from .._io.specs import IOSpec, read_elem_as_dask
from .._types import Read, ReadCallback, StorageType, Write, WriteCallback
from ._dispatch_io import read_dispatched, write_dispatched
from .merge import concat_on_disk
from .multi_files import AnnCollection
from .pytorch import AnnLoader

if TYPE_CHECKING:
    from typing import Any


_DEPRECATED = ["CSRDataset", "CSCDataset", "read_elem", "write_elem", "sparse_dataset"]


def __getattr__(key: str) -> Any:
    if key in _DEPRECATED:
        msg = f"Importing {key} from `anndata.experimental` is deprecated. Import from `anndata` directly."
        warnings.warn(msg, FutureWarning)
        return getattr(anndata, key)
    msg = f"module {__name__!r} has no attribute {key!r}"
    raise AttributeError(msg)


__all__ = [
    "AnnCollection",
    "AnnLoader",
    "read_elem_as_dask",
    "read_dispatched",
    "write_dispatched",
    "IOSpec",
    "concat_on_disk",
    "Read",
    "Write",
    "ReadCallback",
    "WriteCallback",
    "StorageType",
]
