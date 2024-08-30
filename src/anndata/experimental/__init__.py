from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

from anndata._core.sparse_dataset import (  # noqa: F401
    CSCDataset,
    CSRDataset,
    sparse_dataset,
)
from anndata._io.specs import (  # noqa: F401
    IOSpec,
    read_elem,
    read_elem_as_dask,
    write_elem,
)

from .._types import Read, ReadCallback, StorageType, Write, WriteCallback
from ._dispatch_io import read_dispatched, write_dispatched
from .merge import concat_on_disk
from .multi_files import AnnCollection
from .pytorch import AnnLoader

if TYPE_CHECKING:
    from typing import Any

__deprecated__ = ["CSRDataset", "CSCDataset", "read_elem", "write_elem"]


def __getattr__(key: str) -> Any:
    if key in __deprecated__:
        msg = f"Importing {key} from `anndata.experimental` is deprecated. Import from `anndata` directly."
        warnings.warn(msg, FutureWarning)
        import anndata

        return getattr(anndata, key)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "AnnCollection",
    "AnnLoader",
    "read_elem_as_dask",
    "read_dispatched",
    "write_dispatched",
    "IOSpec",
    "concat_on_disk",
    "sparse_dataset",
    "Read",
    "Write",
    "ReadCallback",
    "WriteCallback",
    "StorageType",
]


def __dir__():
    return sorted(__all__ + __deprecated__)
