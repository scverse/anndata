from __future__ import annotations

import warnings
from types import MappingProxyType
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


# Map old name in `anndata.experimental` to new name in `anndata`
_DEPRECATED = MappingProxyType(
    dict(
        (kv if isinstance(kv, tuple) else (kv, kv))
        for kv in (
            "CSRDataset",
            "CSCDataset",
            "sparse_dataset",
            "read_elem",
            "write_elem",
            ("RWAble", "AxisStorable"),
            ("InMemoryElem", "RWAble"),
        )
    )
)


def __getattr__(attr_name: str) -> Any:
    if new_name := _DEPRECATED.get(attr_name):
        rename = f"using its new name {new_name} " if new_name != attr_name else ""
        msg = (
            f"Importing {attr_name} from `anndata.experimental` is deprecated. "
            f"Import {rename}from `anndata` directly."
        )
        warnings.warn(msg, FutureWarning)
        return getattr(anndata, new_name)
    msg = f"module {__name__!r} has no attribute {attr_name!r}"
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
