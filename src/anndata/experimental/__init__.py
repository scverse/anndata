from __future__ import annotations

import warnings
from types import MappingProxyType
from typing import TYPE_CHECKING

import anndata

from .._io.specs import IOSpec, read_elem_lazy
from .._types import Read, ReadCallback, StorageType, Write, WriteCallback
from ._dispatch_io import read_dispatched, write_dispatched
from .backed import read_lazy
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
            ("CSRDataset", "abc.CSRDataset"),
            ("CSCDataset", "abc.CSCDataset"),
            "sparse_dataset",
            "read_elem",
            "write_elem",
            ("RWAble", "typing.AxisStorable"),
            ("InMemoryElem", "typing.RWAble"),
        )
    )
)


def __getattr__(attr_name: str) -> Any:
    if new_path := _DEPRECATED.get(attr_name):
        msg = (
            f"Importing {attr_name} from `anndata.experimental` is deprecated. "
            f"Import anndata.{new_path} instead."
        )
        warnings.warn(msg, FutureWarning)
        # hacky import_object_by_name, but we test all these
        mod = anndata
        while "." in new_path:
            mod_name, new_path = new_path.split(".", 1)
            mod = getattr(mod, mod_name)
        return getattr(mod, new_path)
    msg = f"module {__name__!r} has no attribute {attr_name!r}"
    raise AttributeError(msg)


__all__ = [
    "AnnCollection",
    "AnnLoader",
    "read_elem_lazy",
    "read_dispatched",
    "write_dispatched",
    "IOSpec",
    "concat_on_disk",
    "Read",
    "read_lazy",
    "Write",
    "ReadCallback",
    "WriteCallback",
    "StorageType",
]
