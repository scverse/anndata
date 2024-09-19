from __future__ import annotations

from types import MappingProxyType
from typing import TYPE_CHECKING

from .._types import Read, ReadCallback, StorageType, Write, WriteCallback
from ..io.specs import IOSpec, read_elem_as_dask
from ..utils import module_get_attr_redirect
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
    return module_get_attr_redirect(
        attr_name, old_module_path="experimental", deprecated_mapping=_DEPRECATED
    )


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
