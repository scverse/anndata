from __future__ import annotations

from .multi_files import AnnCollection
from .pytorch import AnnLoader

from anndata._io.specs import read_elem, write_elem, IOSpec
from ._dispatch_io import read_dispatched, write_dispatched
from .merge import concat_on_disk
from .backed._io import read_backed

__all__ = [
    "AnnCollection",
    "AnnLoader",
    "read_elem",
    "write_elem",
    "read_dispatched",
    "write_dispatched",
    "IOSpec",
    "concat_on_disk",
    "read_backed",
]
