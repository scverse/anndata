from .multi_files import AnnCollection
from .pytorch import AnnLoader

from anndata._io.specs import read_elem, write_elem

__all__ = [
    "AnnCollection",
    "AnnLoader",
    "read_elem",
    "write_elem",
    "read_dispatched",
]


def read_dispatched(store, callback):
    from anndata._io.specs import Reader, _REGISTRY

    reader = Reader(_REGISTRY, callback=callback)

    return reader.read_elem(store)
