from types import MappingProxyType

from .multi_files import AnnCollection
from .pytorch import AnnLoader

from anndata._io.specs import read_elem, write_elem

__all__ = [
    "AnnCollection",
    "AnnLoader",
    "read_elem",
    "write_elem",
    "read_dispatched",
    "write_dispatched",
]


def read_dispatched(store, callback):
    from anndata._io.specs import Reader, _REGISTRY

    reader = Reader(_REGISTRY, callback=callback)

    return reader.read_elem(store)


def write_dispatched(store, key, elem, callback, dataset_kwargs=MappingProxyType({})):
    from anndata._io.specs import Writer, _REGISTRY

    writer = Writer(_REGISTRY, callback=callback)

    writer.write_elem(store, key, elem, dataset_kwargs=dataset_kwargs)
