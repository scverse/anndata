from __future__ import annotations

from types import MappingProxyType
from typing import Callable, Any, Union

from .multi_files import AnnCollection
from .pytorch import AnnLoader

from anndata._io.specs import read_elem, write_elem, IOSpec
from anndata.compat import H5Array, H5Group, ZarrArray, ZarrGroup

StorageType = Union[H5Array, H5Group, ZarrArray, ZarrGroup]
StorageGroupType = Union[H5Group, ZarrGroup]

__all__ = [
    "AnnCollection",
    "AnnLoader",
    "read_elem",
    "write_elem",
    "read_dispatched",
    "write_dispatched",
    "IOSpec",
]


def read_dispatched(
    store: StorageType,
    callback: Callable[[Callable[[StorageType], Any], str, StorageType, IOSpec], Any],
) -> Any:
    """
    Read store, calling the callback at each element.

    Params
    ------
    store
        Storage container (e.g. `h5py.Group`, `zarr.Group`). This must have anndata
        element specifications.
    callback
        Function to call at each anndata encoded element. See details below for
        signature.


    The callback has the following signature:

    * `read_func` (`Callable`): A callable which takes the encoded element and returns it's decoded value.
      This is the default decoding function, and what to call if you don't want to modify the decoding.
      It will call this callback again at the next element encoding it sees.
    * `key` (`str`): They absolute key of the element in the store. This will be an absolute key.
    * `elem` (`StorageType`): The encoded element.
    * `spec` (`IOSpec`): The specification of the element.
    """
    from anndata._io.specs import Reader, _REGISTRY

    reader = Reader(_REGISTRY, callback=callback)

    return reader.read_elem(store)


def write_dispatched(
    store: StorageGroupType,
    key: str,
    elem: Any,
    callback: Callable[
        [Callable[[StorageType, str, Any], None], StorageGroupType, str, Any, dict],
        None,
    ],
    *,
    dataset_kwargs=MappingProxyType({}),
) -> None:
    """
    Write elem to store, recusively calling callback at each sub-element.

    Params
    ------
    store
        Storage container to be written to.
    key
        Key to write element to. To write to the root group, use "/".
    elem
        The element to write. Probably an AnnData.
    callback
        Function called when writing each element. See below for signature.
    dataset_kwargs
        Keyword arguments to pass to the dataset creation function.


    The callback has the following signature:

    * `write_func` (`Callable`): A callable which takes the in memory element and writes it to the store.
      This is the default encoding function, and what to call if you don't want to change behaviour at this level.
    * `store` (`StorageGroupType`): The store to write to.
    * `key` (`str`): The key to write elem into store at. This will be an absolute key.
    * `elem` (`Any`): The element to write.
    * `dataset_kwargs` (`dict`): Keyword arguments to pass to the dataset creation function.
    """
    from anndata._io.specs import Writer, _REGISTRY

    writer = Writer(_REGISTRY, callback=callback)

    writer.write_elem(store, key, elem, dataset_kwargs=dataset_kwargs)
