from __future__ import annotations

from types import MappingProxyType
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any

    from anndata._types import (
        GroupStorageType,
        ReadCallback,
        StorageType,
        WriteCallback,
    )
    from anndata.typing import RWAble


def read_dispatched(
    elem: StorageType,
    callback: ReadCallback,
) -> RWAble:
    """
    Read elem, calling the callback at each sub-element.

    Params
    ------
    elem
        Storage container (e.g. `h5py.Group`, `zarr.Group`).
        This must have anndata element specifications.
    callback
        Function to call at each anndata encoded element.

    See Also
    --------
    :doc:`/tutorials/notebooks/{read,write}_dispatched`
    """
    from anndata._io.specs import _REGISTRY, Reader

    reader = Reader(_REGISTRY, callback=callback)

    return reader.read_elem(elem)


def write_dispatched(
    store: GroupStorageType,
    key: str,
    elem: RWAble,
    callback: WriteCallback,
    *,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
) -> None:
    """
    Write elem to store, recursively calling callback at each sub-element.

    Params
    ------
    store
        Storage container to be written to.
    key
        Key to write element to. To write to the root group, use "/".
    elem
        The element to write. Probably an AnnData.
    callback
        Function called when writing each element.
    dataset_kwargs
        Keyword arguments to pass to the dataset creation function.

    See Also
    --------
    :doc:`/tutorials/notebooks/{read,write}_dispatched`
    """
    from anndata._io.specs import _REGISTRY, Writer

    writer = Writer(_REGISTRY, callback=callback)

    writer.write_elem(store, key, elem, dataset_kwargs=dataset_kwargs)
