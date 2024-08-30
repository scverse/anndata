"""
Defines some useful types for this library. Should probably be cleaned up before thinking about exporting.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, TypeVar, Union

import numpy as np
import pandas as pd

from anndata._core.anndata import AnnData

from ._core.storage import ArrayDataStructureType
from .compat import (
    H5Array,
    H5Group,
    ZarrArray,
    ZarrGroup,
)

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any, TypeAlias

    from anndata._io.specs.registry import LazyDataStructures, LazyReader

    from ._io.specs.registry import IOSpec, Reader, Writer

__all__ = [
    "ArrayStorageType",
    "GroupStorageType",
    "StorageType",
]

InMemoryArrayOrScalarType: TypeAlias = Union[
    pd.DataFrame, np.number, str, ArrayDataStructureType
]
RWAble: TypeAlias = Union[
    InMemoryArrayOrScalarType, dict[str, "RWAble"], list["RWAble"]
]  # noqa: TCH010
InMemoryElem: TypeAlias = Union[
    RWAble,
    AnnData,
    pd.Categorical,
    pd.api.extensions.ExtensionArray,
]

ArrayStorageType: TypeAlias = Union[ZarrArray, H5Array]
GroupStorageType: TypeAlias = Union[ZarrGroup, H5Group]
StorageType: TypeAlias = Union[ArrayStorageType, GroupStorageType]

# NOTE: If you change these, be sure to update `autodoc_type_aliases` in docs/conf.py!
ContravariantInMemoryType = TypeVar(
    "ContravariantInMemoryType", bound="InMemoryElem", contravariant=True
)
CovariantInMemoryType = TypeVar(
    "CovariantInMemoryType", bound="InMemoryElem", covariant=True
)
InvariantInMemoryType = TypeVar("InvariantInMemoryType", bound="InMemoryElem")

SCo = TypeVar("SCo", covariant=True, bound=StorageType)
SCon = TypeVar("SCon", contravariant=True, bound=StorageType)


class _ReadInternal(Protocol[SCon, CovariantInMemoryType]):
    def __call__(self, elem: SCon, *, _reader: Reader) -> CovariantInMemoryType: ...


class _ReadLazyInternal(Protocol[SCon]):
    def __call__(
        self, elem: SCon, *, _reader: LazyReader, chunks: tuple[int, ...] | None = None
    ) -> LazyDataStructures: ...


class Read(Protocol[SCon, CovariantInMemoryType]):
    def __call__(self, elem: SCon) -> CovariantInMemoryType:
        """Low-level reading function for an element.

        Parameters
        ----------
        elem
            The element to read from.
        Returns
        -------
            The element read from the store.
        """
        ...


class ReadLazy(Protocol[SCon]):
    def __call__(
        self, elem: SCon, *, chunks: tuple[int, ...] | None = None
    ) -> LazyDataStructures:
        """Low-level reading function for a lazy element.

        Parameters
        ----------
        elem
            The element to read from.
        chunks
            The chunk size to be used.
        Returns
        -------
            The lazy element read from the store.
        """
        ...


class _WriteInternal(Protocol[ContravariantInMemoryType]):
    def __call__(
        self,
        f: StorageType,
        k: str,
        v: ContravariantInMemoryType,
        *,
        _writer: Writer,
        dataset_kwargs: Mapping[str, Any],
    ) -> None: ...


class Write(Protocol[ContravariantInMemoryType]):
    def __call__(
        self,
        f: StorageType,
        k: str,
        v: ContravariantInMemoryType,
        *,
        dataset_kwargs: Mapping[str, Any],
    ) -> None:
        """Low-level writing function for an element.

        Parameters
        ----------
        f
            The store to which `elem` should be written.
        k
            The key to read in from the group.
        v
            The element to write out.
        dataset_kwargs
            Keyword arguments to be passed to a library-level io function, like `chunks` for :doc:`zarr:index`.
        """
        ...


class ReadCallback(Protocol[SCo, InvariantInMemoryType]):
    def __call__(
        self,
        /,
        read_func: Read[SCo, InvariantInMemoryType],
        elem_name: str,
        elem: StorageType,
        *,
        iospec: IOSpec,
    ) -> InvariantInMemoryType:
        """
        Callback used in :func:`anndata.experimental.read_dispatched` to customize reading an element from a store.

        Params
        ------
        read_func
            :func:`anndata.experimental.read_elem` function to call to read the current element given the ``iospec``.
        elem_name
            The key to read in from the group.
        elem
            The element to read from.
        iospec
            Internal AnnData encoding specification for the element.

        Returns
        -------
            The element read from the store.
        """
        ...


class WriteCallback(Protocol[InvariantInMemoryType]):
    def __call__(
        self,
        /,
        write_func: Write[InvariantInMemoryType],
        store: StorageType,
        elem_name: str,
        elem: InvariantInMemoryType,
        *,
        iospec: IOSpec,
        dataset_kwargs: Mapping[str, Any],
    ) -> None:
        """
        Callback used in :func:`anndata.experimental.write_dispatched` to customize writing an element to a store.

        Params
        ------
        write_func
            :func:`anndata.experimental.write_elem` function to call to read the current element given the ``iospec``.
        store
            The store to which `elem` should be written.
        elem_name
            The key to read in from the group.
        elem
            The element to write out.
        iospec
            Internal AnnData encoding specification for the element.
        dataset_kwargs
            Keyword arguments to be passed to a library-level io function, like `chunks` for :doc:`zarr:index`.
        """
        ...
