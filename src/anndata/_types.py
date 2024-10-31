"""
Defines some useful types for this library. Should probably be cleaned up before thinking about exporting.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal, Protocol, TypeVar

from .compat import (
    H5Array,
    H5Group,
    ZarrArray,
    ZarrGroup,
)
from .typing import RWAble

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any, TypeAlias

    from ._io.specs.registry import (
        IOSpec,
        LazyDataStructures,
        LazyReader,
        Reader,
        Writer,
    )

__all__ = [
    "ArrayStorageType",
    "GroupStorageType",
    "StorageType",
]

ArrayStorageType: TypeAlias = ZarrArray | H5Array
GroupStorageType: TypeAlias = ZarrGroup | H5Group
StorageType: TypeAlias = ArrayStorageType | GroupStorageType

# NOTE: If you change these, be sure to update `autodoc_type_aliases` in docs/conf.py!
ContravariantRWAble = TypeVar("ContravariantRWAble", bound=RWAble, contravariant=True)
CovariantRWAble = TypeVar("CovariantRWAble", bound=RWAble, covariant=True)
InvariantRWAble = TypeVar("InvariantRWAble", bound=RWAble)

SCo = TypeVar("SCo", covariant=True, bound=StorageType)
SCon = TypeVar("SCon", contravariant=True, bound=StorageType)


class _ReadInternal(Protocol[SCon, CovariantRWAble]):
    def __call__(self, elem: SCon, *, _reader: Reader) -> CovariantRWAble: ...


class _ReadLazyInternal(Protocol[SCon]):
    def __call__(
        self, elem: SCon, *, _reader: LazyReader, chunks: tuple[int, ...] | None = None
    ) -> LazyDataStructures: ...


class Read(Protocol[SCon, CovariantRWAble]):
    def __call__(self, elem: SCon) -> CovariantRWAble:
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


class _WriteInternal(Protocol[ContravariantRWAble]):
    def __call__(
        self,
        f: StorageType,
        k: str,
        v: ContravariantRWAble,
        *,
        _writer: Writer,
        dataset_kwargs: Mapping[str, Any],
    ) -> None: ...


class Write(Protocol[ContravariantRWAble]):
    def __call__(
        self,
        f: StorageType,
        k: str,
        v: ContravariantRWAble,
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


class ReadCallback(Protocol[SCo, InvariantRWAble]):
    def __call__(
        self,
        /,
        read_func: Read[SCo, InvariantRWAble],
        elem_name: str,
        elem: StorageType,
        *,
        iospec: IOSpec,
    ) -> InvariantRWAble:
        """
        Callback used in :func:`anndata.experimental.read_dispatched` to customize reading an element from a store.

        Params
        ------
        read_func
            :func:`anndata.io.read_elem` function to call to read the current element given the ``iospec``.
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


class WriteCallback(Protocol[InvariantRWAble]):
    def __call__(
        self,
        /,
        write_func: Write[InvariantRWAble],
        store: StorageType,
        elem_name: str,
        elem: InvariantRWAble,
        *,
        iospec: IOSpec,
        dataset_kwargs: Mapping[str, Any],
    ) -> None:
        """
        Callback used in :func:`anndata.experimental.write_dispatched` to customize writing an element to a store.

        Params
        ------
        write_func
            :func:`anndata.io.write_elem` function to call to read the current element given the ``iospec``.
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


AnnDataElem = Literal[
    "obs",
    "var",
    "obsm",
    "varm",
    "obsp",
    "varp",
    "layers",
    "X",
    "raw",
    "uns",
]

Join_T = Literal["inner", "outer"]
