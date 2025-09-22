"""
Defines some useful types for this library. Should probably be cleaned up before thinking about exporting.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal, Protocol, TypeVar

from . import typing
from .compat import H5Array, H5Group, ZarrArray, ZarrGroup

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any, TypeAlias

    from anndata._core.xarray import Dataset2D

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
    "_ReadInternal",
    "_ReadLazyInternal",
    "_WriteInternal",
]

ArrayStorageType: TypeAlias = ZarrArray | H5Array
GroupStorageType: TypeAlias = ZarrGroup | H5Group
StorageType: TypeAlias = ArrayStorageType | GroupStorageType

# NOTE: If you change these, be sure to update `autodoc_type_aliases` in docs/conf.py!
RWAble_contra = TypeVar("RWAble_contra", bound=typing.RWAble, contravariant=True)
RWAble_co = TypeVar("RWAble_co", bound=typing.RWAble, covariant=True)
RWAble = TypeVar("RWAble", bound=typing.RWAble)

S_co = TypeVar("S_co", covariant=True, bound=StorageType)
S_contra = TypeVar("S_contra", contravariant=True, bound=StorageType)


class Dataset2DIlocIndexer(Protocol):
    def __getitem__(self, idx: Any) -> Dataset2D: ...


class _ReadInternal(Protocol[S_contra, RWAble_co]):
    def __call__(self, elem: S_contra, *, _reader: Reader) -> RWAble_co: ...


class _ReadLazyInternal(Protocol[S_contra]):
    def __call__(
        self,
        elem: S_contra,
        *,
        _reader: LazyReader,
        chunks: tuple[int, ...] | None = None,
    ) -> LazyDataStructures: ...


class Read(Protocol[S_contra, RWAble_co]):
    def __call__(self, elem: S_contra) -> RWAble_co:
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


class ReadLazy(Protocol[S_contra]):
    def __call__(
        self, elem: S_contra, *, chunks: tuple[int, ...] | None = None
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


class _WriteInternal(Protocol[RWAble_contra]):
    def __call__(
        self,
        f: StorageType,
        k: str,
        v: RWAble_contra,
        *,
        _writer: Writer,
        dataset_kwargs: Mapping[str, Any],
    ) -> None: ...


class Write(Protocol[RWAble_contra]):
    def __call__(
        self,
        f: StorageType,
        k: str,
        v: RWAble_contra,
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


class ReadCallback(Protocol[S_co, RWAble]):
    def __call__(
        self,
        /,
        read_func: Read[S_co, RWAble],
        elem_name: str,
        elem: StorageType,
        *,
        iospec: IOSpec,
    ) -> RWAble:
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


class WriteCallback(Protocol[RWAble]):
    def __call__(
        self,
        /,
        write_func: Write[RWAble],
        store: StorageType,
        elem_name: str,
        elem: RWAble,
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
