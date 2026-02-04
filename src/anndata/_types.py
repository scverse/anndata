"""
Defines some useful types for this library. Should probably be cleaned up before thinking about exporting.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal, Protocol

from .compat import H5Array, H5Group, ZarrArray, ZarrGroup
from .utils import set_module

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any, TypeAlias

    from anndata._core.xarray import Dataset2D

    from . import typing
    from ._io.specs.registry import (
        IOSpec,
        LazyDataStructures,
        LazyReader,
        Reader,
        Writer,
    )
else:  # https://github.com/tox-dev/sphinx-autodoc-typehints/issues/580
    type S = StorageType
    type RWAble = typing.RWAble


__all__ = [
    "StorageType",
    "_ArrayStorageType",
    "_GroupStorageType",
    "_ReadInternal",
    "_ReadLazyInternal",
    "_WriteInternal",
]

# These two are not public, so we don’t make them `type`s
_ArrayStorageType: TypeAlias = ZarrArray | H5Array  # noqa: UP040
_GroupStorageType: TypeAlias = ZarrGroup | H5Group  # noqa: UP040
type StorageType = _ArrayStorageType | _GroupStorageType


@set_module("anndata.experimental")
class Dataset2DIlocIndexer(Protocol):
    def __getitem__(self, idx: Any) -> Dataset2D: ...


class _ReadInternal[S: StorageType, RWAble: typing.RWAble](Protocol):
    def __call__(self, elem: S, *, _reader: Reader) -> RWAble: ...


class _ReadLazyInternal[S: StorageType](Protocol):
    def __call__(
        self,
        elem: S,
        *,
        _reader: LazyReader,
        chunks: tuple[int, ...] | None = None,
    ) -> LazyDataStructures: ...


@set_module("anndata.experimental")
class Read[S: StorageType, RWAble: typing.RWAble](Protocol):
    def __call__(self, elem: S) -> RWAble:
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


class ReadLazy[S](Protocol):
    def __call__(
        self, elem: S, *, chunks: tuple[int, ...] | None = None
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


class _WriteInternal[RWAble: typing.RWAble](Protocol):
    def __call__(
        self,
        f: StorageType,
        k: str,
        v: RWAble,
        *,
        _writer: Writer,
        dataset_kwargs: Mapping[str, Any],
    ) -> None: ...


@set_module("anndata.experimental")
class Write[RWAble: typing.RWAble](Protocol):
    def __call__(
        self,
        f: StorageType,
        k: str,
        v: RWAble,
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
            Keyword arguments to be passed to a library-level io function, like `chunks` for :mod:`zarr`.
        """
        ...


@set_module("anndata.experimental")
class ReadCallback[S: StorageType, RWAble: typing.RWAble](Protocol):
    def __call__(
        self,
        /,
        read_func: Read[S, RWAble],
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


@set_module("anndata.experimental")
class WriteCallback[RWAble: typing.RWAble](Protocol):
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
            Keyword arguments to be passed to a library-level io function, like `chunks` for :mod:`zarr`.
        """
        ...


type AnnDataElem = Literal[
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

type Join_T = Literal["inner", "outer"]
