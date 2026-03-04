"""
Defines some useful types for this library. Should probably be cleaned up before thinking about exporting.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal, Protocol

from . import typing
from .compat import H5Array, H5Group, ZarrArray, ZarrGroup
from .utils import set_module

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any, TypeAlias

    from anndata._core.xarray import Dataset2D
    from anndata.abc import CSCDataset, CSRDataset

    #: Objects returned by :class:`DaskReader` — always include a :class:`~dask.array.Array`.
    from anndata.compat import DaskArray
    from anndata.experimental.backed._lazy_arrays import CategoricalArray, MaskedArray

    from ._io.specs.registry import BackedReader, DaskReader, IOSpec, Reader, Writer

    type DaskDataStructures = DaskArray | Dataset2D | CategoricalArray | MaskedArray

    #: Objects returned by :class:`BackedReader` — file-backed, never dask.
    type BackedDataStructures = (
        ZarrArray
        | H5Array
        | CSRDataset
        | CSCDataset
        | Dataset2D
        | CategoricalArray
        | MaskedArray
    )

    # Legacy umbrella alias kept for external code that imported it
    type LazyDataStructures = DaskDataStructures | BackedDataStructures

else:  # https://github.com/tox-dev/sphinx-autodoc-typehints/issues/580
    type S = StorageType
    type RWAble = typing.RWAble


__all__ = [
    "StorageType",
    "_ArrayStorageType",
    "_GroupStorageType",
    "_ReadBackedInternal",
    "_ReadDaskInternal",
    "_ReadInternal",
    "_ReadLazyInternal",
    "_WriteInternal",
]

# These two are not public, so we don't make them `type`s
_ArrayStorageType: TypeAlias = ZarrArray | H5Array  # noqa: UP040
_GroupStorageType: TypeAlias = ZarrGroup | H5Group  # noqa: UP040

type StorageType = _ArrayStorageType | _GroupStorageType


@set_module("anndata.experimental")
class Dataset2DIlocIndexer(Protocol):
    def __getitem__(self, idx: Any) -> Dataset2D: ...


class _ReadInternal[S: StorageType, RWAble: typing.RWAble](Protocol):
    def __call__(self, elem: S, *, _reader: Reader) -> RWAble: ...


class _ReadDaskInternal[S: StorageType](Protocol):
    """Internal protocol for functions registered on :data:`~anndata._io.specs.registry._DASK_REGISTRY`.

    The ``_reader`` is always a :class:`~anndata._io.specs.registry.DaskReader` and
    the optional ``chunks`` kwarg controls dask chunking.
    Return type is one of :data:`DaskDataStructures`.
    """

    def __call__(
        self,
        elem: S,
        *,
        _reader: DaskReader,
        chunks: tuple[int, ...] | None = None,
    ) -> DaskDataStructures: ...


class _ReadBackedInternal[S: StorageType](Protocol):
    """Internal protocol for functions registered on :data:`~anndata._io.specs.registry._BACKED_REGISTRY`.

    The ``_reader`` is always a :class:`~anndata._io.specs.registry.BackedReader`.
    There is no ``chunks`` kwarg — backed reads are zero-copy and unchunked.
    Return type is one of :data:`BackedDataStructures`.
    """

    def __call__(
        self,
        elem: S,
        *,
        _reader: BackedReader,
    ) -> BackedDataStructures: ...


# Legacy alias: the old single protocol covered both modes.
_ReadLazyInternal = _ReadDaskInternal


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


class ReadDask[S](Protocol):
    """Public callable type produced by :meth:`~anndata._io.specs.registry.DaskReader.read_elem`
    after the internal ``_reader`` argument is bound via :func:`functools.partial`.
    Returns one of :data:`DaskDataStructures`.
    """

    def __call__(
        self, elem: S, *, chunks: tuple[int, ...] | None = None
    ) -> DaskDataStructures:
        """Low-level reading function for a dask-backed element.

        Parameters
        ----------
        elem
            The element to read from.
        chunks
            The dask chunk size to be used.
        Returns
        -------
        A dask-backed lazy element (see :data:`DaskDataStructures`).
        """
        ...


class ReadBacked[S](Protocol):
    """Public callable type produced by :meth:`~anndata._io.specs.registry.BackedReader.read_elem`
    after the internal ``_reader`` argument is bound via :func:`functools.partial`.
    Returns one of :data:`BackedDataStructures`.
    """

    def __call__(self, elem: S) -> BackedDataStructures:
        """Low-level reading function for a file-backed (non-dask) element.

        Parameters
        ----------
        elem
            The element to read from.
        Returns
        -------
        A file-backed lazy element (see :data:`BackedDataStructures`).
        """
        ...


# Legacy alias
ReadLazy = ReadDask


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
