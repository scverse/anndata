"""
Defines some useful types for this library. Should probably be cleaned up before thinking about exporting.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, TypeVar, runtime_checkable

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

    from anndata._core.anndata import AnnData

    from ._io.specs.registry import DaskReader, IOSpec, Reader, Writer
    from .compat import DaskArray

__all__ = [
    "ArrayStorageType",
    "GroupStorageType",
    "StorageType",
    "_ReadInternal",
    "_ReadDaskInternal",
    "_WriteInternal",
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


class _ReadDaskInternal(Protocol[SCon]):
    def __call__(
        self, elem: SCon, *, _reader: DaskReader, chunks: tuple[int, ...] | None = None
    ) -> DaskArray: ...


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


class ReadDask(Protocol[SCon]):
    def __call__(
        self, elem: SCon, *, chunks: tuple[int, ...] | None = None
    ) -> DaskArray:
        """Low-level reading function for a dask element.

        Parameters
        ----------
        elem
            The element to read from.
        chunks
            The chunk size to be used.
        Returns
        -------
            The dask element read from the store.
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


NS = TypeVar("NS", covariant=True)


@runtime_checkable
class ExtensionNamespace(Protocol[NS]):
    """Protocol for extension namespaces.

    Enforces that the namespace initializer accepts a class with the proper `__init__` method.
    Protocol's can't enforce that the `__init__` accepts the correct types. See
    `_check_namespace_signature` for that. This is mainly useful for static type
    checking with mypy and IDEs.
    """

    def __init__(self, adata: AnnData) -> None:
        """
        Used to enforce the correct signature for extension namespaces.
        """
        ...
