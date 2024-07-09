"""
Defines some useful types for this library. Should probably be cleaned up before thinking about exporting.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, TypeVar, Union

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy import sparse

from ._core.sparse_dataset import BaseCompressedSparseDataset
from .compat import (
    AwkArray,
    CupyArray,
    CupySparseMatrix,
    DaskArray,
    H5Array,
    H5Group,
    SpArray,
    ZappyArray,
    ZarrArray,
    ZarrGroup,
)

if TYPE_CHECKING:
    from types import MappingProxyType
    from typing import TypeAlias

    from ._io.specs.registry import IOSpec, Reader, Writer

__all__ = [
    "ArrayStorageType",
    "GroupStorageType",
    "StorageType",
]

InMemoryArrayOrScalarType: TypeAlias = Union[
    NDArray,
    np.ma.MaskedArray,
    sparse.spmatrix,
    SpArray,
    H5Array,
    ZarrArray,
    ZappyArray,
    BaseCompressedSparseDataset,
    DaskArray,
    CupyArray,
    CupySparseMatrix,
    AwkArray,
    pd.DataFrame,
    np.number,
    str,
]

ArrayStorageType = Union[ZarrArray, H5Array]
GroupStorageType = Union[ZarrGroup, H5Group]
StorageType = Union[ArrayStorageType, GroupStorageType]

ContravariantInMemoryType = TypeVar(
    "ContravariantInMemoryType",
    bound="InMemoryReadElem",  # noqa: F821
    contravariant=True,
)
CovariantInMemoryType = TypeVar(
    "CovariantInMemoryType",
    bound="InMemoryReadElem",  # noqa: F821
    covariant=True,
)
InvariantInMemoryType = TypeVar("InvariantInMemoryType", bound="InMemoryReadElem")  # noqa: F821


class Read(Protocol[CovariantInMemoryType]):
    def __call__(
        self,
        elem: StorageType,
        _reader: Reader,
        *,
        dataset_kwargs: MappingProxyType,
    ) -> CovariantInMemoryType:
        """Low-level reading function for an element.

        Parameters
        ----------
        elem
            The element to read from.
        _reader
            The :class:`anndata.experimental.Reader` instance.
        dataset_kwargs
            Keyword arguments to be passed to a library-level io function, like `chunks` for :doc:`dask:index`.

        Returns
        -------
            The element read from the store.
        """
        ...


class Write(Protocol[ContravariantInMemoryType]):
    def __call__(
        self,
        f: GroupStorageType,
        k: str,
        v: ContravariantInMemoryType,
        _writer: Writer,
        dataset_kwargs: MappingProxyType,
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
        _writer
            The :class:`anndata.experimental.Writer` instance.
        dataset_kwargs
            Keyword arguments to be passed to a library-level io function, like `chunks` for :doc:`dask:index`.
        """
        ...


class ReadCallback(Protocol[InvariantInMemoryType]):
    def __call__(
        self,
        /,
        read_func: Read[InvariantInMemoryType],
        elem_name: str,
        elem: StorageType,
        iospec: IOSpec,
        *,
        dataset_kwargs: MappingProxyType,
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
        dataset_kwargs
            Keyword arguments to be passed to a library-level io function, like `chunks` for :doc:`dask:index`.

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
        store: GroupStorageType,
        elem_name: str,
        elem: InvariantInMemoryType,
        *,
        iospec: IOSpec,
        dataset_kwargs: MappingProxyType,
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
