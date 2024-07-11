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
    from collections.abc import Mapping
    from typing import Any, TypeAlias

    from ._io.specs.registry import IOSpec, Reader, Writer
    from .compat import H5File

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
    "ContravariantInMemoryType", bound="InMemoryReadElem", contravariant=True
)
CovariantInMemoryType = TypeVar(
    "CovariantInMemoryType", bound="InMemoryReadElem", covariant=True
)
InvariantInMemoryType = TypeVar("InvariantInMemoryType", bound="InMemoryReadElem")


class _ReadInternal(Protocol[CovariantInMemoryType]):
    def __call__(
        self,
        elem: StorageType | H5File,
        *,
        _reader: Reader,
    ) -> CovariantInMemoryType: ...


class Read(Protocol[CovariantInMemoryType]):
    def __call__(
        self,
        elem: StorageType | H5File,
    ) -> CovariantInMemoryType:
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


class ReadCallback(Protocol[InvariantInMemoryType]):
    def __call__(
        self,
        /,
        read_func: Read[InvariantInMemoryType],
        elem_name: str,
        elem: StorageType,
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


if TYPE_CHECKING:
    # Needs to be at the end because Sphinx’s type import suffers from circular imports
    from ._io.specs.registry import InMemoryReadElem
