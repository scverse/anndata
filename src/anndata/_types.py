"""
Defines some useful types for this library. Should probably be cleaned up before thinking about exporting.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal, Protocol, runtime_checkable

from .compat import H5Array, H5Group, ZarrArray, ZarrGroup

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any, Self

    import pandas as pd

    from anndata._core.xarray import Dataset2D

    from . import typing
    from ._io.specs.registry import (
        IOSpec,
        LazyDataStructures,
        LazyReader,
        Reader,
        Writer,
    )


__all__ = [
    "ArrayStorageType",
    "DataFrameLike",
    "DataFrameLikeIlocIndexer",
    "GroupStorageType",
    "StorageType",
    "_ReadInternal",
    "_ReadLazyInternal",
    "_WriteInternal",
]

type ArrayStorageType = ZarrArray | H5Array
type GroupStorageType = ZarrGroup | H5Group
type StorageType = ArrayStorageType | GroupStorageType

# circumvent https://github.com/tox-dev/sphinx-autodoc-typehints/issues/580
type S = StorageType
type RWAble = "typing.RWAble"


@runtime_checkable
class DataFrameLikeIlocIndexer(Protocol):
    """Protocol for iloc-style indexers on DataFrame-like objects.

    This protocol defines the minimal interface for positional-based indexing
    that AnnData requires. Both :class:`pandas.DataFrame` and
    :class:`~anndata.experimental.backed.Dataset2D` provide compatible
    ``iloc`` accessors.

    Examples
    --------
    >>> import pandas as pd
    >>> from anndata._types import DataFrameLikeIlocIndexer
    >>> df = pd.DataFrame({"a": [1, 2, 3]})
    >>> isinstance(df.iloc, DataFrameLikeIlocIndexer)
    True
    """

    def __getitem__(self, idx: Any) -> Self: ...


@runtime_checkable
class DataFrameLike(Protocol):
    """Protocol for DataFrame-like objects usable in AnnData.

    This runtime-checkable protocol defines the minimal DataFrame API that
    AnnData uses internally for ``obs``, ``var``, and similar dataframe-like
    data containers. Any class implementing this protocol can be used as a
    drop-in replacement for :class:`pandas.DataFrame` in these contexts.

    The required interface includes:

    - :attr:`index`: Row labels as a :class:`pandas.Index`
    - :attr:`columns`: Column labels as a :class:`pandas.Index`
    - :attr:`shape`: Tuple of (n_rows, n_columns)
    - :attr:`iloc`: Positional indexer returning a :class:`DataFrameLikeIlocIndexer`
    - :meth:`reindex`: Method to reindex rows

    Examples
    --------
    >>> import pandas as pd
    >>> from anndata._types import DataFrameLike
    >>> df = pd.DataFrame({"a": [1, 2, 3]})
    >>> isinstance(df, DataFrameLike)
    True

    See Also
    --------
    :class:`~anndata.experimental.backed.Dataset2D`
        An xarray-based implementation of this protocol.
    """

    @property
    def index(self) -> pd.Index:
        """Row labels of the DataFrame-like object."""
        ...

    @property
    def columns(self) -> pd.Index:
        """Column labels of the DataFrame-like object."""
        ...

    @property
    def shape(self) -> tuple[int, int]:
        """Shape of the DataFrame-like object as (n_rows, n_columns)."""
        ...

    @property
    def iloc(self) -> DataFrameLikeIlocIndexer:
        """Positional indexer for the DataFrame-like object."""
        ...

    def reindex(
        self,
        index: pd.Index | None = None,
        axis: Literal[0] = 0,
        fill_value: Any = ...,
    ) -> Self:
        """Reindex the DataFrame-like object to match a new index.

        Parameters
        ----------
        index
            New index to conform to.
        axis
            Axis to reindex along (only 0 is supported).
        fill_value
            Value to use for missing values.

        Returns
        -------
        Reindexed DataFrame-like object.
        """
        ...


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
            Keyword arguments to be passed to a library-level io function, like `chunks` for :doc:`zarr:index`.
        """
        ...


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
            Keyword arguments to be passed to a library-level io function, like `chunks` for :doc:`zarr:index`.
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
