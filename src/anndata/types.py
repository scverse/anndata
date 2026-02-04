from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, runtime_checkable

if TYPE_CHECKING:
    from typing import Any, Literal, Self

    from pandas import Index

    from ._core.anndata import AnnData


@runtime_checkable
class ExtensionNamespace(Protocol):
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


@runtime_checkable
class DataFrameLikeIlocIndexer(Protocol):
    """Protocol for iloc-style indexers on DataFrame-like objects.

    Only requires `__getitem__`.
    """

    def __getitem__(self, idx: Any) -> Any: ...


@runtime_checkable
class DataFrameLike(Protocol):
    """Protocol for DataFrame-like objects usable in AnnData.

    See Also
    --------
    :class:`~anndata.experimental.backed.Dataset2D`
        An xarray-based implementation of this protocol.
    """

    @property
    def index(self) -> Index:
        """Row labels of the DataFrame-like object."""
        ...

    @property
    def columns(self) -> Index:
        """Column labels of the DataFrame-like object."""
        ...

    @columns.setter
    def columns(self, v: Any) -> None:
        """Setter for columns"""
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
        *,
        index: Index | None = None,
        axis: Literal[0, 1] | None = 0,
        fill_value: Any = ...,
        **kwargs,
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
