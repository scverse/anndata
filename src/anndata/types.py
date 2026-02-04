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

    def __getitem__(self, idx: Any) -> Any: ...


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
