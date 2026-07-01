from __future__ import annotations

import re
from typing import TYPE_CHECKING, Protocol, cast, runtime_checkable

import narwhals as nw
from narwhals._utils import Implementation

if TYPE_CHECKING:
    from typing import Any, Self

    import pandas as pd
    from narwhals._pandas_like.dataframe import PandasLikeDataFrame
    from narwhals.typing import EagerAllowed, IntoBackend
    from narwhals.utils import Version

    from .xarray import Dataset2D


@runtime_checkable
class DataFrameLikeIlocIndexer(Protocol):
    """Positional indexer, as in ``df.iloc[...]``."""

    def __getitem__(self, idx: Any) -> Any: ...


@runtime_checkable
class DataFrameLike(Protocol):
    """Structural contract an AnnData ``obs``/``var`` must satisfy to be stored as-is.

    :class:`pandas.DataFrame` and :class:`~anndata._core.xarray.Dataset2D` conform. Index-less
    frames (polars, pyarrow, cuDF, ...) do not — bring them in with :func:`from_backend`.
    """

    def __len__(self) -> int: ...
    @property
    def index(self) -> pd.Index: ...
    @property
    def columns(self) -> pd.Index: ...
    @columns.setter
    def columns(self, value: Any) -> None: ...
    @property
    def shape(self) -> tuple[int, int]: ...
    @property
    def iloc(self) -> DataFrameLikeIlocIndexer: ...
    def reindex(self, *, index: Any = None, axis: Any = 0, **kwargs: Any) -> Self: ...


NATIVE_PACKAGE = "anndata"


def is_native(native_object: object, /) -> bool:
    """Return whether ``native_object`` is a :class:`Dataset2D`."""
    from .xarray import Dataset2D

    return isinstance(native_object, Dataset2D)


class Dataset2DNamespace:
    """Routes a :class:`Dataset2D` to its pandas realisation as a compliant frame.

    A ``Dataset2D`` realises to pandas (:meth:`~Dataset2D.to_memory`), which narwhals already
    supports, so we route through its ``PandasLikeDataFrame`` instead of reimplementing the
    compliant protocol. Wrapping realises the frame in memory (eager only). Row labels stay the
    index, recoverable via :func:`narwhals.maybe_get_index`.

    ``PandasLikeDataFrame.from_native`` reads ``_implementation`` and ``_version`` off its
    ``context``, so this namespace doubles as that context.
    """

    _implementation: Implementation = Implementation.PANDAS

    def __init__(self, *, version: Version) -> None:
        self._version = version

    def from_native(self, native_object: Dataset2D, /) -> PandasLikeDataFrame:
        from narwhals._pandas_like.dataframe import PandasLikeDataFrame

        return PandasLikeDataFrame.from_native(native_object.to_memory(), context=self)


def __narwhals_namespace__(version: Version) -> Dataset2DNamespace:
    """Return the compliant namespace narwhals uses to wrap a :class:`Dataset2D`."""
    return Dataset2DNamespace(version=version)


# pandas serialises an unnamed index as an "__index_level_<n>__" placeholder column.
_INDEX_PLACEHOLDER = re.compile(r"__index_level_\d+__")

_KNOWN_BACKENDS = ("pandas", "polars", "pyarrow", "modin", "cudf")
_EAGER_BACKENDS = frozenset(Implementation.from_backend(b) for b in _KNOWN_BACKENDS)


def to_backend(
    frame: Any, backend: str | IntoBackend, *, index_name: str = "index"
) -> Any:
    """Outgest: return ``frame`` as a native DataFrame of ``backend``.

    ``frame`` is a ``DataFrameLike`` (:class:`pandas.DataFrame` or :class:`Dataset2D`); any other
    narwhals-wrappable frame also works. ``backend`` is an eager backend: a string (``"pandas"``,
    ``"polars"``, ``"pyarrow"``, ``"modin"``, ``"cudf"``), the backend module, or a
    :class:`narwhals.Implementation`. When the index is unnamed it rides along as an ``index_name``
    column (non-pandas backends only). The result is not ``DataFrameLike`` for index-less backends.
    """
    impl = Implementation.from_backend(backend)
    if impl not in _EAGER_BACKENDS:
        msg = f"Unsupported DataFrame backend {backend!r}; expected one of {_KNOWN_BACKENDS}."
        raise ValueError(msg)

    nwf = nw.from_native(frame)
    if impl is Implementation.PANDAS:
        return nwf.to_pandas()
    table = _name_unnamed_index(nwf.to_arrow(), index_name)
    return nw.from_arrow(
        table, backend=cast("IntoBackend[EagerAllowed]", impl)
    ).to_native()


def from_backend(frame: Any, *, index_name: str | None = None) -> DataFrameLike:
    """Ingest: normalize ``frame`` to a stored ``DataFrameLike``.

    Frames that already conform (:class:`pandas.DataFrame`, :class:`Dataset2D`) pass through.
    A frame from another backend (polars, pyarrow, cuDF, ...) is converted to pandas; if
    ``index_name`` names one of its columns, that column becomes the index (the inverse of
    :func:`to_backend`).
    """
    if isinstance(frame, DataFrameLike):
        return frame
    df = nw.from_native(frame).to_pandas()
    if index_name is not None and index_name in df.columns:
        df = df.set_index(index_name)
    return df


def try_from_backend(
    frame: Any, *, index_name: str | None = None
) -> DataFrameLike | None:
    """Like :func:`from_backend`, but return ``None`` instead of raising when ``frame`` isn't a frame.

    Used where the input might not be a dataframe at all (obs/var construction, the writer's
    normalize-on-dispatch-miss path): an eager frame from another backend (polars, pyarrow,
    cuDF, ...) becomes a pandas ``DataFrameLike``; a Series, ndarray, mapping, etc. returns
    ``None`` so the caller can handle it.
    """
    wrapped = nw.from_native(frame, pass_through=True)
    if wrapped is frame or not isinstance(wrapped, nw.DataFrame):
        return None
    return from_backend(frame, index_name=index_name)


def _name_unnamed_index(table, index_name: str):
    """Rename a pandas ``__index_level_<n>__`` placeholder column to ``index_name``."""
    renamed = [
        index_name if _INDEX_PLACEHOLDER.fullmatch(c) else c for c in table.column_names
    ]
    return table.rename_columns(renamed) if renamed != table.column_names else table
