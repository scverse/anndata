from __future__ import annotations

import re
from typing import TYPE_CHECKING, Protocol, cast, runtime_checkable

import narwhals as nw
import pandas as pd

if TYPE_CHECKING:
    from typing import Any, Self

    from narwhals._pandas_like.dataframe import PandasLikeDataFrame
    from narwhals.typing import EagerAllowed, IntoBackend
    from narwhals.utils import Version

    from .xarray import Dataset2D


class DataFrameLikeIlocIndexer(Protocol):
    """Positional indexer, as in ``df.iloc[...]``."""

    def __getitem__(self, idx: Any) -> Any: ...


@runtime_checkable
class DataFrameLike(Protocol):
    """Common structural contract for a native frame stored as AnnData ``obs``/``var``.

    Axis identity is exposed either through an ``index`` (see :class:`IndexedDataFrameLike`) or
    through the canonical ``obs_names``/``var_names`` column for index-less backends.
    """

    @property
    def columns(self) -> Any: ...
    @property
    def shape(self) -> tuple[int, int]: ...


@runtime_checkable
class IndexedDataFrameLike(DataFrameLike, Protocol):
    """A stored frame whose axis identity is represented by a mutable pandas-style index."""

    @property
    def index(self) -> pd.Index: ...
    @index.setter
    def index(self, value: Any) -> None: ...
    @property
    def columns(self) -> Any: ...
    @columns.setter
    def columns(self, value: Any) -> None: ...
    @property
    def iloc(self) -> DataFrameLikeIlocIndexer: ...
    def reindex(self, *, index: Any = None, axis: Any = 0, **kwargs: Any) -> Self: ...
    def equals(self, other: object) -> bool: ...


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

    _implementation: nw.Implementation = nw.Implementation.PANDAS

    def __init__(self, *, version: Version) -> None:
        self._version = version

    def from_native(self, native_object: Dataset2D, /) -> PandasLikeDataFrame:
        from narwhals._pandas_like.dataframe import PandasLikeDataFrame

        return PandasLikeDataFrame.from_native(native_object.to_memory(), context=self)


def __narwhals_namespace__(version: Version) -> Dataset2DNamespace:
    """Return the compliant namespace narwhals uses to wrap a :class:`Dataset2D`."""
    return Dataset2DNamespace(version=version)


_INDEX_PLACEHOLDER = re.compile(r"__index_level_\d+__")

_KNOWN_BACKENDS = ("pandas", "polars", "pyarrow", "modin", "cudf")
_EAGER_BACKENDS = frozenset(nw.Implementation.from_backend(b) for b in _KNOWN_BACKENDS)
_CANONICAL_AXIS_NAMES = ("obs_names", "var_names")


def to_backend(
    frame: Any, backend: str | IntoBackend, *, index_name: str | None = None
) -> Any:
    """Outgest: return ``frame`` as a native DataFrame of ``backend``.

    ``frame`` is a ``DataFrameLike`` (an indexed pandas/Dataset2D frame or a supported
    index-less native frame); any other narwhals-wrappable frame also works. ``backend`` is an
    eager backend: a string (``"pandas"``,
    ``"polars"``, ``"pyarrow"``, ``"modin"``, ``"cudf"``), the backend module, or a
    :class:`narwhals.Implementation`. When the index is unnamed it rides along as an ``index_name``
    column (non-pandas backends only).
    """
    impl = nw.Implementation.from_backend(backend)
    if impl not in _EAGER_BACKENDS:
        msg = f"Unsupported DataFrame backend {backend!r}; expected one of {_KNOWN_BACKENDS}."
        raise ValueError(msg)

    nwf = nw.from_native(frame)
    if impl is nw.Implementation.PANDAS:
        df = nwf.to_pandas()
        if (
            not isinstance(frame, IndexedDataFrameLike)
            and index_name is not None
            and index_name in df.columns
        ):
            df = df.set_index(index_name)
        return df
    table = nwf.to_arrow()
    if isinstance(frame, IndexedDataFrameLike):
        source_name = frame.index.name
        source_name = source_name if isinstance(source_name, str) else None
        target_name = index_name or source_name or "index"
        if source_name != target_name and target_name in frame.columns:
            msg = (
                f"Cannot represent index as reserved column {target_name!r}: "
                "an annotation column already uses that name."
            )
            raise ValueError(msg)
        table = _name_index_column(
            table, source_name=source_name, index_name=target_name
        )
    return nw.from_arrow(
        table, backend=cast("IntoBackend[EagerAllowed]", impl)
    ).to_native()


def from_backend(frame: Any, *, index_name: str | None = None) -> DataFrameLike:
    """Ingest: normalize ``frame`` to a stored ``DataFrameLike``.

    Indexed frames that AnnData can safely mutate (:class:`pandas.DataFrame`,
    :class:`Dataset2D`) pass through. A frame from another backend is converted to pandas; if
    ``index_name`` names one of its columns, that column becomes the index (the inverse of
    :func:`to_backend`).
    """
    from .xarray import Dataset2D

    if isinstance(frame, pd.DataFrame | Dataset2D):
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
    if _try_from_native(frame) is None:
        return None
    return from_backend(frame, index_name=index_name)


def try_ensure_axis_frame(frame: Any, *, index_name: str) -> DataFrameLike | None:
    """Return a native 2-D frame with axis identity available as an index or canonical column."""
    from .xarray import Dataset2D

    if isinstance(frame, pd.DataFrame | Dataset2D):
        return frame
    wrapped = _try_from_native(frame)
    if wrapped is None:
        return None
    if isinstance(frame, IndexedDataFrameLike):
        return from_backend(frame, index_name=index_name)
    if index_name not in wrapped.columns:
        index = pd.RangeIndex(wrapped.shape[0]).astype(str)
        wrapped = wrapped.with_columns(
            nw.new_series(index_name, index, backend=wrapped.implementation)
        )
    return cast("DataFrameLike", wrapped.to_native())


def axis_index(frame: DataFrameLike, *, index_name: str) -> pd.Index:
    """Return logical axis identity from a native index or canonical identity column."""
    if isinstance(frame, IndexedDataFrameLike):
        return frame.index
    wrapped = _from_native(frame)
    if index_name not in wrapped.columns:
        msg = f"Index-less frame is missing required identity column {index_name!r}."
        raise ValueError(msg)
    return pd.Index(wrapped[index_name].to_list(), name=index_name)


def set_axis_index(
    frame: DataFrameLike, index: pd.Index, *, index_name: str
) -> DataFrameLike:
    """Replace logical axis identity, mutating indexed frames and rebuilding immutable frames."""
    if isinstance(frame, IndexedDataFrameLike):
        frame.index = index
        return frame
    wrapped = _from_native(frame)
    return cast(
        "DataFrameLike",
        wrapped.with_columns(
            nw.new_series(index_name, index, backend=wrapped.implementation)
        ).to_native(),
    )


def relabel_axis_identity(
    frame: DataFrameLike, *, source_name: str, target_name: str
) -> DataFrameLike:
    """Relabel column-backed identity when a frame moves between AnnData axes."""
    if isinstance(frame, IndexedDataFrameLike) or source_name == target_name:
        return frame
    wrapped = _from_native(frame)
    if source_name not in wrapped.columns:
        msg = f"Index-less frame is missing required identity column {source_name!r}."
        raise ValueError(msg)
    if target_name in wrapped.columns:
        msg = (
            f"Cannot relabel axis identity to reserved column {target_name!r}: "
            "an annotation column already uses that name."
        )
        raise ValueError(msg)
    return cast(
        "DataFrameLike",
        wrapped.rename({source_name: target_name}).to_native(),
    )


def subset_frame(frame: DataFrameLike, index: Any) -> DataFrameLike:
    """Positionally subset rows/columns of either indexed or index-less native frames."""
    if isinstance(frame, IndexedDataFrameLike):
        return frame.iloc[index]
    return cast("DataFrameLike", _from_native(frame)[index].to_native())


def frame_equal(left: DataFrameLike, right: object) -> bool:
    """Compare native frames while preserving like-representation semantics for indexed frames."""
    if isinstance(left, IndexedDataFrameLike):
        return left.equals(right)
    wrapped_right = _try_from_native(right)
    if wrapped_right is None:
        return False
    return _from_native(left).to_arrow().equals(wrapped_right.to_arrow())


def copy_frame(frame: DataFrameLike) -> DataFrameLike:
    """Copy a native frame without materializing lazy storage."""
    from .xarray import Dataset2D

    if isinstance(frame, Dataset2D):
        return frame.copy()
    return cast("DataFrameLike", _from_native(frame).clone().to_native())


def frame_annotation_columns(
    frame: DataFrameLike, *, index_name: str | None = None
) -> list[str]:
    """Return logical annotation columns, excluding column-backed axis identity."""
    columns = (
        frame.columns
        if isinstance(frame, IndexedDataFrameLike)
        else _from_native(frame).columns
    )
    return [column for column in columns if column != index_name]


def column_backed_axis_name(frame: DataFrameLike) -> str | None:
    """Return a canonical identity field when ``frame`` represents identity as a column."""
    if isinstance(frame, IndexedDataFrameLike):
        return None
    columns = _from_native(frame).columns
    return next((name for name in _CANONICAL_AXIS_NAMES if name in columns), None)


def _name_index_column(table, *, source_name: str | None, index_name: str):
    """Give an indexed frame's Arrow identity field its canonical column name."""
    renamed = [
        index_name
        if (source_name is not None and c == source_name)
        or (source_name is None and _INDEX_PLACEHOLDER.fullmatch(c))
        else c
        for c in table.column_names
    ]
    return table.rename_columns(renamed) if renamed != table.column_names else table


def _from_native(frame: DataFrameLike) -> nw.DataFrame[Any]:
    return cast("nw.DataFrame[Any]", nw.from_native(cast("Any", frame)))


def _try_from_native(frame: Any) -> nw.DataFrame[Any] | None:
    wrapped = nw.from_native(frame, pass_through=True)
    return (
        wrapped if wrapped is not frame and isinstance(wrapped, nw.DataFrame) else None
    )
