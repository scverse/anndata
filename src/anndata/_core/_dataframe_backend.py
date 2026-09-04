"""
Backend-agnostic access to the dataframe stored as AnnData ``obs``/``var``.

A stored frame keeps its axis names either in a mutable pandas-style ``.index``
(:class:`pandas.DataFrame`, :class:`~anndata._core.xarray.Dataset2D`) or, for backends
without an index (polars, pyarrow, …), in an ``obs_names``/``var_names`` column.
:func:`_adapt` is the only place that distinction is made; everything below goes through
the adapter it returns.

Frames that defer evaluation are neither, since narwhals exposes no row order or index for
a :class:`narwhals.LazyFrame`. :func:`_ingest_lazy` converts them on the way in, so
:func:`_adapt` never sees one.
"""

from __future__ import annotations

import re
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Protocol, cast, runtime_checkable

import narwhals as nw
import pandas as pd

from .xarray import Dataset2D

if TYPE_CHECKING:
    from typing import Any, Literal, NoReturn, Self, TypeGuard

    import pyarrow as pa
    from dask.dataframe import DataFrame as DaskDataFrame
    from narwhals.typing import EagerAllowed, IntoBackend

    type Dim = Literal["obs", "var"]


class DataFrameLikeIlocIndexer(Protocol):
    """Positional indexer, as in ``df.iloc[…]``."""

    def __getitem__(self, idx: Any) -> Any: ...


_conforms: dict[tuple[type, type], bool] = {}


class _CachedProtocolMeta(type(Protocol)):
    """Metaclass caching structural checks by type.

    ``isinstance`` against a runtime-checkable protocol re-inspects every member on each
    call, which costs ~1.8 µs for a non-frame and is a quarter of an ``obsm``/``layers``
    assignment. Frames conform as whole types, so the answer is cached per type.

    Never annotate an attribute on this metaclass: under ``from __future__ import
    annotations`` it lands in ``__protocol_attrs__`` and joins the structural check.
    """

    def __instancecheck__(cls, instance: object) -> bool:
        key = (cls, type(instance))
        if (conforms := _conforms.get(key)) is None:
            conforms = _conforms[key] = super().__instancecheck__(instance)
        return conforms


@runtime_checkable
class DataFrameLike(Protocol, metaclass=_CachedProtocolMeta):
    """A native dataframe AnnData can store as ``obs``/``var``.

    Axis names live either in an index (see :class:`IndexedDataFrameLike`) or, for
    index-less backends, in the ``obs_names``/``var_names`` column.
    """

    @property
    def columns(self) -> Any: ...
    @property
    def shape(self) -> tuple[int, int]: ...


@runtime_checkable
class IndexedDataFrameLike(DataFrameLike, Protocol):
    """A stored frame keeping its axis names in a mutable pandas-style index."""

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


_INDEX_PLACEHOLDER = re.compile(r"__index_level_\d+__")

_EAGER_BACKEND_NAMES = ("pandas", "polars", "pyarrow", "modin", "cudf")
_KNOWN_BACKENDS = (*_EAGER_BACKEND_NAMES, "dask")
_EAGER_BACKENDS = frozenset(
    nw.Implementation.from_backend(b) for b in _EAGER_BACKEND_NAMES
)
_DIMS: tuple[Dim, ...] = ("obs", "var")


def _dim_column(dim: Dim) -> str:
    """Name of the column holding the axis names of ``dim``."""
    return f"{dim}_names"


def other_dim(dim: Dim) -> Dim:
    """The dimension a frame moves to when the object is transposed."""
    return "var" if dim == "obs" else "obs"


class _FrameAdapter(ABC):
    """The frame operations depending on where a stored frame keeps its axis names."""

    __slots__ = ("_frame",)

    def __init__(self, frame: DataFrameLike) -> None:
        self._frame = frame

    @abstractmethod
    def index(self, dim: Dim) -> pd.Index:
        """Axis names as a pandas index."""

    @abstractmethod
    def set_index(self, index: pd.Index, dim: Dim) -> DataFrameLike:
        """Replace the axis names."""

    @abstractmethod
    def relabel(self, source: Dim, target: Dim) -> DataFrameLike:
        """Move the axis names from the ``source`` to the ``target`` dimension."""

    @abstractmethod
    def all_columns(self) -> list[str]:
        """Every column, including one holding the axis names."""

    @abstractmethod
    def subset(self, idx: Any) -> DataFrameLike:
        """Positionally subset rows, and columns too when ``idx`` is a 2-tuple."""

    @abstractmethod
    def equals(self, other: object) -> bool:
        """Compare to ``other`` in the terms this representation compares in."""

    @abstractmethod
    def column_backed_dim(self) -> Dim | None:
        """The dimension whose names are held in a column, if any."""

    @abstractmethod
    def native_backend(self) -> nw.Implementation | None:
        """The backend holding the axis names column, or ``None`` when indexed."""

    @abstractmethod
    def to_pandas(self, dim: Dim | None) -> pd.DataFrame:
        """Convert to pandas, moving the axis names back to the index."""

    @abstractmethod
    def to_arrow(self, dim: Dim | None) -> pa.Table:
        """Convert to Arrow, exposing the axis names as a named column."""

    def true_index(self, dim: Dim) -> pd.Index:
        """The axis names AnnData reports, which are :meth:`index` unless deferred."""
        return self.index(dim)

    def column(self, key: str) -> Any:
        """One column, in the array form this representation keeps it in."""
        return self.wrapped()[key].to_numpy()

    def copy(self) -> DataFrameLike:
        return cast("DataFrameLike", self.wrapped().clone().to_native())

    def to_dask(self, dim: Dim | None) -> DaskDataFrame:
        """Convert to a dask DataFrame, keeping the axis names in the index."""
        import dask.dataframe as dd

        return dd.from_pandas(self.to_pandas(dim), npartitions=1)

    def wrapped(self) -> nw.DataFrame[Any]:
        return cast("nw.DataFrame[Any]", nw.from_native(cast("Any", self._frame)))


class _IndexedFrame(_FrameAdapter):
    """Names live in a mutable pandas-style index, so most operations act in place."""

    _frame: IndexedDataFrameLike

    def index(self, dim: Dim) -> pd.Index:
        return self._frame.index

    def set_index(self, index: pd.Index, dim: Dim) -> DataFrameLike:
        self._frame.index = index
        return self._frame

    def relabel(self, source: Dim, target: Dim) -> DataFrameLike:
        return self._frame

    def all_columns(self) -> list[str]:
        return list(self._frame.columns)

    def subset(self, idx: Any) -> DataFrameLike:
        return self._frame.iloc[idx]

    def equals(self, other: object) -> bool:
        return self._frame.equals(other)

    def column_backed_dim(self) -> Dim | None:
        return None

    def native_backend(self) -> nw.Implementation | None:
        return None

    def column(self, key: str) -> Any:
        return self._frame[key].array

    def to_pandas(self, dim: Dim | None) -> pd.DataFrame:
        return self.wrapped().to_pandas()

    def to_arrow(self, dim: Dim | None) -> pa.Table:
        source_name = self._frame.index.name
        source_name = source_name if isinstance(source_name, str) else None
        target_name = _dim_column(dim) if dim is not None else (source_name or "index")
        if source_name != target_name and target_name in self._frame.columns:
            msg = (
                f"Cannot represent index as reserved column {target_name!r}: "
                "an annotation column already uses that name."
            )
            raise ValueError(msg)
        return _name_index_column(
            self.wrapped().to_arrow(), source_name=source_name, index_name=target_name
        )


class _Dataset2DFrame(_IndexedFrame):
    """A :class:`Dataset2D`, whose storage may still be lazy.

    Names, columns and subsetting read the dataset directly, so they stay lazy on the
    inherited implementations. The rest reach the frame through narwhals, which wraps a
    lazy ``Dataset2D`` as a :class:`narwhals.LazyFrame` that those eager operations cannot
    use, so :meth:`wrapped` reads it into memory first. Copying and conversion to dask
    override that to stay lazy.
    """

    _frame: Dataset2D

    def true_index(self, dim: Dim) -> pd.Index:
        return self._frame.true_index

    def column(self, key: str) -> Any:
        return self._frame[key].variable

    def copy(self) -> DataFrameLike:
        return self._frame.copy()

    def to_arrow(self, dim: Dim | None) -> pa.Table:
        # `to_memory` can rename the index (`true_index_dim`, the dummy concat/merge
        # names), so the inherited version would look for the wrong one.
        return self._in_memory().to_arrow(dim)

    def to_dask(self, dim: Dim | None) -> DaskDataFrame:
        return self._frame.to_dask_dataframe()

    def wrapped(self) -> nw.DataFrame[Any]:
        return self._in_memory().wrapped()

    def _in_memory(self) -> _IndexedFrame:
        return _IndexedFrame(self._frame.to_memory())


class _ColumnBackedFrame(_FrameAdapter):
    """Names live in the ``{dim}_names`` column, so operations rebuild the frame."""

    def index(self, dim: Dim) -> pd.Index:
        name = _dim_column(dim)
        wrapped = self.wrapped()
        if name not in wrapped.columns:
            msg = f"Index-less frame is missing required column {name!r}."
            raise ValueError(msg)
        return pd.Index(wrapped[name].to_list(), name=name)

    def set_index(self, index: pd.Index, dim: Dim) -> DataFrameLike:
        wrapped = self.wrapped()
        return cast(
            "DataFrameLike",
            wrapped.with_columns(
                nw.new_series(_dim_column(dim), index, backend=wrapped.implementation)
            ).to_native(),
        )

    def relabel(self, source: Dim, target: Dim) -> DataFrameLike:
        if source == target:
            return self._frame
        source_name, target_name = _dim_column(source), _dim_column(target)
        wrapped = self.wrapped()
        if source_name not in wrapped.columns:
            msg = f"Index-less frame is missing required column {source_name!r}."
            raise ValueError(msg)
        if target_name in wrapped.columns:
            msg = (
                f"Cannot relabel axis names to reserved column {target_name!r}: "
                "an annotation column already uses that name."
            )
            raise ValueError(msg)
        return cast(
            "DataFrameLike", wrapped.rename({source_name: target_name}).to_native()
        )

    def all_columns(self) -> list[str]:
        return list(self.wrapped().columns)

    def subset(self, idx: Any) -> DataFrameLike:
        return cast("DataFrameLike", self.wrapped()[idx].to_native())

    def equals(self, other: object) -> bool:
        wrapped_other = _try_wrap(other)
        if wrapped_other is None:
            return False
        return self.wrapped().to_arrow().equals(wrapped_other.to_arrow())

    def column_backed_dim(self) -> Dim | None:
        columns = self.wrapped().columns
        return next((d for d in _DIMS if _dim_column(d) in columns), None)

    def native_backend(self) -> nw.Implementation | None:
        return self.wrapped().implementation

    def to_pandas(self, dim: Dim | None) -> pd.DataFrame:
        df = self.wrapped().to_pandas()
        if dim is not None and (name := _dim_column(dim)) in df.columns:
            df = df.set_index(name)
        return df

    def to_arrow(self, dim: Dim | None) -> pa.Table:
        return self.wrapped().to_arrow()


def _adapt(frame: DataFrameLike) -> _FrameAdapter:
    """Pick the adapter for where ``frame`` keeps its axis names."""
    if isinstance(frame, Dataset2D):
        return _Dataset2DFrame(frame)
    if isinstance(frame, IndexedDataFrameLike):
        return _IndexedFrame(frame)
    if _try_wrap_lazy(frame) is not None:
        msg = (
            f"{type(frame).__name__} defers evaluation, so AnnData stores it as a Dataset2D "
            "rather than directly. Reaching this means it was stored without being ingested."
        )
        raise TypeError(msg)
    return _ColumnBackedFrame(frame)


def _is_stored_frame(frame: Any) -> TypeGuard[pd.DataFrame | Dataset2D]:
    """Whether ``frame`` is already in a form AnnData stores and can safely mutate."""
    return isinstance(frame, pd.DataFrame | Dataset2D)


def to_backend(
    frame: Any, backend: str | IntoBackend, *, dim: Dim | None = None
) -> Any:
    """Return ``frame`` as a native dataframe of ``backend``.

    ``frame`` is a :class:`DataFrameLike`, or any other frame narwhals can wrap.
    ``backend`` is a name (``"pandas"``, ``"polars"``, ``"pyarrow"``, ``"modin"``,
    ``"cudf"``, ``"dask"``, ``"narwhals"``), the backend module, or a
    :class:`narwhals.Implementation`. ``"narwhals"`` wraps the frame as it is stored,
    giving a :class:`narwhals.LazyFrame` when the storage is lazy.
    Given ``dim``, the axis names are kept: in the index for the pandas-like backends
    (``"pandas"``, ``"dask"``), in the ``{dim}_names`` column otherwise. Only ``"dask"``
    stays lazy; every other backend is eager, so a lazy frame is read into memory.
    """
    if backend == "narwhals":
        return nw.from_native(cast("Any", frame))
    impl = nw.Implementation.from_backend(backend)
    if impl is nw.Implementation.DASK:
        return _adapt(frame).to_dask(dim)
    if impl not in _EAGER_BACKENDS:
        msg = f"Unsupported DataFrame backend {backend!r}; expected one of {(*_KNOWN_BACKENDS, 'narwhals')}."
        raise ValueError(msg)
    adapter = _adapt(frame)
    if impl is nw.Implementation.PANDAS:
        return adapter.to_pandas(dim)
    return nw.from_arrow(
        adapter.to_arrow(dim), backend=cast("IntoBackend[EagerAllowed]", impl)
    ).to_native()


def unwrap_narwhals(frame: Any) -> Any:
    """Return the native dataframe behind a narwhals wrapper, or ``frame`` unchanged.

    AnnData stores native dataframes, so a :class:`narwhals.DataFrame` or
    :class:`narwhals.LazyFrame` handed to it is unwrapped on the way in. Both satisfy
    :class:`DataFrameLike` structurally, so without this they would be stored as the
    wrapper, which nothing downstream, IO least of all, knows what to do with.
    """
    return (
        frame.to_native() if isinstance(frame, nw.DataFrame | nw.LazyFrame) else frame
    )


def _from_backend(frame: Any, *, dim: Dim | None = None) -> DataFrameLike:
    """Convert ``frame`` to one AnnData can store and safely mutate.

    Frames that already are one (:class:`pandas.DataFrame`, :class:`Dataset2D`) pass
    through; any other becomes pandas, with a ``{dim}_names`` column moved to the index.
    """
    if _is_stored_frame(frame):
        return frame
    return _ColumnBackedFrame(frame).to_pandas(dim)


def _ingest_as_pandas(frame: Any, *, dim: Dim | None = None) -> DataFrameLike | None:
    """Convert any dataframe to pandas, giving ``None`` when ``frame`` is not one.

    For the writer's fallback on an unregistered type, which needs something pandas
    understands and reads every value anyway, so a lazy frame is collected.
    """
    frame = unwrap_narwhals(frame)
    if _is_stored_frame(frame):
        return frame
    if (lazy := _try_wrap_lazy(frame)) is not None:
        frame = lazy.collect().to_native()
    if _try_wrap(frame) is None:
        return None
    return _from_backend(frame, dim=dim)


def _ingest_lazy(frame: Any, *, dim: Dim) -> DataFrameLike | None:
    """Convert a lazy dataframe to its stored form, giving ``None`` if it is not lazy.

    The operations here need row order and an index, which narwhals does not offer over a
    :class:`narwhals.LazyFrame`. Dask frames become a :class:`Dataset2D`, which is lazy and
    does offer both; other lazy backends have no array form to keep, so they are collected.
    """
    # A stored `Dataset2D` also wraps as a `LazyFrame`, but it is the destination here.
    if _is_stored_frame(frame):
        return None
    lazy = _try_wrap_lazy(frame)
    if lazy is None:
        return None
    if lazy.implementation is nw.Implementation.DASK:
        return Dataset2D.from_dask_dataframe(frame, dim_name=_dim_column(dim))
    return cast("DataFrameLike", lazy.collect().to_native())


def _ingest_axis_frame(frame: Any, *, dim: Dim) -> DataFrameLike | None:
    """Ingest the dataframe *declaring* the ``dim`` axis, giving ``None`` if it is not one.

    Unlike :func:`_ingest_as_pandas`, an index-less frame stays in its own backend, and
    gains a positional ``{dim}_names`` column when it has none, since declaring an axis is
    what gets to name it.
    """
    if _is_stored_frame(frame):
        return frame
    if (ingested := _ingest_lazy(frame, dim=dim)) is not None:
        if _is_stored_frame(ingested):
            return ingested
        frame = ingested
    wrapped = _try_wrap(frame)
    if wrapped is None:
        return None
    if isinstance(frame, IndexedDataFrameLike):
        return _from_backend(frame, dim=dim)
    name = _dim_column(dim)
    if name not in wrapped.columns:
        index = pd.RangeIndex(wrapped.shape[0]).astype(str)
        wrapped = wrapped.with_columns(
            nw.new_series(name, index, backend=wrapped.implementation)
        )
    return cast("DataFrameLike", wrapped.to_native())


def axis_index(frame: DataFrameLike, *, dim: Dim) -> pd.Index:
    """Axis names, from a native index or the ``{dim}_names`` column."""
    return _adapt(frame).index(dim)


def true_axis_index(frame: DataFrameLike, *, dim: Dim) -> pd.Index:
    """The axis names AnnData reports to users.

    The same as :func:`axis_index`, except that a lazily read frame answers with the names
    it keeps out of memory (cell names, say) rather than whatever in-memory index stands in
    for them. Use this when handing names back to a caller and :func:`axis_index` when
    operating on the frame itself.
    """
    return _adapt(frame).true_index(dim)


def set_axis_index(frame: DataFrameLike, index: pd.Index, *, dim: Dim) -> DataFrameLike:
    """Replace the axis names, in place for indexed frames and by rebuilding otherwise."""
    return _adapt(frame).set_index(index, dim)


def relabel_axis(frame: DataFrameLike, *, source: Dim, target: Dim) -> DataFrameLike:
    """Relabel a ``{dim}_names`` column when a frame moves between AnnData axes."""
    return _adapt(frame).relabel(source, target)


def subset_frame(frame: DataFrameLike, index: Any) -> DataFrameLike:
    """Positionally subset rows, and columns too when ``index`` is a 2-tuple."""
    return _adapt(frame).subset(index)


def frame_equal(left: DataFrameLike, right: object) -> bool:
    """Compare two native frames in the terms their representation compares in."""
    return _adapt(left).equals(right)


def copy_frame(frame: DataFrameLike) -> DataFrameLike:
    """Copy a native frame without reading lazy storage."""
    return _adapt(frame).copy()


def frame_annotation_columns(
    frame: DataFrameLike, *, dim: Dim | None = None
) -> list[str]:
    """The annotation columns, excluding one holding the axis names."""
    columns = _adapt(frame).all_columns()
    if dim is None:
        return columns
    name = _dim_column(dim)
    return [column for column in columns if column != name]


def column_backed_dim(frame: DataFrameLike) -> Dim | None:
    """The dimension whose names ``frame`` holds in a column, if it does."""
    return _adapt(frame).column_backed_dim()


def native_backend(frame: DataFrameLike) -> nw.Implementation | None:
    """The backend of an index-less frame, or ``None`` when the names are an index.

    A backend here is one an index-less representation can be rebuilt in, so indexed frames
    answer ``None`` even though pandas is a backend in its own right.
    """
    return _adapt(frame).native_backend()


def frame_column(frame: DataFrameLike, key: str) -> Any:
    """One annotation column, in the array form its representation keeps it in.

    Indexed frames answer with the pandas array behind the column, a lazy frame with its
    xarray variable, and index-less frames with NumPy, since their columns have no array
    form AnnData can otherwise hold on to.
    """
    return _adapt(frame).column(key)


def supports_categorical_ops(frame: object) -> TypeGuard[pd.DataFrame]:
    """Whether pandas' ``.cat`` accessor operations apply to ``frame``.

    Categorical upkeep is pandas-only, as narwhals exposes no cross-backend ``.cat``
    equivalent. Implicit maintenance (dropping unused categories when building a view,
    converting strings on write) skips frames that fail this check, so non-pandas
    annotations keep their backend dtypes; explicit user-invoked operations go through
    :func:`reject_categorical_op` rather than silently doing nothing.
    """
    return isinstance(frame, pd.DataFrame)


def reject_categorical_op(operation: str) -> NoReturn:
    """Reject ``operation`` when it is invoked on a non-pandas annotation."""
    msg = (
        f"{operation} is currently supported only for pandas-backed annotations. "
        "Convert and reassign the annotation as pandas first."
    )
    raise NotImplementedError(msg)


def _name_index_column(
    table: pa.Table, *, source_name: str | None, index_name: str
) -> pa.Table:
    """Rename the Arrow field holding an indexed frame's axis names to ``index_name``."""
    renamed = [
        index_name
        if (source_name is not None and c == source_name)
        or (source_name is None and _INDEX_PLACEHOLDER.fullmatch(c))
        else c
        for c in table.column_names
    ]
    return table.rename_columns(renamed) if renamed != table.column_names else table


def _try_wrap(frame: Any) -> nw.DataFrame[Any] | None:
    """Wrap ``frame`` when it is a native eager dataframe.

    The callers ask whether AnnData can store ``frame``, and only the native frames
    narwhals wraps qualify, so a narwhals object gives ``None``. Those reach the ingest
    functions already unwrapped by :func:`unwrap_narwhals`.
    """
    wrapped = nw.from_native(frame, pass_through=True)
    return (
        wrapped if wrapped is not frame and isinstance(wrapped, nw.DataFrame) else None
    )


def _try_wrap_lazy(frame: Any) -> nw.LazyFrame[Any] | None:
    """Wrap ``frame`` when its backend defers evaluation (dask, polars lazy)."""
    wrapped = nw.from_native(frame, pass_through=True)
    return wrapped if isinstance(wrapped, nw.LazyFrame) else None
