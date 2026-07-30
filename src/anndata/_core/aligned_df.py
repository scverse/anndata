from __future__ import annotations

from collections.abc import Mapping
from functools import singledispatch
from typing import TYPE_CHECKING, cast

import narwhals as nw
import pandas as pd
from pandas.api.types import is_string_dtype

from .._settings import settings
from .._warnings import ImplicitModificationWarning
from ..compat import XDataset, pandas_as_str
from ..utils import warn
from ._dataframe_backend import _ingest_axis_frame, axis_index, set_axis_index
from .xarray import Dataset2D

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Any, Literal

    from ._dataframe_backend import DataFrameLike


@singledispatch
def _gen_dataframe(
    anno: Any,
    index_names: Iterable[str],
    *,
    source: Literal["X", "shape"],
    attr: Literal["obs", "var"],
    length: int | None = None,
) -> DataFrameLike:
    """Coerce `anno` to a frame that can be stored as `obs`/`var`.

    This fallback handles any dataframe narwhals supports: index-less frames such as
    polars and pyarrow stay as they are, dask frames become a lazy `Dataset2D`, and other
    pandas-like frames are converted to pandas.
    """
    if (frame := _ingest_axis_frame(anno, dim=attr)) is None:
        msg = f"Cannot convert {type(anno)} to {attr} DataFrame"
        raise ValueError(msg)
    _check_length(frame, source=source, attr=attr, length=length)
    return _coerce_str_index(frame, attr=attr)


@_gen_dataframe.register(Mapping)
@_gen_dataframe.register(type(None))
def _gen_dataframe_mapping(
    anno: Mapping[str, Any] | None,
    index_names: Iterable[str],
    *,
    source: Literal["X", "shape"],
    attr: Literal["obs", "var"],
    length: int | None = None,
) -> pd.DataFrame:
    df = _dataframe_from_mapping(anno, list(index_names), length=length)
    _check_length(df, source=source, attr=attr, length=length)
    return _coerce_str_index(df, attr=attr)


@_gen_dataframe.register(pd.DataFrame)
def _gen_dataframe_df(
    anno: pd.DataFrame,
    index_names: Iterable[str],
    *,
    source: Literal["X", "shape"],
    attr: Literal["obs", "var"],
    length: int | None = None,
) -> pd.DataFrame:
    _reject_pandas_multiindex(anno)
    _check_length(anno, source=source, attr=attr, length=length)
    return _coerce_str_index(_coerce_pandas_columns(anno), attr=attr)


@_gen_dataframe.register(nw.DataFrame)
@_gen_dataframe.register(nw.LazyFrame)
def _gen_dataframe_narwhals(
    anno: nw.DataFrame[Any] | nw.LazyFrame[Any],
    index_names: Iterable[str],
    *,
    source: Literal["X", "shape"],
    attr: Literal["obs", "var"],
    length: int | None = None,
) -> DataFrameLike:
    # Dispatch again on the native frame, so a pandas-backed one keeps its index.
    return _gen_dataframe(
        anno.to_native(), index_names, source=source, attr=attr, length=length
    )


@_gen_dataframe.register(Dataset2D)
def _gen_dataframe_dataset2d(
    anno: Dataset2D,
    index_names: Iterable[str],
    *,
    source: Literal["X", "shape"],
    attr: Literal["obs", "var"],
    length: int | None = None,
) -> Dataset2D:
    # Left as it is: the index came off disk and is already what it should be.
    _check_length(anno, source=source, attr=attr, length=length)
    return anno


@_gen_dataframe.register(XDataset)
def _gen_dataframe_xdataset(
    anno: XDataset,
    index_names: Iterable[str],
    *,
    source: Literal["X", "shape"],
    attr: Literal["obs", "var"],
    length: int | None = None,
) -> Dataset2D:
    return _gen_dataframe_dataset2d(
        Dataset2D(anno), index_names, source=source, attr=attr, length=length
    )


def _check_length(
    anno: DataFrameLike,
    *,
    source: Literal["X", "shape"],
    attr: Literal["obs", "var"],
    length: int | None,
) -> None:
    if length is not None and length != anno.shape[0]:
        raise _mk_df_error(source, attr, length, anno.shape[0])


def _coerce_str_index[T: DataFrameLike](anno: T, *, attr: Literal["obs", "var"]) -> T:
    """Stringify the axis names where declaring them as obs/var names requires it."""
    index = axis_index(anno, dim=attr)
    if not _needs_str_index(index):
        return anno
    warn("Transforming to str index.", ImplicitModificationWarning)
    return cast("T", set_axis_index(anno, pandas_as_str(index), dim=attr))


def _dataframe_from_mapping(
    anno: Mapping[str, Any] | None, index_names: list[str], *, length: int | None
) -> pd.DataFrame:
    """Build a pandas dataframe from a mapping (or `None`), mirroring the constructor."""
    if anno is None or len(anno) == 0:
        anno = {}

    def mk_index(length: int) -> pd.Index:
        return pd.RangeIndex(0, length, name=None).astype(str)

    for index_name in index_names:
        if index_name not in anno:
            continue
        df = pd.DataFrame(
            anno,
            index=anno[index_name],
            columns=[k for k in anno if k != index_name],
        )
        break
    else:
        df = pd.DataFrame(
            anno,
            index=None if length is None else mk_index(length),
            columns=None if anno else pd.array([], dtype="str"),
        )

    if length is None:
        df.index = mk_index(len(df))
    return df


def _reject_pandas_multiindex(anno: pd.DataFrame) -> None:
    """Disallow a `pd.MultiIndex` row index on obs/var declaration, unless opted out."""
    if isinstance(anno.index, pd.MultiIndex) and settings.restrict_index_types:
        msg = (
            "pandas.MultiIndex not supported as index for obs or var on declaration.\n"
            "You can set `obs_names` manually although most operations after will error or convert to str.\n"
            "You can also opt out of `settings.restrict_index_types` which will allow pandas.MultiIndex."
        )
        raise ValueError(msg)


def _needs_str_index(index: pd.Index) -> bool:
    """Whether declaring `index` as obs/var names requires coercing it to strings."""
    return (
        settings.restrict_index_types and not is_string_dtype(index[~index.isna()])
    ) or pd.api.types.is_integer_dtype(index)


def _coerce_pandas_columns(anno: pd.DataFrame) -> pd.DataFrame:
    """Shallow-copy so the caller’s frame is untouched, and stringify empty columns."""
    anno = anno.copy(deep=False)
    if not len(anno.columns):
        anno.columns = pandas_as_str(anno.columns)
    return anno


def _mk_df_error(
    source: Literal["X", "shape"],
    attr: Literal["obs", "var"],
    expected: int,
    actual: int,
):
    what = "row" if attr == "obs" else "column"
    if source == "X":
        msg = (
            f"Observations annot. `{attr}` must have as many rows as `X` has {what}s "
            f"({expected}), but has {actual} rows."
        )
    else:
        msg = (
            f"`shape` is inconsistent with `{attr}` "
            f"({actual} {what}s instead of {expected})"
        )
    return ValueError(msg)
