from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING

import pandas as pd
from pandas.api.types import is_string_dtype

from .._settings import settings
from .._warnings import ImplicitModificationWarning
from ..compat import XDataset, pandas_as_str
from ..utils import warn
from ._dataframe_backend import axis_index, set_axis_index, try_ensure_axis_frame
from .xarray import Dataset2D

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Any, Literal

    from ._dataframe_backend import DataFrameLike


def _gen_dataframe(
    anno: Any,
    index_names: Iterable[str],
    *,
    source: Literal["X", "shape"],
    attr: Literal["obs", "var"],
    length: int | None = None,
) -> DataFrameLike:
    """Coerce ``anno`` to a stored ``obs``/``var`` frame (a :class:`DataFrameLike`).

    Accepts ``None``/mappings (built into a pandas frame), an :class:`xarray.Dataset` /
    :class:`Dataset2D`, or a Narwhals-native eager frame. Supported index-less frames such as
    Polars and PyArrow remain native; mutable pandas-compatible foreign frames normalize to
    pandas.
    """
    index_names = list(index_names)
    if not isinstance(anno, pd.DataFrame | Dataset2D):
        if isinstance(anno, XDataset):
            anno = Dataset2D(anno)
        elif anno is None or isinstance(anno, Mapping):
            anno = _dataframe_from_mapping(anno, index_names, length=length)
        else:
            coerced = try_ensure_axis_frame(anno, index_name=index_names[0])
            if coerced is None:
                msg = f"Cannot convert {type(anno)} to {attr} DataFrame"
                raise ValueError(msg)
            anno = coerced
    if isinstance(anno, pd.DataFrame):
        _reject_pandas_multiindex(anno)
    if length is not None and length != anno.shape[0]:
        raise _mk_df_error(source, attr, length, anno.shape[0])
    if isinstance(anno, pd.DataFrame):
        anno = _coerce_pandas_df(anno)
    elif not isinstance(anno, Dataset2D):
        index_name = index_names[0]
        index = axis_index(anno, index_name=index_name)
        if (
            settings.restrict_index_types and not is_string_dtype(index[~index.isna()])
        ) or pd.api.types.is_integer_dtype(index):
            warn("Transforming to str index.", ImplicitModificationWarning)
            anno = set_axis_index(anno, pandas_as_str(index), index_name=index_name)
    return anno


def _dataframe_from_mapping(
    anno: Mapping[str, Any] | None, index_names: list[str], *, length: int | None
) -> pd.DataFrame:
    """Build a pandas DataFrame from a mapping (or ``None``), mirroring the constructor."""
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
    """Disallow a pandas ``MultiIndex`` row index on obs/var declaration (unless opted out)."""
    if isinstance(anno.index, pd.MultiIndex) and settings.restrict_index_types:
        msg = (
            "pandas.MultiIndex not supported as index for obs or var on declaration.\n"
            "You can set `obs_names` manually although most operations after will error or convert to str.\n"
            "You can also opt out of `settings.restrict_index_types` which will allow pandas.MultiIndex."
        )
        raise ValueError(msg)


def _coerce_pandas_df(anno: pd.DataFrame) -> pd.DataFrame:
    """pandas-only index/column hygiene applied on obs/var declaration."""
    anno = anno.copy(deep=False)
    if (
        settings.restrict_index_types
        and not is_string_dtype(anno.index[~anno.index.isna()])
    ) or pd.api.types.is_integer_dtype(anno.index):
        msg = "Transforming to str index."
        warn(msg, ImplicitModificationWarning)
        anno.index = pandas_as_str(anno.index)
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
