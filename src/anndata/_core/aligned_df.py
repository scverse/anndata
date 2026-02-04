from __future__ import annotations

from collections.abc import Mapping
from functools import singledispatch
from typing import TYPE_CHECKING

import pandas as pd
from pandas.api.types import is_string_dtype

from .._warnings import ImplicitModificationWarning
from ..compat import XDataset, pandas_as_str
from ..types import DataFrameLike
from ..utils import warn
from .xarray import Dataset2D

if TYPE_CHECKING:
    from collections.abc import Iterable
    from typing import Any, Literal


@singledispatch
def _gen_dataframe(
    anno: Any,
    index_names: Iterable[str],
    *,
    source: Literal["X", "shape"],
    attr: Literal["obs", "var"],
    length: int | None = None,
) -> DataFrameLike:  # pragma: no cover
    # Check if anno satisfies the DataFrameLike protocol
    # This allows any DataFrameLike-compliant object to be used as obs/var
    if isinstance(anno, DataFrameLike):
        if length is not None and anno.shape[0] != length:
            raise _mk_df_error(source, attr, length, anno.shape[0])
        return anno
    msg = f"Cannot convert {type(anno)} to {attr} DataFrame"
    raise ValueError(msg)


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
    if anno is None or len(anno) == 0:
        anno = {}

    def mk_index(l: int) -> pd.Index:
        return pd.RangeIndex(0, l, name=None).astype(str)

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
    elif length != len(df):
        raise _mk_df_error(source, attr, length, len(df))
    return df


@_gen_dataframe.register(pd.DataFrame)
def _gen_dataframe_df(
    anno: pd.DataFrame,
    index_names: Iterable[str],
    *,
    source: Literal["X", "shape"],
    attr: Literal["obs", "var"],
    length: int | None = None,
):
    if isinstance(anno.index, pd.MultiIndex):
        msg = (
            "pandas.MultiIndex not supported as index for obs or var on declaration.\n\
            You can set `obs_names` manually although most operations after will error or convert to str.\n\
            This behavior will likely be clarified in a future breaking release."
        )
        raise ValueError(msg)
    if length is not None and length != len(anno):
        raise _mk_df_error(source, attr, length, len(anno))
    anno = anno.copy(deep=False)
    if not is_string_dtype(anno.index[~anno.index.isna()]):
        msg = "Transforming to str index."
        warn(msg, ImplicitModificationWarning)
        anno.index = pandas_as_str(anno.index)
    if not len(anno.columns):
        anno.columns = pandas_as_str(anno.columns)
    return anno


@_gen_dataframe.register(pd.Series)
@_gen_dataframe.register(pd.Index)
def _gen_dataframe_1d(
    anno: pd.Series | pd.Index,
    index_names: Iterable[str],
    *,
    source: Literal["X", "shape"],
    attr: Literal["obs", "var"],
    length: int | None = None,
):
    msg = f"Cannot convert {type(anno)} to {attr} DataFrame"
    raise ValueError(msg)


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


@_gen_dataframe.register(Dataset2D)
def _gen_dataframe_xr(
    anno: Dataset2D,
    index_names: Iterable[str],
    *,
    source: Literal["X", "shape"],
    attr: Literal["obs", "var"],
    length: int | None = None,
):
    return anno


@_gen_dataframe.register(XDataset)
def _gen_dataframe_xdataset(
    anno: XDataset,
    index_names: Iterable[str],
    *,
    source: Literal["X", "shape"],
    attr: Literal["obs", "var"],
    length: int | None = None,
):
    return Dataset2D(anno)
