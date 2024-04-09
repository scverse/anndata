from __future__ import annotations

import warnings
from functools import singledispatch
from typing import TYPE_CHECKING, Any, Literal

import pandas as pd
from pandas.api.types import is_string_dtype

from .._warnings import ImplicitModificationWarning

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping


@singledispatch
def _gen_dataframe(
    anno: Mapping[str, Any],
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
            columns=[k for k in anno.keys() if k != index_name],
        )
        break
    else:
        df = pd.DataFrame(
            anno,
            index=None if length is None else mk_index(length),
            columns=None if len(anno) else [],
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
    if length is not None and length != len(anno):
        raise _mk_df_error(source, attr, length, len(anno))
    anno = anno.copy(deep=False)
    if not is_string_dtype(anno.index):
        warnings.warn("Transforming to str index.", ImplicitModificationWarning)
        anno.index = anno.index.astype(str)
    if not len(anno.columns):
        anno.columns = anno.columns.astype(str)
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
    raise ValueError(f"Cannot convert {type(anno)} to {attr} DataFrame")


def _mk_df_error(
    source: Literal["X", "shape"],
    attr: Literal["obs", "var"],
    expected: int,
    actual: int,
):
    if source == "X":
        what = "row" if attr == "obs" else "column"
        msg = (
            f"Observations annot. `{attr}` must have as many rows as `X` has {what}s "
            f"({expected}), but has {actual} rows."
        )
    else:
        msg = (
            f"`shape` is inconsistent with `{attr}` "
            "({actual} {what}s instead of {expected})"
        )
    return ValueError(msg)
