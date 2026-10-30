"""Map in-memory dataframe columns onto the on-disk types, see :ref:`type-mapping`."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from pandas.api.types import infer_dtype

from ...compat import PANDAS_3, PANDAS_STRING_ARRAY_TYPES

if TYPE_CHECKING:
    from pandas.api.extensions import ExtensionArray


__all__ = [
    "UnsupportedColumnTypeError",
    "add_column_note",
    "normalize_column",
]


class UnsupportedColumnTypeError(TypeError):
    """An object-dtype dataframe column cannot be converted for writing."""


_CONVERT_HINT = (
    "Convert it to a supported dtype before writing, e.g. with `.astype(...)`. "
    "For missing values, use a nullable dtype such as 'Int64', 'boolean', or 'string'"
    + (
        "."
        if PANDAS_3
        else " (with `anndata.settings.allow_write_nullable_strings = True`)."
    )
)


def _unsupported(name: str, kind: str, reason: str) -> UnsupportedColumnTypeError:
    return UnsupportedColumnTypeError(
        f"Cannot write column {name!r} (dtype object, inferred as {kind!r}): "
        f"{reason}. {_CONVERT_HINT}"
    )


def normalize_column(
    values: np.ndarray | ExtensionArray, *, name: str
) -> np.ndarray | ExtensionArray:
    """Make an object-dtype column writable, return any other column unchanged."""
    if isinstance(values, np.ndarray) and values.dtype == object:
        return _object_column_to_writable(values, name=name)
    return values


_DATETIME_HINT = (
    "datetime/timedelta columns are not supported yet, see scverse/anndata#2238. "
    "Convert the column before writing, e.g. with `.astype(str)` or `.astype('int64')`."
)
_NULLABLE_FLOAT_HINT = (
    "There is no on-disk encoding for nullable floats yet. "
    "Call `.to_numpy(np.float64, na_value=np.nan)` on the column "
    "if collapsing NA into NaN is acceptable."
)
_ARROW_HINT = (
    "pyarrow-backed columns are not supported. Convert the column before writing, "
    "e.g. with `.convert_dtypes(dtype_backend='numpy_nullable')` or `.astype(...)`."
)
_GENERIC_HINT = (
    "Convert the column to a numpy, nullable integer/boolean, string, "
    "or categorical column before writing."
)


def _column_hint(values: np.ndarray | ExtensionArray) -> str:
    if isinstance(values, pd.arrays.DatetimeArray | pd.arrays.TimedeltaArray):
        return _DATETIME_HINT
    if isinstance(values, pd.arrays.FloatingArray):
        return _NULLABLE_FLOAT_HINT
    if isinstance(values, pd.arrays.ArrowExtensionArray) and not isinstance(
        values, tuple(PANDAS_STRING_ARRAY_TYPES)
    ):
        return _ARROW_HINT
    return _GENERIC_HINT


def add_column_note(
    e: BaseException, values: np.ndarray | ExtensionArray, name: str
) -> None:
    """Explain a registry miss for a dataframe column: which column, which dtype, what to do."""
    e.add_note(
        f"Column {name!r} (dtype {values.dtype}) has no on-disk encoding. "
        f"{_column_hint(values)}"
    )


def _to_integer_or_boolean(
    values: np.ndarray, kind: str, *, name: str
) -> np.ndarray | ExtensionArray:
    """`values` as `int64`/`bool`, or as `Int64`/`boolean` if any are missing."""
    nullable_dtype: pd.Int64Dtype | pd.BooleanDtype
    numpy_dtype: np.dtype
    if kind == "integer":
        nullable_dtype, numpy_dtype = pd.Int64Dtype(), np.dtype(np.int64)
    else:
        nullable_dtype, numpy_dtype = pd.BooleanDtype(), np.dtype(np.bool_)
    missing = pd.isna(values)
    try:
        if not missing.any():
            return values.astype(numpy_dtype)
        filled = np.where(missing, numpy_dtype.type(0), values).astype(numpy_dtype)
    except (TypeError, ValueError, OverflowError) as e:
        msg = f"converting it to {nullable_dtype.name!r} failed"
        raise _unsupported(name, kind, msg) from e
    return nullable_dtype.construct_array_type()(filled, missing)


def _to_float64(values: np.ndarray, kind: str, *, name: str) -> np.ndarray:
    """`values` as `float64` with NaN for missing values, rounding integers above 2**53."""
    try:
        floats = pd.array(values, dtype=pd.Float64Dtype())
    except (TypeError, ValueError, OverflowError) as e:
        raise _unsupported(name, kind, "converting it to 'Float64' failed") from e
    return floats.to_numpy(np.float64, na_value=np.nan)


def _object_column_to_writable(
    values: np.ndarray, *, name: str
) -> np.ndarray | ExtensionArray:
    """Convert an object-dtype column for writing, see :ref:`type-mapping`."""
    kind = infer_dtype(values, skipna=True)
    if kind in {"integer", "boolean"}:
        return _to_integer_or_boolean(values, kind, name=name)
    if kind in {"floating", "mixed-integer-float"}:
        return _to_float64(values, kind, name=name)
    if (
        len(values) == 0
        or (kind in {"string", "bytes"} and not pd.isna(values).any())
        or (kind == "mixed" and all(isinstance(v, str | bytes) for v in values))
    ):
        return values
    raise _unsupported(name, kind, "it has no supported on-disk representation")
