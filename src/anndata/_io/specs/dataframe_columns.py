"""\
Map in-memory dataframe columns onto the on-disk type system.

anndata's on-disk format supports a small, closed set of column types
(see :doc:`/fileformat-prose`), while pandas columns are open-ended.
:func:`normalize_column` is the single place where the two meet:
a column is converted only when that is lossless and unambiguous,
and otherwise rejected with an error naming the column and a way out.
"""

from __future__ import annotations

from types import MappingProxyType
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ...compat import PANDAS_STRING_ARRAY_TYPES

if TYPE_CHECKING:
    from collections.abc import Mapping

    from pandas.api.extensions import ExtensionArray


__all__ = ["UnsupportedColumnTypeError", "normalize_column"]


class UnsupportedColumnTypeError(TypeError):
    """A dataframe column has no lossless representation in the on-disk format.

    Attributes
    ----------
    path
        Name or path of the column, as given to :func:`normalize_column`.
    dtype
        The column's in-memory dtype.
    kind
        What :func:`pandas.api.types.infer_dtype` made of an object column, if applicable.
    reason
        Why the column cannot be written.
    hint
        What the user can do about it.
    """

    def __init__(
        self,
        path: str,
        dtype: object,
        reason: str,
        *,
        hint: str = "",
        kind: str | None = None,
    ) -> None:
        self.path = path
        self.dtype = dtype
        self.kind = kind
        self.reason = reason
        self.hint = hint
        inferred = "" if kind is None else f", inferred as {kind!r}"
        msg = f"Cannot write column {path!r} (dtype {dtype}{inferred}): {reason}."
        super().__init__(f"{msg} {hint}" if hint else msg)


_CONVERT_HINT = (
    "Convert it to a supported dtype before writing, e.g. with `.astype(...)`. "
    "For missing values, use a nullable dtype such as 'Int64', 'boolean', or 'string'."
)
_ARROW_HINT = (
    "Convert it before writing, e.g. with "
    "`.convert_dtypes(dtype_backend='numpy_nullable')` or `.astype(...)`."
)
_EXTENSION_HINT = (
    "Convert it to a numpy, nullable ('Int64', 'boolean', 'string'), "
    "or categorical dtype before writing."
)

# Column families written as they are. The registry decides whether a writer exists
# (e.g. datetime arrays currently have none, see scverse/anndata#2238).
_PASSTHROUGH_ARRAY_TYPES: tuple[type[ExtensionArray], ...] = (
    pd.Categorical,
    pd.arrays.IntegerArray,
    pd.arrays.BooleanArray,
    *PANDAS_STRING_ARRAY_TYPES,
    pd.arrays.DatetimeArray,
    pd.arrays.TimedeltaArray,
)


def normalize_column(
    values: np.ndarray | ExtensionArray, *, path: str
) -> np.ndarray | ExtensionArray:
    """\
    Return a writable version of a dataframe column's backing array.

    The conversion is lossless or does not happen:

    * numpy arrays, categoricals, nullable integer/boolean/string arrays,
      and datetime/timedelta arrays are returned unchanged;
    * nullable floating arrays (``Float64``/``Float32``) become numpy floats
      of the same width with ``NaN`` for missing values;
    * object arrays of integers, booleans, or floats (with or without missing values)
      become the matching nullable or numpy array, see :func:`_object_column_to_writable`;
    * pyarrow-backed arrays and every other extension array raise
      :exc:`UnsupportedColumnTypeError`.

    Parameters
    ----------
    values
        The column's backing array, i.e. ``series._values``.
    path
        How to refer to the column in error messages.
    """
    if isinstance(values, np.ndarray):
        if values.dtype == object:
            return _object_column_to_writable(values, path)
        return values
    if isinstance(values, pd.arrays.NumpyExtensionArray):
        return normalize_column(values.to_numpy(), path=path)
    if isinstance(values, _PASSTHROUGH_ARRAY_TYPES):
        return values
    if isinstance(values, pd.arrays.FloatingArray):
        # Nullable floating point has no on-disk encoding, but numpy floats of the
        # same width with NaN hold every value (pandas does the same for non-object input).
        return values.to_numpy(dtype=np.dtype(values.dtype.type), na_value=np.nan)
    if isinstance(values, pd.arrays.ArrowExtensionArray):
        raise UnsupportedColumnTypeError(
            path,
            values.dtype,
            "pyarrow-backed columns are not supported",
            hint=_ARROW_HINT,
        )
    if isinstance(values, pd.api.extensions.ExtensionArray):
        raise UnsupportedColumnTypeError(
            path,
            values.dtype,
            f"{type(values).__name__} columns are not supported",
            hint=_EXTENSION_HINT,
        )
    return values


# Object-dtype columns that pandas infers as integer, boolean, or floating
# (e.g. `np.array([1, 2, None, 3])`) have an unambiguous on-disk representation.
_OBJECT_COLUMN_CONVERSIONS: Mapping[
    str, tuple[pd.Int64Dtype | pd.BooleanDtype | pd.Float64Dtype, np.dtype | None]
] = MappingProxyType({
    "integer": (pd.Int64Dtype(), np.dtype(np.int64)),
    "boolean": (pd.BooleanDtype(), np.dtype(np.bool_)),
    # Nullable floating point is not a supported on-disk type,
    # so always fall back to numpy float64 with NaN, like pandas would for non-object input.
    "floating": (pd.Float64Dtype(), None),
    "mixed-integer-float": (pd.Float64Dtype(), None),
})
# On x86 `np.longdouble` is 80-bit, so even small values may not fit float64 exactly
_LONGDOUBLE_IS_WIDER = np.finfo(np.longdouble).nmant > np.finfo(np.float64).nmant


def _object_column_to_writable(
    values: np.ndarray, path: str
) -> np.ndarray | ExtensionArray:
    """\
    Convert an object-dtype column into something we can write.

    Columns inferred as integer or boolean become the matching nullable pandas array
    if they contain missing values, and the plain numpy dtype otherwise.
    Columns inferred as floating become numpy float64 (missing values become NaN)
    if that loses no precision.
    String/bytes columns without missing values and empty columns are returned unchanged.
    Everything else raises an :exc:`UnsupportedColumnTypeError` naming the column.
    """
    kind = pd.api.types.infer_dtype(values, skipna=True)
    if (conversion := _OBJECT_COLUMN_CONVERSIONS.get(kind)) is None:
        # The string writers handle (mixtures of) str and bytes without missing values
        if (
            len(values) == 0
            or pd.api.types.infer_dtype(values, skipna=False) in {"string", "bytes"}
            or (kind == "mixed" and all(isinstance(v, str | bytes) for v in values))
        ):
            return values
        raise UnsupportedColumnTypeError(
            path,
            "object",
            "it has no supported on-disk representation",
            hint=_CONVERT_HINT,
            kind=kind,
        )
    nullable_dtype, numpy_dtype = conversion
    try:
        if (
            kind == "integer"
            and pd.api.types.infer_dtype(values, skipna=False) == "integer"
        ):
            # No missing values at all, so a direct cast is cheapest and lossless.
            # (Only done for ints: a bool cast would turn None into False, and
            # numpy casts some missing markers, e.g. complex NaN, to an int.)
            assert numpy_dtype is not None
            return values.astype(numpy_dtype)
        converted = pd.array(values, dtype=nullable_dtype)
    except (TypeError, ValueError, OverflowError) as e:
        raise UnsupportedColumnTypeError(
            path,
            "object",
            f"converting it to {nullable_dtype.name!r} failed",
            hint=_CONVERT_HINT,
            kind=kind,
        ) from e
    if numpy_dtype is not None:
        return converted if converted.isna().any() else converted.to_numpy(numpy_dtype)
    floats = converted.to_numpy(dtype=np.float64, na_value=np.nan)
    # float64 holds integers below 2**53 and all float64/float32 values exactly, so only
    # larger magnitudes need checking (missing values are NaN here and fail the cutoff,
    # ±inf pass the comparison below). Extended-precision scalars can be rounded at any
    # magnitude, so if any are present every value is checked.
    suspect = np.abs(floats) >= 2.0**53
    if _LONGDOUBLE_IS_WIDER and any(
        issubclass(t, np.longdouble) for t in {type(v) for v in values}
    ):
        suspect = ~np.isnan(floats)
    if suspect.any():
        # Python compares int and float exactly, numpy scalars do not (they promote to
        # float64), so compare Python scalars
        original = (
            x.item() if isinstance(x, np.generic) else x for x in values[suspect]
        )
        if any(x != f for x, f in zip(original, floats[suspect].tolist(), strict=True)):
            raise UnsupportedColumnTypeError(
                path,
                "object",
                "its values cannot be represented exactly as float64",
                hint=_CONVERT_HINT,
                kind=kind,
            )
    return floats
