"""Map in-memory dataframe columns onto the on-disk types, see :ref:`type-mapping`."""

from __future__ import annotations

from types import MappingProxyType
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ...compat import PANDAS_STRING_ARRAY_TYPES

if TYPE_CHECKING:
    from collections.abc import Mapping

    from pandas.api.extensions import ExtensionArray


__all__ = [
    "UnsupportedColumnTypeError",
    "add_column_note",
    "normalize_column",
]


class UnsupportedColumnTypeError(TypeError):
    """An object-dtype dataframe column cannot be converted for writing.

    `path` names the column, `kind` is what :func:`~pandas.api.types.infer_dtype`
    made of it, `reason` says why it cannot be written and `hint` what to do about it.
    """

    def __init__(
        self, path: str, reason: str, hint: str = "", kind: str | None = None
    ) -> None:
        self.path, self.reason, self.hint, self.kind = path, reason, hint, kind
        inferred = "" if kind is None else f", inferred as {kind!r}"
        msg = f"Cannot write column {path!r} (dtype object{inferred}): {reason}."
        super().__init__(f"{msg} {hint}" if hint else msg)

    def __reduce__(self) -> tuple[type, tuple[str, str, str, str | None], dict]:
        # for pickling and copying: the default passes only the message to `__init__`
        return type(self), (self.path, self.reason, self.hint, self.kind), self.__dict__


_CONVERT_HINT = (
    "Convert it to a supported dtype before writing, e.g. with `.astype(...)`. "
    "For missing values, use a nullable dtype such as 'Int64', 'boolean', or 'string'."
)


def normalize_column(
    values: np.ndarray | ExtensionArray, *, path: str
) -> np.ndarray | ExtensionArray:
    """\
    Return a writable version of a dataframe column's backing array (``series._values``).

    Only object-dtype columns are converted, see :func:`_object_column_to_writable`;
    everything else is returned unchanged for the IO registry to judge.
    Idempotent, and the input is never modified.
    `path` names the column in error messages.
    """
    if isinstance(values, np.ndarray) and values.dtype == object:
        return _object_column_to_writable(values, path=path)
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
    if isinstance(values, pd.arrays.DatetimeArray | pd.arrays.TimedeltaArray) or (
        isinstance(values, np.ndarray) and values.dtype.kind in "Mm"
    ):
        return _DATETIME_HINT
    if isinstance(values, pd.arrays.FloatingArray):
        return _NULLABLE_FLOAT_HINT
    # `ArrowStringArray` has a writer, other arrow arrays do not
    if isinstance(values, pd.arrays.ArrowExtensionArray) and not isinstance(
        values, tuple(PANDAS_STRING_ARRAY_TYPES)
    ):
        return _ARROW_HINT
    return _GENERIC_HINT


def add_column_note(
    e: BaseException, values: np.ndarray | ExtensionArray, path: str
) -> None:
    """Explain a registry miss for a dataframe column: which column, which dtype, what to do."""
    e.add_note(
        f"Column {path!r} (dtype {values.dtype}) has no on-disk encoding. "
        f"{_column_hint(values)}"
    )


# object columns with one unambiguous on-disk type, by `infer_dtype` kind
_OBJECT_COLUMN_CONVERSIONS: Mapping[
    str, tuple[pd.Int64Dtype | pd.BooleanDtype | pd.Float64Dtype, np.dtype | None]
] = MappingProxyType({
    "integer": (pd.Int64Dtype(), np.dtype(np.int64)),
    "boolean": (pd.BooleanDtype(), np.dtype(np.bool_)),
    # no nullable float on disk: numpy float64 with NaN, as pandas does
    "floating": (pd.Float64Dtype(), None),
    "mixed-integer-float": (pd.Float64Dtype(), None),
})
# 80-bit on x86: even small values may not fit float64
_LONGDOUBLE_IS_WIDER = np.finfo(np.longdouble).nmant > np.finfo(np.float64).nmant


def _object_column_to_writable(
    values: np.ndarray, *, path: str
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
        # the string writers take str/bytes mixtures without missing values
        if (
            len(values) == 0
            or pd.api.types.infer_dtype(values, skipna=False) in {"string", "bytes"}
            or (kind == "mixed" and all(isinstance(v, str | bytes) for v in values))
        ):
            return values
        raise UnsupportedColumnTypeError(
            path,
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
            # ints only: a bool cast turns None into False
            assert numpy_dtype is not None
            return values.astype(numpy_dtype)
        converted = pd.array(values, dtype=nullable_dtype)
    except (TypeError, ValueError, OverflowError) as e:
        raise UnsupportedColumnTypeError(
            path,
            f"converting it to {nullable_dtype.name!r} failed",
            hint=_CONVERT_HINT,
            kind=kind,
        ) from e
    if numpy_dtype is not None:
        return converted if converted.isna().any() else converted.to_numpy(numpy_dtype)
    floats = converted.to_numpy(dtype=np.float64, na_value=np.nan)
    # only magnitudes >= 2**53 can lose precision in float64
    suspect = np.abs(floats) >= 2.0**53
    if _LONGDOUBLE_IS_WIDER and any(
        issubclass(t, np.longdouble) for t in {type(v) for v in values}
    ):
        # extended precision can be rounded at any magnitude
        suspect = ~np.isnan(floats)
    if suspect.any():
        # numpy scalars promote to float64 on comparison, Python ints do not
        original = (
            x.item() if isinstance(x, np.generic) else x for x in values[suspect]
        )
        if any(x != f for x, f in zip(original, floats[suspect].tolist(), strict=True)):
            raise UnsupportedColumnTypeError(
                path,
                "its values cannot be represented exactly as float64",
                hint=_CONVERT_HINT,
                kind=kind,
            )
    return floats
