"""
Utility functions for HTML representation.

This module provides:
- Serialization checking using the anndata IO registry
- String-to-category warning detection
- Color list detection and validation
- HTML escaping and sanitization
- Memory size formatting
"""

from __future__ import annotations

import html
import re
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .._repr_constants import (
    DICT_PREVIEW_KEYS,
    DICT_PREVIEW_KEYS_LARGE,
    LIST_PREVIEW_ITEMS,
    STRING_INLINE_LIMIT,
)
from .._settings import settings

if TYPE_CHECKING:
    from collections.abc import Iterable

    import pandas as pd

    from anndata import AnnData

    from .registry import FormatterContext


def _check_serializable_single(obj: object) -> tuple[bool, str]:
    """Check if a single (non-container) object is serializable."""
    # Handle None
    if obj is None:
        return True, ""

    # Use the actual IO registry
    try:
        from .._io.specs.registry import _REGISTRY

        _REGISTRY.get_spec(obj)  # type: ignore[arg-type]
        return True, ""
    except (KeyError, TypeError):
        pass

    # Check for basic Python types that are serializable
    if isinstance(obj, (bool, int, float, str, bytes)):
        return True, ""

    # Check numpy scalar types
    if isinstance(obj, np.generic):
        return True, ""

    return (
        False,
        f"Type '{type(obj).__module__}.{type(obj).__name__}' has no registered writer",
    )


_SERIALIZABLE_FAST_PATH_MIN_LEN = 1000


def _is_plain_array(seq: list | tuple) -> bool:
    """Whether ``seq`` converts to a numeric/bool/string numpy array (no objects)."""
    try:
        return np.asarray(seq).dtype.kind in "biufcUS"
    except Exception:  # noqa: BLE001
        # e.g. ragged nested lists
        return False


def is_serializable(
    obj: object,
    *,
    _depth: int = 0,
    _max_depth: int = 10,
) -> tuple[bool, str]:
    """
    Check if an object can be serialized to H5AD/Zarr.

    Uses the actual anndata IO registry to check if a type has a registered writer.
    For containers (dict, list), recursively checks all elements.

    Parameters
    ----------
    obj
        Object to check
    _depth
        Current recursion depth (internal)
    _max_depth
        Maximum recursion depth to prevent infinite loops

    Returns
    -------
    tuple of (is_serializable, reason_if_not)
    """
    if _depth > _max_depth:
        return False, "Maximum nesting depth exceeded"

    # Check containers recursively
    if isinstance(obj, dict):
        for k, v in obj.items():
            ok, reason = is_serializable(v, _depth=_depth + 1, _max_depth=_max_depth)
            if not ok:
                return False, f"Key '{k}': {reason}"
        return True, ""

    if isinstance(obj, (list, tuple)):
        # Fast path: long lists of scalars become a plain numpy array on write,
        # so one vectorized dtype check replaces walking millions of elements.
        if len(obj) <= _SERIALIZABLE_FAST_PATH_MIN_LEN or not _is_plain_array(obj):
            for i, v in enumerate(obj):
                ok, reason = is_serializable(
                    v, _depth=_depth + 1, _max_depth=_max_depth
                )
                if not ok:
                    return False, f"Index {i}: {reason}"
        return True, ""

    return _check_serializable_single(obj)


def should_warn_string_column(
    series: pd.Series, n_unique: int | None
) -> tuple[bool, str]:
    """
    Check if a string column will be auto-converted to categorical on save.

    This replicates the logic from AnnData.strings_to_categoricals()
    (see _core/anndata.py:1249-1259):
    - Column must be string type (infer_dtype == "string")
    - Number of unique values must be less than total values

    Parameters
    ----------
    series
        Pandas Series to check
    n_unique
        Pre-computed nunique value (None if skipped due to unique_limit or lazy)

    Returns
    -------
    tuple of (should_warn, warning_message)
    """
    # Can't check if n_unique wasn't computed
    if n_unique is None:
        return False, ""

    from pandas.api.types import infer_dtype

    # Same check as AnnData.strings_to_categoricals()
    dtype_str = infer_dtype(series)
    if dtype_str != "string":
        return False, ""

    n_total = len(series)
    if n_unique < n_total:
        return (
            True,
            (
                f"String column ({n_unique} unique). "
                f"Will be converted to categorical on save."
            ),
        )

    return False, ""


def _is_color_string(s: str) -> bool:
    """Check if a string looks like a color value."""
    if s.startswith("#"):
        return True
    s_lower = s.lower()
    if s_lower in _NAMED_COLORS:
        return True
    return s_lower.startswith(("rgb(", "rgba("))


def sanitize_css_color(color: str) -> str | None:  # noqa: PLR0911
    """
    Sanitize a color string for safe use in CSS style attributes.

    Returns the sanitized color if valid, or None if the color is invalid
    or potentially dangerous (contains CSS injection attempts).

    This is critical for security - color values go into style attributes
    and must not allow CSS injection (e.g., "red; background-image: url(...)").

    Note: Multiple returns are intentional for clarity in validating different
    color formats (hex, named, rgb/rgba).

    Parameters
    ----------
    color
        The color string to sanitize

    Returns
    -------
    The sanitized color string, or None if invalid/unsafe
    """
    if not isinstance(color, str):
        return None

    color = color.strip()
    if not color:
        return None

    # Length limit to prevent DoS via very long strings
    if len(color) > 50:
        return None

    # Hex colors: #RGB, #RRGGBB, or #RRGGBBAA (strict whitelist)
    if color.startswith("#"):
        hex_part = color[1:]
        if len(hex_part) in (3, 4, 6, 8) and all(
            c in "0123456789abcdefABCDEF" for c in hex_part
        ):
            return color
        return None

    # Named colors - must exactly match a known CSS color name (whitelist)
    color_lower = color.lower()
    if color_lower in _NAMED_COLORS:
        return color_lower

    # rgb() and rgba() - WHITELIST approach: only allow safe characters
    if color_lower.startswith("rgb"):
        # Only these characters can appear in valid rgb/rgba colors
        safe_chars = set("rgbaRGBA0123456789(),. %")
        if not all(c in safe_chars for c in color):
            return None
        # Validate rgb/rgba format strictly with regex
        rgb_pattern = r"^rgba?\(\s*\d{1,3}%?\s*,\s*\d{1,3}%?\s*,\s*\d{1,3}%?\s*(,\s*(0|1|0?\.\d+))?\s*\)$"
        if re.match(rgb_pattern, color_lower):
            return color
        return None

    # Reject everything else - no hsl(), var(), url(), expression(), etc.
    return None


def is_color_list(key: str, value: object) -> bool:
    """
    Check if a value is a color list following the *_colors convention.

    Parameters
    ----------
    key
        The key name (should end with '_colors')
    value
        The value to check

    Returns
    -------
    True if this appears to be a color list
    """
    if not isinstance(key, str) or not key.endswith("_colors"):
        return False
    if not isinstance(value, (list, np.ndarray, tuple)):
        return False
    # Empty list is valid
    if len(value) == 0:
        return True
    # Check first element
    first = value[0]
    return isinstance(first, str) and _is_color_string(first)


def _get_categories_from_column(col: object) -> list:
    """
    Get categories from a categorical column.

    Works for both pandas Series (.cat.categories) and xarray DataArray
    (dtype.categories). Returns empty list if categories cannot be extracted.
    """
    try:
        # Pandas Series
        if hasattr(col, "cat"):
            return list(col.cat.categories)

        # xarray DataArray or other objects with CategoricalDtype
        if hasattr(col, "dtype") and hasattr(col.dtype, "categories"):
            return list(col.dtype.categories)
    except Exception as e:  # noqa: BLE001
        from .._warnings import warn

        warn(
            f"Failed to extract categories from column: {type(e).__name__}: {e}",
            UserWarning,
        )

    return []


def get_categories_for_display(
    col: object,
    context: FormatterContext,
    *,
    is_lazy: bool,
) -> tuple[list, bool, int | None]:
    """
    Get categories for a column, handling lazy loading appropriately.

    Parameters
    ----------
    col
        The column to get categories from
    context
        FormatterContext with display settings
    is_lazy
        Whether this is a lazy column (from read_lazy())

    Returns
    -------
    tuple of (categories_list, was_truncated, n_categories)
        categories_list: List of category values
        was_truncated: True if categories were truncated for lazy columns
        n_categories: Total number of categories (if known)
    """
    if is_lazy:
        from .lazy import get_lazy_categories

        return get_lazy_categories(col, context)

    # Non-lazy categorical - use unified accessor
    categories = _get_categories_from_column(col)
    return categories, False, len(categories) if categories else None


@dataclass(frozen=True)
class ColumnColors:
    """Colors stored for a categorical column in ``uns["{column}_colors"]``."""

    n_total: int
    """Number of colors stored."""

    head: list[str]
    """The first (up to ``limit``) colors, as strings."""


def get_column_colors(
    adata: AnnData, column_name: str, *, limit: int
) -> ColumnColors | None:
    """
    Read the colors for a categorical column from ``uns``, if there are any.

    Only the first ``limit`` colors are materialized: for lazy AnnData the
    color array is a dask array, so this reads just the displayed part.

    Parameters
    ----------
    adata
        AnnData object (or object with ``.uns``; objects without return None)
    column_name
        Name of the column (colors key will be ``"{column_name}_colors"``)
    limit
        Maximum number of colors to load

    Returns
    -------
    The color count and the first ``limit`` colors, or None if there are no
    colors or they are not a sequence.
    """
    try:
        uns = adata.uns
        color_key = f"{column_name}_colors"
        if color_key not in uns:
            return None
        colors = uns[color_key]
        if isinstance(colors, str | bytes) or not hasattr(colors, "__len__"):
            return None
        if limit == 0:
            return ColumnColors(n_total=len(colors), head=[])
        head = colors[:limit]
        if hasattr(head, "compute"):  # dask (lazy AnnData)
            head = head.compute()
        return ColumnColors(n_total=len(colors), head=[str(c) for c in head])
    except Exception:  # noqa: BLE001
        # Missing/broken uns or an unsliceable value: treat as "no colors"
        return None


def count_invalid_colors(colors: Iterable[object]) -> int:
    """
    Count colors that fail sanitization.

    Parameters
    ----------
    colors
        Color values to check

    Returns
    -------
    Number of colors that fail sanitize_css_color validation
    """
    return sum(1 for c in colors if sanitize_css_color(str(c)) is None)


def format_invalid_colors_warning(invalid_count: int, *, has_more: bool = False) -> str:
    """
    Format a warning message for invalid colors.

    Parameters
    ----------
    invalid_count
        Number of invalid colors found
    has_more
        If True, adds "+" suffix to indicate more unchecked colors

    Returns
    -------
    Formatted warning message like "2 invalid colors" or "2+ invalid colors"
    """
    suffix = "+" if has_more else ""
    s = "s" if invalid_count > 1 else ""
    return f"{invalid_count}{suffix} invalid color{s}"


def format_index_preview(index: pd.Index, preview_n: int = 5) -> str:
    """Format a preview of a pandas Index.

    Shows first and last items with ellipsis in between for long indices.
    Handles bytes index values (from older h5ad files) by decoding them.

    Parameters
    ----------
    index
        The pandas Index to preview
    preview_n
        Number of items to show at the start and end

    Returns
    -------
    Comma-separated preview string, or ``<em>empty</em>`` for empty indices.
    """
    n = len(index)
    if n == 0:
        return "<em>empty</em>"

    def _format_value(x: object) -> str:
        """Format a single index value, decoding bytes if needed."""
        if isinstance(x, bytes):
            try:
                return x.decode("utf-8")
            except UnicodeDecodeError:
                return x.decode("latin-1")
        return str(x)

    if n <= preview_n * 2:
        items = [escape_html(_format_value(x)) for x in index]
    else:
        first = [escape_html(_format_value(x)) for x in index[:preview_n]]
        last = [escape_html(_format_value(x)) for x in index[-preview_n:]]
        items = [*first, "...", *last]

    return ", ".join(items)


def escape_html(text: str) -> str:
    """Escape HTML special characters and replace null bytes.

    Null bytes in user data (e.g., column names like ``"null\\x00byte"``)
    break HTML parsers and cause truncated rendering. They are replaced
    with the Unicode replacement character U+FFFD.
    """
    return html.escape(str(text).replace("\x00", "\ufffd"))


def format_memory_size(size_bytes: float) -> str:
    """Format memory size in human-readable form."""
    if size_bytes < 0:
        return "Unknown"

    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(size_bytes) < 1024:
            if unit == "B":
                return f"{int(size_bytes)} {unit}"
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024

    return f"{size_bytes:.1f} PB"


def format_number(n: object) -> str:
    """Format a number with thousand separators.

    Accepts int, float, or anything else (e.g. fallback values like "?"),
    which is returned as its string form.
    """
    if not isinstance(n, int | float | np.integer | np.floating):
        return str(n)
    if isinstance(n, float | np.floating):
        if not np.isfinite(n):
            return str(n)
        if n == int(n):
            n = int(n)
        else:
            return f"{n:,.2f}"
    return f"{n:,}"


def get_anndata_version() -> str:
    """Get the anndata version string."""
    try:
        from importlib.metadata import PackageNotFoundError, version

        return version("anndata")
    except PackageNotFoundError:
        return "unknown"


def is_view(obj: object) -> bool:
    """Check if an object is a view (for AnnData-like objects)."""
    try:
        return getattr(obj, "is_view", False)
    except Exception:  # noqa: BLE001
        return False


def is_backed(obj: object) -> bool:
    """Check if an object is backed (for AnnData-like objects)."""
    try:
        return getattr(obj, "isbacked", False)
    except Exception:  # noqa: BLE001
        return False


def get_backing_info(obj: object) -> dict[str, bool | str | None]:
    """Get information about backing for an AnnData-like object."""
    try:
        if not is_backed(obj):
            return {"backed": False}

        filename = str(getattr(obj, "filename", None) or "")
        info: dict[str, bool | str | None] = {
            "backed": True,
            "filename": filename,
        }

        # Try to get file status
        file_obj = getattr(obj, "file", None)
        if file_obj is not None:
            info["is_open"] = getattr(file_obj, "is_open", None)

        # Detect format from filename
        if filename:
            if filename.endswith(".h5ad"):
                info["format"] = "H5AD"
            elif ".zarr" in filename:
                info["format"] = "Zarr"
            else:
                info["format"] = "Unknown"

        return info
    except Exception:  # noqa: BLE001
        return {"backed": False}


def _load_css_colors() -> frozenset[str]:
    """Load CSS named colors from static file.

    The colors are loaded from static/css_colors.txt which contains the
    147 CSS3 named colors. This file can be easily updated if needed.

    Returns
    -------
    frozenset of lowercase color names
    """
    from functools import cache
    from importlib.resources import files

    @cache
    def _load() -> frozenset[str]:
        content = (
            files("anndata._repr.static")
            .joinpath("css_colors.txt")
            .read_text(encoding="utf-8")
        )
        colors = set()
        for line in content.splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                colors.add(line.lower())
        return frozenset(colors)

    return _load()


# CSS named colors for color detection in _is_color_string().
# Loaded from static/css_colors.txt - see that file for the full list.
# Colors can also be specified as hex (#RGB, #RRGGBB), rgb(), or rgba().
_NAMED_COLORS = _load_css_colors()


# -----------------------------------------------------------------------------
# Value preview functions
# -----------------------------------------------------------------------------


def preview_string(value: str, max_len: int) -> str:
    """Preview a string value."""
    if len(value) <= max_len:
        return f'"{value}"'
    return f'"{value[:max_len]}..."'


def preview_number(value: float | np.integer | np.floating) -> str:
    """Preview a numeric value."""
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, (int, np.integer)):
        return str(value)
    # Float - format nicely (nan/inf have no integer form)
    if np.isfinite(value) and value == int(value):
        return str(int(value))
    return f"{value:.6g}"


def preview_dict(value: dict) -> str:
    """Preview a dict value."""
    n_keys = len(value)
    if n_keys == 0:
        return "{}"
    if n_keys <= DICT_PREVIEW_KEYS:
        keys_preview = ", ".join(str(k) for k in list(value.keys())[:DICT_PREVIEW_KEYS])
        return f"{{{keys_preview}}}"
    keys_preview = ", ".join(
        str(k) for k in list(value.keys())[:DICT_PREVIEW_KEYS_LARGE]
    )
    return f"{{{keys_preview}, ...}} ({n_keys} keys)"


def preview_sequence(value: list | tuple) -> str:
    """Preview a list or tuple value."""
    n_items = len(value)
    bracket = "[]" if isinstance(value, list) else "()"
    if n_items == 0:
        return bracket
    if n_items <= LIST_PREVIEW_ITEMS:
        try:
            items = [preview_item(v) for v in value[:LIST_PREVIEW_ITEMS]]
            if all(items):
                return f"{bracket[0]}{', '.join(items)}{bracket[1]}"
        except Exception:  # noqa: BLE001
            # Intentional broad catch: preview generation is best-effort
            pass
    return f"({n_items} items)"


def preview_item(value: object) -> str:
    """Generate a short preview for a single item (for list/tuple previews)."""
    if isinstance(value, str):
        if len(value) <= STRING_INLINE_LIMIT:
            return f'"{value}"'
        truncate_at = STRING_INLINE_LIMIT - 3  # Leave room for "..."
        return f'"{value[:truncate_at]}..."'
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        return str(value)
    if value is None:
        return "None"
    return ""  # Empty string means skip


def validate_key(key: str) -> tuple[bool, str, bool]:
    """Check if a key name is valid for HDF5/Zarr serialization.

    Key names (column names, uns keys, etc.) are validated because certain
    characters cause issues with the underlying storage formats (HDF5 and Zarr).

    Parameters
    ----------
    key
        Key name to validate

    Returns
    -------
    tuple of (is_valid, reason, is_hard_error)
        is_valid: False if there's an issue
        reason: Description of the issue
        is_hard_error: True means write fails NOW, False means deprecation warning
    """
    if not isinstance(key, str):
        return False, f"Non-string key ({type(key).__name__})", True
    if "/" in key:
        # Whether a slash aborts an h5ad write is user-configurable; mirror the
        # setting so the repr's severity matches what a write would actually do.
        disallowed = settings.disallow_forward_slash_in_h5ad
        reason = "Contains '/'" if disallowed else "Contains '/' (deprecated)"
        return False, reason, disallowed
    return True, "", False
