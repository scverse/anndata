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
from typing import TYPE_CHECKING

import numpy as np

from .._repr_constants import (
    DICT_PREVIEW_KEYS,
    DICT_PREVIEW_KEYS_LARGE,
    LIST_PREVIEW_ITEMS,
    STRING_INLINE_LIMIT,
)

if TYPE_CHECKING:
    from typing import Any

    import pandas as pd

    from anndata import AnnData

    from .registry import FormatterContext


def _check_serializable_single(obj: Any) -> tuple[bool, str]:
    """Check if a single (non-container) object is serializable."""
    # Handle None
    if obj is None:
        return True, ""

    # Use the actual IO registry
    try:
        from .._io.specs.registry import _REGISTRY

        _REGISTRY.get_spec(obj)
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


def is_serializable(
    obj: Any,
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
        for i, v in enumerate(obj):
            ok, reason = is_serializable(v, _depth=_depth + 1, _max_depth=_max_depth)
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
            f"String column ({n_unique} unique). "
            f"Will be converted to categorical on save.",
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


def is_color_list(key: str, value: Any) -> bool:
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


def _get_categories_from_column(col: Any) -> list:
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
    col: Any,
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


def _compute_if_dask(obj: Any) -> Any:
    """
    Compute a dask array/object if it is one, otherwise return as-is.

    For lazy AnnData, uns values may be dask arrays that need to be
    computed to get the actual values.
    """
    if hasattr(obj, "compute"):
        return obj.compute()
    return obj


def get_matching_column_colors(
    adata: AnnData,
    column_name: str,
    *,
    limit: int | None = None,
) -> list[str] | None:
    """
    Get colors for a column from uns if they exist.

    This function is called by CategoricalFormatter which already verified
    the column is categorical. It just looks up and returns the colors.
    Color count validation is done separately by check_color_category_mismatch.

    Parameters
    ----------
    adata
        AnnData object
    column_name
        Name of the column to get colors for
    limit
        If provided, only load the first `limit` colors. This avoids loading
        all colors from disk when only displaying partial categories.

    Returns
    -------
    List of color strings if colors exist, None otherwise
    """
    colors = _get_colors_from_uns(adata, column_name, limit=limit)
    return list(colors) if colors is not None else None


def check_color_category_mismatch(
    adata: AnnData,
    column_name: str,
    n_categories: int,
) -> str | None:
    """
    Check if colors exist but don't match category count.

    Called by _render_dataframe_entry for categorical columns. The caller
    already knows this is categorical and has the category count.

    Parameters
    ----------
    adata
        AnnData object (or object with .uns attribute)
    column_name
        Name of the column to check
    n_categories
        Number of categories in the column

    Returns
    -------
    Warning message if mismatch, None otherwise
    """
    colors = _get_colors_from_uns(adata, column_name)
    if colors is None:
        return None

    if len(colors) != n_categories:
        return f"Color mismatch: {len(colors)} colors for {n_categories} categories"

    return None


def _get_colors_from_uns(
    adata: AnnData,
    column_name: str,
    limit: int | None = None,
) -> Any | None:
    """Get colors from uns for a column, handling lazy loading.

    Parameters
    ----------
    adata
        AnnData object (or object with .uns attribute)
    column_name
        Name of the column (colors key will be "{column_name}_colors")
    limit
        If provided, only load the first `limit` colors (for dask arrays)

    Returns
    -------
    Colors array/list if found, None otherwise
    """
    # Handle objects without .uns (e.g., Raw)
    if not hasattr(adata, "uns"):
        return None

    color_key = f"{column_name}_colors"
    if color_key not in adata.uns:
        return None

    colors = adata.uns[color_key]

    # For lazy AnnData with dask arrays, slice before computing
    if limit is not None and hasattr(colors, "compute"):
        return colors[:limit].compute()

    # Compute if dask array (for lazy AnnData)
    return _compute_if_dask(colors)


def escape_html(text: str) -> str:
    """Escape HTML special characters."""
    return html.escape(str(text))


def sanitize_for_id(text: str) -> str:
    """Sanitize a string for use as an HTML id attribute."""
    # Replace non-alphanumeric chars with underscore
    sanitized = re.sub(r"[^a-zA-Z0-9_-]", "_", str(text))
    # Ensure it starts with a letter
    if sanitized and not sanitized[0].isalpha():
        sanitized = "id_" + sanitized
    return sanitized


def truncate_string(text: str, max_length: int = 100) -> str:
    """Truncate a string and add ellipsis if needed."""
    text = str(text)
    if len(text) <= max_length:
        return text
    return text[: max_length - 3] + "..."


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


def format_number(n: float | str) -> str:
    """Format a number with thousand separators.

    Accepts int, float, or str (for fallback values like "?").
    """
    if isinstance(n, str):
        return n
    if isinstance(n, float):
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


def is_view(obj: Any) -> bool:
    """Check if an object is a view (for AnnData-like objects)."""
    try:
        return getattr(obj, "is_view", False)
    except Exception:  # noqa: BLE001
        return False


def is_backed(obj: Any) -> bool:
    """Check if an object is backed (for AnnData-like objects)."""
    try:
        return getattr(obj, "isbacked", False)
    except Exception:  # noqa: BLE001
        return False


def get_backing_info(obj: Any) -> dict[str, Any]:
    """Get information about backing for an AnnData-like object."""
    try:
        if not is_backed(obj):
            return {"backed": False}

        filename = str(getattr(obj, "filename", None) or "")
        info: dict[str, Any] = {
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
    from functools import lru_cache
    from importlib.resources import files

    @lru_cache(maxsize=1)
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
    # Float - format nicely
    if value == int(value):
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


def preview_item(value: Any) -> str:
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


def generate_value_preview(value: Any, max_len: int = 100) -> str:
    """Generate a human-readable preview of a value.

    Returns empty string if no meaningful preview can be generated.
    """
    if value is None:
        return "None"
    if isinstance(value, str):
        return preview_string(value, max_len)
    if isinstance(value, (bool, int, float, np.integer, np.floating)):
        return preview_number(value)
    if isinstance(value, dict):
        return preview_dict(value)
    if isinstance(value, (list, tuple)):
        return preview_sequence(value)
    # No preview for complex types
    return ""


def get_setting(name: str, *, default: Any) -> Any:
    """Get a setting value from anndata.settings, falling back to default.

    Parameters
    ----------
    name
        The setting name (e.g., "repr_html_max_items")
    default
        Default value if setting is not available

    Returns
    -------
    The setting value or default
    """
    try:
        from anndata import settings

        return getattr(settings, name, default)
    except (ImportError, AttributeError):
        return default


def check_column_name(name: str) -> tuple[bool, str, bool]:
    """Check if a column name is valid for HDF5/Zarr serialization.

    Column/key names are validated because certain characters cause issues
    with the underlying storage formats (HDF5 and Zarr).

    Parameters
    ----------
    name
        Column or key name to validate

    Returns
    -------
    tuple of (is_valid, reason, is_hard_error)
        is_valid: False if there's an issue
        reason: Description of the issue
        is_hard_error: True means write fails NOW, False means deprecation warning
    """
    if not isinstance(name, str):
        return False, f"Non-string name ({type(name).__name__})", True
    # Slashes will be disallowed in h5 stores (FutureWarning)
    if "/" in name:
        return False, "Contains '/' (deprecated)", False
    return True, "", False
