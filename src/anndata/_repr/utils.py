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

if TYPE_CHECKING:
    from typing import Any

    import pandas as pd

    from anndata import AnnData


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


def get_matching_column_colors(
    adata: AnnData,
    column_name: str,
) -> list[str] | None:
    """
    Get colors for a categorical column if they exist and match.

    Parameters
    ----------
    adata
        AnnData object
    column_name
        Name of the column to get colors for

    Returns
    -------
    List of color strings if colors exist and match, None otherwise
    """
    # Handle objects without .uns (e.g., Raw)
    if not hasattr(adata, "uns"):
        return None

    color_key = f"{column_name}_colors"
    if color_key not in adata.uns:
        return None

    colors = adata.uns[color_key]

    # Find the column in obs or var
    col = None
    if hasattr(adata, "obs") and column_name in adata.obs.columns:
        col = adata.obs[column_name]
    elif hasattr(adata, "var") and column_name in adata.var.columns:
        col = adata.var[column_name]

    if col is None:
        return None

    # Must be categorical
    if not hasattr(col, "cat"):
        return None

    n_categories = len(col.cat.categories)
    if len(colors) != n_categories:
        return None  # Mismatch

    return list(colors)


def check_color_category_mismatch(
    adata: AnnData,
    column_name: str,
) -> str | None:
    """
    Check if colors exist but don't match category count.

    Parameters
    ----------
    adata
        AnnData object (or object with .uns attribute)
    column_name
        Name of the column to check

    Returns
    -------
    Warning message if mismatch, None otherwise
    """
    # Handle objects without .uns (e.g., Raw)
    if not hasattr(adata, "uns"):
        return None

    color_key = f"{column_name}_colors"
    if color_key not in adata.uns:
        return None

    colors = adata.uns[color_key]

    for df in (adata.obs, adata.var):
        if column_name in df.columns and hasattr(df[column_name], "cat"):
            n_cats = len(df[column_name].cat.categories)
            if len(colors) != n_cats:
                return f"Color mismatch: {len(colors)} colors for {n_cats} categories"

    return None


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


def format_number(n: float) -> str:
    """Format a number with thousand separators."""
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
    return getattr(obj, "is_view", False)


def is_backed(obj: Any) -> bool:
    """Check if an object is backed (for AnnData-like objects)."""
    return getattr(obj, "isbacked", False)


def is_lazy_series(series: Any) -> bool:
    """
    Check if a Series-like object is lazy (backed by remote/lazy storage).

    This detects Series from Dataset2D (xarray-backed DataFrames used in
    lazy AnnData) to prevent operations that would trigger data loading.
    """
    # Check if it's from a Dataset2D (lazy DataFrame)
    # Dataset2D columns have a .data attribute that's an xarray DataArray
    if hasattr(series, "data") and hasattr(series.data, "dims"):
        return True
    # Check for xarray Variable backing
    return hasattr(series, "_variable")


def get_backing_info(obj: Any) -> dict[str, Any]:
    """Get information about backing for an AnnData-like object."""
    if not is_backed(obj):
        return {"backed": False}

    info = {
        "backed": True,
        "filename": str(getattr(obj, "filename", None)),
    }

    # Try to get file status
    file_obj = getattr(obj, "file", None)
    if file_obj is not None:
        info["is_open"] = getattr(file_obj, "is_open", None)

    # Detect format from filename
    filename = info["filename"]
    if filename:
        if filename.endswith(".h5ad"):
            info["format"] = "H5AD"
        elif ".zarr" in filename:
            info["format"] = "Zarr"
        else:
            info["format"] = "Unknown"

    return info


# Basic named colors for color detection
_NAMED_COLORS = frozenset({
    "red",
    "green",
    "blue",
    "yellow",
    "cyan",
    "magenta",
    "black",
    "white",
    "gray",
    "grey",
    "orange",
    "pink",
    "purple",
    "brown",
    "navy",
    "teal",
    "olive",
    "maroon",
    "lime",
    "aqua",
    "silver",
    "fuchsia",
    # CSS color names (partial list)
    "aliceblue",
    "antiquewhite",
    "aquamarine",
    "azure",
    "beige",
    "bisque",
    "blanchedalmond",
    "blueviolet",
    "burlywood",
    "cadetblue",
    "chartreuse",
    "chocolate",
    "coral",
    "cornflowerblue",
    "cornsilk",
    "crimson",
    "darkblue",
    "darkcyan",
    "darkgoldenrod",
    "darkgray",
    "darkgreen",
    "darkkhaki",
    "darkmagenta",
    "darkolivegreen",
    "darkorange",
    "darkorchid",
    "darkred",
    "darksalmon",
    "darkseagreen",
    "darkslateblue",
    "darkslategray",
    "darkturquoise",
    "darkviolet",
    "deeppink",
    "deepskyblue",
    "dimgray",
    "dodgerblue",
    "firebrick",
    "floralwhite",
    "forestgreen",
    "gainsboro",
    "ghostwhite",
    "gold",
    "goldenrod",
    "greenyellow",
    "honeydew",
    "hotpink",
    "indianred",
    "indigo",
    "ivory",
    "khaki",
    "lavender",
    "lavenderblush",
    "lawngreen",
    "lemonchiffon",
    "lightblue",
    "lightcoral",
    "lightcyan",
    "lightgoldenrodyellow",
    "lightgray",
    "lightgreen",
    "lightpink",
    "lightsalmon",
    "lightseagreen",
    "lightskyblue",
    "lightslategray",
    "lightsteelblue",
    "lightyellow",
    "limegreen",
    "linen",
    "mediumaquamarine",
    "mediumblue",
    "mediumorchid",
    "mediumpurple",
    "mediumseagreen",
    "mediumslateblue",
    "mediumspringgreen",
    "mediumturquoise",
    "mediumvioletred",
    "midnightblue",
    "mintcream",
    "mistyrose",
    "moccasin",
    "navajowhite",
    "oldlace",
    "olivedrab",
    "orangered",
    "orchid",
    "palegoldenrod",
    "palegreen",
    "paleturquoise",
    "palevioletred",
    "papayawhip",
    "peachpuff",
    "peru",
    "plum",
    "powderblue",
    "rosybrown",
    "royalblue",
    "saddlebrown",
    "salmon",
    "sandybrown",
    "seagreen",
    "seashell",
    "sienna",
    "skyblue",
    "slateblue",
    "slategray",
    "snow",
    "springgreen",
    "steelblue",
    "tan",
    "thistle",
    "tomato",
    "turquoise",
    "violet",
    "wheat",
    "whitesmoke",
    "yellowgreen",
})


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
    if n_keys <= 3:
        keys_preview = ", ".join(str(k) for k in list(value.keys())[:3])
        return f"{{{keys_preview}}}"
    keys_preview = ", ".join(str(k) for k in list(value.keys())[:2])
    return f"{{{keys_preview}, ...}} ({n_keys} keys)"


def preview_sequence(value: list | tuple) -> str:
    """Preview a list or tuple value."""
    n_items = len(value)
    bracket = "[]" if isinstance(value, list) else "()"
    if n_items == 0:
        return bracket
    if n_items <= 3:
        try:
            items = [preview_item(v) for v in value[:3]]
            if all(items):
                return f"{bracket[0]}{', '.join(items)}{bracket[1]}"
        except Exception:  # noqa: BLE001
            # Intentional broad catch: preview generation is best-effort
            pass
    return f"({n_items} items)"


def preview_item(value: Any) -> str:
    """Generate a short preview for a single item (for list/tuple previews)."""
    if isinstance(value, str):
        if len(value) <= 20:
            return f'"{value}"'
        return f'"{value[:17]}..."'
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
