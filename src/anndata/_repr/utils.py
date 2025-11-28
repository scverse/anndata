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
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from anndata import AnnData


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

    # Handle None
    if obj is None:
        return True, ""

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

    # Use the actual IO registry
    try:
        from anndata._io.specs.registry import _REGISTRY

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

    return False, f"Type '{type(obj).__module__}.{type(obj).__name__}' has no registered writer"


def should_warn_string_column(series: pd.Series) -> tuple[bool, str]:
    """
    Check if a string column will be auto-converted to categorical on save.

    This replicates the logic from AnnData.strings_to_categoricals():
    - Column must be string type (infer_dtype == "string")
    - Number of unique values must be less than total values

    Parameters
    ----------
    series
        Pandas Series to check

    Returns
    -------
    tuple of (should_warn, warning_message)
    """
    from pandas.api.types import infer_dtype

    dtype_str = infer_dtype(series)
    if dtype_str != "string":
        return False, ""

    try:
        n_unique = series.nunique()
        n_total = len(series)
    except Exception:
        return False, ""

    if n_unique < n_total:
        return (
            True,
            f"String column ({n_unique} unique values). "
            f"Will be converted to categorical on save.",
        )

    return False, ""


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
    if not key.endswith("_colors"):
        return False

    if not isinstance(value, (list, np.ndarray, tuple)):
        return False

    # Empty list is valid
    if len(value) == 0:
        return True

    # Check first element
    first = value[0] if len(value) > 0 else None
    if first is None:
        return False

    if isinstance(first, str):
        # Hex color
        if first.startswith("#"):
            return True
        # Named color (basic check)
        if first.lower() in _NAMED_COLORS:
            return True
        # RGB/RGBA string like "rgb(255, 0, 0)"
        if first.lower().startswith(("rgb(", "rgba(")):
            return True

    return False


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
    color_key = f"{column_name}_colors"
    if color_key not in adata.uns:
        return None

    colors = adata.uns[color_key]

    # Find the column in obs or var
    col = None
    if column_name in adata.obs.columns:
        col = adata.obs[column_name]
    elif column_name in adata.var.columns:
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
        AnnData object
    column_name
        Name of the column to check

    Returns
    -------
    Warning message if mismatch, None otherwise
    """
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


def format_memory_size(size_bytes: int | float) -> str:
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


def format_number(n: int | float) -> str:
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
        from importlib.metadata import version

        return version("anndata")
    except Exception:
        return "unknown"


def is_view(obj: Any) -> bool:
    """Check if an object is a view (for AnnData-like objects)."""
    return getattr(obj, "is_view", False)


def is_backed(obj: Any) -> bool:
    """Check if an object is backed (for AnnData-like objects)."""
    return getattr(obj, "isbacked", False)


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
_NAMED_COLORS = frozenset(
    {
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
    }
)
