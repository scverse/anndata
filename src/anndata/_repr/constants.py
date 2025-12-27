"""
Constants for HTML representation.

This module contains default values for repr_html settings.
It has no dependencies to avoid circular imports, allowing both
_settings.py and _repr modules to import from here.
"""

from __future__ import annotations

# Display behavior
DEFAULT_FOLD_THRESHOLD = 5  # Auto-fold sections with more than N entries
DEFAULT_MAX_DEPTH = 3  # Maximum recursion depth for nested objects
DEFAULT_MAX_ITEMS = 200  # Maximum items to show per section
DEFAULT_MAX_STRING_LENGTH = 100  # Truncate strings longer than this
DEFAULT_PREVIEW_ITEMS = 5  # Number of items to show in previews (first/last)
DEFAULT_MAX_CATEGORIES = 100  # Max category values to display inline
DEFAULT_UNIQUE_LIMIT = 1_000_000  # Max rows to compute unique counts (0 to disable)

# Column widths (pixels)
DEFAULT_MAX_FIELD_WIDTH = 400  # Max width for field name column
DEFAULT_TYPE_WIDTH = 220  # Width for type column

# Inline styles for graceful degradation (hidden until JS enables)
STYLE_HIDDEN = "display:none;"
STYLE_SECTION_CONTENT = "padding:0;overflow:hidden;"
STYLE_SECTION_TABLE = "width:100%;border-collapse:collapse;"

# Warning messages
NOT_SERIALIZABLE_MSG = "Not serializable to H5AD/Zarr"
