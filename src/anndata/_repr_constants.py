"""
Constants for HTML representation.

This module contains default values for repr_html settings.
It is located outside the _repr/ package to avoid loading the full
_repr module when _settings.py imports these constants at anndata
import time. Python loads parent packages before submodules, so
importing from _repr.constants would trigger _repr/__init__.py.
"""

from __future__ import annotations

# Display behavior
DEFAULT_FOLD_THRESHOLD = 5  # Auto-fold sections with more than N entries
DEFAULT_MAX_DEPTH = 3  # Maximum recursion depth for nested objects
DEFAULT_MAX_ITEMS = 200  # Maximum items to show per section
DEFAULT_MAX_STRING_LENGTH = 100  # Truncate strings longer than this
DEFAULT_PREVIEW_ITEMS = 5  # Number of items to show in previews (first/last)
# Max category values to display inline (used by render_category_list in components.py)
# Note: DataFrame columns in obsm/varm have no limit - see DF_COLS_PREVIEW_LIMIT comment
DEFAULT_MAX_CATEGORIES = 100
DEFAULT_MAX_LAZY_CATEGORIES = (
    100  # Max categories to load for lazy categoricals (0 to skip)
)
DEFAULT_UNIQUE_LIMIT = 1_000_000  # Max rows to compute unique counts (0 to disable)

# Column widths (pixels)
DEFAULT_MAX_FIELD_WIDTH = 400  # Max width for field name column
DEFAULT_TYPE_WIDTH = 220  # Width for type column

# Field name column width calculation constants
# These values are empirically tuned for the default 13px monospace font
CHAR_WIDTH_PX = 8  # Average character width for monospace at 13px font-size
COPY_BUTTON_PADDING_PX = 30  # Extra space for copy button and cell padding
MIN_FIELD_WIDTH_PX = 80  # Minimum column width to avoid cramped display
DEFAULT_FIELD_WIDTH_PX = 100  # Default when no field names exist

# Inline styles for graceful degradation (hidden until JS enables)
STYLE_HIDDEN = "display:none;"
STYLE_SECTION_CONTENT = "padding:0;overflow:hidden;"
STYLE_SECTION_TABLE = "width:100%;border-collapse:collapse;"
STYLE_CAT_DOT = "width:8px;height:8px;border-radius:50%;display:inline-block;"

# Warning messages
NOT_SERIALIZABLE_MSG = "Not serializable to H5AD/Zarr"

# Preview truncation limits
TOOLTIP_TRUNCATE_LENGTH = 500  # Max chars for tooltip full text
ERROR_TRUNCATE_LENGTH = 200  # Max chars for error messages
COLOR_PREVIEW_LIMIT = 15  # Max color swatches to show
DICT_PREVIEW_KEYS = 3  # Keys to show in dict preview (small dicts)
DICT_PREVIEW_KEYS_LARGE = 2  # Keys to show in dict preview (large dicts)
LIST_PREVIEW_ITEMS = 3  # Items to show in list preview
STRING_INLINE_LIMIT = 20  # Max string length before truncating inline

# DataFrame column preview limits (for DataFrames in uns only)
# These control the compact inline preview shown in the type cell, e.g. "[col1, col2, ...]"
# Used by DataFrameFormatter in formatters.py for uns entries.
#
# Note: DataFrames in obsm/varm render ALL columns with CSS truncation + wrap button
# (render_entry_preview_cell() in components.py). No constant limits this - all columns
# are in the HTML, CSS truncates the display, and the wrap button expands to show all.
# This differs from categories which are limited by DEFAULT_MAX_CATEGORIES below.
DF_COLS_PREVIEW_LIMIT = 5  # Max columns to show in compact preview
DF_COLS_PREVIEW_MAX_LEN = 40  # Max total chars for column list string

# Table structure
ENTRY_TABLE_COLSPAN = 3  # Number of columns in entry table (name, type, preview)

# CSS class names for entry rows
CSS_ENTRY = "adata-entry"
CSS_ENTRY_NAME = "adata-entry-name"
CSS_ENTRY_TYPE = "adata-entry-type"
CSS_ENTRY_PREVIEW = "adata-entry-preview"
CSS_TEXT_MUTED = "adata-text-muted"
CSS_NESTED_ROW = "adata-nested-row"
CSS_NESTED_CONTENT = "adata-nested-content"
CSS_NESTED_ANNDATA = "adata-nested-anndata"
