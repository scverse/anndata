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
DEFAULT_MAX_README_SIZE = 100_000  # Max README size in chars (100KB, 0 to disable)

# Column widths (pixels)
DEFAULT_MAX_FIELD_WIDTH = 400  # Max width for field name column
DEFAULT_TYPE_WIDTH = 220  # Width for type column

# Field name column width calculation constants
# These values are empirically tuned for the default 13px monospace font
CHAR_WIDTH_PX = 8  # Average character width for monospace at 13px font-size
COPY_BUTTON_PADDING_PX = 30  # Extra space for copy button and cell padding
MIN_FIELD_WIDTH_PX = 80  # Minimum column width to avoid cramped display
DEFAULT_FIELD_WIDTH_PX = 100  # Default when no field names exist

# Inline style for graceful degradation (hidden until JS enables).
# JS sets different display values per element, so this must stay inline.
STYLE_HIDDEN = "display:none;"

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

# CSS class names for entry rows (BEM: anndata-entry block)
CSS_ENTRY = "anndata-entry"
CSS_ENTRY_NAME = "anndata-entry__name"
CSS_ENTRY_TYPE = "anndata-entry__type"
CSS_ENTRY_PREVIEW = "anndata-entry__preview"
CSS_TEXT_MUTED = "anndata-text--muted"
CSS_TEXT_ERROR = "anndata-text--error"
CSS_TEXT_WARNING = "anndata-text--warning"
CSS_NESTED_ROW = "anndata-entry--nested"
CSS_NESTED_CONTENT = "anndata-entry__nested-content"
CSS_NESTED_ANNDATA = "anndata-entry__nested-anndata"

# CSS class names for dtype spans (BEM: anndata-dtype block with modifiers)
# These provide visual differentiation for different data types
# Basic types
CSS_DTYPE_INT = "anndata-dtype--int"
CSS_DTYPE_FLOAT = "anndata-dtype--float"
CSS_DTYPE_BOOL = "anndata-dtype--bool"
CSS_DTYPE_STRING = "anndata-dtype--string"
CSS_DTYPE_OBJECT = "anndata-dtype--object"
# Container/structured types
CSS_DTYPE_CATEGORY = "anndata-dtype--category"
CSS_DTYPE_DATAFRAME = "anndata-dtype--dataframe"
CSS_DTYPE_ANNDATA = "anndata-dtype--anndata"
# Array types
CSS_DTYPE_NDARRAY = "anndata-dtype--ndarray"
CSS_DTYPE_SPARSE = "anndata-dtype--sparse"
# Specialized array types
CSS_DTYPE_DASK = "anndata-dtype--dask"
CSS_DTYPE_GPU = "anndata-dtype--gpu"
CSS_DTYPE_TPU = "anndata-dtype--tpu"
CSS_DTYPE_AWKWARD = "anndata-dtype--awkward"
CSS_DTYPE_ARRAY_API = "anndata-dtype--array-api"
# Extension/unknown
CSS_DTYPE_EXTENSION = "anndata-dtype--extension"
CSS_DTYPE_UNKNOWN = "anndata-dtype--unknown"

# CSS class names for badges (BEM: anndata-badge block with modifiers)
CSS_BADGE = "anndata-badge"
CSS_BADGE_VIEW = "anndata-badge--view"
CSS_BADGE_BACKED = "anndata-badge--backed"
CSS_BADGE_LAZY = "anndata-badge--lazy"
CSS_BADGE_EXTENSION = "anndata-badge--extension"

# CSS class names for color swatches (BEM: anndata-colors block)
CSS_COLORS = "anndata-colors"
CSS_COLORS_SWATCH = "anndata-colors__swatch"
CSS_COLORS_SWATCH_INVALID = "anndata-colors__swatch--invalid"

# Section names (canonical strings used for data-section attributes and keys)
SECTION_X = "X"
SECTION_OBS = "obs"
SECTION_VAR = "var"
SECTION_UNS = "uns"
SECTION_OBSM = "obsm"
SECTION_VARM = "varm"
SECTION_LAYERS = "layers"
SECTION_OBSP = "obsp"
SECTION_VARP = "varp"
SECTION_RAW = "raw"

# Internal AnnData attributes to skip when detecting unknown sections.
#
# Why explicit list instead of introspection?
# -------------------------------------------
# We intentionally maintain this list manually rather than using introspection
# (e.g., inspect.getmembers(AnnData, property)) because:
#
# 1. NEW data slots added to AnnData should appear in "unknown sections" until
#    a proper renderer is implemented. Introspection would hide them.
#
# 2. This list contains only internal/meta attributes that are NOT data slots:
#    - Shape/size metadata (shape, n_obs, n_vars)
#    - Index accessors (obs_names, var_names)
#    - File/backing info (filename, file, isbacked)
#    - View status (is_view, isview)
#    - Transpose accessor (T)
#
# What else is skipped (not in this list)?
# ----------------------------------------
# - Data sections (obs, var, uns, obsm, etc.) via SECTION_ORDER
# - Custom sections registered via SectionFormatter (e.g., obst, vart from TreeData)
# - Methods are filtered by `not callable()` at runtime
#
# When to update this list:
# - Add entries when AnnData gains new internal/meta properties (not data slots)
# - Do NOT add new data slots here - they should appear in "unknown" until rendered
INTERNAL_ANNDATA_ATTRS = frozenset({
    # Shape/size metadata
    "shape",
    "n_obs",
    "n_vars",
    # Index accessors
    "obs_names",
    "var_names",
    # File/backing info
    "filename",
    "file",
    "isbacked",
    # View status
    "is_view",
    "isview",  # Deprecated alias, triggers warning if accessed
    # Transpose accessor
    "T",
})
