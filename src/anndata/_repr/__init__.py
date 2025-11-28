"""
Rich HTML representation for AnnData objects in Jupyter notebooks.

This module provides an extensible HTML representation system with:
- Foldable sections with auto-collapse
- Search/filter functionality
- Color visualization for categorical data
- Serialization warnings
- Support for nested AnnData objects
- Graceful handling of unknown types

The system is designed to be extensible via a registry pattern,
allowing new data types (e.g., TreeData, MuData, SpatialData) to
register custom formatters without modifying core code.
"""

from __future__ import annotations

# Constants - these can be overridden via settings
DEFAULT_FOLD_THRESHOLD = 5  # Auto-fold sections with more than N entries
DEFAULT_MAX_DEPTH = 3  # Maximum recursion depth for nested objects
DEFAULT_MAX_ITEMS = 200  # Maximum items to show per section
DEFAULT_MAX_STRING_LENGTH = 100  # Truncate strings longer than this
DEFAULT_PREVIEW_ITEMS = 5  # Number of items to show in previews (first/last)
DEFAULT_MAX_CATEGORIES = 5  # Max category values to display inline

# Documentation base URL
DOCS_BASE_URL = "https://anndata.readthedocs.io/en/latest/"

# Section order for display
SECTION_ORDER = ("X", "obs", "var", "uns", "obsm", "varm", "layers", "obsp", "varp", "raw")

# Import main functionality
from anndata._repr.html import generate_repr_html
from anndata._repr.registry import (
    FormatterRegistry,
    formatter_registry,
    register_formatter,
    SectionFormatter,
    TypeFormatter,
)

__all__ = [
    # Constants
    "DEFAULT_FOLD_THRESHOLD",
    "DEFAULT_MAX_DEPTH",
    "DEFAULT_MAX_ITEMS",
    "DEFAULT_MAX_STRING_LENGTH",
    "DEFAULT_PREVIEW_ITEMS",
    "DEFAULT_MAX_CATEGORIES",
    "DOCS_BASE_URL",
    "SECTION_ORDER",
    # Main function
    "generate_repr_html",
    # Registry for extensibility
    "FormatterRegistry",
    "formatter_registry",
    "register_formatter",
    "SectionFormatter",
    "TypeFormatter",
]
