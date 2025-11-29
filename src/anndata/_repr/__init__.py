"""
Rich HTML representation for AnnData objects in Jupyter notebooks.

This module provides an extensible HTML representation system with:
- Foldable sections with auto-collapse
- Search/filter functionality
- Color visualization for categorical data
- Value previews for simple types in uns (strings, numbers, dicts, lists)
- Serialization warnings
- Support for nested AnnData objects
- Graceful handling of unknown types

Extensibility
-------------
The system is designed to be extensible via registry patterns:

**Type Formatters** (for Python object types):
    New data types (e.g., TreeData, MuData, SpatialData) can register
    custom formatters without modifying core code.

**Uns Renderers** (for serialized data in uns):
    Packages can register custom HTML renderers for data stored in uns
    that contains a type hint. This is useful for packages that store
    complex data as JSON strings for H5AD/Zarr compatibility.

    Note: Since arbitrary Python objects cannot be serialized to H5AD/Zarr,
    packages should store their data as serializable types (strings, dicts,
    lists, arrays) with a type hint that their renderer can interpret.

    Security: Data in uns NEVER triggers code execution. Packages must
    be explicitly imported by the user, and only then can their registered
    renderers process data with matching type hints. Unrecognized hints
    fall back to safe JSON/string preview.

    Example for package authors::

        # In mypackage/__init__.py
        try:
            from anndata._repr import register_uns_renderer, UnsRendererOutput

            def render_my_config(value, context):
                # Parse and render the stored data
                return UnsRendererOutput(
                    html='<span>Custom preview</span>',
                    type_label="my config",
                )

            register_uns_renderer("mypackage.config", render_my_config)
        except ImportError:
            pass  # anndata not available

    Example data structure in uns::

        adata.uns["my_config"] = {
            "__anndata_repr__": "mypackage.config",
            "data": '{"setting": "value"}'
        }

    If mypackage is imported, its renderer will be used. Otherwise,
    a fallback preview shows: "[mypackage.config] (import mypackage to enable)".
"""

from __future__ import annotations

# Constants - these can be overridden via settings
DEFAULT_FOLD_THRESHOLD = 5  # Auto-fold sections with more than N entries
DEFAULT_MAX_DEPTH = 3  # Maximum recursion depth for nested objects
DEFAULT_MAX_ITEMS = 200  # Maximum items to show per section
DEFAULT_MAX_STRING_LENGTH = 100  # Truncate strings longer than this
DEFAULT_PREVIEW_ITEMS = 5  # Number of items to show in previews (first/last)
DEFAULT_MAX_CATEGORIES = 20  # Max category values to display inline
DEFAULT_UNIQUE_LIMIT = 1_000_000  # Max rows to compute unique counts (0 to disable)

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
    # Uns renderer registry (for custom serialized data visualization)
    UnsRendererOutput,
    UnsRendererRegistry,
    uns_renderer_registry,
    register_uns_renderer,
    extract_uns_type_hint,
    UNS_TYPE_HINT_KEY,
)

__all__ = [
    # Constants
    "DEFAULT_FOLD_THRESHOLD",
    "DEFAULT_MAX_DEPTH",
    "DEFAULT_MAX_ITEMS",
    "DEFAULT_MAX_STRING_LENGTH",
    "DEFAULT_PREVIEW_ITEMS",
    "DEFAULT_MAX_CATEGORIES",
    "DEFAULT_UNIQUE_LIMIT",
    "DOCS_BASE_URL",
    "SECTION_ORDER",
    # Main function
    "generate_repr_html",
    # Registry for extensibility (type formatters)
    "FormatterRegistry",
    "formatter_registry",
    "register_formatter",
    "SectionFormatter",
    "TypeFormatter",
    # Uns renderer registry (for custom serialized data visualization)
    "UnsRendererOutput",
    "UnsRendererRegistry",
    "uns_renderer_registry",
    "register_uns_renderer",
    "extract_uns_type_hint",
    "UNS_TYPE_HINT_KEY",
]
