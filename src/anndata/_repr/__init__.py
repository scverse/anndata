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
The system is designed to be extensible via two registry patterns:

**TypeFormatter** (for custom visualization of values):
    Register a formatter to customize how specific types are displayed.
    Can match by Python type OR by embedded type hints in data.

    Attributes:
        - ``priority``: Higher priority formatters are checked first (default: 0)
        - ``sections``: Tuple of section names to restrict formatter to (default: None = all)

    Example - format by Python type::

        from anndata._repr import register_formatter, TypeFormatter, FormattedOutput


        @register_formatter
        class MyArrayFormatter(TypeFormatter):
            sections = ("obsm", "varm")  # Only apply to obsm/varm

            def can_format(self, obj):
                return isinstance(obj, MyArrayType)

            def format(self, obj, context):
                return FormattedOutput(
                    type_name=f"MyArray {obj.shape}",
                    css_class="dtype-myarray",
                )

    Example - format by embedded type hint (for tagged data in uns)::

        from anndata._repr import register_formatter, TypeFormatter, FormattedOutput
        from anndata._repr import extract_uns_type_hint


        @register_formatter
        class MyConfigFormatter(TypeFormatter):
            priority = 100  # Check before fallback
            sections = ("uns",)  # Only apply to uns

            def can_format(self, obj):
                hint, _ = extract_uns_type_hint(obj)
                return hint == "mypackage.config"

            def format(self, obj, context):
                hint, data = extract_uns_type_hint(obj)
                return FormattedOutput(
                    type_name="config",
                    html_content="<span>Custom config preview</span>",
                )

    Data structure for type hints (works in any section)::

        adata.uns["my_config"] = {
            "__anndata_repr__": "mypackage.config",
            "data": '{"setting": "value"}',
        }

    When a package registers a formatter and the user imports that package,
    the formatter will automatically handle matching tagged data. Without
    the import, a fallback shows: "[mypackage.config] (import mypackage)".

    See :func:`extract_uns_type_hint` for full documentation on this pattern.

**SectionFormatter** (for adding new sections):
    Register a formatter to add entirely new sections (like TreeData's obst/vart).

    Example::

        from anndata._repr import register_formatter, SectionFormatter
        from anndata._repr import FormattedEntry, FormattedOutput


        @register_formatter
        class ObstSectionFormatter(SectionFormatter):
            section_name = "obst"
            after_section = "obsm"  # Position after obsm

            def should_show(self, obj):
                return hasattr(obj, "obst") and len(obj.obst) > 0

            def get_entries(self, obj, context):
                return [
                    FormattedEntry(
                        key=k,
                        output=FormattedOutput(type_name=f"Tree ({v.n_nodes} nodes)"),
                    )
                    for k, v in obj.obst.items()
                ]
"""

from __future__ import annotations

# Constants - these can be overridden via settings
DEFAULT_FOLD_THRESHOLD = 5  # Auto-fold sections with more than N entries
DEFAULT_MAX_DEPTH = 3  # Maximum recursion depth for nested objects
DEFAULT_MAX_ITEMS = 200  # Maximum items to show per section
DEFAULT_MAX_STRING_LENGTH = 100  # Truncate strings longer than this
DEFAULT_PREVIEW_ITEMS = 5  # Number of items to show in previews (first/last)
DEFAULT_MAX_CATEGORIES = 100  # Max category values to display inline
DEFAULT_UNIQUE_LIMIT = 1_000_000  # Max rows to compute unique counts (0 to disable)

# Documentation base URL
DOCS_BASE_URL = "https://anndata.readthedocs.io/en/latest/"

# Section order for display
SECTION_ORDER = (
    "X",
    "obs",
    "var",
    "uns",
    "obsm",
    "varm",
    "layers",
    "obsp",
    "varp",
    "raw",
)

# Import main functionality
from anndata._repr.html import generate_repr_html  # noqa: E402
from anndata._repr.registry import (  # noqa: E402
    UNS_TYPE_HINT_KEY,
    FormattedEntry,
    FormattedOutput,
    FormatterContext,
    # Type formatter registry
    FormatterRegistry,
    SectionFormatter,
    TypeFormatter,
    # Type hint extraction (for tagged data in uns)
    extract_uns_type_hint,
    formatter_registry,
    register_formatter,
)

__all__ = [  # noqa: RUF022  # organized by category, not alphabetically
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
    # Registry for extensibility
    "FormatterRegistry",
    "formatter_registry",
    "register_formatter",
    "SectionFormatter",
    "TypeFormatter",
    "FormattedOutput",
    "FormattedEntry",
    "FormatterContext",
    # Type hint extraction (for tagged data in uns)
    "extract_uns_type_hint",
    "UNS_TYPE_HINT_KEY",
]
