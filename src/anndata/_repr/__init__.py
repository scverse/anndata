"""
Rich HTML representation for AnnData objects in Jupyter notebooks.

Module Architecture
-------------------
This package uses a layered import hierarchy to avoid circular imports.
When modifying imports, maintain this order:

.. code-block:: text

    _repr_constants.py (outside _repr/, no internal imports)
    │   Constants only. Imported by _settings.py at anndata import time.
    │   Must not import anything from anndata.
    │
    └─► utils.py (depends only on external: numpy, pandas)
        │   HTML escaping, formatting, serialization checks.
        │
        └─► components.py (depends on: utils)
            │   UI building blocks: badges, buttons, icons.
            │
            └─► registry.py (depends on: _repr_constants)
                │   Formatter registry, TypeFormatter, SectionFormatter.
                │   NOTE: formatters.py imports registry for registration.
                │
                └─► core.py (depends on: components, registry, utils)
                    │   Shared rendering primitives: render_section().
                    │
                    ├─► sections.py (depends on: core, components, registry, utils)
                    │   │   Section-specific renderers (obs, var, uns, etc.).
                    │   │   Uses late import of html.render_formatted_entry.
                    │   │
                    │   └─► html.py (depends on: core, sections, components, registry, utils)
                    │           Main orchestrator: generate_repr_html().
                    │           Side-effect import of formatters.py for registration.
                    │
                    └─► formatters.py (depends on: registry, components, utils)
                            Built-in type formatters. Auto-registers on import.
                            Late import of html.generate_repr_html for AnnDataFormatter.

    __init__.py (imports from all modules for public API)

Key patterns for avoiding circular imports:
1. ``_repr_constants.py`` is outside ``_repr/`` - safe to import anywhere
2. Side-effect imports use ``from . import formatters as _formatters  # noqa: F401``
3. Late imports inside functions for cross-module dependencies
4. Type hints use ``if TYPE_CHECKING:`` blocks

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
                    css_class="anndata-dtype--myarray",
                    # preview_html provides HTML for the preview column (rightmost)
                    preview_html=f'<span class="anndata-text--muted">({obj.n_items} items)</span>',
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
                    preview_html="<span>Custom config preview</span>",
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

Building Custom _repr_html_
---------------------------
For packages with AnnData-adjacent objects (like SpatialData, MuData) that need
their own ``_repr_html_``, you can reuse anndata's CSS, JavaScript, and helpers.

**Basic structure**::

    from anndata._repr import get_css, get_javascript


    class MyData:
        def _repr_html_(self):
            container_id = f"mydata-{id(self)}"
            return f'''
                {get_css()}
                <div class="anndata-repr" id="{container_id}">
                    <div class="anndata-header">
                        <span class="anndata-header__type">MyData</span>
                        <span class="anndata-header__shape">100 items</span>
                    </div>
                    <div class="anndata-repr__sections">
                        <!-- sections go here -->
                    </div>
                </div>
                {get_javascript(container_id)}
            '''

**CSS classes** (stable, can be used directly):

    Classes follow `BEM naming convention <http://getbem.com/>`_:
    ``anndata-{block}__{element}--{modifier}``

    **Blocks** (top-level components):
        - ``anndata-repr``: Main container (required for JS and styling)
        - ``anndata-header``: Header row (flexbox, contains type/shape/badges)
        - ``anndata-footer``: Footer row (version, memory info)
        - ``anndata-section``: Individual section wrapper
        - ``anndata-entry``: Table row for data entries
        - ``anndata-badge``: Status badges (View, Backed, Lazy)
        - ``anndata-dtype``: Data type indicators

    **Elements** (parts of blocks, use ``__``):
        - ``anndata-header__type``: Type name span in header
        - ``anndata-header__shape``: Shape/dimensions span in header
        - ``anndata-section__header``: Section header (clickable to fold)
        - ``anndata-section__content``: Section content (rows)
        - ``anndata-entry__name``: Entry name cell
        - ``anndata-entry__type``: Entry type cell
        - ``anndata-entry__preview``: Entry preview cell

    **Modifiers** (variants, use ``--``):
        - ``anndata-section--collapsed``: Collapsed section state
        - ``anndata-badge--view``: View badge variant
        - ``anndata-dtype--category``: Categorical dtype styling

**CSS variables** (set on ``.anndata-repr`` element):

    - ``--anndata-name-col-width``: Width of name column (default: 150px)
    - ``--anndata-type-col-width``: Width of type column (default: 220px)

**Using render helpers** for consistent section rendering::

    from anndata._repr import (
        get_css,
        get_javascript,
        render_section,
        render_formatted_entry,
        render_badge,
        render_search_box,
        FormattedEntry,
        FormattedOutput,
    )


    def _repr_html_(self):
        container_id = f"mydata-{id(self)}"
        parts = [get_css()]

        # Header
        parts.append(f'''
            <div class="anndata-repr" id="{container_id}">
            <div class="anndata-header">
                <span class="anndata-header__type">MyData</span>
                {render_badge("Zarr", "anndata-badge--backed")}
                <span style="flex-grow:1;"></span>
                {render_search_box(container_id)}
            </div>
            <div class="anndata-repr__sections">
        ''')

        # Build section entries
        entries = []
        for key, value in self.items.items():
            entry = FormattedEntry(
                key=key,
                output=FormattedOutput(
                    type_name=f"array {value.shape}",
                    css_class="anndata-dtype--ndarray",
                ),
            )
            entries.append(render_formatted_entry(entry))

        # Render section
        parts.append(
            render_section(
                "items",
                "\\n".join(entries),
                n_items=len(self.items),
            )
        )

        parts.append("</div></div>")
        parts.append(get_javascript(container_id))
        return "\\n".join(parts)

**Embedding nested AnnData** with full interactivity::

    from anndata._repr import generate_repr_html, FormattedEntry, FormattedOutput

    nested_html = generate_repr_html(adata, depth=1, max_depth=3)
    entry = FormattedEntry(
        key="table",
        output=FormattedOutput(
            type_name=f"AnnData ({adata.n_obs} x {adata.n_vars})",
            expanded_html=nested_html,  # Collapsible content below the row
        ),
    )

**Complete example**: See ``MockSpatialData`` in ``tests/visual_inspect_repr_html.py``
for a full implementation with images, labels, points, shapes, and nested tables.
"""

from __future__ import annotations

# Import constants from dedicated module (single source of truth)
# Note: _repr_constants is outside _repr/ to avoid loading the full _repr
# package when _settings.py imports constants at anndata import time.
from .._repr_constants import (
    DEFAULT_FOLD_THRESHOLD,
    DEFAULT_MAX_CATEGORIES,
    DEFAULT_MAX_DEPTH,
    DEFAULT_MAX_FIELD_WIDTH,
    DEFAULT_MAX_ITEMS,
    DEFAULT_MAX_LAZY_CATEGORIES,
    DEFAULT_MAX_STRING_LENGTH,
    DEFAULT_PREVIEW_ITEMS,
    DEFAULT_TYPE_WIDTH,
    DEFAULT_UNIQUE_LIMIT,
    NOT_SERIALIZABLE_MSG,
)

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
# Inline styles for graceful degradation (from single source of truth)
from .._repr_constants import STYLE_HIDDEN  # noqa: E402

# Building blocks for packages that want to create their own _repr_html_
# These allow reusing anndata's styling while building custom representations
from .components import (  # noqa: E402
    TypeCellConfig,
    render_badge,
    render_copy_button,
    render_fold_icon,
    render_header_badges,
    render_search_box,
    render_warning_icon,
)
from .css import get_css  # noqa: E402
from .html import (  # noqa: E402
    generate_repr_html,
    render_formatted_entry,
    render_section,
)
from .javascript import get_javascript  # noqa: E402
from .registry import (  # noqa: E402
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

# HTML rendering helpers for building custom sections
from .utils import (  # noqa: E402
    check_column_name,
    escape_html,
    format_memory_size,
    format_number,
)

__all__ = [  # noqa: RUF022  # organized by category, not alphabetically
    # Constants
    "DEFAULT_FOLD_THRESHOLD",
    "DEFAULT_MAX_DEPTH",
    "DEFAULT_MAX_ITEMS",
    "DEFAULT_MAX_STRING_LENGTH",
    "DEFAULT_PREVIEW_ITEMS",
    "DEFAULT_MAX_CATEGORIES",
    "DEFAULT_MAX_LAZY_CATEGORIES",
    "DEFAULT_UNIQUE_LIMIT",
    "DEFAULT_MAX_FIELD_WIDTH",
    "DEFAULT_TYPE_WIDTH",
    "DOCS_BASE_URL",
    "SECTION_ORDER",
    "NOT_SERIALIZABLE_MSG",
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
    # Building blocks for custom _repr_html_ implementations
    "get_css",
    "get_javascript",
    "escape_html",
    "format_number",
    "format_memory_size",
    "render_section",
    "render_formatted_entry",
    "STYLE_HIDDEN",
    # UI component helpers
    "render_search_box",
    "render_fold_icon",
    "render_copy_button",
    "render_badge",
    "render_header_badges",
    "render_warning_icon",
    "TypeCellConfig",
    # Validation helpers
    "check_column_name",
]
