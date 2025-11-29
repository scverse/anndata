"""
Main HTML generator for AnnData representation.

This module generates the complete HTML representation by:
1. Building the header with badges
2. Rendering the search box
3. Generating metadata (version, memory)
4. Rendering each section (X, obs, var, uns, etc.)
5. Handling nested objects recursively
"""

from __future__ import annotations

import uuid
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from anndata._repr import (
    DEFAULT_FOLD_THRESHOLD,
    DEFAULT_MAX_CATEGORIES,
    DEFAULT_MAX_DEPTH,
    DEFAULT_MAX_ITEMS,
    DEFAULT_MAX_STRING_LENGTH,
    DEFAULT_PREVIEW_ITEMS,
    DEFAULT_UNIQUE_LIMIT,
    DOCS_BASE_URL,
    SECTION_ORDER,
)
from anndata._repr.css import get_css
from anndata._repr.javascript import get_javascript
from anndata._repr.registry import (
    FormattedEntry,
    FormattedOutput,
    FormatterContext,
    extract_uns_type_hint,
    formatter_registry,
)
from anndata._repr.utils import (
    check_color_category_mismatch,
    escape_html,
    format_memory_size,
    format_number,
    get_anndata_version,
    get_backing_info,
    get_matching_column_colors,
    is_backed,
    is_color_list,
    is_view,
    should_warn_string_column,
)

if TYPE_CHECKING:
    from anndata import AnnData

# Import formatters to register them (side-effect import)
import anndata._repr.formatters  # noqa: F401


def generate_repr_html(
    adata: AnnData,
    *,
    depth: int = 0,
    max_depth: int | None = None,
    fold_threshold: int | None = None,
    max_items: int | None = None,
    show_header: bool = True,
    show_search: bool = True,
    _container_id: str | None = None,
) -> str:
    """
    Generate HTML representation for an AnnData object.

    Parameters
    ----------
    adata
        The AnnData object to represent
    depth
        Current recursion depth (for nested AnnData in .uns)
    max_depth
        Maximum recursion depth. Uses settings/default if None.
    fold_threshold
        Auto-fold sections with more entries than this. Uses settings/default if None.
    max_items
        Maximum items to show per section. Uses settings/default if None.
    show_header
        Whether to show the header (for nested display)
    show_search
        Whether to show the search box (only at top level)
    _container_id
        Internal: container ID for scoping

    Returns
    -------
    HTML string
    """
    # Get settings with defaults
    if max_depth is None:
        max_depth = _get_setting("repr_html_max_depth", DEFAULT_MAX_DEPTH)
    if fold_threshold is None:
        fold_threshold = _get_setting(
            "repr_html_fold_threshold", DEFAULT_FOLD_THRESHOLD
        )
    if max_items is None:
        max_items = _get_setting("repr_html_max_items", DEFAULT_MAX_ITEMS)

    # Check if HTML repr is enabled
    if not _get_setting("repr_html_enabled", True):
        # Fallback to text repr
        return f"<pre>{escape_html(repr(adata))}</pre>"

    # Check max depth
    if depth >= max_depth:
        return _render_max_depth_indicator(adata)

    # Generate unique container ID
    container_id = _container_id or f"anndata-repr-{uuid.uuid4().hex[:8]}"

    # Create formatter context
    context = FormatterContext(
        depth=depth,
        max_depth=max_depth,
        adata_ref=adata,
    )

    # Build HTML parts
    parts = []

    # CSS and JS only at top level
    if depth == 0:
        parts.append(get_css())

    # Container
    parts.append(f'<div class="anndata-repr" id="{container_id}" data-depth="{depth}">')

    # Header (with search box integrated on the right)
    if show_header:
        parts.append(
            _render_header(
                adata, show_search=show_search and depth == 0, container_id=container_id
            )
        )

    # Index preview (only at top level)
    if depth == 0:
        parts.append(_render_index_preview(adata))

    # Sections container
    parts.append('<div class="adata-sections">')

    # X as a simple entry (like layers)
    parts.append(_render_x_entry(adata, context))

    # Get custom sections grouped by their position
    custom_sections_after = _get_custom_sections_by_position(adata)

    # Standard sections with custom sections inserted at their positions
    for section in SECTION_ORDER:
        if section == "X":
            continue  # Already rendered
        if section == "raw":
            parts.append(_render_raw_section(adata, context, fold_threshold, max_items))
        elif section in ("obs", "var"):
            parts.append(
                _render_dataframe_section(
                    adata, section, context, fold_threshold, max_items
                )
            )
        elif section == "uns":
            parts.append(
                _render_uns_section(
                    adata, context, fold_threshold, max_items, max_depth
                )
            )
        else:
            parts.append(
                _render_mapping_section(
                    adata, section, context, fold_threshold, max_items
                )
            )

        # Render any custom sections that should appear after this section
        if section in custom_sections_after:
            for section_formatter in custom_sections_after[section]:
                parts.append(
                    _render_custom_section(
                        adata, section_formatter, context, fold_threshold, max_items
                    )
                )

    # Render any custom sections that don't have a specific position (appear at end)
    if None in custom_sections_after:
        for section_formatter in custom_sections_after[None]:
            parts.append(
                _render_custom_section(
                    adata, section_formatter, context, fold_threshold, max_items
                )
            )

    parts.append("</div>")  # adata-sections

    # Footer with metadata (only at top level)
    if depth == 0:
        parts.append(_render_footer(adata))

    parts.append("</div>")  # anndata-repr

    # JavaScript (only at top level)
    if depth == 0:
        parts.append(get_javascript(container_id))

    return "\n".join(parts)


# =============================================================================
# Custom Section Support
# =============================================================================


def _get_custom_sections_by_position(adata: Any) -> dict[str | None, list]:
    """
    Get registered custom section formatters grouped by their position.

    Returns a dict mapping after_section -> list of formatters.
    None key contains formatters that should appear at the end.
    """
    from collections import defaultdict

    result = defaultdict(list)

    for section_name in formatter_registry.get_registered_sections():
        formatter = formatter_registry.get_section_formatter(section_name)
        if formatter is None:
            continue

        # Skip standard sections (they're handled separately)
        if section_name in SECTION_ORDER:
            continue

        # Check if this section should be shown for this object
        try:
            if not formatter.should_show(adata):
                continue
        except Exception:  # noqa: BLE001
            # Intentional broad catch: custom formatters shouldn't break the repr
            continue

        # Group by position
        after = getattr(formatter, "after_section", None)
        result[after].append(formatter)

    return dict(result)


def _render_custom_section(
    adata: Any,
    formatter: Any,  # SectionFormatter
    context: FormatterContext,
    fold_threshold: int,
    max_items: int,
) -> str:
    """Render a custom section using its registered formatter."""
    try:
        entries = formatter.get_entries(adata, context)
    except Exception as e:  # noqa: BLE001
        # Intentional broad catch: custom formatters shouldn't crash the entire repr
        from anndata._warnings import warn

        warn(
            f"Custom section formatter '{formatter.section_name}' failed: {e}",
            UserWarning,
        )
        return ""

    if not entries:
        return ""

    n_items = len(entries)
    should_collapse = n_items > fold_threshold

    section_name = formatter.section_name
    display_name = getattr(formatter, "display_name", section_name)
    doc_url = getattr(formatter, "doc_url", None)
    tooltip = getattr(formatter, "tooltip", "")

    parts = [
        f'<div class="anndata-sec" data-section="{escape_html(section_name)}" '
        f'data-should-collapse="{str(should_collapse).lower()}">'
    ]

    # Header
    parts.append(
        _render_section_header(display_name, f"({n_items} items)", doc_url, tooltip)
    )

    # Content
    content_style = "padding:0;overflow:hidden;"
    parts.append(f'<div class="anndata-seccontent" style="{content_style}">')
    table_style = "width:100%;border-collapse:collapse;font-size:12px;"
    parts.append(f'<table class="adata-table" style="{table_style}">')

    for i, entry in enumerate(entries):
        if i >= max_items:
            parts.append(_render_truncation_indicator(n_items - max_items))
            break
        parts.append(_render_formatted_entry(entry, section_name))

    parts.append("</table>")
    parts.append("</div>")  # anndata-seccontent
    parts.append("</div>")  # anndata-sec

    return "\n".join(parts)


def _render_formatted_entry(entry: FormattedEntry, section: str) -> str:
    """Render a FormattedEntry from a custom section formatter."""
    output = entry.output

    # Button styles (hidden by default - JS shows them)
    btn_style = "display:none;border:none;background:transparent;cursor:pointer;font-size:11px;padding:2px;"
    expand_btn_style = "display:none;padding:2px 8px;font-size:11px;border-radius:4px;cursor:pointer;margin-left:8px;"

    entry_class = "adata-entry"
    if output.warnings:
        entry_class += " warning"
    if not output.is_serializable:
        entry_class += " error"

    has_expandable_content = output.html_content and output.is_expandable

    parts = [
        f'<tr class="{entry_class}" data-key="{escape_html(entry.key)}" '
        f'data-dtype="{escape_html(output.type_name)}">'
    ]

    # Name
    parts.append('<td class="adata-entry-name">')
    parts.append(escape_html(entry.key))
    parts.append(
        f'<button class="adata-copy-btn" style="{btn_style}" data-copy="{escape_html(entry.key)}" title="Copy name">📋</button>'
    )
    parts.append("</td>")

    # Type
    parts.append('<td class="adata-entry-type">')
    if output.warnings or not output.is_serializable:
        warnings_list = output.warnings.copy()
        if not output.is_serializable:
            warnings_list.insert(0, "Not serializable to H5AD/Zarr")
        title = escape_html("; ".join(warnings_list))
        parts.append(f'<span class="{output.css_class} dtype-warning" title="{title}">')
        parts.append(f"{escape_html(output.type_name)} ⚠️")
        parts.append("</span>")
    else:
        parts.append(
            f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>'
        )

    if has_expandable_content:
        parts.append(
            f'<button class="adata-expand-btn" style="{expand_btn_style}" aria-expanded="false">Expand ▼</button>'
        )

    if output.html_content and not output.is_expandable:
        parts.append(
            f'<div class="adata-custom-content" style="margin-top:4px;">{output.html_content}</div>'
        )

    parts.append("</td>")

    # Meta (empty for custom sections, or could be customized)
    parts.append('<td class="adata-entry-meta"></td>')

    parts.append("</tr>")

    # Expandable content
    if has_expandable_content:
        parts.append('<tr class="adata-nested-row">')
        parts.append('<td colspan="3" class="adata-nested-content">')
        parts.append(f'<div class="adata-custom-expanded">{output.html_content}</div>')
        parts.append("</td>")
        parts.append("</tr>")

    return "\n".join(parts)


# =============================================================================
# Section Renderers
# =============================================================================


def _render_header(
    adata: AnnData, *, show_search: bool = False, container_id: str = ""
) -> str:
    """Render the header with type, shape, badges, and optional search box."""
    # Use inline styles for layout only - colors handled by CSS for dark mode support
    header_style = (
        "display:flex;flex-wrap:wrap;align-items:center;gap:8px;"
        "padding:8px 12px;"
        "font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;"
    )
    parts = [f'<div class="anndata-hdr" style="{header_style}">']

    # Type name - allow for extension types
    type_name = type(adata).__name__
    type_style = "font-weight:600;font-size:14px;"
    parts.append(
        f'<span class="adata-type" style="{type_style}">{escape_html(type_name)}</span>'
    )

    # Shape
    shape_str = f"{format_number(adata.n_obs)} obs × {format_number(adata.n_vars)} vars"
    shape_style = "font-family:ui-monospace,monospace;font-size:12px;"
    parts.append(f'<span class="adata-shape" style="{shape_style}">{shape_str}</span>')

    # Badges
    if is_view(adata):
        parts.append('<span class="adata-badge adata-badge-view">View</span>')

    if is_backed(adata):
        backing = get_backing_info(adata)
        filename = backing.get("filename", "")
        format_str = backing.get("format", "")
        status = "Open" if backing.get("is_open") else "Closed"
        parts.append(
            f'<span class="adata-badge adata-badge-backed">'
            f"📁 {format_str} ({status})</span>"
        )
        # Inline file path (full path, no truncation)
        if filename:
            path_style = (
                "font-family:ui-monospace,monospace;font-size:11px;"
                "color:var(--anndata-text-secondary, #6c757d);"
            )
            parts.append(
                f'<span class="adata-file-path" style="{path_style}">'
                f"{escape_html(filename)}"
                f"</span>"
            )

    # Check for extension type (not standard AnnData)
    if type_name != "AnnData":
        parts.append(
            f'<span class="adata-badge adata-badge-extension">{type_name}</span>'
        )

    # Search box on the right (spacer pushes it right)
    if show_search:
        spacer_style = "flex-grow:1;"
        parts.append(f'<span style="{spacer_style}"></span>')
        # Search input hidden by default (JS shows it) - filter indicator uses CSS .active class
        search_style = "display:none;padding:4px 8px;font-size:11px;border-radius:4px;outline:none;width:150px;"
        search_id = f"{container_id}-search" if container_id else "anndata-search"
        parts.append(
            f'<input type="text" id="{search_id}" name="{search_id}" class="adata-search-input" style="{search_style}" '
            f'placeholder="Search..." aria-label="Search fields">'
        )
        # Filter indicator visibility controlled by CSS (.adata-filter-indicator.active) - no inline style
        parts.append('<span class="adata-filter-indicator"></span>')

    parts.append("</div>")
    return "\n".join(parts)


def _render_footer(adata: AnnData) -> str:
    """Render the footer with version and memory info."""
    footer_style = (
        "display:flex;justify-content:space-between;padding:4px 12px;"
        "font-size:10px;"
        "font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;"
    )
    parts = [f'<div class="anndata-ftr" style="{footer_style}">']

    # Version
    version = get_anndata_version()
    parts.append(f"<span>anndata v{version}</span>")

    # Memory usage
    try:
        mem_bytes = adata.__sizeof__()
        mem_str = format_memory_size(mem_bytes)
        parts.append(f'<span title="Estimated memory usage">~{mem_str}</span>')
    except Exception:  # noqa: BLE001
        # Intentional broad catch: __sizeof__ may not be implemented or may fail
        pass

    parts.append("</div>")
    return "\n".join(parts)


def _render_index_preview(adata: AnnData) -> str:
    """Render preview of obs_names and var_names."""
    parts = ['<div class="adata-index-preview">']

    # obs_names preview
    obs_preview = _format_index_preview(adata.obs_names, "obs_names")
    parts.append(f"<div><strong>obs_names:</strong> {obs_preview}</div>")

    # var_names preview
    var_preview = _format_index_preview(adata.var_names, "var_names")
    parts.append(f"<div><strong>var_names:</strong> {var_preview}</div>")

    parts.append("</div>")
    return "\n".join(parts)


def _format_index_preview(index: pd.Index, name: str) -> str:
    """Format a preview of an index."""
    n = len(index)
    if n == 0:
        return "<em>empty</em>"

    preview_n = DEFAULT_PREVIEW_ITEMS
    if n <= preview_n * 2:
        # Show all
        items = [escape_html(str(x)) for x in index]
    else:
        # Show first and last
        first = [escape_html(str(x)) for x in index[:preview_n]]
        last = [escape_html(str(x)) for x in index[-preview_n:]]
        items = [*first, "...", *last]

    return ", ".join(items)


def _render_x_entry(adata: AnnData, context: FormatterContext) -> str:
    """Render X as a single compact entry row."""
    X = adata.X

    # Inline styles for layout only - colors handled by CSS
    row_style = "display:flex;align-items:center;gap:12px;padding:6px 12px;"
    name_style = "font-family:ui-monospace,monospace;font-weight:600;min-width:60px;"
    type_style = "font-family:ui-monospace,monospace;font-size:11px;"

    parts = [f'<div class="adata-x-entry" style="{row_style}">']
    parts.append(f'<span style="{name_style}">X</span>')

    if X is None:
        parts.append(f'<span style="{type_style}"><em>None</em></span>')
    else:
        # Format the X matrix
        output = formatter_registry.format_value(X, context)

        # Build compact type string
        type_parts = [output.type_name]

        # Add sparsity info inline for sparse matrices
        if "sparsity" in output.details and output.details["sparsity"] is not None:
            sparsity = output.details["sparsity"]
            nnz = output.details.get("nnz", "?")
            type_parts.append(f"{sparsity:.1%} sparse ({format_number(nnz)} stored)")

        # Chunk info for Dask
        if "chunks" in output.details:
            type_parts.append(f"chunks={output.details['chunks']}")

        # Backed info
        if is_backed(adata):
            type_parts.append("📁 on disk")

        type_str = " · ".join(type_parts)
        parts.append(
            f'<span class="{output.css_class}" style="{type_style}">{escape_html(type_str)}</span>'
        )

    parts.append("</div>")
    return "\n".join(parts)


def _render_dataframe_section(
    adata: AnnData,
    section: str,
    context: FormatterContext,
    fold_threshold: int,
    max_items: int,
) -> str:
    """Render obs or var section."""
    df: pd.DataFrame = getattr(adata, section)
    n_cols = len(df.columns)

    # Doc URL and tooltip for this section
    doc_url = f"{DOCS_BASE_URL}generated/anndata.AnnData.{section}.html"
    tooltip = "Observation annotations" if section == "obs" else "Variable annotations"

    if n_cols == 0:
        return _render_empty_section(section, doc_url, tooltip)

    # Should this section be collapsed? (only via JS, default is expanded)
    should_collapse = n_cols > fold_threshold

    # Section - no inline colors to allow dark mode CSS to work
    parts = [
        f'<div class="anndata-sec" data-section="{section}" '
        f'data-should-collapse="{str(should_collapse).lower()}">'
    ]

    # Header
    parts.append(
        _render_section_header(section, f"({n_cols} columns)", doc_url, tooltip)
    )

    # Content - always visible by default (JS can hide it)
    content_style = "padding:0;overflow:hidden;"
    parts.append(f'<div class="anndata-seccontent" style="{content_style}">')
    table_style = "width:100%;border-collapse:collapse;font-size:12px;"
    parts.append(f'<table class="adata-table" style="{table_style}">')

    # Render each column
    for i, col_name in enumerate(df.columns):
        if i >= max_items:
            parts.append(_render_truncation_indicator(n_cols - max_items))
            break

        col = df[col_name]
        parts.append(_render_dataframe_entry(adata, section, col_name, col, context))

    parts.append("</table>")
    parts.append("</div>")  # anndata-seccontent
    parts.append("</div>")  # anndata-sec

    return "\n".join(parts)


def _render_dataframe_entry(
    adata: AnnData,
    section: str,
    col_name: str,
    col: pd.Series,
    context: FormatterContext,
) -> str:
    """Render a single DataFrame column entry."""
    # Format the column
    output = formatter_registry.format_value(col, context)

    # Check for string->category warning (skip for large columns)
    entry_warnings = list(output.warnings)
    unique_limit = _get_setting("repr_html_unique_limit", DEFAULT_UNIQUE_LIMIT)
    should_warn, warn_msg = should_warn_string_column(col, unique_limit)
    if should_warn:
        entry_warnings.append(warn_msg)

    # Check for color mismatch
    color_warning = check_color_category_mismatch(adata, col_name)
    if color_warning:
        entry_warnings.append(color_warning)

    # Get colors if categorical
    colors = get_matching_column_colors(adata, col_name)

    # Copy button hidden by default - JS shows it
    btn_style = "display:none;border:none;background:transparent;cursor:pointer;font-size:11px;padding:2px;"

    # Add warning class if needed (CSS handles color)
    entry_class = "adata-entry"
    if entry_warnings:
        entry_class += " warning"

    # Build row
    parts = [
        f'<tr class="{entry_class}" data-key="{escape_html(col_name)}" '
        f'data-dtype="{escape_html(output.type_name)}">'
    ]

    # Name cell
    parts.append('<td class="adata-entry-name">')
    parts.append(escape_html(col_name))
    parts.append(
        f'<button class="adata-copy-btn" style="{btn_style}" data-copy="{escape_html(col_name)}" title="Copy name">📋</button>'
    )
    parts.append("</td>")

    # Check if this is a categorical column (for wrap button)
    is_categorical = hasattr(col, "cat")
    categories = list(col.cat.categories) if is_categorical else []
    max_cats = _get_setting("repr_html_max_categories", DEFAULT_MAX_CATEGORIES)
    n_cats = min(len(categories), max_cats) if is_categorical else 0

    # Type cell
    parts.append('<td class="adata-entry-type">')
    if entry_warnings:
        title = escape_html("; ".join(entry_warnings))
        parts.append(f'<span class="{output.css_class} dtype-warning" title="{title}">')
        parts.append(f"{escape_html(output.type_name)} ⚠️")
        parts.append("</span>")
    else:
        parts.append(
            f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>'
        )

    # Add wrap button for categories in the type column
    if is_categorical and n_cats > 0:
        parts.append(
            '<button class="adata-cats-wrap-btn" title="Toggle multi-line view">⋯</button>'
        )
    parts.append("</td>")

    # Meta cell - show category values with colors
    parts.append('<td class="adata-entry-meta">')

    if is_categorical:
        cat_style = "display:inline-flex;align-items:center;gap:3px;margin-right:8px;"
        dot_style = "width:8px;height:8px;border-radius:50%;display:inline-block;"

        # Category list container (can be toggled to multi-line with CSS class)
        parts.append('<span class="adata-cats-list">')
        for i, cat in enumerate(categories[:max_cats]):
            cat_name = escape_html(str(cat))
            color = colors[i] if colors and i < len(colors) else None
            parts.append(f'<span class="adata-cat-item" style="{cat_style}">')
            if color:
                parts.append(
                    f'<span style="{dot_style}background:{escape_html(color)};"></span>'
                )
            parts.append(f"<span>{cat_name}</span>")
            parts.append("</span>")

        if len(categories) > max_cats:
            remaining = len(categories) - max_cats
            parts.append(f'<span class="adata-text-muted">...+{remaining}</span>')
        parts.append("</span>")

    elif hasattr(col, "nunique"):
        # Skip nunique() for very large columns to avoid performance issues
        unique_limit = _get_setting("repr_html_unique_limit", DEFAULT_UNIQUE_LIMIT)
        if unique_limit > 0 and len(col) <= unique_limit:
            try:
                n_unique = col.nunique()
                parts.append(
                    f'<span class="adata-text-muted">({n_unique} unique)</span>'
                )
            except Exception:  # noqa: BLE001
                # Intentional broad catch: nunique() can fail on unhashable types
                pass

    parts.append("</td>")

    parts.append("</tr>")
    return "\n".join(parts)


def _render_mapping_section(
    adata: AnnData,
    section: str,
    context: FormatterContext,
    fold_threshold: int,
    max_items: int,
) -> str:
    """Render obsm, varm, layers, obsp, varp sections."""
    mapping = getattr(adata, section, None)
    if mapping is None:
        return ""

    keys = list(mapping.keys())
    n_items = len(keys)

    # Doc URL and tooltip for this section
    doc_url = f"{DOCS_BASE_URL}generated/anndata.AnnData.{section}.html"
    tooltip = _get_section_tooltip(section)

    if n_items == 0:
        return _render_empty_section(section, doc_url, tooltip)

    should_collapse = n_items > fold_threshold

    # Section - no inline colors to allow dark mode CSS to work
    parts = [
        f'<div class="anndata-sec" data-section="{section}" '
        f'data-should-collapse="{str(should_collapse).lower()}">'
    ]

    # Header
    parts.append(
        _render_section_header(section, f"({n_items} items)", doc_url, tooltip)
    )

    # Content - always visible by default
    content_style = "padding:0;overflow:hidden;"
    parts.append(f'<div class="anndata-seccontent" style="{content_style}">')
    table_style = "width:100%;border-collapse:collapse;font-size:12px;"
    parts.append(f'<table class="adata-table" style="{table_style}">')

    for i, key in enumerate(keys):
        if i >= max_items:
            parts.append(_render_truncation_indicator(n_items - max_items))
            break

        value = mapping[key]
        parts.append(_render_mapping_entry(key, value, context, section))

    parts.append("</table>")
    parts.append("</div>")
    parts.append("</div>")

    return "\n".join(parts)


def _render_mapping_entry(
    key: str,
    value: Any,
    context: FormatterContext,
    section: str,
) -> str:
    """Render a single mapping entry."""
    output = formatter_registry.format_value(value, context)

    # Button styles (hidden by default - JS shows them)
    btn_style = "display:none;border:none;background:transparent;cursor:pointer;font-size:11px;padding:2px;"
    expand_btn_style = "display:none;padding:2px 8px;font-size:11px;border-radius:4px;cursor:pointer;margin-left:8px;"

    # Build class list for CSS styling
    entry_class = "adata-entry"
    if output.warnings:
        entry_class += " warning"
    if not output.is_serializable:
        entry_class += " error"

    # Check if we have expandable custom HTML content
    has_expandable_content = output.html_content and output.is_expandable

    parts = [
        f'<tr class="{entry_class}" data-key="{escape_html(key)}" '
        f'data-dtype="{escape_html(output.type_name)}">'
    ]

    # Name
    parts.append('<td class="adata-entry-name">')
    parts.append(escape_html(key))
    parts.append(
        f'<button class="adata-copy-btn" style="{btn_style}" data-copy="{escape_html(key)}" title="Copy name">📋</button>'
    )
    parts.append("</td>")

    # Type
    parts.append('<td class="adata-entry-type">')
    if output.warnings or not output.is_serializable:
        warnings = output.warnings.copy()
        if not output.is_serializable:
            warnings.insert(0, "Not serializable to H5AD/Zarr")
        title = escape_html("; ".join(warnings))
        parts.append(f'<span class="{output.css_class} dtype-warning" title="{title}">')
        parts.append(f"{escape_html(output.type_name)} ⚠️")
        parts.append("</span>")
    else:
        parts.append(
            f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>'
        )

    # Add expand button for custom HTML content
    if has_expandable_content:
        parts.append(
            f'<button class="adata-expand-btn" style="{expand_btn_style}" aria-expanded="false">Expand ▼</button>'
        )

    # Add wrap button for DataFrame columns list
    has_columns_list = output.details.get("has_columns_list", False)
    if has_columns_list:
        parts.append(
            '<button class="adata-cols-wrap-btn" title="Toggle multi-line view">⋯</button>'
        )

    # Inline (non-expandable) custom HTML content
    if output.html_content and not output.is_expandable:
        parts.append(
            f'<div class="adata-custom-content" style="margin-top:4px;">{output.html_content}</div>'
        )

    parts.append("</td>")

    # Meta - show shape/cols for obsm/varm, or custom meta_preview (DataFrame columns)
    parts.append('<td class="adata-entry-meta">')
    if has_columns_list and "columns" in output.details:
        # Render DataFrame columns as a wrappable list
        columns = output.details["columns"]
        col_str = ", ".join(escape_html(str(c)) for c in columns)
        parts.append(f'<span class="adata-cols-list">[{col_str}]</span>')
    elif "meta_preview" in output.details:
        # Other meta preview (non-DataFrame)
        parts.append(
            f'<span title="{escape_html(output.details.get("meta_preview_full", output.details["meta_preview"]))}">{escape_html(output.details["meta_preview"])}</span>'
        )
    elif "shape" in output.details and section in ("obsm", "varm"):
        shape = output.details["shape"]
        if len(shape) >= 2:
            parts.append(f"({format_number(shape[1])} cols)")
    parts.append("</td>")

    parts.append("</tr>")

    # Expandable custom HTML content (hidden by default, shown on expand)
    if has_expandable_content:
        parts.append('<tr class="adata-nested-row">')
        parts.append('<td colspan="3" class="adata-nested-content">')
        parts.append(f'<div class="adata-custom-expanded">{output.html_content}</div>')
        parts.append("</td>")
        parts.append("</tr>")

    return "\n".join(parts)


def _render_uns_section(
    adata: AnnData,
    context: FormatterContext,
    fold_threshold: int,
    max_items: int,
    max_depth: int,
) -> str:
    """Render the uns section with special handling."""
    uns = adata.uns
    keys = list(uns.keys())
    n_items = len(keys)

    # Doc URL and tooltip
    doc_url = f"{DOCS_BASE_URL}generated/anndata.AnnData.uns.html"
    tooltip = "Unstructured annotation"

    if n_items == 0:
        return _render_empty_section("uns", doc_url, tooltip)

    should_collapse = n_items > fold_threshold

    # Section - use data-should-collapse for JS to handle (consistent with other sections)
    parts = [
        f'<div class="anndata-sec" data-section="uns" '
        f'data-should-collapse="{str(should_collapse).lower()}">'
    ]

    # Header
    parts.append(_render_section_header("uns", f"({n_items} items)", doc_url, tooltip))

    # Content - with inline styles for consistency
    content_style = "padding:0;overflow:hidden;"
    parts.append(f'<div class="anndata-seccontent" style="{content_style}">')
    table_style = "width:100%;border-collapse:collapse;font-size:12px;"
    parts.append(f'<table class="adata-table" style="{table_style}">')

    for i, key in enumerate(keys):
        if i >= max_items:
            parts.append(_render_truncation_indicator(n_items - max_items))
            break

        value = uns[key]
        parts.append(_render_uns_entry(adata, key, value, context, max_depth))

    parts.append("</table>")
    parts.append("</div>")
    parts.append("</div>")

    return "\n".join(parts)


def _render_uns_entry(
    adata: AnnData,
    key: str,
    value: Any,
    context: FormatterContext,
    max_depth: int,
) -> str:
    """Render a single uns entry with special type handling.

    Rendering priority:
    1. Custom TypeFormatter (checks type hint via can_format)
    2. Color list (keys ending in _colors)
    3. Nested AnnData objects
    4. Generic preview for simple types (str, int, float, bool, small dict/list)
    5. Default type info only
    """
    # 1. Check for custom TypeFormatter (may check type hint in can_format)
    # First try formatting - if a TypeFormatter matches and provides html_content, use it
    output = formatter_registry.format_value(value, context)
    if output.html_content:
        # A TypeFormatter provided custom HTML - render it
        return _render_uns_entry_with_custom_html(key, output)

    # Check if there's an unhandled type hint (no formatter matched but hint exists)
    type_hint, cleaned_value = extract_uns_type_hint(value)
    if type_hint is not None:
        # Type hint present but no formatter registered - show helpful message
        package_name = type_hint.split(".")[0] if "." in type_hint else type_hint
        return _render_uns_entry_with_preview(
            key,
            cleaned_value,
            context,
            preview_note=f"[{type_hint}] (import {package_name} to enable)",
        )

    # 2. Check for color list
    if is_color_list(key, value):
        return _render_color_list_entry(key, value)

    # 3. Check for nested AnnData
    if type(value).__name__ == "AnnData" and hasattr(value, "n_obs"):
        return _render_nested_anndata_entry(key, value, context, max_depth)

    # 4. Generic preview for simple/small types
    return _render_uns_entry_with_preview(key, value, context)


def _render_uns_entry_with_preview(
    key: str,
    value: Any,
    context: FormatterContext,
    preview_note: str | None = None,
) -> str:
    """Render an uns entry with a value preview in the meta column.

    Generates previews for:
    - Strings: truncated to max length
    - Numbers (int, float): shown directly
    - Booleans: shown directly
    - Small dicts: key count and first few keys
    - Small lists: length and first few items
    - Other types: type info only
    """
    output = formatter_registry.format_value(value, context)
    max_str_len = _get_setting("repr_html_max_string_length", DEFAULT_MAX_STRING_LENGTH)

    # Generate preview based on type
    preview = _generate_value_preview(value, max_str_len)
    if preview_note:
        preview = f"{preview_note} {preview}" if preview else preview_note

    # Inline styles for layout - colors via CSS
    name_style = "padding:6px 12px;font-family:ui-monospace,monospace;font-weight:500;"
    type_style = "padding:6px 12px;font-family:ui-monospace,monospace;font-size:11px;"
    meta_style = "padding:6px 12px;font-size:11px;text-align:left;max-width:300px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;"
    btn_style = "display:none;border:none;background:transparent;cursor:pointer;font-size:11px;padding:2px;"

    entry_class = "adata-entry"
    if output.warnings:
        entry_class += " warning"
    if not output.is_serializable:
        entry_class += " error"

    parts = [
        f'<tr class="{entry_class}" data-key="{escape_html(key)}" data-dtype="{escape_html(output.type_name)}">'
    ]

    # Name
    parts.append(f'<td class="adata-entry-name" style="{name_style}">')
    parts.append(escape_html(key))
    parts.append(
        f'<button class="adata-copy-btn" style="{btn_style}" data-copy="{escape_html(key)}" title="Copy name">📋</button>'
    )
    parts.append("</td>")

    # Type
    parts.append(f'<td class="adata-entry-type" style="{type_style}">')
    if output.warnings or not output.is_serializable:
        warn_list = output.warnings.copy()
        if not output.is_serializable:
            warn_list.insert(0, "Not serializable to H5AD/Zarr")
        title = escape_html("; ".join(warn_list))
        parts.append(f'<span class="{output.css_class} dtype-warning" title="{title}">')
        parts.append(f"{escape_html(output.type_name)} ⚠️")
        parts.append("</span>")
    else:
        parts.append(
            f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>'
        )
    parts.append("</td>")

    # Meta - value preview
    parts.append(f'<td class="adata-entry-meta" style="{meta_style}">')
    if preview:
        parts.append(
            f'<span class="adata-text-muted" title="{escape_html(str(value)[:500])}">{escape_html(preview)}</span>'
        )
    parts.append("</td>")
    parts.append("</tr>")

    return "\n".join(parts)


def _generate_value_preview(value: Any, max_len: int = 100) -> str:
    """Generate a human-readable preview of a value.

    Returns empty string if no meaningful preview can be generated.
    """
    # Strings
    if isinstance(value, str):
        if len(value) <= max_len:
            return f'"{value}"'
        return f'"{value[:max_len]}..."'

    # Numbers
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, (int, np.integer)):
        return str(value)
    if isinstance(value, (float, np.floating)):
        # Format floats nicely
        if value == int(value):
            return str(int(value))
        return f"{value:.6g}"

    # Small dicts
    if isinstance(value, dict):
        n_keys = len(value)
        if n_keys == 0:
            return "{}"
        if n_keys <= 3:
            keys_preview = ", ".join(str(k) for k in list(value.keys())[:3])
            return f"{{{keys_preview}}}"
        keys_preview = ", ".join(str(k) for k in list(value.keys())[:2])
        return f"{{{keys_preview}, ...}} ({n_keys} keys)"

    # Small lists/tuples
    if isinstance(value, (list, tuple)):
        n_items = len(value)
        if n_items == 0:
            return "[]" if isinstance(value, list) else "()"
        if n_items <= 3:
            # Show items if they're simple
            try:
                items = [_preview_item(v) for v in value[:3]]
                if all(items):
                    bracket = "[]" if isinstance(value, list) else "()"
                    return f"{bracket[0]}{', '.join(items)}{bracket[1]}"
            except Exception:  # noqa: BLE001
                # Intentional broad catch: preview generation is best-effort
                pass
        return f"({n_items} items)"

    # None
    if value is None:
        return "None"

    # No preview for complex types
    return ""


def _preview_item(value: Any) -> str:
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


def _render_uns_entry_with_custom_html(key: str, output: FormattedOutput) -> str:
    """Render an uns entry with custom HTML from a TypeFormatter.

    The output should have html_content set.
    """
    # Inline styles for layout
    name_style = "padding:6px 12px;font-family:ui-monospace,monospace;font-weight:500;"
    type_style = "padding:6px 12px;font-family:ui-monospace,monospace;font-size:11px;"
    meta_style = "padding:6px 12px;font-size:11px;text-align:left;"
    btn_style = "display:none;border:none;background:transparent;cursor:pointer;font-size:11px;padding:2px;"

    type_label = output.type_name

    parts = [
        f'<tr class="adata-entry" data-key="{escape_html(key)}" data-dtype="{escape_html(type_label)}">'
    ]

    # Name
    parts.append(f'<td class="adata-entry-name" style="{name_style}">')
    parts.append(escape_html(key))
    parts.append(
        f'<button class="adata-copy-btn" style="{btn_style}" data-copy="{escape_html(key)}" title="Copy name">📋</button>'
    )
    parts.append("</td>")

    # Type
    parts.append(f'<td class="adata-entry-type" style="{type_style}">')
    parts.append(f'<span class="{output.css_class}">{escape_html(type_label)}</span>')
    parts.append("</td>")

    # Meta - custom HTML content
    parts.append(f'<td class="adata-entry-meta" style="{meta_style}">')
    # Note: HTML is trusted from registered TypeFormatter (package must be imported)
    parts.append(output.html_content)
    parts.append("</td>")
    parts.append("</tr>")

    return "\n".join(parts)


def _render_color_list_entry(key: str, value: Any) -> str:
    """Render a color list entry with swatches."""
    colors = list(value) if hasattr(value, "__iter__") else []
    n_colors = len(colors)

    # Inline styles for layout - colors via CSS
    name_style = "padding:6px 12px;font-family:ui-monospace,monospace;font-weight:500;"
    type_style = "padding:6px 12px;font-family:ui-monospace,monospace;font-size:11px;"
    meta_style = "padding:6px 12px;font-size:11px;text-align:right;"
    btn_style = "display:none;border:none;background:transparent;cursor:pointer;font-size:11px;padding:2px;"

    parts = [
        f'<tr class="adata-entry" data-key="{escape_html(key)}" data-dtype="colors">'
    ]

    # Name
    parts.append(f'<td class="adata-entry-name" style="{name_style}">')
    parts.append(escape_html(key))
    parts.append(
        f'<button class="adata-copy-btn" style="{btn_style}" data-copy="{escape_html(key)}" title="Copy name">📋</button>'
    )
    parts.append("</td>")

    # Type
    parts.append(f'<td class="adata-entry-type" style="{type_style}">')
    parts.append(f'<span class="dtype-object">colors ({n_colors})</span>')
    parts.append("</td>")

    # Meta - color swatches
    parts.append(f'<td class="adata-entry-meta" style="{meta_style}">')
    parts.append('<span class="adata-color-swatches">')
    for color in colors[:15]:  # Limit preview
        parts.append(
            f'<span class="adata-color-swatch" style="background:{escape_html(str(color))}" title="{escape_html(str(color))}"></span>'
        )
    if n_colors > 15:
        parts.append(f'<span class="adata-text-muted">+{n_colors - 15}</span>')
    parts.append("</span>")
    parts.append("</td>")
    parts.append("</tr>")

    return "\n".join(parts)


def _render_nested_anndata_entry(
    key: str,
    value: Any,
    context: FormatterContext,
    max_depth: int,
) -> str:
    """Render a nested AnnData entry."""
    n_obs = getattr(value, "n_obs", "?")
    n_vars = getattr(value, "n_vars", "?")

    can_expand = context.depth < max_depth - 1

    # Inline styles for layout - colors via CSS
    name_style = "padding:6px 12px;font-family:ui-monospace,monospace;font-weight:500;"
    type_style = "padding:6px 12px;font-family:ui-monospace,monospace;font-size:11px;"
    meta_style = "padding:6px 12px;font-size:11px;text-align:right;"
    btn_style = "display:none;border:none;background:transparent;cursor:pointer;font-size:11px;padding:2px;"
    # Expand button hidden by default - JS shows it
    expand_btn_style = "display:none;padding:2px 8px;font-size:11px;border-radius:4px;cursor:pointer;margin-left:8px;"

    parts = [
        f'<tr class="adata-entry" data-key="{escape_html(key)}" data-dtype="AnnData">'
    ]

    # Name
    parts.append(f'<td class="adata-entry-name" style="{name_style}">')
    parts.append(escape_html(key))
    parts.append(
        f'<button class="adata-copy-btn" style="{btn_style}" data-copy="{escape_html(key)}" title="Copy name">📋</button>'
    )
    parts.append("</td>")

    # Type
    parts.append(f'<td class="adata-entry-type" style="{type_style}">')
    parts.append(
        f'<span class="dtype-anndata">AnnData ({format_number(n_obs)} × {format_number(n_vars)})</span>'
    )
    if can_expand:
        parts.append(
            f'<button class="adata-expand-btn" style="{expand_btn_style}" aria-expanded="false">Expand ▼</button>'
        )
    parts.append("</td>")

    parts.append(f'<td class="adata-entry-meta" style="{meta_style}"></td>')
    parts.append("</tr>")

    # Nested content (hidden by default)
    if can_expand:
        parts.append('<tr class="adata-nested-row">')
        parts.append('<td colspan="3" class="adata-nested-content">')
        parts.append('<div class="adata-nested-anndata">')
        # Recursive call
        nested_html = generate_repr_html(
            value,
            depth=context.depth + 1,
            max_depth=max_depth,
            show_header=True,
            show_search=False,
        )
        parts.append(nested_html)
        parts.append("</div>")
        parts.append("</td>")
        parts.append("</tr>")

    return "\n".join(parts)


def _render_raw_section(
    adata: AnnData,
    context: FormatterContext,
    fold_threshold: int,
    max_items: int,
) -> str:
    """Render the raw section."""
    raw = getattr(adata, "raw", None)
    if raw is None:
        return ""

    # Raw section always starts collapsed (use data-should-collapse for consistency)
    parts = ['<div class="anndata-sec" data-section="raw" data-should-collapse="true">']

    # Header
    doc_url = f"{DOCS_BASE_URL}generated/anndata.AnnData.raw.html"
    n_vars = getattr(raw, "n_vars", "?")
    parts.append(
        _render_section_header(
            "raw",
            f"(n_vars = {format_number(n_vars)})",
            doc_url,
            "Raw data (original unprocessed)",
        )
    )

    # Content with inline styles
    content_style = "padding:0;overflow:hidden;"
    parts.append(f'<div class="anndata-seccontent" style="{content_style}">')

    # Info items with inline styles
    info_style = "padding:6px 12px;font-size:12px;"

    # raw.X info
    if hasattr(raw, "X") and raw.X is not None:
        output = formatter_registry.format_value(raw.X, context)
        parts.append(
            f'<div style="{info_style}"><strong>raw.X:</strong> <span class="{output.css_class}">{escape_html(output.type_name)}</span></div>'
        )

    # raw.var columns
    if hasattr(raw, "var") and len(raw.var.columns) > 0:
        parts.append(
            f'<div style="{info_style}"><strong>raw.var:</strong> {len(raw.var.columns)} columns</div>'
        )

    # raw.varm
    if hasattr(raw, "varm") and len(raw.varm) > 0:
        parts.append(
            f'<div style="{info_style}"><strong>raw.varm:</strong> {len(raw.varm)} items</div>'
        )

    parts.append("</div>")
    parts.append("</div>")

    return "\n".join(parts)


# =============================================================================
# Helper Functions
# =============================================================================


def _render_section_header(
    name: str,
    count_str: str,
    doc_url: str | None,
    tooltip: str,
) -> str:
    """Render a section header - colors handled by CSS for dark mode support."""
    # Layout-only inline styles - colors handled by CSS
    header_style = (
        "display:flex;align-items:center;gap:8px;padding:8px 12px;cursor:pointer;"
        "font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;"
    )
    name_style = "font-weight:600;"
    count_style = "font-size:11px;"
    # Fold icon hidden by default, shown via JS - centered for proper rotation
    fold_style = (
        "display:none;width:16px;height:16px;font-size:10px;"
        "align-items:center;justify-content:center;transform-origin:center;flex-shrink:0;"
    )
    link_style = "margin-left:auto;padding:2px 6px;font-size:11px;text-decoration:none;"

    parts = [f'<div class="anndata-sechdr" style="{header_style}">']
    parts.append(f'<span class="adata-fold-icon" style="{fold_style}">▼</span>')
    parts.append(
        f'<span class="anndata-sec-name" style="{name_style}">{escape_html(name)}</span>'
    )
    parts.append(
        f'<span class="anndata-sec-count" style="{count_style}">{escape_html(count_str)}</span>'
    )
    if doc_url:
        parts.append(
            f'<a class="adata-help-link" style="{link_style}" href="{escape_html(doc_url)}" target="_blank" title="{escape_html(tooltip)}">?</a>'
        )
    parts.append("</div>")
    return "\n".join(parts)


def _render_empty_section(
    name: str,
    doc_url: str | None = None,
    tooltip: str = "",
) -> str:
    """Render an empty section indicator."""
    # Use data-should-collapse for consistency - empty sections always collapsed
    header_style = (
        "display:flex;align-items:center;gap:8px;padding:8px 12px;cursor:pointer;"
        "font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;"
    )
    # Fold icon hidden by default, shown via JS - centered for proper rotation
    fold_style = (
        "display:none;width:16px;height:16px;font-size:10px;"
        "align-items:center;justify-content:center;transform-origin:center;flex-shrink:0;"
    )
    name_style = "font-weight:600;"
    count_style = "font-size:11px;"
    link_style = "margin-left:auto;padding:2px 6px;font-size:11px;text-decoration:none;"
    content_style = "padding:0;overflow:hidden;"
    empty_style = "padding:8px 12px;font-size:11px;font-style:italic;"

    # Build help link if doc_url provided
    help_link = ""
    if doc_url:
        help_link = f'<a class="adata-help-link" style="{link_style}" href="{escape_html(doc_url)}" target="_blank" title="{escape_html(tooltip)}">?</a>'

    return f"""
<div class="anndata-sec" data-section="{escape_html(name)}" data-should-collapse="true">
    <div class="anndata-sechdr" style="{header_style}">
        <span class="adata-fold-icon" style="{fold_style}">▼</span>
        <span class="anndata-sec-name" style="{name_style}">{escape_html(name)}</span>
        <span class="anndata-sec-count" style="{count_style}">(empty)</span>
        {help_link}
    </div>
    <div class="anndata-seccontent" style="{content_style}">
        <div class="adata-empty" style="{empty_style}">No entries</div>
    </div>
</div>
"""


def _render_truncation_indicator(remaining: int) -> str:
    """Render a truncation indicator."""
    return f'<tr><td colspan="3" class="adata-truncated">... and {format_number(remaining)} more</td></tr>'


def _render_max_depth_indicator(adata: AnnData) -> str:
    """Render indicator when max depth is reached."""
    n_obs = getattr(adata, "n_obs", "?")
    n_vars = getattr(adata, "n_vars", "?")
    return f'<div class="adata-max-depth">AnnData ({format_number(n_obs)} × {format_number(n_vars)}) - max depth reached</div>'


def _get_section_tooltip(section: str) -> str:
    """Get tooltip text for a section."""
    tooltips = {
        "obs": "Observation (cell) annotations",
        "var": "Variable (gene) annotations",
        "uns": "Unstructured annotation",
        "obsm": "Multi-dimensional observation annotations",
        "varm": "Multi-dimensional variable annotations",
        "layers": "Additional data layers (same shape as X)",
        "obsp": "Pairwise observation annotations",
        "varp": "Pairwise variable annotations",
        "raw": "Raw data (original unprocessed)",
    }
    return tooltips.get(section, "")


def _get_setting(name: str, default: Any) -> Any:
    """Get a setting value, falling back to default."""
    try:
        from anndata import settings

        return getattr(settings, name, default)
    except (ImportError, AttributeError):
        return default
