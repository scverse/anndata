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
import warnings
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
    formatter_registry,
    # Uns renderer registry for custom serialized data
    uns_renderer_registry,
    extract_uns_type_hint,
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
    is_serializable,
    is_view,
    sanitize_for_id,
    should_warn_string_column,
    truncate_string,
)

if TYPE_CHECKING:
    from anndata import AnnData

# Import formatters to register them
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
        fold_threshold = _get_setting("repr_html_fold_threshold", DEFAULT_FOLD_THRESHOLD)
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
        parts.append(_render_header(adata, show_search=show_search and depth == 0, container_id=container_id))

    # Index preview (only at top level)
    if depth == 0:
        parts.append(_render_index_preview(adata))

    # Sections container
    parts.append('<div class="adata-sections">')

    # X as a simple entry (like layers)
    parts.append(_render_x_entry(adata, context))

    # Standard sections
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
                _render_uns_section(adata, context, fold_threshold, max_items, max_depth)
            )
        else:
            parts.append(
                _render_mapping_section(
                    adata, section, context, fold_threshold, max_items
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
# Section Renderers
# =============================================================================


def _render_header(adata: AnnData, *, show_search: bool = False, container_id: str = "") -> str:
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
    parts.append(f'<span class="adata-type" style="{type_style}">{escape_html(type_name)}</span>')

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
            f'📁 {format_str} ({status})</span>'
        )
        # Inline file path (truncated with full path on hover)
        if filename:
            path_style = (
                "font-family:ui-monospace,monospace;font-size:11px;"
                "color:var(--anndata-text-secondary, #6c757d);max-width:300px;"
                "overflow:hidden;text-overflow:ellipsis;white-space:nowrap;"
            )
            parts.append(
                f'<span class="adata-file-path" style="{path_style}" title="{escape_html(filename)}">'
                f'{escape_html(filename)}'
                f'</span>'
            )

    # Check for extension type (not standard AnnData)
    if type_name != "AnnData":
        parts.append(f'<span class="adata-badge adata-badge-extension">{type_name}</span>')

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
    except Exception:
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
        items = first + ["..."] + last

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
        parts.append(f'<span class="{output.css_class}" style="{type_style}">{escape_html(type_str)}</span>')

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

    if n_cols == 0:
        return _render_empty_section(section)

    # Should this section be collapsed? (only via JS, default is expanded)
    should_collapse = n_cols > fold_threshold

    # Section - no inline colors to allow dark mode CSS to work
    parts = [
        f'<div class="anndata-sec" data-section="{section}" '
        f'data-should-collapse="{str(should_collapse).lower()}">'
    ]

    # Header
    doc_url = f"{DOCS_BASE_URL}generated/anndata.AnnData.{section}.html"
    tooltip = "Observation annotations" if section == "obs" else "Variable annotations"
    parts.append(_render_section_header(section, f"({n_cols} columns)", doc_url, tooltip))

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

    # Check for string->category warning
    warnings = list(output.warnings)
    should_warn, warn_msg = should_warn_string_column(col)
    if should_warn:
        warnings.append(warn_msg)

    # Check for color mismatch
    color_warning = check_color_category_mismatch(adata, col_name)
    if color_warning:
        warnings.append(color_warning)

    # Get colors if categorical
    colors = get_matching_column_colors(adata, col_name)

    # Inline styles for layout only - colors handled by CSS
    row_style = ""
    name_style = "padding:6px 12px;font-family:ui-monospace,monospace;font-weight:500;"
    type_style = "padding:6px 12px;font-family:ui-monospace,monospace;font-size:11px;"
    # Copy button hidden by default - JS shows it
    btn_style = "display:none;border:none;background:transparent;cursor:pointer;font-size:11px;padding:2px;"

    # Add warning class if needed (CSS handles color)
    entry_class = "adata-entry"
    if warnings:
        entry_class += " warning"

    # Build row
    parts = [
        f'<tr class="{entry_class}" data-key="{escape_html(col_name)}" '
        f'data-dtype="{escape_html(output.type_name)}">'
    ]

    # Name cell
    parts.append(f'<td class="adata-entry-name" style="{name_style}">')
    parts.append(escape_html(col_name))
    parts.append(f'<button class="adata-copy-btn" style="{btn_style}" data-copy="{escape_html(col_name)}" title="Copy name">📋</button>')
    parts.append("</td>")

    # Type cell (colors now shown with category names in meta cell)
    parts.append(f'<td class="adata-entry-type" style="{type_style}">')
    if warnings:
        title = escape_html("; ".join(warnings))
        parts.append(f'<span class="{output.css_class} dtype-warning" title="{title}">')
        parts.append(f"{escape_html(output.type_name)} ⚠️")
        parts.append("</span>")
    else:
        parts.append(f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>')
    parts.append("</td>")

    # Meta cell - show category values with colors
    meta_style = "padding:6px 12px;font-size:11px;text-align:left;"
    parts.append(f'<td class="adata-entry-meta" style="{meta_style}">')

    if hasattr(col, "cat"):
        categories = list(col.cat.categories)
        max_cats = _get_setting("repr_html_max_categories", DEFAULT_MAX_CATEGORIES)
        cat_style = "display:inline-flex;align-items:center;gap:3px;margin-right:8px;"
        dot_style = "width:8px;height:8px;border-radius:50%;display:inline-block;"

        for i, cat in enumerate(categories[:max_cats]):
            cat_name = escape_html(str(cat))
            # Get color for this category if available
            color = colors[i] if colors and i < len(colors) else None

            parts.append(f'<span style="{cat_style}">')
            if color:
                parts.append(f'<span style="{dot_style}background:{escape_html(color)};"></span>')
            parts.append(f'<span>{cat_name}</span>')
            parts.append('</span>')

        if len(categories) > max_cats:
            remaining = len(categories) - max_cats
            parts.append(f'<span class="adata-text-muted">...+{remaining}</span>')

    elif hasattr(col, "nunique"):
        # Skip nunique() for very large columns to avoid performance issues
        unique_limit = _get_setting("repr_html_unique_limit", DEFAULT_UNIQUE_LIMIT)
        if unique_limit > 0 and len(col) <= unique_limit:
            try:
                n_unique = col.nunique()
                parts.append(f'<span class="adata-text-muted">({n_unique} unique)</span>')
            except Exception:
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

    if n_items == 0:
        return _render_empty_section(section)

    should_collapse = n_items > fold_threshold

    # Section - no inline colors to allow dark mode CSS to work
    parts = [
        f'<div class="anndata-sec" data-section="{section}" '
        f'data-should-collapse="{str(should_collapse).lower()}">'
    ]

    # Header
    doc_url = f"{DOCS_BASE_URL}generated/anndata.AnnData.{section}.html"
    tooltip = _get_section_tooltip(section)
    parts.append(_render_section_header(section, f"({n_items} items)", doc_url, tooltip))

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

    # Inline styles for layout only - colors handled by CSS
    name_style = "padding:6px 12px;font-family:ui-monospace,monospace;font-weight:500;"
    type_style = "padding:6px 12px;font-family:ui-monospace,monospace;font-size:11px;"
    btn_style = "display:none;border:none;background:transparent;cursor:pointer;font-size:11px;padding:2px;"

    # Build class list for CSS styling
    entry_class = "adata-entry"
    if output.warnings:
        entry_class += " warning"
    if not output.is_serializable:
        entry_class += " error"

    parts = [
        f'<tr class="{entry_class}" data-key="{escape_html(key)}" '
        f'data-dtype="{escape_html(output.type_name)}">'
    ]

    # Name
    parts.append(f'<td class="adata-entry-name" style="{name_style}">')
    parts.append(escape_html(key))
    parts.append(f'<button class="adata-copy-btn" style="{btn_style}" data-copy="{escape_html(key)}" title="Copy name">📋</button>')
    parts.append("</td>")

    # Type
    parts.append(f'<td class="adata-entry-type" style="{type_style}">')
    if output.warnings or not output.is_serializable:
        warnings = output.warnings.copy()
        if not output.is_serializable:
            warnings.insert(0, "Not serializable to H5AD/Zarr")
        title = escape_html("; ".join(warnings))
        parts.append(f'<span class="{output.css_class} dtype-warning" title="{title}">')
        parts.append(f"{escape_html(output.type_name)} ⚠️")
        parts.append("</span>")
    else:
        parts.append(f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>')
    parts.append("</td>")

    # Meta - show shape/cols for obsm/varm only (layers are always n_vars cols)
    meta_style = "padding:6px 12px;font-size:11px;text-align:right;"
    parts.append(f'<td class="adata-entry-meta" style="{meta_style}">')
    if "shape" in output.details and section in ("obsm", "varm"):
        shape = output.details["shape"]
        if len(shape) >= 2:
            parts.append(f"({format_number(shape[1])} cols)")
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

    if n_items == 0:
        return _render_empty_section("uns")

    should_collapse = n_items > fold_threshold

    # Section - use data-should-collapse for JS to handle (consistent with other sections)
    parts = [
        f'<div class="anndata-sec" data-section="uns" '
        f'data-should-collapse="{str(should_collapse).lower()}">'
    ]

    # Header
    doc_url = f"{DOCS_BASE_URL}generated/anndata.AnnData.uns.html"
    parts.append(_render_section_header("uns", f"({n_items} items)", doc_url, "Unstructured annotation"))

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
        parts.append(
            _render_uns_entry(adata, key, value, context, max_depth)
        )

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
    1. Custom renderer via type hint (__anndata_repr__)
    2. Color list (keys ending in _colors)
    3. Nested AnnData objects
    4. Generic preview for simple types (str, int, float, bool, small dict/list)
    5. Default type info only
    """
    # 1. Check for custom renderer via type hint
    type_hint, cleaned_value = extract_uns_type_hint(value)
    if type_hint is not None:
        renderer = uns_renderer_registry.get_renderer(type_hint)
        if renderer is not None:
            try:
                result = renderer(cleaned_value, context)
                return _render_custom_uns_entry(key, type_hint, result)
            except Exception as e:
                # Fall through to default rendering with warning
                warnings.warn(
                    f"Custom renderer for '{type_hint}' failed: {e}",
                    stacklevel=2,
                )
        # Type hint present but no renderer registered - show helpful message
        # Extract package name from type hint (e.g., "mypackage.config" -> "mypackage")
        package_name = type_hint.split(".")[0] if "." in type_hint else type_hint
        return _render_uns_entry_with_preview(
            key, cleaned_value, context,
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

    parts = [f'<tr class="{entry_class}" data-key="{escape_html(key)}" data-dtype="{escape_html(output.type_name)}">']

    # Name
    parts.append(f'<td class="adata-entry-name" style="{name_style}">')
    parts.append(escape_html(key))
    parts.append(f'<button class="adata-copy-btn" style="{btn_style}" data-copy="{escape_html(key)}" title="Copy name">📋</button>')
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
        parts.append(f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>')
    parts.append("</td>")

    # Meta - value preview
    parts.append(f'<td class="adata-entry-meta" style="{meta_style}">')
    if preview:
        parts.append(f'<span class="adata-text-muted" title="{escape_html(str(value)[:500])}">{escape_html(preview)}</span>')
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
            except Exception:
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


def _render_custom_uns_entry(key: str, type_hint: str, result: Any) -> str:
    """Render an uns entry with custom HTML from a registered renderer.

    The result should be an UnsRendererOutput with sanitized HTML.
    """
    from anndata._repr.registry import UnsRendererOutput

    if not isinstance(result, UnsRendererOutput):
        # Fallback if renderer didn't return proper type
        return _render_uns_entry_with_preview(
            key, result, FormatterContext(),
            preview_note=f"[{type_hint}]",
        )

    # Inline styles for layout
    name_style = "padding:6px 12px;font-family:ui-monospace,monospace;font-weight:500;"
    type_style = "padding:6px 12px;font-family:ui-monospace,monospace;font-size:11px;"
    meta_style = "padding:6px 12px;font-size:11px;text-align:left;"
    btn_style = "display:none;border:none;background:transparent;cursor:pointer;font-size:11px;padding:2px;"

    type_label = result.type_label or type_hint

    parts = [f'<tr class="adata-entry" data-key="{escape_html(key)}" data-dtype="{escape_html(type_label)}">']

    # Name
    parts.append(f'<td class="adata-entry-name" style="{name_style}">')
    parts.append(escape_html(key))
    parts.append(f'<button class="adata-copy-btn" style="{btn_style}" data-copy="{escape_html(key)}" title="Copy name">📋</button>')
    parts.append("</td>")

    # Type
    parts.append(f'<td class="adata-entry-type" style="{type_style}">')
    parts.append(f'<span class="dtype-extension">{escape_html(type_label)}</span>')
    parts.append("</td>")

    # Meta - custom HTML content (already sanitized by renderer)
    parts.append(f'<td class="adata-entry-meta" style="{meta_style}">')
    # Note: HTML is trusted from registered renderer (package must be imported)
    parts.append(result.html)
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

    parts = [f'<tr class="adata-entry" data-key="{escape_html(key)}" data-dtype="colors">']

    # Name
    parts.append(f'<td class="adata-entry-name" style="{name_style}">')
    parts.append(escape_html(key))
    parts.append(f'<button class="adata-copy-btn" style="{btn_style}" data-copy="{escape_html(key)}" title="Copy name">📋</button>')
    parts.append("</td>")

    # Type with color swatches
    parts.append(f'<td class="adata-entry-type" style="{type_style}">')
    parts.append(f'<span class="dtype-object">colors ({n_colors})</span>')
    parts.append('<span class="adata-color-swatches">')
    for color in colors[:15]:  # Limit preview
        parts.append(f'<span class="adata-color-swatch" style="background:{escape_html(str(color))}" title="{escape_html(str(color))}"></span>')
    if n_colors > 15:
        parts.append(f'<span class="adata-text-muted">+{n_colors - 15}</span>')
    parts.append("</span>")
    parts.append("</td>")

    parts.append(f'<td class="adata-entry-meta" style="{meta_style}"></td>')
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

    parts = [f'<tr class="adata-entry" data-key="{escape_html(key)}" data-dtype="AnnData">']

    # Name
    parts.append(f'<td class="adata-entry-name" style="{name_style}">')
    parts.append(escape_html(key))
    parts.append(f'<button class="adata-copy-btn" style="{btn_style}" data-copy="{escape_html(key)}" title="Copy name">📋</button>')
    parts.append("</td>")

    # Type
    parts.append(f'<td class="adata-entry-type" style="{type_style}">')
    parts.append(f'<span class="dtype-anndata">AnnData ({format_number(n_obs)} × {format_number(n_vars)})</span>')
    if can_expand:
        parts.append(f'<button class="adata-expand-btn" style="{expand_btn_style}" aria-expanded="false">Expand ▼</button>')
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
    parts.append(_render_section_header("raw", f"(n_vars = {format_number(n_vars)})", doc_url, "Raw data (original unprocessed)"))

    # Content with inline styles
    content_style = "padding:0;overflow:hidden;"
    parts.append(f'<div class="anndata-seccontent" style="{content_style}">')

    # Info items with inline styles
    info_style = "padding:6px 12px;font-size:12px;"

    # raw.X info
    if hasattr(raw, "X") and raw.X is not None:
        output = formatter_registry.format_value(raw.X, context)
        parts.append(f'<div style="{info_style}"><strong>raw.X:</strong> <span class="{output.css_class}">{escape_html(output.type_name)}</span></div>')

    # raw.var columns
    if hasattr(raw, "var") and len(raw.var.columns) > 0:
        parts.append(f'<div style="{info_style}"><strong>raw.var:</strong> {len(raw.var.columns)} columns</div>')

    # raw.varm
    if hasattr(raw, "varm") and len(raw.varm) > 0:
        parts.append(f'<div style="{info_style}"><strong>raw.varm:</strong> {len(raw.varm)} items</div>')

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
    parts.append(f'<span class="anndata-sec-name" style="{name_style}">{escape_html(name)}</span>')
    parts.append(f'<span class="anndata-sec-count" style="{count_style}">{escape_html(count_str)}</span>')
    if doc_url:
        parts.append(f'<a class="adata-help-link" style="{link_style}" href="{escape_html(doc_url)}" target="_blank" title="{escape_html(tooltip)}">?</a>')
    parts.append("</div>")
    return "\n".join(parts)


def _render_empty_section(name: str) -> str:
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
    content_style = "padding:0;overflow:hidden;"
    empty_style = "padding:8px 12px;font-size:11px;font-style:italic;"

    return f"""
<div class="anndata-sec" data-section="{escape_html(name)}" data-should-collapse="true">
    <div class="anndata-sechdr" style="{header_style}">
        <span class="adata-fold-icon" style="{fold_style}">▼</span>
        <span class="anndata-sec-name" style="{name_style}">{escape_html(name)}</span>
        <span class="anndata-sec-count" style="{count_style}">(empty)</span>
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
