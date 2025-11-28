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
    DEFAULT_MAX_DEPTH,
    DEFAULT_MAX_ITEMS,
    DEFAULT_PREVIEW_ITEMS,
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

    # Header
    if show_header:
        parts.append(_render_header(adata))

    # Search box (only at top level)
    if show_search and depth == 0:
        parts.append(_render_search_box())

    # Metadata bar
    if depth == 0:
        parts.append(_render_metadata(adata))

    # Index preview (only at top level)
    if depth == 0:
        parts.append(_render_index_preview(adata))

    # Sections container
    parts.append('<div class="ad-sections">')

    # X section
    parts.append(_render_x_section(adata, context))

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

    parts.append("</div>")  # ad-sections
    parts.append("</div>")  # anndata-repr

    # JavaScript (only at top level)
    if depth == 0:
        parts.append(get_javascript(container_id))

    return "\n".join(parts)


# =============================================================================
# Section Renderers
# =============================================================================


def _render_header(adata: AnnData) -> str:
    """Render the header with type, shape, and badges."""
    parts = ['<div class="ad-header">']

    # Type name - allow for extension types
    type_name = type(adata).__name__
    parts.append(f'<span class="ad-type">{escape_html(type_name)}</span>')

    # Shape
    shape_str = f"n_obs × n_vars = {format_number(adata.n_obs)} × {format_number(adata.n_vars)}"
    parts.append(f'<span class="ad-shape">{shape_str}</span>')

    # Badges
    if is_view(adata):
        parts.append('<span class="ad-badge ad-badge-view">View</span>')

    if is_backed(adata):
        backing = get_backing_info(adata)
        title = escape_html(backing.get("filename", ""))
        format_str = backing.get("format", "")
        status = "Open" if backing.get("is_open") else "Closed"
        parts.append(
            f'<span class="ad-badge ad-badge-backed" title="{title}">'
            f'📁 {format_str} ({status})</span>'
        )

    # Check for extension type (not standard AnnData)
    if type_name != "AnnData":
        parts.append(f'<span class="ad-badge ad-badge-extension">{type_name}</span>')

    parts.append("</div>")
    return "\n".join(parts)


def _render_search_box() -> str:
    """Render the search/filter input."""
    return """
<div class="ad-search">
    <input type="text" class="ad-search-input" placeholder="Search fields..." aria-label="Search fields">
    <span class="ad-filter-indicator"></span>
</div>
"""


def _render_metadata(adata: AnnData) -> str:
    """Render the metadata bar with version and memory info."""
    parts = ['<div class="ad-metadata">']

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

    # Creation date if available
    if hasattr(adata, "uns") and "created_date" in adata.uns:
        parts.append(f'<span>Created: {escape_html(str(adata.uns["created_date"]))}</span>')

    parts.append("</div>")
    return "\n".join(parts)


def _render_index_preview(adata: AnnData) -> str:
    """Render preview of obs_names and var_names."""
    parts = ['<div class="ad-index-preview">']

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


def _render_x_section(adata: AnnData, context: FormatterContext) -> str:
    """Render the X matrix section."""
    parts = ['<div class="ad-x-section">']
    parts.append('<dl class="ad-x-info">')

    X = adata.X

    if X is None:
        parts.append("<dt>X</dt><dd><em>None</em></dd>")
    else:
        # Format the X matrix
        output = formatter_registry.format_value(X, context)

        parts.append(f"<dt>Type:</dt><dd><span class=\"{output.css_class}\">{escape_html(output.type_name)}</span></dd>")

        # Shape
        if "shape" in output.details:
            shape = output.details["shape"]
            shape_str = " × ".join(format_number(s) for s in shape)
            parts.append(f"<dt>Shape:</dt><dd>{shape_str}</dd>")

        # Dtype
        if "dtype" in output.details:
            parts.append(f"<dt>Dtype:</dt><dd>{escape_html(str(output.details['dtype']))}</dd>")

        # Sparsity for sparse matrices
        if "sparsity" in output.details and output.details["sparsity"] is not None:
            sparsity = output.details["sparsity"]
            nnz = output.details.get("nnz", "?")
            parts.append(f"<dt>Sparsity:</dt><dd>{sparsity:.1%} sparse ({format_number(nnz)} stored)</dd>")

        # Chunk info for Dask
        if "chunks" in output.details:
            chunks = output.details["chunks"]
            parts.append(f"<dt>Chunks:</dt><dd>{chunks}</dd>")

        # Backed info
        if is_backed(adata):
            parts.append("<dt>Storage:</dt><dd>📁 On disk</dd>")

    parts.append("</dl>")
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

    # Should this section be collapsed?
    collapsed = n_cols > fold_threshold

    parts = [f'<div class="ad-section{"" if not collapsed else " collapsed"}" data-section="{section}">']

    # Header
    doc_url = f"{DOCS_BASE_URL}generated/anndata.AnnData.{section}.html"
    tooltip = "Observation annotations" if section == "obs" else "Variable annotations"
    parts.append(_render_section_header(section, f"({n_cols} columns)", doc_url, tooltip))

    # Content
    parts.append('<div class="ad-section-content">')
    parts.append('<table class="ad-table">')

    # Render each column
    for i, col_name in enumerate(df.columns):
        if i >= max_items:
            parts.append(_render_truncation_indicator(n_cols - max_items))
            break

        col = df[col_name]
        parts.append(_render_dataframe_entry(adata, section, col_name, col, context))

    parts.append("</table>")
    parts.append("</div>")  # ad-section-content
    parts.append("</div>")  # ad-section

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

    # Build entry class
    entry_class = "ad-entry"
    if warnings:
        entry_class += " warning"

    # Build row
    parts = [f'<tr class="{entry_class}" data-key="{escape_html(col_name)}" data-dtype="{escape_html(output.type_name)}">']

    # Name cell
    parts.append('<td class="ad-entry-name">')
    parts.append(escape_html(col_name))
    parts.append(f'<button class="ad-copy-btn" data-copy="{escape_html(col_name)}" title="Copy name">📋</button>')
    parts.append("</td>")

    # Type cell
    parts.append('<td class="ad-entry-type">')
    if warnings:
        title = escape_html("; ".join(warnings))
        parts.append(f'<span class="{output.css_class} dtype-warning" title="{title}">')
        parts.append(f"{escape_html(output.type_name)} ⚠️")
        parts.append("</span>")
    else:
        parts.append(f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>')

    # Color swatches
    if colors:
        parts.append('<span class="ad-color-swatches">')
        for color in colors[:10]:  # Limit to 10 swatches
            parts.append(f'<span class="ad-color-swatch" style="background:{escape_html(color)}" title="{escape_html(color)}"></span>')
        if len(colors) > 10:
            parts.append(f"<span>+{len(colors) - 10}</span>")
        parts.append("</span>")

    parts.append("</td>")

    # Meta cell
    parts.append('<td class="ad-entry-meta">')
    if hasattr(col, "cat"):
        parts.append(f"({len(col.cat.categories)} categories)")
    elif hasattr(col, "nunique"):
        try:
            parts.append(f"({col.nunique()} unique)")
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

    collapsed = n_items > fold_threshold

    parts = [f'<div class="ad-section{"" if not collapsed else " collapsed"}" data-section="{section}">']

    # Header
    doc_url = f"{DOCS_BASE_URL}generated/anndata.AnnData.{section}.html"
    tooltip = _get_section_tooltip(section)
    parts.append(_render_section_header(section, f"({n_items} items)", doc_url, tooltip))

    # Content
    parts.append('<div class="ad-section-content">')
    parts.append('<table class="ad-table">')

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

    entry_class = "ad-entry"
    if output.warnings:
        entry_class += " warning"
    if not output.is_serializable:
        entry_class += " error"

    parts = [f'<tr class="{entry_class}" data-key="{escape_html(key)}" data-dtype="{escape_html(output.type_name)}">']

    # Name
    parts.append('<td class="ad-entry-name">')
    parts.append(escape_html(key))
    parts.append(f'<button class="ad-copy-btn" data-copy="{escape_html(key)}" title="Copy name">📋</button>')
    parts.append("</td>")

    # Type
    parts.append('<td class="ad-entry-type">')
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

    # Meta - show shape/cols for obsm/varm
    parts.append('<td class="ad-entry-meta">')
    if "shape" in output.details and section in ("obsm", "varm", "layers"):
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

    collapsed = n_items > fold_threshold

    parts = [f'<div class="ad-section{"" if not collapsed else " collapsed"}" data-section="uns">']

    # Header
    doc_url = f"{DOCS_BASE_URL}generated/anndata.AnnData.uns.html"
    parts.append(_render_section_header("uns", f"({n_items} items)", doc_url, "Unstructured annotation"))

    # Content
    parts.append('<div class="ad-section-content">')
    parts.append('<table class="ad-table">')

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
    """Render a single uns entry with special type handling."""
    parts = []

    # Check for color list
    if is_color_list(key, value):
        return _render_color_list_entry(key, value)

    # Check for nested AnnData
    if type(value).__name__ == "AnnData" and hasattr(value, "n_obs"):
        return _render_nested_anndata_entry(key, value, context, max_depth)

    # Regular entry
    output = formatter_registry.format_value(value, context)

    entry_class = "ad-entry"
    if output.warnings:
        entry_class += " warning"
    if not output.is_serializable:
        entry_class += " error"

    parts.append(f'<tr class="{entry_class}" data-key="{escape_html(key)}" data-dtype="{escape_html(output.type_name)}">')

    # Name
    parts.append('<td class="ad-entry-name">')
    parts.append(escape_html(key))
    parts.append(f'<button class="ad-copy-btn" data-copy="{escape_html(key)}" title="Copy name">📋</button>')
    parts.append("</td>")

    # Type
    parts.append('<td class="ad-entry-type">')
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

    # Meta
    parts.append('<td class="ad-entry-meta"></td>')
    parts.append("</tr>")

    return "\n".join(parts)


def _render_color_list_entry(key: str, value: Any) -> str:
    """Render a color list entry with swatches."""
    colors = list(value) if hasattr(value, "__iter__") else []
    n_colors = len(colors)

    parts = [f'<tr class="ad-entry" data-key="{escape_html(key)}" data-dtype="colors">']

    # Name
    parts.append('<td class="ad-entry-name">')
    parts.append(escape_html(key))
    parts.append(f'<button class="ad-copy-btn" data-copy="{escape_html(key)}" title="Copy name">📋</button>')
    parts.append("</td>")

    # Type with color swatches
    parts.append('<td class="ad-entry-type">')
    parts.append(f'<span class="dtype-object">colors ({n_colors})</span>')
    parts.append('<span class="ad-color-swatches">')
    for color in colors[:15]:  # Limit preview
        parts.append(f'<span class="ad-color-swatch" style="background:{escape_html(str(color))}" title="{escape_html(str(color))}"></span>')
    if n_colors > 15:
        parts.append(f"<span>+{n_colors - 15}</span>")
    parts.append("</span>")
    parts.append("</td>")

    parts.append('<td class="ad-entry-meta"></td>')
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

    parts = [f'<tr class="ad-entry" data-key="{escape_html(key)}" data-dtype="AnnData">']

    # Name
    parts.append('<td class="ad-entry-name">')
    parts.append(escape_html(key))
    parts.append(f'<button class="ad-copy-btn" data-copy="{escape_html(key)}" title="Copy name">📋</button>')
    parts.append("</td>")

    # Type
    parts.append('<td class="ad-entry-type">')
    parts.append(f'<span class="dtype-anndata">AnnData ({format_number(n_obs)} × {format_number(n_vars)})</span>')
    if can_expand:
        parts.append('<button class="ad-expand-btn" aria-expanded="false">Expand ▼</button>')
    parts.append("</td>")

    parts.append('<td class="ad-entry-meta"></td>')
    parts.append("</tr>")

    # Nested content (hidden by default)
    if can_expand:
        parts.append('<tr class="ad-nested-row">')
        parts.append('<td colspan="3" class="ad-nested-content">')
        parts.append('<div class="ad-nested-anndata">')
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

    parts = ['<div class="ad-section collapsed" data-section="raw">']

    # Header
    doc_url = f"{DOCS_BASE_URL}generated/anndata.AnnData.raw.html"
    n_vars = getattr(raw, "n_vars", "?")
    parts.append(_render_section_header("raw", f"(n_vars = {format_number(n_vars)})", doc_url, "Raw data (original unprocessed)"))

    # Content
    parts.append('<div class="ad-section-content">')

    # raw.X info
    if hasattr(raw, "X") and raw.X is not None:
        output = formatter_registry.format_value(raw.X, context)
        parts.append(f'<div class="ad-raw-info"><strong>raw.X:</strong> <span class="{output.css_class}">{escape_html(output.type_name)}</span></div>')

    # raw.var columns
    if hasattr(raw, "var") and len(raw.var.columns) > 0:
        parts.append(f'<div class="ad-raw-info"><strong>raw.var:</strong> {len(raw.var.columns)} columns</div>')

    # raw.varm
    if hasattr(raw, "varm") and len(raw.varm) > 0:
        parts.append(f'<div class="ad-raw-info"><strong>raw.varm:</strong> {len(raw.varm)} items</div>')

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
    """Render a section header."""
    parts = ['<div class="ad-section-header">']
    parts.append('<span class="ad-fold-icon">▼</span>')
    parts.append(f'<span class="ad-section-name">{escape_html(name)}</span>')
    parts.append(f'<span class="ad-section-count">{escape_html(count_str)}</span>')
    if doc_url:
        parts.append(f'<a class="ad-help-link" href="{escape_html(doc_url)}" target="_blank" title="{escape_html(tooltip)}">?</a>')
    parts.append("</div>")
    return "\n".join(parts)


def _render_empty_section(name: str) -> str:
    """Render an empty section indicator."""
    return f"""
<div class="ad-section collapsed" data-section="{escape_html(name)}">
    <div class="ad-section-header">
        <span class="ad-fold-icon">▼</span>
        <span class="ad-section-name">{escape_html(name)}</span>
        <span class="ad-section-count">(empty)</span>
    </div>
    <div class="ad-section-content">
        <div class="ad-empty">No entries</div>
    </div>
</div>
"""


def _render_truncation_indicator(remaining: int) -> str:
    """Render a truncation indicator."""
    return f'<tr><td colspan="3" class="ad-truncated">... and {format_number(remaining)} more</td></tr>'


def _render_max_depth_indicator(adata: AnnData) -> str:
    """Render indicator when max depth is reached."""
    n_obs = getattr(adata, "n_obs", "?")
    n_vars = getattr(adata, "n_vars", "?")
    return f'<div class="ad-max-depth">AnnData ({format_number(n_obs)} × {format_number(n_vars)}) - max depth reached</div>'


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
