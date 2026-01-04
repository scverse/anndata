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
from typing import TYPE_CHECKING

from .._repr_constants import (
    STYLE_HIDDEN,
    STYLE_SECTION_CONTENT,
    STYLE_SECTION_TABLE,
)
from . import (
    DEFAULT_FOLD_THRESHOLD,
    DEFAULT_MAX_DEPTH,
    DEFAULT_MAX_FIELD_WIDTH,
    DEFAULT_MAX_ITEMS,
    DEFAULT_PREVIEW_ITEMS,
    DEFAULT_TYPE_WIDTH,
    SECTION_ORDER,
)
from .components import (
    render_badge,
    render_copy_button,
    render_fold_icon,
    render_search_box,
    render_warning_icon,
)
from .css import get_css

# Section renderers (imported here to avoid circular imports at module load time)
# These are used by _render_all_sections and _render_section
from .javascript import get_javascript
from .registry import (
    FormatterContext,
    formatter_registry,
)
from .utils import (
    escape_html,
    format_memory_size,
    format_number,
    get_anndata_version,
    get_backing_info,
    is_backed,
    is_view,
)

if TYPE_CHECKING:
    from typing import Any

    import pandas as pd

    from anndata import AnnData

    from .registry import FormattedEntry, FormattedOutput

# Import formatters to register them (side-effect import)
from . import formatters as _formatters  # noqa: F401

# Approximate character width in pixels for monospace font at 13px
CHAR_WIDTH_PX = 8

# CSS classes in css.py provide full styling with dark mode support.
# These inline styles are minimal fallbacks that ensure basic readability
# if CSS fails to load (e.g., email clients, restrictive embeds).
# Colors and theming are intentionally CSS-only to support dark mode.
# STYLE_HIDDEN, STYLE_SECTION_CONTENT, STYLE_SECTION_TABLE are imported
# from constants.py (single source of truth).

STYLE_SECTION_INFO = "padding:6px 12px;"

# Category color dot - needs inline for dynamic background color
STYLE_CAT_DOT = "width:8px;height:8px;border-radius:50%;display:inline-block;"


def _calculate_field_name_width(adata: AnnData, max_width: int) -> int:
    """
    Calculate the optimal field name column width based on longest field name.

    Collects field names from obs, var, uns, obsm, varm, layers, obsp, varp
    and returns a pixel width that fits the longest name (up to max_width).
    """
    all_names: list[str] = []

    # obs/var column names
    if hasattr(adata, "obs") and adata.obs is not None:
        all_names.extend(adata.obs.columns.tolist())
    if hasattr(adata, "var") and adata.var is not None:
        all_names.extend(adata.var.columns.tolist())

    # Mapping sections (uns, obsm, varm, layers, obsp, varp)
    for attr in ("uns", "obsm", "varm", "layers", "obsp", "varp"):
        try:
            mapping = getattr(adata, attr, None)
            if mapping is not None:
                all_names.extend(mapping.keys())
        except Exception:  # noqa: BLE001
            # Skip sections that fail to access (will show error during rendering)
            pass

    if not all_names:
        return 100  # Minimum default

    # Find longest name
    max_len = max(len(name) for name in all_names)

    # Convert to pixels (with some padding for copy button)
    # Add ~30px for padding and copy button
    width_px = (max_len * CHAR_WIDTH_PX) + 30

    # Clamp to reasonable range
    return max(80, min(width_px, max_width))


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
        max_depth = _get_setting("repr_html_max_depth", default=DEFAULT_MAX_DEPTH)
    if fold_threshold is None:
        fold_threshold = _get_setting(
            "repr_html_fold_threshold", default=DEFAULT_FOLD_THRESHOLD
        )
    if max_items is None:
        max_items = _get_setting("repr_html_max_items", default=DEFAULT_MAX_ITEMS)

    # Check if HTML repr is enabled
    if not _get_setting("repr_html_enabled", default=True):
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

    # Calculate field name column width based on content
    max_field_width = _get_setting(
        "repr_html_max_field_width", default=DEFAULT_MAX_FIELD_WIDTH
    )
    field_width = _calculate_field_name_width(adata, max_field_width)

    # Get type column width from settings
    type_width = _get_setting("repr_html_type_width", default=DEFAULT_TYPE_WIDTH)

    # Container with computed column widths as CSS variables
    style = f"--anndata-name-col-width: {field_width}px; --anndata-type-col-width: {type_width}px;"
    parts.append(
        f'<div class="anndata-repr" id="{container_id}" data-depth="{depth}" style="{style}">'
    )

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
    parts.append(_render_x_entry(adata, context))
    parts.extend(
        _render_all_sections(adata, context, fold_threshold, max_items, max_depth)
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


def _render_all_sections(
    adata: AnnData,
    context: FormatterContext,
    fold_threshold: int,
    max_items: int,
    max_depth: int,
) -> list[str]:
    """Render all standard and custom sections."""
    parts = []
    custom_sections_after = _get_custom_sections_by_position(adata)

    for section in SECTION_ORDER:
        if section == "X":
            # X is already rendered, but check for custom sections after X
            if "X" in custom_sections_after:
                parts.extend(
                    _render_custom_section(
                        adata, section_formatter, context, fold_threshold, max_items
                    )
                    for section_formatter in custom_sections_after["X"]
                )
            continue
        parts.append(
            _render_section(
                adata,
                section,
                context,
                fold_threshold=fold_threshold,
                max_items=max_items,
                max_depth=max_depth,
            )
        )

        # Render custom sections after this section
        if section in custom_sections_after:
            parts.extend(
                _render_custom_section(
                    adata, section_formatter, context, fold_threshold, max_items
                )
                for section_formatter in custom_sections_after[section]
            )

    # Custom sections at end (no specific position)
    if None in custom_sections_after:
        parts.extend(
            _render_custom_section(
                adata, section_formatter, context, fold_threshold, max_items
            )
            for section_formatter in custom_sections_after[None]
        )

    # Detect and show unknown sections (mapping-like attributes not in SECTION_ORDER)
    unknown_sections = _detect_unknown_sections(adata)
    if unknown_sections:
        parts.append(_render_unknown_sections(unknown_sections))

    return parts


def _render_section(
    adata: AnnData,
    section: str,
    context: FormatterContext,
    *,
    fold_threshold: int,
    max_items: int,
    max_depth: int,
) -> str:
    """Render a single standard section."""
    try:
        if section == "X":
            return _render_x_entry(adata, context)
        if section == "raw":
            return _render_raw_section(adata, context, fold_threshold, max_items)
        if section in ("obs", "var"):
            return _render_dataframe_section(
                adata, section, context, fold_threshold, max_items
            )
        if section == "uns":
            return _render_uns_section(
                adata, context, fold_threshold, max_items, max_depth
            )
        return _render_mapping_section(
            adata, section, context, fold_threshold, max_items
        )
    except Exception as e:  # noqa: BLE001
        # Show error instead of hiding the section
        return _render_error_entry(section, str(e))


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
        from .._warnings import warn

        warn(
            f"Custom section formatter '{formatter.section_name}' failed: {e}",
            UserWarning,
        )
        return ""

    if not entries:
        return ""

    n_items = len(entries)
    section_name = formatter.section_name

    # Render entries (with truncation)
    rows = []
    for i, entry in enumerate(entries):
        if i >= max_items:
            rows.append(_render_truncation_indicator(n_items - max_items))
            break
        rows.append(render_formatted_entry(entry, section_name))

    # Use render_section for consistent structure
    return render_section(
        getattr(formatter, "display_name", section_name),
        "\n".join(rows),
        n_items=n_items,
        doc_url=getattr(formatter, "doc_url", None),
        tooltip=getattr(formatter, "tooltip", ""),
        should_collapse=n_items > fold_threshold,
        section_id=section_name,
    )


def _render_entry_meta_content(
    output: FormattedOutput, section: str, *, has_columns_list: bool
) -> str:
    """Render meta column content for an entry.

    Priority: meta_content > columns list > meta_preview > shape
    """
    if output.meta_content:
        return output.meta_content
    if has_columns_list and "columns" in output.details:
        columns = output.details["columns"]
        col_str = ", ".join(escape_html(str(c)) for c in columns)
        return f'<span class="adata-cols-list">[{col_str}]</span>'
    if "meta_preview" in output.details:
        full_preview = output.details.get(
            "meta_preview_full", output.details["meta_preview"]
        )
        return f'<span title="{escape_html(full_preview)}">{escape_html(output.details["meta_preview"])}</span>'
    if "shape" in output.details and section in ("obsm", "varm"):
        shape = output.details["shape"]
        if len(shape) >= 2:
            return f"({format_number(shape[1])} cols)"
    return ""


def render_formatted_entry(
    entry: FormattedEntry,
    section: str = "",
    *,
    extra_warnings: list[str] | None = None,
    is_hard_error: bool = False,
) -> str:
    """
    Render a FormattedEntry as a table row.

    This is a public API for packages building their own _repr_html_.
    It provides the same flexibility as internal code by accepting
    FormattedEntry/FormattedOutput objects.

    Parameters
    ----------
    entry
        A FormattedEntry containing the key and FormattedOutput
    section
        Optional section name (used for meta column rendering)
    extra_warnings
        Additional warnings to display (e.g., key validation warnings)
    is_hard_error
        Whether there's a hard error (e.g., invalid key name) in addition
        to serializability issues

    Returns
    -------
    HTML string for the table row(s)

    Examples
    --------
    ::

        from . import (
            FormattedEntry,
            FormattedOutput,
            render_formatted_entry,
        )

        entry = FormattedEntry(
            key="my_array",
            output=FormattedOutput(
                type_name="ndarray (100, 50) float32",
                css_class="dtype-ndarray",
                tooltip="My custom array",
                warnings=["Some warning"],
            ),
        )
        html = render_formatted_entry(entry)

    With expandable nested content::

        nested_html = generate_repr_html(adata, depth=1)
        entry = FormattedEntry(
            key="cell_table",
            output=FormattedOutput(
                type_name="AnnData (150 × 30)",
                css_class="dtype-anndata",
                html_content=nested_html,
                is_expandable=True,
            ),
        )
        html = render_formatted_entry(entry)

    With key validation warnings::

        entry = FormattedEntry(
            key="bad/key",
            output=FormattedOutput(...),
        )
        html = render_formatted_entry(
            entry, extra_warnings=["Contains '/' (deprecated)"]
        )
    """
    output = entry.output
    extra_warnings = extra_warnings or []

    # Compute entry CSS classes
    all_warnings = extra_warnings + list(output.warnings)
    has_error = not output.is_serializable or is_hard_error

    entry_class = "adata-entry"
    if all_warnings:
        entry_class += " warning"
    if has_error:
        entry_class += " error"

    has_expandable_content = output.html_content and output.is_expandable

    parts = [
        f'<tr class="{entry_class}" data-key="{escape_html(entry.key)}" '
        f'data-dtype="{escape_html(output.type_name)}">'
    ]

    # Name
    parts.append(_render_name_cell(entry.key))

    # Type cell
    parts.append('<td class="adata-entry-type">')

    # Type span with optional tooltip
    if output.tooltip:
        parts.append(
            f'<span class="{output.css_class}" title="{escape_html(output.tooltip)}">'
            f"{escape_html(output.type_name)}</span>"
        )
    else:
        parts.append(
            f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>'
        )

    # Warning icon
    parts.append(render_warning_icon(all_warnings, is_error=has_error))

    # Expand button for expandable content
    if has_expandable_content:
        parts.append(
            f'<button class="adata-expand-btn" style="{STYLE_HIDDEN}" aria-expanded="false">Expand ▼</button>'
        )

    # Columns list toggle button (for DataFrames in obsm/varm)
    has_columns_list = output.details.get("has_columns_list", False)
    if has_columns_list:
        parts.append(
            '<button class="adata-cols-wrap-btn" title="Toggle multi-line view">⋯</button>'
        )

    # Inline custom content (non-expandable)
    if output.html_content and not output.is_expandable:
        parts.append(
            f'<div class="adata-custom-content" style="margin-top:4px;">{output.html_content}</div>'
        )

    parts.append("</td>")

    # Meta column (for data previews, dimensions, etc.)
    parts.append('<td class="adata-entry-meta">')
    parts.append(
        _render_entry_meta_content(output, section, has_columns_list=has_columns_list)
    )
    parts.append("</td>")

    parts.append("</tr>")

    # Expandable content row
    if has_expandable_content:
        parts.append('<tr class="adata-nested-row">')
        parts.append('<td colspan="3" class="adata-nested-content">')
        parts.append(f'<div class="adata-custom-expanded">{output.html_content}</div>')
        parts.append("</td>")
        parts.append("</tr>")

    return "\n".join(parts)


def _render_name_cell(name: str) -> str:
    """Render a name cell with copy button and tooltip for truncated names.

    The structure uses flexbox so the copy button stays visible even when
    the name text overflows and shows ellipsis.
    """
    escaped_name = escape_html(name)
    return (
        f'<td class="adata-entry-name">'
        f'<div class="adata-entry-name-inner">'
        f'<span class="adata-name-text" title="{escaped_name}">{escaped_name}</span>'
        f"{render_copy_button(name, 'Copy name')}"
        f"</div>"
        f"</td>"
    )


def _render_header(
    adata: AnnData, *, show_search: bool = False, container_id: str = ""
) -> str:
    """Render the header with type, shape, badges, and optional search box."""
    parts = ['<div class="anndata-hdr">']

    # Type name - allow for extension types
    type_name = type(adata).__name__
    parts.append(f'<span class="adata-type">{escape_html(type_name)}</span>')

    # Shape
    shape_str = f"{format_number(adata.n_obs)} obs × {format_number(adata.n_vars)} vars"
    parts.append(f'<span class="adata-shape">{shape_str}</span>')

    # Badges - use render_badge() helper
    if is_view(adata):
        parts.append(render_badge("View", "adata-badge-view"))

    if is_backed(adata):
        backing = get_backing_info(adata)
        filename = backing.get("filename", "")
        format_str = backing.get("format", "")
        status = "Open" if backing.get("is_open") else "Closed"
        parts.append(render_badge(f"{format_str} ({status})", "adata-badge-backed"))
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
        parts.append(render_badge(type_name, "adata-badge-extension"))

    # README icon if uns["README"] exists with a string
    readme_content = adata.uns.get("README") if hasattr(adata, "uns") else None
    if isinstance(readme_content, str) and readme_content.strip():
        escaped_readme = escape_html(readme_content)
        # Truncate for no-JS tooltip (first 500 chars)
        tooltip_text = readme_content[:500]
        if len(readme_content) > 500:
            tooltip_text += "..."
        escaped_tooltip = escape_html(tooltip_text)

        parts.append(
            f'<span class="adata-readme-icon" '
            f'data-readme="{escaped_readme}" '
            f'title="{escaped_tooltip}" '
            f'role="button" tabindex="0" aria-label="View README">'
            f"ⓘ"
            f"</span>"
        )

    # Search box on the right (spacer pushes it right) - use render_search_box() helper
    if show_search:
        parts.append('<span style="flex-grow:1;"></span>')
        parts.append(render_search_box(container_id))

    parts.append("</div>")
    return "\n".join(parts)


def _render_footer(adata: AnnData) -> str:
    """Render the footer with version and memory info."""
    parts = ['<div class="anndata-ftr">']

    # Version
    version = get_anndata_version()
    parts.append(f"<span>anndata v{version}</span>")

    # Memory usage
    try:
        mem_bytes = adata.__sizeof__()
        mem_str = format_memory_size(mem_bytes)
        parts.append(f'<span title="Estimated memory usage">~{mem_str}</span>')
    except Exception:  # noqa: BLE001
        # Broad catch: __sizeof__ recursively calls into user data which could raise anything
        pass

    parts.append("</div>")
    return "\n".join(parts)


def _render_index_preview(adata: AnnData) -> str:
    """Render preview of obs_names and var_names."""
    parts = ['<div class="adata-index-preview">']

    # obs_names preview
    obs_preview = _format_index_preview(adata.obs_names)
    parts.append(f"<div><strong>obs_names:</strong> {obs_preview}</div>")

    # var_names preview
    var_preview = _format_index_preview(adata.var_names)
    parts.append(f"<div><strong>var_names:</strong> {var_preview}</div>")

    parts.append("</div>")
    return "\n".join(parts)


def _format_index_preview(index: pd.Index) -> str:
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


def _render_x_entry(obj: Any, context: FormatterContext) -> str:
    """Render X as a single compact entry row.

    Works with both AnnData and Raw objects.
    """
    X = obj.X

    parts = ['<div class="adata-x-entry">']
    parts.append("<span>X</span>")

    if X is None:
        parts.append("<span><em>None</em></span>")
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

        # Backed info (only for AnnData, not Raw)
        if is_backed(obj):
            type_parts.append("on disk")

        type_str = " · ".join(type_parts)
        parts.append(f'<span class="{output.css_class}">{escape_html(type_str)}</span>')

    parts.append("</div>")
    return "\n".join(parts)


# Section renderers are imported from sections.py
from .sections import (  # noqa: E402
    _detect_unknown_sections,
    _render_dataframe_section,
    _render_error_entry,
    _render_mapping_section,
    _render_raw_section,
    _render_unknown_sections,
    _render_uns_section,
)


def _render_section_header(
    name: str,
    count_str: str,
    doc_url: str | None,
    tooltip: str,
) -> str:
    """Render a section header - colors handled by CSS for dark mode support."""
    parts = ['<div class="anndata-sechdr">']
    parts.append(render_fold_icon())  # Use helper for fold icon
    parts.append(f'<span class="anndata-sec-name">{escape_html(name)}</span>')
    parts.append(f'<span class="anndata-sec-count">{escape_html(count_str)}</span>')
    if doc_url:
        parts.append(
            f'<a class="adata-help-link"  href="{escape_html(doc_url)}" target="_blank" title="{escape_html(tooltip)}">?</a>'
        )
    parts.append("</div>")
    return "\n".join(parts)


def _render_empty_section(
    name: str,
    doc_url: str | None = None,
    tooltip: str = "",
) -> str:
    """Render an empty section indicator."""
    # Build help link if doc_url provided
    help_link = ""
    if doc_url:
        help_link = f'<a class="adata-help-link"  href="{escape_html(doc_url)}" target="_blank" title="{escape_html(tooltip)}">?</a>'

    # Use render_fold_icon() helper for consistency
    fold_icon = render_fold_icon()

    return f"""
<div class="anndata-sec" data-section="{escape_html(name)}" data-should-collapse="true">
    <div class="anndata-sechdr">
        {fold_icon}
        <span class="anndata-sec-name">{escape_html(name)}</span>
        <span class="anndata-sec-count">(empty)</span>
        {help_link}
    </div>
    <div class="anndata-seccontent" style="{STYLE_SECTION_CONTENT}">
        <div class="adata-empty">No entries</div>
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


def _get_setting(name: str, *, default: Any) -> Any:
    """Get a setting value, falling back to default."""
    try:
        from anndata import settings

        return getattr(settings, name, default)
    except (ImportError, AttributeError):
        return default


def render_section(  # noqa: PLR0913
    name: str,
    entries_html: str,
    *,
    n_items: int,
    doc_url: str | None = None,
    tooltip: str = "",
    should_collapse: bool = False,
    section_id: str | None = None,
    count_str: str | None = None,
) -> str:
    """
    Render a complete section with header and content.

    This is a public API for packages building their own _repr_html_.
    It is also used internally for consistency.

    Parameters
    ----------
    name
        Display name for the section header (e.g., 'images', 'tables')
    entries_html
        HTML content for the section body (table rows)
    n_items
        Number of items (used for empty check and default count string)
    doc_url
        URL for the help link (? icon)
    tooltip
        Tooltip text for the help link
    should_collapse
        Whether this section should start collapsed
    section_id
        ID for the section in data-section attribute (defaults to name)
    count_str
        Custom count string for header (defaults to "(N items)")

    Returns
    -------
    HTML string for the complete section

    Examples
    --------
    ::

        from . import (
            FormattedEntry,
            FormattedOutput,
            render_formatted_entry,
        )

        rows = []
        for key, info in items.items():
            entry = FormattedEntry(
                key=key,
                output=FormattedOutput(
                    type_name=info["type"], css_class="dtype-ndarray"
                ),
            )
            rows.append(render_formatted_entry(entry))

        html = render_section(
            "images",
            "\\n".join(rows),
            n_items=len(items),
            doc_url="https://docs.example.com/images",
            tooltip="Image data",
        )
    """
    if section_id is None:
        section_id = name

    if n_items == 0:
        return _render_empty_section(name, doc_url, tooltip)

    if count_str is None:
        count_str = f"({n_items} items)"

    parts = [
        f'<div class="anndata-sec" data-section="{escape_html(section_id)}" '
        f'data-should-collapse="{str(should_collapse).lower()}">'
    ]

    # Header
    parts.append(_render_section_header(name, count_str, doc_url, tooltip))

    # Content
    parts.append(f'<div class="anndata-seccontent" style="{STYLE_SECTION_CONTENT}">')
    parts.append(f'<table class="adata-table" style="{STYLE_SECTION_TABLE}">')
    parts.append(entries_html)
    parts.append("</table></div></div>")

    return "\n".join(parts)
