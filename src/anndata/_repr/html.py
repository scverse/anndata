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

import numpy as np

from anndata._repr import (
    DEFAULT_FOLD_THRESHOLD,
    DEFAULT_MAX_CATEGORIES,
    DEFAULT_MAX_DEPTH,
    DEFAULT_MAX_FIELD_WIDTH,
    DEFAULT_MAX_ITEMS,
    DEFAULT_MAX_STRING_LENGTH,
    DEFAULT_PREVIEW_ITEMS,
    DEFAULT_TYPE_WIDTH,
    DEFAULT_UNIQUE_LIMIT,
    DOCS_BASE_URL,
    SECTION_ORDER,
)
from anndata._repr.constants import (
    NOT_SERIALIZABLE_MSG,
    STYLE_HIDDEN,
    STYLE_SECTION_CONTENT,
    STYLE_SECTION_TABLE,
)
from anndata._repr.css import get_css
from anndata._repr.javascript import get_javascript
from anndata._repr.registry import (
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
    is_lazy_series,
    is_view,
    should_warn_string_column,
)

if TYPE_CHECKING:
    from typing import Any

    import pandas as pd

    from anndata import AnnData
    from anndata._repr.registry import FormattedEntry, FormattedOutput

# Import formatters to register them (side-effect import)
import anndata._repr.formatters  # noqa: F401
from anndata._repr.formatters import check_column_name

# Approximate character width in pixels for monospace font at 13px
CHAR_WIDTH_PX = 8

# =============================================================================
# Inline styles for graceful degradation
# =============================================================================
# CSS classes in css.py provide full styling with dark mode support.
# These inline styles are minimal fallbacks that ensure basic readability
# if CSS fails to load (e.g., email clients, restrictive embeds).
# Colors and theming are intentionally CSS-only to support dark mode.
# STYLE_HIDDEN, STYLE_SECTION_CONTENT, STYLE_SECTION_TABLE are imported
# from constants.py (single source of truth).

STYLE_SECTION_INFO = "padding:6px 12px;"

# Category color dot - needs inline for dynamic background color
STYLE_CAT_DOT = "width:8px;height:8px;border-radius:50%;display:inline-block;"


def render_warning_icon(warnings: list[str], is_error: bool = False) -> str:
    """Render warning icon with tooltip if there are warnings or errors.

    Parameters
    ----------
    warnings
        List of warning messages to show in tooltip.
    is_error
        If True, prepends "Not serializable to H5AD/Zarr" to warnings.

    Returns
    -------
    HTML string for warning icon, or empty string if no warnings/error.
    """
    if not warnings and not is_error:
        return ""
    all_warnings = list(warnings)
    if is_error:
        all_warnings.insert(0, NOT_SERIALIZABLE_MSG)
    title = escape_html("; ".join(all_warnings))
    return f'<span class="adata-warning-icon" title="{title}">(!)</span>'


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


# =============================================================================
# Custom Section Support
# =============================================================================


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
        from anndata._warnings import warn

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


def render_formatted_entry(entry: FormattedEntry, section: str = "") -> str:
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
        Optional section name (for future use)

    Returns
    -------
    HTML string for the table row(s)

    Examples
    --------
    ::

        from anndata._repr import (
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
    """
    output = entry.output

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
    parts.append(_render_name_cell(entry.key))

    # Type
    parts.append('<td class="adata-entry-type">')
    parts.append(
        f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>'
    )
    parts.append(render_warning_icon(list(output.warnings), not output.is_serializable))

    if has_expandable_content:
        parts.append(
            f'<button class="adata-expand-btn" style="{STYLE_HIDDEN}" aria-expanded="false">Expand ▼</button>'
        )

    if output.html_content and not output.is_expandable:
        parts.append(
            f'<div class="adata-custom-content" style="margin-top:4px;">{output.html_content}</div>'
        )

    parts.append("</td>")

    # Meta column (for data previews, dimensions, etc.)
    parts.append('<td class="adata-entry-meta">')
    if output.meta_content:
        parts.append(output.meta_content)
    parts.append("</td>")

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
# UI Component Helpers (for external packages)
# =============================================================================
# These functions generate HTML for common interactive UI elements.
# External packages (SpatialData, MuData, etc.) can use these instead of
# hardcoding CSS classes and inline styles.


def render_search_box(container_id: str = "") -> str:
    """
    Render a search box with filter indicator and search mode toggles.

    The search box is hidden by default and shown when JavaScript is enabled.
    It filters entries across all sections by key, type, or content.
    Includes toggle buttons for case-sensitive search and regex mode.

    Parameters
    ----------
    container_id
        Unique ID for the container (used for label association)

    Returns
    -------
    HTML string for the search box

    Example
    -------
    >>> parts = ['<div class="anndata-hdr">']
    >>> parts.append('<span class="adata-type">SpatialData</span>')
    >>> parts.append('<span style="flex-grow:1;"></span>')  # Spacer
    >>> parts.append(render_search_box(container_id))
    >>> parts.append("</div>")
    """
    search_id = f"{container_id}-search" if container_id else "anndata-search"
    return (
        f'<span class="adata-search-box" style="{STYLE_HIDDEN}">'
        f'<input type="text" id="{search_id}" name="{search_id}" '
        f'class="adata-search-input" '
        f'placeholder="Search..." aria-label="Search fields">'
        f'<span class="adata-search-toggles">'
        f'<button type="button" class="adata-search-toggle adata-toggle-case" '
        f'title="Match case" aria-label="Match case" aria-pressed="false">Aa</button>'
        f'<button type="button" class="adata-search-toggle adata-toggle-regex" '
        f'title="Use regular expression" aria-label="Use regular expression" aria-pressed="false">.*</button>'
        f"</span>"
        f"</span>"
        f'<span class="adata-filter-indicator"></span>'
    )


def render_fold_icon() -> str:
    """
    Render a fold/expand icon for section headers.

    The icon is hidden by default and shown when JavaScript is enabled.
    It rotates when the section is collapsed.

    Returns
    -------
    HTML string for the fold icon

    Example
    -------
    >>> parts = ['<div class="anndata-sechdr">']
    >>> parts.append(render_fold_icon())
    >>> parts.append('<span class="anndata-sec-name">images</span>')
    >>> parts.append("</div>")
    """
    return f'<span class="adata-fold-icon" style="{STYLE_HIDDEN}">▼</span>'


def render_copy_button(text: str, tooltip: str = "Copy") -> str:
    """
    Render a copy-to-clipboard button.

    The button is hidden by default and shown when JavaScript is enabled.
    When clicked, it copies the specified text to the clipboard.

    Parameters
    ----------
    text
        The text to copy when clicked
    tooltip
        Tooltip text (default: "Copy")

    Returns
    -------
    HTML string for the copy button

    Example
    -------
    >>> html = f"<span>{name}</span>{render_copy_button(name, 'Copy name')}"
    """
    escaped_text = escape_html(text)
    escaped_tooltip = escape_html(tooltip)
    return (
        f'<button class="adata-copy-btn" style="{STYLE_HIDDEN}" '
        f'data-copy="{escaped_text}" title="{escaped_tooltip}" '
        f'aria-label="{escaped_tooltip}"></button>'
    )


def render_badge(
    text: str,
    variant: str = "",
    tooltip: str = "",
) -> str:
    """
    Render a badge (pill-shaped label).

    Parameters
    ----------
    text
        Badge text
    variant
        Variant class for styling. Built-in variants:
        - "" (default gray)
        - "adata-badge-view" (blue, for views)
        - "adata-badge-backed" (orange, for backed mode)
        - "adata-badge-sparse" (green, for sparse matrices)
        - "adata-badge-dask" (purple, for Dask arrays)
        - "adata-badge-extension" (for extension types)
    tooltip
        Tooltip text on hover

    Returns
    -------
    HTML string for the badge

    Example
    -------
    >>> render_badge("Zarr", "adata-badge-backed", "Backed by Zarr store")
    """
    escaped_text = escape_html(text)
    title_attr = f' title="{escape_html(tooltip)}"' if tooltip else ""
    # Always include base class, optionally add variant
    css_class = f"adata-badge {variant}".strip() if variant else "adata-badge"
    return f'<span class="{css_class}"{title_attr}>{escaped_text}</span>'


def render_header_badges(
    *,
    is_view: bool = False,
    is_backed: bool = False,
    backing_path: str | None = None,
    backing_format: str | None = None,
) -> str:
    """
    Render standard header badges for view/backed status.

    Parameters
    ----------
    is_view
        Whether this is a view
    is_backed
        Whether this is backed by a file
    backing_path
        Path to the backing file (for tooltip)
    backing_format
        Format of the backing file ("H5AD", "Zarr", etc.)

    Returns
    -------
    HTML string with badges

    Example
    -------
    >>> badges = render_header_badges(
    ...     is_backed=True,
    ...     backing_path="/data/sample.zarr",
    ...     backing_format="Zarr",
    ... )
    """
    parts = []
    if is_view:
        parts.append(
            render_badge("View", "adata-badge-view", "This is a view of another object")
        )
    if is_backed:
        tooltip = f"Backed by {backing_path}" if backing_path else "Backed mode"
        label = backing_format or "Backed"
        parts.append(render_badge(label, "adata-badge-backed", tooltip))
    return "".join(parts)


# =============================================================================
# Name Cell Renderer
# =============================================================================


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


# =============================================================================
# Section Renderers
# =============================================================================


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

    # Render entries (with truncation)
    rows = []
    for i, col_name in enumerate(df.columns):
        if i >= max_items:
            rows.append(_render_truncation_indicator(n_cols - max_items))
            break
        col = df[col_name]
        rows.append(_render_dataframe_entry(adata, section, col_name, col, context))

    return render_section(
        section,
        "\n".join(rows),
        n_items=n_cols,
        doc_url=doc_url,
        tooltip=tooltip,
        should_collapse=n_cols > fold_threshold,
        count_str=f"({n_cols} columns)",
    )


def _render_category_list(
    categories: list,
    colors: list[str] | None,
    max_cats: int,
) -> str:
    """Render a list of category values with optional color dots."""
    parts = ['<span class="adata-cats-list">']
    for i, cat in enumerate(categories[:max_cats]):
        cat_name = escape_html(str(cat))
        color = colors[i] if colors and i < len(colors) else None
        parts.append('<span class="adata-cat-item">')
        if color:
            parts.append(
                f'<span style="{STYLE_CAT_DOT}background:{escape_html(color)};"></span>'
            )
        parts.append(f"<span>{cat_name}</span>")
        parts.append("</span>")

    if len(categories) > max_cats:
        remaining = len(categories) - max_cats
        parts.append(f'<span class="adata-text-muted">...+{remaining}</span>')
    parts.append("</span>")
    return "".join(parts)


def _render_unique_count(col: pd.Series) -> str:
    """Render unique count for a non-categorical column."""
    # Show "(lazy)" for lazy series to avoid triggering data loading
    if is_lazy_series(col):
        return '<span class="adata-text-muted">(lazy)</span>'

    unique_limit = _get_setting("repr_html_unique_limit", default=DEFAULT_UNIQUE_LIMIT)
    if unique_limit > 0 and len(col) <= unique_limit:
        try:
            n_unique = col.nunique()
            return f'<span class="adata-text-muted">({n_unique} unique)</span>'
        except TypeError:
            # nunique() fails on unhashable types (e.g., object columns with lists/dicts)
            pass
    return ""


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
    unique_limit = _get_setting("repr_html_unique_limit", default=DEFAULT_UNIQUE_LIMIT)
    should_warn, warn_msg = should_warn_string_column(col, unique_limit)
    if should_warn:
        entry_warnings.append(warn_msg)

    # Check column name validity (issue #321)
    name_valid, name_reason, name_hard_error = check_column_name(col_name)
    name_error = False
    if not name_valid:
        entry_warnings.append(name_reason)
        name_error = name_hard_error

    # Check for color mismatch
    color_warning = check_color_category_mismatch(adata, col_name)
    if color_warning:
        entry_warnings.append(color_warning)

    # Get colors if categorical
    colors = get_matching_column_colors(adata, col_name)

    # Copy button hidden by default - JS shows it

    # Add warning/error class (CSS handles color - error is red like in uns)
    entry_class = "adata-entry"
    if entry_warnings:
        entry_class += " warning"
    if not output.is_serializable or name_error:
        entry_class += " error"

    # Build row
    parts = [
        f'<tr class="{entry_class}" data-key="{escape_html(col_name)}" '
        f'data-dtype="{escape_html(output.type_name)}">'
    ]

    # Name cell
    parts.append(_render_name_cell(col_name))

    # Check if this is a categorical column (for wrap button)
    is_categorical = hasattr(col, "cat")
    categories = list(col.cat.categories) if is_categorical else []
    max_cats = _get_setting("repr_html_max_categories", default=DEFAULT_MAX_CATEGORIES)
    n_cats = min(len(categories), max_cats) if is_categorical else 0

    # Type cell
    parts.append('<td class="adata-entry-type">')
    parts.append(
        f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>'
    )
    is_error = not output.is_serializable or name_error
    parts.append(render_warning_icon(entry_warnings, is_error))

    # Add wrap button for categories in the type column
    if is_categorical and n_cats > 0:
        parts.append(
            '<button class="adata-cats-wrap-btn" title="Toggle multi-line view">⋯</button>'
        )
    parts.append("</td>")

    # Meta cell - show category values with colors or unique count
    parts.append('<td class="adata-entry-meta">')
    if is_categorical:
        parts.append(_render_category_list(categories, colors, max_cats))
    elif hasattr(col, "nunique"):
        parts.append(_render_unique_count(col))
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

    # Render entries (with truncation)
    rows = []
    for i, key in enumerate(keys):
        if i >= max_items:
            rows.append(_render_truncation_indicator(n_items - max_items))
            break
        value = mapping[key]
        rows.append(_render_mapping_entry(key, value, context, section))

    return render_section(
        section,
        "\n".join(rows),
        n_items=n_items,
        doc_url=doc_url,
        tooltip=tooltip,
        should_collapse=n_items > fold_threshold,
    )


def _render_type_cell(
    output: FormattedOutput,
    *,
    has_expandable_content: bool,
    extra_warnings: list[str] | None = None,
    key_hard_error: bool = False,
) -> list[str]:
    """Render the type cell for a mapping entry."""
    parts = ['<td class="adata-entry-type">']

    parts.append(
        f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>'
    )
    all_warnings = (extra_warnings or []) + list(output.warnings)
    has_error = not output.is_serializable or key_hard_error
    parts.append(render_warning_icon(all_warnings, has_error))

    if has_expandable_content:
        parts.append(
            f'<button class="adata-expand-btn" style="{STYLE_HIDDEN}" aria-expanded="false">Expand ▼</button>'
        )

    has_columns_list = output.details.get("has_columns_list", False)
    if has_columns_list:
        parts.append(
            '<button class="adata-cols-wrap-btn" title="Toggle multi-line view">⋯</button>'
        )

    if output.html_content and not output.is_expandable:
        parts.append(
            f'<div class="adata-custom-content" style="margin-top:4px;">{output.html_content}</div>'
        )

    parts.append("</td>")
    return parts


def _render_mapping_meta_cell(output: FormattedOutput, section: str) -> list[str]:
    """Render the meta cell for a mapping entry."""
    parts = ['<td class="adata-entry-meta">']

    has_columns_list = output.details.get("has_columns_list", False)
    if has_columns_list and "columns" in output.details:
        columns = output.details["columns"]
        col_str = ", ".join(escape_html(str(c)) for c in columns)
        parts.append(f'<span class="adata-cols-list">[{col_str}]</span>')
    elif "meta_preview" in output.details:
        full_preview = output.details.get(
            "meta_preview_full", output.details["meta_preview"]
        )
        parts.append(
            f'<span title="{escape_html(full_preview)}">{escape_html(output.details["meta_preview"])}</span>'
        )
    elif "shape" in output.details and section in ("obsm", "varm"):
        shape = output.details["shape"]
        if len(shape) >= 2:
            parts.append(f"({format_number(shape[1])} cols)")

    parts.append("</td>")
    return parts


def _render_mapping_entry(
    key: str,
    value: Any,
    context: FormatterContext,
    section: str,
) -> str:
    """Render a single mapping entry."""
    output = formatter_registry.format_value(value, context)

    # Check key name validity (issue #321)
    key_valid, key_reason, key_hard_error = check_column_name(key)
    key_warnings = []
    if not key_valid:
        key_warnings.append(key_reason)

    # Build class list for CSS styling
    entry_class = "adata-entry"
    if output.warnings or key_warnings:
        entry_class += " warning"
    if not output.is_serializable or key_hard_error:
        entry_class += " error"

    has_expandable_content = output.html_content and output.is_expandable

    parts = [
        f'<tr class="{entry_class}" data-key="{escape_html(key)}" '
        f'data-dtype="{escape_html(output.type_name)}">'
    ]

    # Name cell
    parts.append(_render_name_cell(key))

    # Type cell
    parts.extend(
        _render_type_cell(
            output,
            has_expandable_content=has_expandable_content,
            extra_warnings=key_warnings,
            key_hard_error=key_hard_error,
        )
    )

    # Meta cell
    parts.extend(_render_mapping_meta_cell(output, section))

    parts.append("</tr>")

    # Expandable content row
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

    # Render entries (with truncation)
    rows = []
    for i, key in enumerate(keys):
        if i >= max_items:
            rows.append(_render_truncation_indicator(n_items - max_items))
            break
        value = uns[key]
        rows.append(_render_uns_entry(adata, key, value, context, max_depth))

    return render_section(
        "uns",
        "\n".join(rows),
        n_items=n_items,
        doc_url=doc_url,
        tooltip=tooltip,
        should_collapse=n_items > fold_threshold,
    )


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
    max_str_len = _get_setting(
        "repr_html_max_string_length", default=DEFAULT_MAX_STRING_LENGTH
    )

    # Generate preview based on type
    preview = _generate_value_preview(value, max_str_len)
    if preview_note:
        preview = f"{preview_note} {preview}" if preview else preview_note

    # Check key name validity
    key_valid, key_reason, key_hard_error = check_column_name(key)
    all_warnings = list(output.warnings)
    if not key_valid:
        all_warnings.append(key_reason)

    entry_class = "adata-entry"
    if all_warnings:
        entry_class += " warning"
    if not output.is_serializable or key_hard_error:
        entry_class += " error"

    parts = [
        f'<tr class="{entry_class}" data-key="{escape_html(key)}" data-dtype="{escape_html(output.type_name)}">'
    ]

    # Name
    parts.append(_render_name_cell(key))

    # Type
    parts.append('<td class="adata-entry-type">')
    parts.append(
        f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>'
    )
    is_error = not output.is_serializable or key_hard_error
    parts.append(render_warning_icon(all_warnings, is_error))
    parts.append("</td>")

    # Meta - value preview
    parts.append('<td class="adata-entry-meta">')
    if preview:
        parts.append(
            f'<span class="adata-text-muted" title="{escape_html(str(value)[:500])}">{escape_html(preview)}</span>'
        )
    parts.append("</td>")
    parts.append("</tr>")

    return "\n".join(parts)


def _preview_string(value: str, max_len: int) -> str:
    """Preview a string value."""
    if len(value) <= max_len:
        return f'"{value}"'
    return f'"{value[:max_len]}..."'


def _preview_number(value: float | np.integer | np.floating) -> str:
    """Preview a numeric value."""
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, (int, np.integer)):
        return str(value)
    # Float - format nicely
    if value == int(value):
        return str(int(value))
    return f"{value:.6g}"


def _preview_dict(value: dict) -> str:
    """Preview a dict value."""
    n_keys = len(value)
    if n_keys == 0:
        return "{}"
    if n_keys <= 3:
        keys_preview = ", ".join(str(k) for k in list(value.keys())[:3])
        return f"{{{keys_preview}}}"
    keys_preview = ", ".join(str(k) for k in list(value.keys())[:2])
    return f"{{{keys_preview}, ...}} ({n_keys} keys)"


def _preview_sequence(value: list | tuple) -> str:
    """Preview a list or tuple value."""
    n_items = len(value)
    bracket = "[]" if isinstance(value, list) else "()"
    if n_items == 0:
        return bracket
    if n_items <= 3:
        try:
            items = [_preview_item(v) for v in value[:3]]
            if all(items):
                return f"{bracket[0]}{', '.join(items)}{bracket[1]}"
        except Exception:  # noqa: BLE001
            # Intentional broad catch: preview generation is best-effort
            pass
    return f"({n_items} items)"


def _generate_value_preview(value: Any, max_len: int = 100) -> str:
    """Generate a human-readable preview of a value.

    Returns empty string if no meaningful preview can be generated.
    """
    if value is None:
        return "None"
    if isinstance(value, str):
        return _preview_string(value, max_len)
    if isinstance(value, (bool, int, float, np.integer, np.floating)):
        return _preview_number(value)
    if isinstance(value, dict):
        return _preview_dict(value)
    if isinstance(value, (list, tuple)):
        return _preview_sequence(value)
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
    type_label = output.type_name

    # Check key name validity
    key_valid, key_reason, key_hard_error = check_column_name(key)
    all_warnings = list(output.warnings)
    if not key_valid:
        all_warnings.append(key_reason)

    # Build class list for CSS styling
    entry_class = "adata-entry"
    if all_warnings:
        entry_class += " warning"
    if not output.is_serializable or key_hard_error:
        entry_class += " error"

    parts = [
        f'<tr class="{entry_class}" data-key="{escape_html(key)}" data-dtype="{escape_html(type_label)}">'
    ]

    # Name
    parts.append(_render_name_cell(key))

    # Type
    parts.append('<td class="adata-entry-type">')
    parts.append(f'<span class="{output.css_class}">{escape_html(type_label)}</span>')
    is_error = not output.is_serializable or key_hard_error
    parts.append(render_warning_icon(all_warnings, is_error))
    parts.append("</td>")

    # Meta - custom HTML content
    parts.append('<td class="adata-entry-meta">')
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

    parts = [
        f'<tr class="adata-entry" data-key="{escape_html(key)}" data-dtype="colors">'
    ]

    # Name
    parts.append(_render_name_cell(key))

    # Type
    parts.append('<td class="adata-entry-type">')
    parts.append(f'<span class="dtype-object">colors ({n_colors})</span>')
    parts.append("</td>")

    # Meta - color swatches
    parts.append('<td class="adata-entry-meta">')
    parts.append('<span class="adata-color-swatches">')
    parts.extend(
        f'<span class="adata-color-swatch" style="background:{escape_html(str(color))}" title="{escape_html(str(color))}"></span>'
        for color in colors[:15]  # Limit preview
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
    # Expand button hidden by default - JS shows it

    parts = [
        f'<tr class="adata-entry" data-key="{escape_html(key)}" data-dtype="AnnData">'
    ]

    # Name
    parts.append(_render_name_cell(key))

    # Type
    parts.append('<td class="adata-entry-type">')
    parts.append(
        f'<span class="dtype-anndata">AnnData ({format_number(n_obs)} × {format_number(n_vars)})</span>'
    )
    if can_expand:
        parts.append(
            f'<button class="adata-expand-btn" style="{STYLE_HIDDEN}" aria-expanded="false">Expand ▼</button>'
        )
    parts.append("</td>")

    parts.append('<td class="adata-entry-meta"></td>')
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


def _detect_unknown_sections(adata) -> list[tuple[str, str]]:
    """Detect mapping-like attributes that aren't in SECTION_ORDER.

    Returns list of (attr_name, type_description) tuples for unknown sections.
    """
    from collections.abc import Mapping

    # Known sections and internal attributes to skip
    known = set(SECTION_ORDER) | {
        # Internal/meta attributes
        "shape",
        "n_obs",
        "n_vars",
        "obs_names",
        "var_names",
        "filename",
        "file",
        "is_view",
        "isbacked",
        "isview",
        "T",
        # Methods (not data)
        "obs_keys",
        "var_keys",
        "uns_keys",
        "obsm_keys",
        "varm_keys",
    }

    # Also exclude sections that have registered custom formatters
    # (including those with should_show=False that suppress display)
    known |= set(formatter_registry.get_registered_sections())

    unknown = []
    for attr in dir(adata):
        # Skip private, known, and callable attributes
        if attr.startswith("_") or attr in known:
            continue

        try:
            val = getattr(adata, attr)
            # Check if it's a data container (mapping-like or has keys())
            if isinstance(val, Mapping) or (
                hasattr(val, "keys")
                and hasattr(val, "__getitem__")
                and not callable(val)
            ):
                # Get type description
                type_name = type(val).__name__
                try:
                    n_items = len(val)
                    type_desc = f"{type_name} ({n_items} items)"
                except Exception:  # noqa: BLE001
                    type_desc = type_name
                unknown.append((attr, type_desc))
        except Exception:  # noqa: BLE001
            # If we can't even access the attribute, note it as inaccessible
            unknown.append((attr, "inaccessible"))

    return unknown


def _render_unknown_sections(unknown_sections: list[tuple[str, str]]) -> str:
    """Render a section showing unknown/unrecognized attributes."""
    parts = [
        '<div class="anndata-sec anndata-sec-unknown" data-section="unknown" '
        'data-should-collapse="true">'
    ]
    parts.append('<div class="anndata-sechdr">')
    parts.append(render_fold_icon())
    parts.append('<span class="anndata-sec-name">other</span>')
    parts.append(f'<span class="anndata-sec-count">({len(unknown_sections)})</span>')
    parts.append("</div>")

    parts.append(f'<div class="anndata-seccontent" style="{STYLE_SECTION_CONTENT}">')
    parts.append(f'<table style="{STYLE_SECTION_TABLE}">')

    for attr_name, type_desc in unknown_sections:
        parts.append(f'<tr class="adata-entry" data-key="{escape_html(attr_name)}">')
        parts.append(_render_name_cell(attr_name))
        parts.append('<td class="adata-entry-type">')
        parts.append(
            f'<span class="dtype-unknown" title="Unrecognized attribute">'
            f"{escape_html(type_desc)}</span>"
        )
        parts.append("</td>")
        parts.append('<td class="adata-entry-meta"></td>')
        parts.append("</tr>")

    parts.append("</table>")
    parts.append("</div>")
    parts.append("</div>")

    return "\n".join(parts)


def _render_error_entry(section: str, error: str) -> str:
    """Render an error indicator for a section that failed to render."""
    error_escaped = escape_html(str(error)[:200])  # Truncate long errors
    return f"""
<div class="anndata-sec anndata-sec-error" data-section="{escape_html(section)}">
    <div class="anndata-sechdr">
        {render_fold_icon()}
        <span class="anndata-sec-name">{escape_html(section)}</span>
        <span class="anndata-sec-count adata-error-badge" style="color: var(--anndata-warning, #dc3545);">(error)</span>
    </div>
    <div class="anndata-seccontent" style="{STYLE_SECTION_CONTENT}">
        <div class="adata-error" style="color: var(--anndata-warning, #dc3545); padding: 4px 8px; font-size: 12px;">
            Failed to render: {error_escaped}
        </div>
    </div>
</div>
"""


def _safe_get_attr(obj, attr: str, default="?"):
    """Safely get an attribute with fallback."""
    try:
        val = getattr(obj, attr, None)
        return val if val is not None else default
    except Exception:  # noqa: BLE001
        return default


def _get_raw_meta_parts(raw) -> list[str]:
    """Build meta info parts for raw section."""
    meta_parts = []
    try:
        if hasattr(raw, "var") and raw.var is not None and len(raw.var.columns) > 0:
            meta_parts.append(f"var: {len(raw.var.columns)} cols")
    except Exception:  # noqa: BLE001
        pass
    try:
        if hasattr(raw, "varm") and raw.varm is not None and len(raw.varm) > 0:
            meta_parts.append(f"varm: {len(raw.varm)}")
    except Exception:  # noqa: BLE001
        pass
    return meta_parts


def _render_raw_section(
    adata: AnnData,
    context: FormatterContext,
    fold_threshold: int,
    max_items: int,
) -> str:
    """Render the raw section as a single expandable row.

    The raw section shows unprocessed data that was saved before filtering/normalization.
    It contains raw.X (the matrix), raw.var (variable annotations), and raw.varm
    (multi-dimensional variable annotations).

    Unlike the main AnnData, raw shares obs with the parent but has its own var
    (which may have more variables than the filtered main data).

    Rendered as a single row with an expand button (no section header).
    When expanded, shows a full AnnData-like repr for Raw contents (X, var, varm).
    The depth parameter prevents infinite recursion.
    """
    raw = getattr(adata, "raw", None)
    if raw is None:
        return ""

    # Safely get dimensions with fallbacks
    n_obs = _safe_get_attr(raw, "n_obs", "?")
    n_vars = _safe_get_attr(raw, "n_vars", "?")
    max_depth = context.max_depth

    # Check if we can expand (same logic as nested AnnData)
    can_expand = context.depth < max_depth - 1

    # Build meta info string safely
    meta_parts = _get_raw_meta_parts(raw)
    meta_text = ", ".join(meta_parts) if meta_parts else ""

    # Single row container (like a minimal section with just one entry)
    parts = ['<div class="anndata-sec anndata-sec-raw" data-section="raw">']
    parts.append(f'<table style="{STYLE_SECTION_TABLE}">')

    # Single row with raw info and expand button
    parts.append('<tr class="adata-entry" data-key="raw" data-dtype="Raw">')
    parts.append(_render_name_cell("raw"))
    parts.append('<td class="adata-entry-type">')
    # Show just dimensions - "raw" is already in the name cell
    type_str = f"{format_number(n_obs)} obs × {format_number(n_vars)} var"
    parts.append(f'<span class="dtype-anndata">{escape_html(type_str)}</span>')
    if can_expand:
        parts.append(
            f'<button class="adata-expand-btn" style="{STYLE_HIDDEN}" aria-expanded="false">Expand ▼</button>'
        )
    parts.append("</td>")
    parts.append(f'<td class="adata-entry-meta">{escape_html(meta_text)}</td>')
    parts.append("</tr>")

    # Nested content (hidden by default, shown on expand)
    if can_expand:
        parts.append('<tr class="adata-nested-row">')
        parts.append('<td colspan="3" class="adata-nested-content">')
        parts.append('<div class="adata-nested-anndata">')

        nested_html = _generate_raw_repr_html(
            raw,
            depth=context.depth + 1,
            max_depth=max_depth,
            fold_threshold=fold_threshold,
            max_items=max_items,
        )
        parts.append(nested_html)

        parts.append("</div>")
        parts.append("</td>")
        parts.append("</tr>")

    parts.append("</table>")
    parts.append("</div>")

    return "\n".join(parts)


def _generate_raw_repr_html(
    raw,
    depth: int = 0,
    max_depth: int | None = None,
    fold_threshold: int | None = None,
    max_items: int | None = None,
) -> str:
    """Generate HTML repr for a Raw object.

    This renders X, var, and varm sections similar to AnnData,
    but without obs, obsm, layers, obsp, varp, uns, or raw sections.

    Parameters
    ----------
    raw
        Raw object to render
    depth
        Current nesting depth
    max_depth
        Maximum nesting depth (defaults to settings or DEFAULT_MAX_DEPTH)
    fold_threshold
        Number of items before a section auto-folds (defaults to settings or DEFAULT_FOLD_THRESHOLD)
    max_items
        Maximum items to display per section (defaults to settings or DEFAULT_MAX_ITEMS)
    """
    # Use configured settings with fallback to defaults
    if max_depth is None:
        max_depth = _get_setting("repr_html_max_depth", default=DEFAULT_MAX_DEPTH)
    if fold_threshold is None:
        fold_threshold = _get_setting(
            "repr_html_fold_threshold", default=DEFAULT_FOLD_THRESHOLD
        )
    if max_items is None:
        max_items = _get_setting("repr_html_max_items", default=DEFAULT_MAX_ITEMS)
    from anndata._repr.registry import FormatterContext

    # Safely get dimensions
    n_obs = _safe_get_attr(raw, "n_obs", "?")
    n_vars = _safe_get_attr(raw, "n_vars", "?")

    context = FormatterContext(
        depth=depth,
        max_depth=max_depth,
        adata_ref=None,
    )

    parts = []

    # Container with header showing Raw shape
    container_id = f"raw-repr-{id(raw)}"
    parts.append(f'<div class="anndata-repr" id="{container_id}">')

    # Header for Raw - same structure as AnnData header
    parts.append('<div class="anndata-hdr">')
    parts.append('<span class="adata-type">Raw</span>')
    shape_str = f"{format_number(n_obs)} obs × {format_number(n_vars)} var"
    parts.append(f'<span class="adata-shape">{shape_str}</span>')
    parts.append("</div>")

    # X section - show matrix info (with error handling)
    try:
        if hasattr(raw, "X") and raw.X is not None:
            parts.append(_render_x_entry(raw, context))
    except Exception as e:  # noqa: BLE001
        parts.append(_render_error_entry("X", str(e)))

    # var section (like AnnData's var)
    # _render_dataframe_section expects an object with a .var attribute
    try:
        if hasattr(raw, "var") and raw.var is not None and len(raw.var.columns) > 0:
            var_context = FormatterContext(
                depth=depth,
                max_depth=max_depth,
                adata_ref=None,
                section="var",
            )
            parts.append(
                _render_dataframe_section(
                    raw,  # Pass raw object, not raw.var
                    "var",
                    var_context,
                    fold_threshold=fold_threshold,
                    max_items=max_items,
                )
            )
    except Exception as e:  # noqa: BLE001
        parts.append(_render_error_entry("var", str(e)))

    # varm section (like AnnData's varm)
    try:
        if hasattr(raw, "varm") and raw.varm is not None and len(raw.varm) > 0:
            varm_context = FormatterContext(
                depth=depth,
                max_depth=max_depth,
                adata_ref=None,
                section="varm",
            )
            parts.append(
                _render_mapping_section(
                    raw,  # Pass raw object, not raw.varm
                    "varm",
                    varm_context,
                    fold_threshold=fold_threshold,
                    max_items=max_items,
                )
            )
    except Exception as e:  # noqa: BLE001
        parts.append(_render_error_entry("varm", str(e)))

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


# =============================================================================
# Public API for building custom _repr_html_
# =============================================================================


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

        from anndata._repr import (
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
