"""
Reusable UI components for HTML representation.

This module provides building blocks for creating consistent HTML representations:
- Warning/error icons with tooltips
- Search box with filter toggles
- Fold/expand icons for collapsible sections
- Copy-to-clipboard buttons
- Status badges (view, backed, sparse, etc.)

These components are designed to be used by both anndata's internal repr
and by external packages (MuData, SpatialData, TreeData) that want to
build compatible representations.
"""

from __future__ import annotations

from .._repr_constants import (
    NOT_SERIALIZABLE_MSG,
    STYLE_HIDDEN,
)
from .utils import escape_html


def render_warning_icon(warnings: list[str], *, is_error: bool = False) -> str:
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
    >>> container_id = "spatialdata-123"
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
    >>> name = "gene_expression"
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
    >>> badge = render_badge("Zarr", "adata-badge-backed", "Backed by Zarr store")
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
    is_lazy: bool = False,
    backing_path: str | None = None,
    backing_format: str | None = None,
) -> str:
    """
    Render standard header badges for view/backed/lazy status.

    Parameters
    ----------
    is_view
        Whether this is a view
    is_backed
        Whether this is backed by a file
    is_lazy
        Whether this uses lazy loading (experimental read_lazy)
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
    if is_lazy:
        parts.append(
            render_badge(
                "Lazy", "adata-badge-lazy", "Lazy loading (experimental read_lazy)"
            )
        )
    return "".join(parts)


def render_name_cell(name: str) -> str:
    """Render a name cell with copy button and tooltip for truncated names.

    The structure uses flexbox so the copy button stays visible even when
    the name text overflows and shows ellipsis.

    Parameters
    ----------
    name
        The field name to display

    Returns
    -------
    HTML string for the table cell
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
