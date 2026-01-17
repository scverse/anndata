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

from dataclasses import dataclass, field

from .._repr_constants import (
    CSS_ENTRY,
    CSS_TEXT_MUTED,
    ENTRY_TABLE_COLSPAN,
    NOT_SERIALIZABLE_MSG,
    STYLE_CAT_DOT,
    STYLE_HIDDEN,
)
from .utils import escape_html, sanitize_css_color


def render_entry_row_open(
    key: str,
    dtype: str,
    *,
    has_warnings: bool = False,
    is_error: bool = False,
    extra_classes: str = "",
) -> str:
    """Render the opening <tr> tag for an entry row.

    This consolidates the repeated pattern of building entry row opening tags
    with data attributes for search/filter functionality.

    Parameters
    ----------
    key
        The entry key (column name, field name, etc.)
    dtype
        The data type string (for data-dtype attribute)
    has_warnings
        Whether the entry has warnings
    is_error
        Whether the entry has errors (not serializable, invalid key)
    extra_classes
        Additional CSS classes to include

    Returns
    -------
    Opening <tr> tag with class and data attributes

    Example
    -------
    >>> open_tag = render_entry_row_open("gene_name", "category", has_warnings=True)
    >>> # Returns: '<tr class="anndata-entry warning" data-key="gene_name" data-dtype="category">'
    """
    # Build CSS class string
    classes = [CSS_ENTRY]
    if extra_classes:
        classes.append(extra_classes)
    if has_warnings:
        classes.append("warning")
    if is_error:
        classes.append("error")
    css_class = " ".join(classes)

    escaped_key = escape_html(key)
    escaped_dtype = escape_html(dtype)
    return f'<tr class="{css_class}" data-key="{escaped_key}" data-dtype="{escaped_dtype}">'


def render_warning_icon(
    warnings: list[str], *, is_not_serializable: bool = False
) -> str:
    """Render warning icon with tooltip if there are warnings or serialization issues.

    Parameters
    ----------
    warnings
        List of warning messages to show in tooltip.
    is_not_serializable
        If True, prepends "Not serializable to H5AD/Zarr" to warnings.

    Returns
    -------
    HTML string for warning icon, or empty string if no warnings.
    """
    if not warnings and not is_not_serializable:
        return ""

    # Build the tooltip message
    if is_not_serializable:
        if warnings:
            # "Not serializable: reason1; reason2"
            reasons = "; ".join(warnings)
            title = f"{NOT_SERIALIZABLE_MSG}: {reasons}"
        else:
            # Just "Not serializable to H5AD/Zarr"
            title = NOT_SERIALIZABLE_MSG
    else:
        # Independent warnings joined with ";"
        title = "; ".join(warnings)

    title = escape_html(title)
    return f'<span class="anndata-entry__warning" title="{title}">(!)</span>'


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
    >>> parts = ['<div class="anndata-header">']
    >>> parts.append('<span class="anndata-header__type">SpatialData</span>')
    >>> parts.append('<span style="flex-grow:1;"></span>')  # Spacer
    >>> parts.append(render_search_box(container_id))
    >>> parts.append("</div>")
    """
    search_id = f"{container_id}-search" if container_id else "anndata-search"
    return (
        f'<span class="anndata-search__box" style="{STYLE_HIDDEN}">'
        f'<input type="text" id="{search_id}" name="{search_id}" '
        f'class="anndata-search__input" '
        f'placeholder="Search..." aria-label="Search fields">'
        f'<span class="anndata-search__toggles">'
        f'<button type="button" class="anndata-search__toggle anndata-search__toggle--case" '
        f'title="Match case" aria-label="Match case" aria-pressed="false">Aa</button>'
        f'<button type="button" class="anndata-search__toggle anndata-search__toggle--regex" '
        f'title="Use regular expression" aria-label="Use regular expression" aria-pressed="false">.*</button>'
        f"</span>"
        f"</span>"
        f'<span class="anndata-search__indicator"></span>'
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
    >>> parts = ['<div class="anndata-section__header">']
    >>> parts.append(render_fold_icon())
    >>> parts.append('<span class="anndata-section__name">images</span>')
    >>> parts.append("</div>")
    """
    return f'<span class="anndata-section__fold" style="{STYLE_HIDDEN}">▼</span>'


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
        f'<button class="anndata-entry__copy" style="{STYLE_HIDDEN}" '
        f'data-copy="{escaped_text}" title="{escaped_tooltip}" '
        f'aria-label="{escaped_tooltip}"></button>'
    )


def _render_wrap_button(css_class: str) -> str:
    """Render a wrap toggle button with the specified CSS class.

    Internal helper used by render_categories_wrap_button and render_columns_wrap_button.
    """
    return f'<button class="{css_class}" title="Expand to multi-line view">▼</button>'


def render_categories_wrap_button() -> str:
    """Render a button to toggle category list between single-line and multi-line.

    Returns
    -------
    HTML string for the wrap button (▼ expands, ▲ collapses)
    """
    return _render_wrap_button("anndata-categories__wrap")


def render_columns_wrap_button() -> str:
    """Render a button to toggle column list between single-line and multi-line.

    Returns
    -------
    HTML string for the wrap button (▼ expands, ▲ collapses)
    """
    return _render_wrap_button("anndata-columns__wrap")


def render_muted_span(text: str) -> str:
    """Render text in a muted span (gray color).

    Parameters
    ----------
    text
        Text to render (will be HTML-escaped)

    Returns
    -------
    HTML string with muted styling
    """
    return f'<span class="{CSS_TEXT_MUTED}">{escape_html(text)}</span>'


def render_expand_button() -> str:
    """Render an expand/collapse button for nested content.

    Returns
    -------
    HTML string for the expand button (hidden by default, shown by JS)
    """
    return (
        f'<button class="anndata-entry__expand" style="{STYLE_HIDDEN}" '
        f'aria-expanded="false">Expand ▼</button>'
    )


def render_nested_content_cell(
    html_content: str, colspan: int = ENTRY_TABLE_COLSPAN
) -> str:
    """Render a table cell containing nested/expanded content.

    Parameters
    ----------
    html_content
        The HTML content to display in the nested cell
    colspan
        Number of columns to span (default: 3)

    Returns
    -------
    HTML string for a table row with nested content
    """
    return (
        f'<tr class="anndata-entry--nested">'
        f'<td colspan="{colspan}" class="anndata-entry__nested-content">'
        f'<div class="anndata-entry__expanded">{html_content}</div>'
        f"</td></tr>"
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
        - "anndata-badge--view" (blue, for views)
        - "anndata-badge--backed" (orange, for backed mode)
        - "anndata-badge--sparse" (green, for sparse matrices)
        - "anndata-badge--dask" (purple, for Dask arrays)
        - "anndata-badge--extension" (for extension types)
    tooltip
        Tooltip text on hover

    Returns
    -------
    HTML string for the badge

    Example
    -------
    >>> badge = render_badge("Zarr", "anndata-badge--backed", "Backed by Zarr store")
    """
    escaped_text = escape_html(text)
    title_attr = f' title="{escape_html(tooltip)}"' if tooltip else ""
    # Always include base class, optionally add variant
    css_class = f"anndata-badge {variant}".strip() if variant else "anndata-badge"
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
            render_badge(
                "View", "anndata-badge--view", "This is a view of another object"
            )
        )
    if is_backed:
        tooltip = f"Backed by {backing_path}" if backing_path else "Backed mode"
        label = backing_format or "Backed"
        parts.append(render_badge(label, "anndata-badge--backed", tooltip))
    if is_lazy:
        parts.append(
            render_badge(
                "Lazy", "anndata-badge--lazy", "Lazy loading (experimental read_lazy)"
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
        f'<td class="anndata-entry__name">'
        f'<div class="anndata-entry__name-inner">'
        f'<span class="anndata-entry__name-text" title="{escaped_name}">{escaped_name}</span>'
        f"{render_copy_button(name, 'Copy name')}"
        f"</div>"
        f"</td>"
    )


def render_category_list(
    categories: list,
    colors: list[str] | None,
    max_cats: int,
    *,
    n_hidden: int = 0,
) -> str:
    """Render a list of category values with optional color dots.

    Parameters
    ----------
    categories
        List of category values to display
    colors
        Optional list of colors matching categories
    max_cats
        Maximum number of categories to show
    n_hidden
        Number of additional hidden categories (for lazy truncation).
        These are added to any truncation from max_cats.

    Returns
    -------
    HTML string for the category list
    """
    parts = ['<span class="anndata-categories">']
    for i, cat in enumerate(categories[:max_cats]):
        cat_name = escape_html(str(cat))
        color = colors[i] if colors and i < len(colors) else None
        parts.append('<span class="anndata-categories__item">')
        if color:
            # Sanitize color to prevent CSS injection
            safe_color = sanitize_css_color(str(color))
            if safe_color:
                parts.append(
                    f'<span style="{STYLE_CAT_DOT}background:{safe_color};"></span>'
                )
            # Skip color dot if color is invalid/unsafe
        parts.append(f"<span>{cat_name}</span>")
        parts.append("</span>")

    # Calculate total hidden: from max_cats truncation + lazy truncation
    hidden_from_max_cats = max(0, len(categories) - max_cats)
    total_hidden = hidden_from_max_cats + n_hidden

    if total_hidden > 0:
        parts.append(f'<span class="{CSS_TEXT_MUTED}">...+{total_hidden}</span>')
    parts.append("</span>")
    return "".join(parts)


@dataclass
class TypeCellConfig:
    """Configuration for rendering a type cell.

    Groups the many parameters of render_entry_type_cell into a single object,
    making call sites cleaner and easier to understand.

    Attributes
    ----------
    type_name
        The type name to display (e.g., "ndarray (100, 50) float32")
    css_class
        CSS class for the type span (e.g., "anndata-dtype--ndarray")
    type_html
        Optional custom HTML content for the type cell
    tooltip
        Optional tooltip for the type label
    warnings
        List of warning messages
    is_not_serializable
        Whether the data cannot be serialized to H5AD/Zarr
    has_expandable_content
        Whether to show expand button
    has_columns_list
        Whether to show columns wrap button
    has_categories_list
        Whether to show categories wrap button
    append_type_html
        If True, type_html is appended below type_name instead of replacing it

    Examples
    --------
    >>> config = TypeCellConfig(
    ...     type_name="ndarray (100, 50) float32",
    ...     css_class="anndata-dtype--ndarray",
    ...     tooltip="Dense array",
    ... )
    >>> html = render_entry_type_cell(config)

    With warnings::

        >>> config = TypeCellConfig(
        ...     type_name="object",
        ...     css_class="anndata-dtype--object",
        ...     warnings=["Custom warning"],
        ...     is_not_serializable=True,
        ... )
    """

    type_name: str
    css_class: str
    type_html: str | None = None
    tooltip: str = ""
    warnings: list[str] = field(default_factory=list)
    is_not_serializable: bool = False
    has_expandable_content: bool = False
    has_columns_list: bool = False
    has_categories_list: bool = False
    append_type_html: bool = False


def render_entry_type_cell(config: TypeCellConfig) -> str:
    """Render the type cell for an entry row.

    This is a unified helper that handles all type cell variations:
    - Type label with optional tooltip
    - Custom type_html (as replacement or appended content)
    - Warning icon
    - Expand/wrap buttons

    The type_html and append_type_html config fields control content rendering:

    1. No type_html: Shows type_name in a styled span
       ``<span class="anndata-dtype--X">type_name</span>``

    2. type_html with append_type_html=False (default): type_html REPLACES type_name
       Used for fully custom type content (e.g., category swatches instead of text)

    3. type_html with append_type_html=True: type_html is shown BELOW type_name
       Used to add extra content while keeping the type label
       (e.g., showing category list below "categorical" label)

    Parameters
    ----------
    config
        TypeCellConfig object with all rendering options

    Returns
    -------
    HTML string for the complete type cell

    Examples
    --------
    >>> config = TypeCellConfig(
    ...     type_name="ndarray (100, 50) float32",
    ...     css_class="anndata-dtype--ndarray",
    ...     tooltip="Dense array",
    ... )
    >>> html = render_entry_type_cell(config)
    """
    type_name = config.type_name
    css_class = config.css_class
    type_html = config.type_html
    tooltip = config.tooltip
    warnings = config.warnings
    is_not_serializable = config.is_not_serializable
    has_expandable_content = config.has_expandable_content
    has_columns_list = config.has_columns_list
    has_categories_list = config.has_categories_list
    append_type_html = config.append_type_html

    parts = ['<td class="anndata-entry__type">']

    # Type content: handle different cases
    if type_html and not append_type_html:
        # type_html replaces the type label entirely
        parts.append(type_html)
    elif tooltip:
        parts.append(
            f'<span class="{css_class}" title="{escape_html(tooltip)}">'
            f"{escape_html(type_name)}</span>"
        )
    else:
        parts.append(f'<span class="{css_class}">{escape_html(type_name)}</span>')

    # Warning icon
    parts.append(
        render_warning_icon(warnings or [], is_not_serializable=is_not_serializable)
    )

    # Expand button for expandable content
    if has_expandable_content:
        parts.append(render_expand_button())

    # Wrap buttons
    if has_columns_list:
        parts.append(render_columns_wrap_button())
    if has_categories_list:
        parts.append(render_categories_wrap_button())

    # Appended type_html (for custom inline rendering below the type)
    if type_html and append_type_html:
        parts.append(
            f'<div class="anndata-entry__custom" style="margin-top:4px;">{type_html}</div>'
        )

    parts.append("</td>")
    return "".join(parts)


def render_entry_preview_cell(
    preview_html: str | None = None,
    preview_text: str | None = None,
) -> str:
    """Render the preview cell (third column) for an entry row.

    Formatters are responsible for producing complete preview content.
    This function just wraps it in the appropriate cell element.

    Parameters
    ----------
    preview_html
        Raw HTML content for preview (highest priority)
    preview_text
        Plain text preview (will be escaped and muted)

    Returns
    -------
    HTML string for the preview cell
    """
    parts = ['<td class="anndata-entry__preview">']

    if preview_html:
        parts.append(preview_html)
    elif preview_text:
        parts.append(render_muted_span(preview_text))

    parts.append("</td>")
    return "".join(parts)
