"""
Core rendering primitives for AnnData HTML representation.

This module contains shared rendering functions used by both:
- html.py (main orchestration)
- sections.py (section-specific renderers)

By extracting these to a separate module, we avoid circular imports
between html.py and sections.py.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .._repr_constants import (
    CSS_DTYPE_CATEGORY,
    CSS_DTYPE_DATAFRAME,
    CSS_TEXT_ERROR,
    CSS_TEXT_MUTED,
    ENTRY_TABLE_COLSPAN,
    STYLE_SECTION_CONTENT,
    STYLE_SECTION_TABLE,
)
from .components import (
    TypeCellConfig,
    render_entry_preview_cell,
    render_entry_row_open,
    render_entry_type_cell,
    render_fold_icon,
    render_name_cell,
    render_nested_content_cell,
)
from .registry import formatter_registry
from .utils import escape_html, format_number

if TYPE_CHECKING:
    from typing import Any

    from .registry import FormattedEntry, FormatterContext


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
            CSS_DTYPE_NDARRAY,
            FormattedEntry,
            FormattedOutput,
            render_formatted_entry,
            render_section,
        )

        rows = []
        for key, info in items.items():
            entry = FormattedEntry(
                key=key,
                output=FormattedOutput(
                    type_name=info["type"], css_class=CSS_DTYPE_NDARRAY
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
        return render_empty_section(name, doc_url, tooltip)

    if count_str is None:
        count_str = f"({n_items} items)"

    parts = [
        f'<div class="anndata-section" data-section="{escape_html(section_id)}" '
        f'data-should-collapse="{str(should_collapse).lower()}">'
    ]

    # Header
    parts.append(_render_section_header(name, count_str, doc_url, tooltip))

    # Content
    parts.append(
        f'<div class="anndata-section__content" style="{STYLE_SECTION_CONTENT}">'
    )
    parts.append(
        f'<table class="anndata-section__table" style="{STYLE_SECTION_TABLE}">'
    )
    parts.append(entries_html)
    parts.append("</table></div></div>")

    return "\n".join(parts)


def _render_section_header(
    name: str,
    count_str: str,
    doc_url: str | None,
    tooltip: str,
) -> str:
    """Render a section header - colors handled by CSS for dark mode support."""
    parts = ['<div class="anndata-section__header">']
    parts.append(render_fold_icon())
    parts.append(f'<span class="anndata-section__name">{escape_html(name)}</span>')
    parts.append(
        f'<span class="anndata-section__count">{escape_html(count_str)}</span>'
    )
    if doc_url:
        parts.append(
            f'<a class="anndata-section__help"  href="{escape_html(doc_url)}" target="_blank" title="{escape_html(tooltip)}">?</a>'
        )
    parts.append("</div>")
    return "\n".join(parts)


def render_empty_section(
    name: str,
    doc_url: str | None = None,
    tooltip: str = "",
) -> str:
    """Render an empty section indicator."""
    # Build help link if doc_url provided
    help_link = ""
    if doc_url:
        help_link = f'<a class="anndata-section__help"  href="{escape_html(doc_url)}" target="_blank" title="{escape_html(tooltip)}">?</a>'

    # Use render_fold_icon() helper for consistency
    fold_icon = render_fold_icon()

    return f"""
<div class="anndata-section" data-section="{escape_html(name)}" data-should-collapse="true">
    <div class="anndata-section__header">
        {fold_icon}
        <span class="anndata-section__name">{escape_html(name)}</span>
        <span class="anndata-section__count">(empty)</span>
        {help_link}
    </div>
    <div class="anndata-section__content" style="{STYLE_SECTION_CONTENT}">
        <div class="anndata-section__empty">No entries</div>
    </div>
</div>
"""


def render_truncation_indicator(remaining: int) -> str:
    """Render a truncation indicator."""
    return f'<tr><td colspan="{ENTRY_TABLE_COLSPAN}" class="anndata-section__truncated">... and {format_number(remaining)} more</td></tr>'


def get_section_tooltip(section: str) -> str:
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


def render_x_entry(obj: Any, context: FormatterContext) -> str:
    """Render X as a single compact entry row.

    Works with AnnData, Raw, and any object with an X attribute.
    Handles missing or broken X attributes gracefully.
    """
    parts = ['<div class="anndata-x__entry">']
    parts.append("<span>X</span>")

    try:
        X = obj.X
    except Exception as e:  # noqa: BLE001
        # Handle missing or broken X attribute gracefully
        error_msg = f"error: {type(e).__name__}"
        parts.append(
            f'<span class="{CSS_TEXT_MUTED}"><em>({escape_html(error_msg)})</em></span>'
        )
        parts.append("</div>")
        return "\n".join(parts)

    if X is None:
        parts.append("<span><em>None</em></span>")
    else:
        # Format the X matrix (formatter includes all info like sparsity, on disk, etc.)
        try:
            output = formatter_registry.format_value(X, context)
            parts.append(
                f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>'
            )
        except Exception as e:  # noqa: BLE001
            error_msg = f"error formatting: {type(e).__name__}"
            parts.append(
                f'<span class="{CSS_TEXT_MUTED}"><em>({escape_html(error_msg)})</em></span>'
            )

    parts.append("</div>")
    return "\n".join(parts)


def render_formatted_entry(
    entry: FormattedEntry,
    section: str = "",
    *,
    extra_warnings: list[str] | None = None,
    append_type_html: bool = False,
    preview_note: str | None = None,
) -> str:
    """
    Render a FormattedEntry as a table row.

    This is the unified entry renderer used both internally and as a public API
    for packages building their own _repr_html_.

    Parameters
    ----------
    entry
        A FormattedEntry containing the key and FormattedOutput
    section
        Optional section name (used for meta column rendering)
    extra_warnings
        Additional warnings to display (e.g., key validation warnings)
    append_type_html
        If True, append type_html below type_name instead of replacing it.
        Used for mapping entries (obsm, varm, etc.) to show extra content.
    preview_note
        Optional note to prepend to preview text (for type hints in uns)

    Returns
    -------
    HTML string for the table row(s)

    Examples
    --------
    ::

        from anndata._repr import (
            CSS_DTYPE_ANNDATA,
            CSS_DTYPE_NDARRAY,
            FormattedEntry,
            FormattedOutput,
            render_formatted_entry,
        )

        entry = FormattedEntry(
            key="my_array",
            output=FormattedOutput(
                type_name="ndarray (100, 50) float32",
                css_class=CSS_DTYPE_NDARRAY,
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
                css_class=CSS_DTYPE_ANNDATA,
                expanded_html=nested_html,
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

    With explicit error::

        entry = FormattedEntry(
            key="broken_data",
            output=FormattedOutput(
                type_name="MyType",
                error="Failed to load: file not found",
            ),
        )
        html = render_formatted_entry(entry)
    """
    output = entry.output
    extra_warnings = extra_warnings or []

    # Compute entry CSS classes
    # Both hard errors and serialization issues get red background
    all_warnings = extra_warnings + list(output.warnings)
    has_error = output.error is not None or not output.is_serializable

    has_expandable_content = output.expanded_html is not None
    # Detect wrap button needs from output css_class
    has_categories = output.css_class == CSS_DTYPE_CATEGORY and bool(
        output.preview_html
    )
    has_columns_list = output.css_class == CSS_DTYPE_DATAFRAME and bool(
        output.preview_html
    )

    # Build row using consolidated helper
    parts = [
        render_entry_row_open(
            entry.key,
            output.type_name,
            has_warnings=bool(all_warnings),
            is_error=has_error,
        )
    ]

    # Name cell
    parts.append(render_name_cell(entry.key))

    # Type cell
    type_cell_config = TypeCellConfig(
        type_name=output.type_name,
        css_class=output.css_class,
        type_html=output.type_html if append_type_html else None,
        tooltip=output.tooltip,
        warnings=all_warnings,
        is_not_serializable=not output.is_serializable,
        has_expandable_content=has_expandable_content,
        has_columns_list=has_columns_list,
        has_categories_list=has_categories,
        append_type_html=append_type_html,
    )
    parts.append(render_entry_type_cell(type_cell_config))

    # Preview cell
    # Error takes precedence over preview/preview_html
    preview_html = output.preview_html
    preview_text = output.preview

    if output.error and not preview_html:
        # Generate error preview if error is set but no preview_html provided
        error_text = escape_html(output.error)
        preview_html = f'<span class="{CSS_TEXT_ERROR}">{error_text}</span>'

    if preview_note and preview_text:
        preview_text = f"{preview_note} {preview_text}"
    elif preview_note:
        preview_text = preview_note

    parts.append(
        render_entry_preview_cell(
            preview_html=preview_html,
            preview_text=preview_text,
        )
    )

    parts.append("</tr>")

    # Expandable content row
    if output.expanded_html is not None:
        parts.append(render_nested_content_cell(output.expanded_html))

    return "\n".join(parts)
