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
    DEFAULT_PREVIEW_ITEMS,
    ERROR_TRUNCATE_LENGTH,
)
from .components import (
    TypeCellConfig,
    render_entry_preview_cell,
    render_entry_row_open,
    render_entry_type_cell,
    render_name_cell,
    render_nested_content,
)
from .registry import formatter_registry
from .utils import escape_html, format_index_preview, format_number

if TYPE_CHECKING:
    from anndata import AnnData, Raw

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
        Custom count string for header (defaults to "(N items)"). Escaped.

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
    if n_items == 0:
        return render_empty_section(name, doc_url, tooltip, section_id=section_id)

    return render_details_section(
        section_id or name,
        name,
        escape_html(count_str or f"({pluralize(n_items, 'item')})"),
        f'<div class="anndata-section__entries">{entries_html}</div>',
        is_open=not should_collapse,
        doc_url=doc_url,
        tooltip=tooltip,
    )


def render_details_section(  # noqa: PLR0913
    section_id: str,
    name: str,
    count_html: str,
    content_html: str,
    *,
    is_open: bool,
    doc_url: str | None = None,
    tooltip: str = "",
    extra_classes: str = "",
) -> str:
    """Render a foldable section: a ``<details>`` with a summary header.

    This is the single place that produces section markup; ``render_section``,
    ``render_empty_section`` and the error/unknown-attribute sections use it.

    Parameters
    ----------
    section_id
        Value for the ``data-section`` attribute (escaped)
    name
        Display name in the header (escaped)
    count_html
        Count label HTML next to the name, e.g. ``"(3 items)"`` (caller escapes)
    content_html
        Section body HTML (caller escapes)
    is_open
        Whether the section starts expanded
    doc_url
        URL for the help link (? icon)
    tooltip
        Tooltip text for the help link
    extra_classes
        Additional CSS classes for the ``<details>`` element
    """
    classes = f"anndata-section {extra_classes}".strip()
    open_attr = " open" if is_open else ""
    help_link = (
        f'<a class="anndata-section__help" href="{escape_html(doc_url)}" '
        f'target="_blank" title="{escape_html(tooltip)}">?</a>'
        if doc_url
        else ""
    )
    return (
        f'<details class="{classes}" data-section="{escape_html(section_id)}"{open_attr}>'
        f"<summary>"
        f'<span class="anndata-section__name">{escape_html(name)}</span>'
        f'<span class="anndata-section__count">{count_html}</span>'
        f"{help_link}"
        f"</summary>"
        f'<div class="anndata-section__content">{content_html}</div>'
        f"</details>"
    )


def pluralize(n: int, noun: str) -> str:
    """Format a count with a correctly pluralized noun, e.g. ``"1 item"``, ``"2 items"``."""
    return f"{format_number(n)} {noun}{'' if n == 1 else 's'}"


def render_empty_section(
    name: str,
    doc_url: str | None = None,
    tooltip: str = "",
    *,
    section_id: str | None = None,
) -> str:
    """Render an empty (collapsed) section indicator."""
    return render_details_section(
        section_id or name,
        name,
        "(empty)",
        '<div class="anndata-section__empty">No entries</div>',
        is_open=False,
        doc_url=doc_url,
        tooltip=tooltip,
    )


def render_error_section(section: str, error: str) -> str:
    """Render an (expanded) error indicator for a section that failed to render."""
    if len(error) > ERROR_TRUNCATE_LENGTH:
        error = error[:ERROR_TRUNCATE_LENGTH] + "..."
    return render_details_section(
        section,
        section,
        '<span class="anndata-badge--error">(error)</span>',
        f'<div class="anndata-entry--error">Failed to render: {escape_html(error)}</div>',
        is_open=True,
        extra_classes="anndata-sec-error",
    )


def render_index_preview(obj: object) -> str:
    """Render a preview of ``obj.obs_names`` and ``obj.var_names``.

    Works for AnnData, Raw and other objects; missing or broken indices are
    shown as "not available".
    """
    parts = ['<div class="anndata-header__index">']
    for attr in ("obs_names", "var_names"):
        try:
            preview = format_index_preview(getattr(obj, attr), DEFAULT_PREVIEW_ITEMS)
        except Exception:  # noqa: BLE001
            preview = "<em>not available</em>"
        parts.append(f"<div><strong>{attr}:</strong> {preview}</div>")
    parts.append("</div>")
    return "".join(parts)


def render_truncation_indicator(remaining: int) -> str:
    """Render a truncation indicator."""
    return f'<div class="anndata-section__truncated">... and {format_number(remaining)} more</div>'


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


def render_x_entry(obj: AnnData | Raw, context: FormatterContext) -> str:
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
            has_expandable_content=has_expandable_content,
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
        has_columns_list=has_columns_list,
        has_categories_list=has_categories,
        append_type_html=append_type_html,
    )
    parts.append(render_entry_type_cell(type_cell_config))

    # Preview cell: error takes precedence over preview_html, which takes
    # precedence over preview
    preview_html = output.preview_html
    preview_text = output.preview
    if output.error:
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

    # Expandable entries use <details>/<summary>; render_nested_content
    # closes the <summary> and adds the nested content div.
    if output.expanded_html is not None:
        parts.append(render_nested_content(output.expanded_html))
        parts.append("</details>")
    else:
        parts.append("</div>")

    return "\n".join(parts)
