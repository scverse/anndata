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
    STYLE_SECTION_CONTENT,
    STYLE_SECTION_TABLE,
)
from .components import render_fold_icon
from .registry import formatter_registry
from .utils import escape_html, format_number

if TYPE_CHECKING:
    from typing import Any

    from .registry import FormatterContext


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
            render_section,
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
        return render_empty_section(name, doc_url, tooltip)

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


def _render_section_header(
    name: str,
    count_str: str,
    doc_url: str | None,
    tooltip: str,
) -> str:
    """Render a section header - colors handled by CSS for dark mode support."""
    parts = ['<div class="anndata-sechdr">']
    parts.append(render_fold_icon())
    parts.append(f'<span class="anndata-sec-name">{escape_html(name)}</span>')
    parts.append(f'<span class="anndata-sec-count">{escape_html(count_str)}</span>')
    if doc_url:
        parts.append(
            f'<a class="adata-help-link"  href="{escape_html(doc_url)}" target="_blank" title="{escape_html(tooltip)}">?</a>'
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


def render_truncation_indicator(remaining: int) -> str:
    """Render a truncation indicator."""
    return f'<tr><td colspan="3" class="adata-truncated">... and {format_number(remaining)} more</td></tr>'


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

    Works with both AnnData and Raw objects.
    """
    X = obj.X

    parts = ['<div class="adata-x-entry">']
    parts.append("<span>X</span>")

    if X is None:
        parts.append("<span><em>None</em></span>")
    else:
        # Format the X matrix (formatter includes all info like sparsity, on disk, etc.)
        output = formatter_registry.format_value(X, context)
        parts.append(
            f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>'
        )

    parts.append("</div>")
    return "\n".join(parts)
