"""
Section-specific renderers for AnnData HTML representation.

This module contains renderers for each section type:
- DataFrame sections (obs, var)
- Mapping sections (obsm, varm, layers, obsp, varp)
- Uns section (unstructured annotations)
- Raw section (unprocessed data)
- Unknown sections (extension attributes)

Error Handling Policy
---------------------
This module uses broad exception handling (``except Exception``) in several places.
This is intentional - user data may contain arbitrary objects that raise unexpected
exceptions when accessed. The repr should never crash; instead it should:

1. Use ``# noqa: BLE001`` to acknowledge the broad catch
2. Provide a fallback (e.g., show "?" or skip the problematic item)
3. Continue rendering the rest of the representation

This ensures a partially-rendered repr is always better than a crashed cell.
"""

from __future__ import annotations

from dataclasses import replace
from typing import TYPE_CHECKING

from .._repr_constants import (
    CSS_DTYPE_ANNDATA,
    CSS_DTYPE_UNKNOWN,
    ERROR_TRUNCATE_LENGTH,
    INTERNAL_ANNDATA_ATTRS,
    STYLE_SECTION_CONTENT,
    STYLE_SECTION_TABLE,
)
from . import (
    get_section_doc_url,
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
from .core import (
    get_section_tooltip,
    render_empty_section,
    render_section,
    render_truncation_indicator,
    render_x_entry,
)
from .registry import (
    FormattedEntry,
    extract_uns_type_hint,
    formatter_registry,
)
from .utils import (
    check_column_name,
    escape_html,
    format_number,
)

if TYPE_CHECKING:
    from typing import Any

    import pandas as pd

    from anndata import AnnData

    from .registry import FormattedOutput, FormatterContext


def _validate_key_and_collect_warnings(
    key: str,
    output: FormattedOutput,
) -> tuple[list[str], str | None]:
    """Validate key name and return key-specific warnings and errors.

    Checks key validity (issue #321) and returns warnings specific to the key.
    The output.warnings are handled separately by render_formatted_entry.

    Parameters
    ----------
    key
        The entry key to validate
    output
        FormattedOutput (unused, kept for API compatibility)

    Returns
    -------
    tuple of (key_warnings, is_key_not_serializable)
        key_warnings: List of key-specific warning messages
        is_key_not_serializable: True if key prevents serialization
    """
    key_valid, key_reason, is_hard_error = check_column_name(key)

    key_warnings: list[str] = []
    is_key_not_serializable = False

    if not key_valid:
        # All key issues are warnings (shown in tooltip)
        # Hard errors (non-string) also mark as not serializable
        key_warnings.append(key_reason)
        if is_hard_error:
            is_key_not_serializable = True

    return key_warnings, is_key_not_serializable


def _render_entry_row(
    key: str,
    output: FormattedOutput,
    *,
    append_type_html: bool = False,
    preview_note: str | None = None,
) -> str:
    """Render an entry row for DataFrame, mapping, or uns sections.

    This is a thin wrapper around render_formatted_entry that handles
    key validation. The renderer only sees FormattedOutput, never the
    original data object.

    Parameters
    ----------
    key
        Entry key/name to display
    output
        FormattedOutput from a TypeFormatter
    append_type_html
        If True, append type_html below type_name (for mapping entries)
    preview_note
        Optional note to prepend to preview (for type hints in uns)

    Returns
    -------
    HTML string for the entry row (and optional expandable content row)
    """
    # Import here to avoid circular import (html.py imports sections.py)
    from .html import render_formatted_entry

    # Validate key and collect warnings
    extra_warnings, is_key_not_serializable = _validate_key_and_collect_warnings(
        key, output
    )

    # If key has serialization issues, mark output as not serializable
    if is_key_not_serializable and output.is_serializable:
        output = replace(output, is_serializable=False)

    entry = FormattedEntry(key=key, output=output)
    return render_formatted_entry(
        entry,
        extra_warnings=extra_warnings if extra_warnings else None,
        append_type_html=append_type_html,
        preview_note=preview_note,
    )


# -----------------------------------------------------------------------------
# DataFrame Section (obs, var)
# -----------------------------------------------------------------------------


def _render_dataframe_section(
    adata: AnnData,
    section: str,
    context: FormatterContext,
) -> str:
    """Render obs or var section."""
    df: pd.DataFrame = getattr(adata, section)
    n_cols = len(df.columns)

    # Doc URL and tooltip for this section
    doc_url = get_section_doc_url(section)
    tooltip = "Observation annotations" if section == "obs" else "Variable annotations"

    if n_cols == 0:
        return render_empty_section(section, doc_url, tooltip)

    # Set section for section-specific formatters (e.g., LazyColumnFormatter)
    section_context = replace(context, section=section)

    # Render entries (with truncation)
    rows = []
    for i, col_name in enumerate(df.columns):
        if i >= context.max_items:
            rows.append(render_truncation_indicator(n_cols - context.max_items))
            break
        col = df[col_name]
        col_context = replace(section_context, column_name=col_name)
        output = formatter_registry.format_value(col, col_context)
        rows.append(_render_entry_row(col_name, output))

    return render_section(
        section,
        "\n".join(rows),
        n_items=n_cols,
        doc_url=doc_url,
        tooltip=tooltip,
        should_collapse=n_cols > context.fold_threshold,
        count_str=f"({n_cols} columns)",
    )


# -----------------------------------------------------------------------------
# Mapping Section (obsm, varm, layers, obsp, varp)
# -----------------------------------------------------------------------------


def _render_mapping_section(
    adata: AnnData,
    section: str,
    context: FormatterContext,
) -> str:
    """Render obsm, varm, layers, obsp, varp sections."""
    mapping = getattr(adata, section, None)
    if mapping is None:
        return ""

    # Get count without creating full list (O(1) for most mappings)
    n_items = len(mapping)

    # Doc URL and tooltip for this section
    doc_url = get_section_doc_url(section)
    tooltip = get_section_tooltip(section)

    if n_items == 0:
        return render_empty_section(section, doc_url, tooltip)

    # Set section for section-specific formatters (e.g., DaskArrayFormatter)
    section_context = replace(context, section=section)

    # Render entries (with truncation) - iterate lazily, stop at max_items
    rows = []
    for i, key in enumerate(mapping.keys()):
        if i >= context.max_items:
            rows.append(render_truncation_indicator(n_items - context.max_items))
            break
        value = mapping[key]
        output = formatter_registry.format_value(value, section_context)
        rows.append(_render_entry_row(key, output, append_type_html=True))

    return render_section(
        section,
        "\n".join(rows),
        n_items=n_items,
        doc_url=doc_url,
        tooltip=tooltip,
        should_collapse=n_items > context.fold_threshold,
    )


# -----------------------------------------------------------------------------
# Uns Section (unstructured annotations)
# -----------------------------------------------------------------------------


def _render_uns_section(
    adata: AnnData,
    context: FormatterContext,
) -> str:
    """Render the uns section with special handling."""
    uns = adata.uns
    # Get count without creating full list (O(1) for dict)
    n_items = len(uns)

    # Doc URL and tooltip
    doc_url = get_section_doc_url("uns")
    tooltip = "Unstructured annotation"

    if n_items == 0:
        return render_empty_section("uns", doc_url, tooltip)

    # Render entries (with truncation) - iterate lazily, stop at max_items
    rows = []
    for i, key in enumerate(uns.keys()):
        if i >= context.max_items:
            rows.append(render_truncation_indicator(n_items - context.max_items))
            break
        value = uns[key]
        rows.append(_render_uns_entry(adata, key, value, context))

    return render_section(
        "uns",
        "\n".join(rows),
        n_items=n_items,
        doc_url=doc_url,
        tooltip=tooltip,
        should_collapse=n_items > context.fold_threshold,
    )


def _render_uns_entry(
    adata: AnnData,
    key: str,
    value: Any,
    context: FormatterContext,
) -> str:
    """Render a single uns entry with special type handling.

    Rendering priority:
    1. Custom TypeFormatter (may handle type hints, color lists, AnnData)
    2. Unhandled type hint (show import suggestion)
    3. Default formatter
    """
    # Pass key to formatter via column_name for key-based detection (e.g., color lists)
    key_context = replace(context, column_name=key)

    # 1. Try formatter first - handles type hints, color lists, AnnData
    output = formatter_registry.format_value(value, key_context)

    # If a custom formatter produced preview_html, use it directly
    if output.preview_html:
        return _render_entry_row(key, output)

    # 2. Check for unhandled type hint (basic formatter matched, not custom)
    type_hint, cleaned_value = extract_uns_type_hint(value)
    if type_hint is not None:
        # Type hint present but no custom formatter - show import suggestion
        package_name = type_hint.split(".")[0] if "." in type_hint else type_hint
        cleaned_output = formatter_registry.format_value(cleaned_value, key_context)
        return _render_entry_row(
            key,
            cleaned_output,
            preview_note=f"[{type_hint}] (import {package_name} to enable)",
        )

    # 3. Use formatter output
    return _render_entry_row(key, output)


# -----------------------------------------------------------------------------
# Unknown Sections (extension attributes)
# -----------------------------------------------------------------------------


def _detect_unknown_sections(adata: AnnData) -> list[tuple[str, str]]:
    """Detect mapping-like attributes that aren't in SECTION_ORDER.

    Returns list of (attr_name, type_description) tuples for unknown sections.
    """
    from collections.abc import Mapping

    from . import SECTION_ORDER

    # Skip sections we render + internal/meta attributes (not data slots)
    # See INTERNAL_ANNDATA_ATTRS docstring for why this is explicit, not introspected
    known = set(SECTION_ORDER) | INTERNAL_ANNDATA_ATTRS

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
    parts.append('<div class="anndata-section__header">')
    parts.append(render_fold_icon())
    parts.append('<span class="anndata-section__name">other</span>')
    parts.append(
        f'<span class="anndata-section__count">({len(unknown_sections)})</span>'
    )
    parts.append("</div>")

    parts.append(
        f'<div class="anndata-section__content" style="{STYLE_SECTION_CONTENT}">'
    )
    parts.append(f'<table style="{STYLE_SECTION_TABLE}">')

    for attr_name, type_desc in unknown_sections:
        parts.append(render_entry_row_open(attr_name, type_desc))
        parts.append(render_name_cell(attr_name))
        parts.append('<td class="anndata-entry__type">')
        parts.append(
            f'<span class="{CSS_DTYPE_UNKNOWN}" title="Unrecognized attribute">'
            f"{escape_html(type_desc)}</span>"
        )
        parts.append("</td>")
        parts.append('<td class="anndata-entry__preview"></td>')
        parts.append("</tr>")

    parts.append("</table>")
    parts.append("</div>")
    parts.append("</div>")

    return "\n".join(parts)


def _render_error_entry(section: str, error: str) -> str:
    """Render an error indicator for a section that failed to render."""
    error_str = str(error)
    if len(error_str) > ERROR_TRUNCATE_LENGTH:
        error_str = error_str[:ERROR_TRUNCATE_LENGTH] + "..."
    error_escaped = escape_html(error_str)
    # Use --anndata-error-color CSS variable (defined in css.py) with fallback
    error_color = "var(--anndata-error-color, #dc3545)"
    return f"""
<div class="anndata-sec anndata-sec-error" data-section="{escape_html(section)}">
    <div class="anndata-section__header">
        {render_fold_icon()}
        <span class="anndata-section__name">{escape_html(section)}</span>
        <span class="anndata-section__count anndata-badge--error" style="color: {error_color};">(error)</span>
    </div>
    <div class="anndata-section__content" style="{STYLE_SECTION_CONTENT}">
        <div class="anndata-entry--error" style="color: {error_color}; padding: 4px 8px; font-size: 12px;">
            Failed to render: {error_escaped}
        </div>
    </div>
</div>
"""


# -----------------------------------------------------------------------------
# Raw Section
# -----------------------------------------------------------------------------


def _safe_get_attr(obj: Any, attr: str, default: Any = "?") -> Any:
    """Safely get an attribute with fallback.

    Parameters
    ----------
    obj
        Object to get attribute from
    attr
        Attribute name
    default
        Default value if attribute is missing or access raises exception

    Returns
    -------
    Attribute value or default
    """
    try:
        val = getattr(obj, attr, None)
        return val if val is not None else default
    except Exception:  # noqa: BLE001
        return default


def _get_raw_meta_parts(raw: Any) -> list[str]:
    """Build meta info parts for raw section.

    Parameters
    ----------
    raw
        Raw object to extract metadata from

    Returns
    -------
    List of metadata strings like ["var: 5 cols", "varm: 2"]
    """
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

    # Check if we can expand (same logic as nested AnnData)
    can_expand = context.depth < context.max_depth - 1

    # Build meta info string safely
    meta_parts = _get_raw_meta_parts(raw)
    meta_text = ", ".join(meta_parts) if meta_parts else ""

    # Single row container (like a minimal section with just one entry)
    parts = ['<div class="anndata-sec anndata-sec-raw" data-section="raw">']
    parts.append(f'<table style="{STYLE_SECTION_TABLE}">')

    # Single row with raw info and expand button
    type_str = f"{format_number(n_obs)} obs × {format_number(n_vars)} var"
    parts.append(render_entry_row_open("raw", "Raw"))
    parts.append(render_name_cell("raw"))
    type_cell_config = TypeCellConfig(
        type_name=type_str,
        css_class=CSS_DTYPE_ANNDATA,
        has_expandable_content=can_expand,
    )
    parts.append(render_entry_type_cell(type_cell_config))
    parts.append(render_entry_preview_cell(preview_text=meta_text))
    parts.append("</tr>")

    # Nested content (hidden by default, shown on expand)
    if can_expand:
        nested_html = _generate_raw_repr_html(raw, context.child("raw"))
        # Wrap in anndata-entry__nested-anndata for specific styling
        wrapped_html = f'<div class="anndata-entry__nested-anndata">{nested_html}</div>'
        parts.append(render_nested_content_cell(wrapped_html))

    parts.append("</table>")
    parts.append("</div>")

    return "\n".join(parts)


def _generate_raw_repr_html(
    raw,
    context: FormatterContext,
) -> str:
    """Generate HTML repr for a Raw object.

    This renders X, var, and varm sections similar to AnnData,
    but without obs, obsm, layers, obsp, varp, uns, or raw sections.

    Parameters
    ----------
    raw
        Raw object to render
    context
        FormatterContext with depth, max_depth, fold_threshold, max_items
    """
    # Safely get dimensions
    n_obs = _safe_get_attr(raw, "n_obs", "?")
    n_vars = _safe_get_attr(raw, "n_vars", "?")

    parts = []

    # Container with header showing Raw shape
    container_id = f"raw-repr-{id(raw)}"
    parts.append(f'<div class="anndata-repr" id="{container_id}">')

    # Header for Raw - same structure as AnnData header
    parts.append('<div class="anndata-header">')
    parts.append('<span class="anndata-header__type">Raw</span>')
    shape_str = f"{format_number(n_obs)} obs × {format_number(n_vars)} var"
    parts.append(f'<span class="anndata-header__shape">{shape_str}</span>')
    parts.append("</div>")

    # X section - show matrix info (with error handling)
    try:
        if hasattr(raw, "X") and raw.X is not None:
            parts.append(render_x_entry(raw, context))
    except Exception as e:  # noqa: BLE001
        parts.append(_render_error_entry("X", str(e)))

    # var section (like AnnData's var)
    # _render_dataframe_section expects an object with a .var attribute
    try:
        if hasattr(raw, "var") and raw.var is not None and len(raw.var.columns) > 0:
            # Raw doesn't have the same structure as AnnData, so clear adata_ref
            var_context = replace(context, adata_ref=None, section="var")
            parts.append(_render_dataframe_section(raw, "var", var_context))
    except Exception as e:  # noqa: BLE001
        parts.append(_render_error_entry("var", str(e)))

    # varm section (like AnnData's varm)
    try:
        if hasattr(raw, "varm") and raw.varm is not None and len(raw.varm) > 0:
            varm_context = replace(context, adata_ref=None, section="varm")
            parts.append(_render_mapping_section(raw, "varm", varm_context))
    except Exception as e:  # noqa: BLE001
        parts.append(_render_error_entry("varm", str(e)))

    parts.append("</div>")

    return "\n".join(parts)
