"""
Section-specific renderers for AnnData HTML representation.

This module contains renderers for each section type:
- DataFrame sections (obs, var)
- Mapping sections (obsm, varm, layers, obsp, varp)
- Uns section (unstructured annotations)
- Raw section (unprocessed data)
- Unknown sections (extension attributes)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .._repr_constants import (
    STYLE_HIDDEN,
    STYLE_SECTION_CONTENT,
    STYLE_SECTION_TABLE,
)
from . import (
    DOCS_BASE_URL,
    SECTION_ORDER,
)
from .components import render_fold_icon, render_name_cell, render_warning_icon
from .core import (
    get_section_tooltip,
    render_empty_section,
    render_section,
    render_truncation_indicator,
    render_x_entry,
)
from .formatters import check_column_name
from .registry import (
    FormatterContext,
    extract_uns_type_hint,
    formatter_registry,
)
from .utils import (
    check_color_category_mismatch,
    escape_html,
    format_number,
    generate_value_preview,
    get_matching_column_colors,
    is_color_list,
    is_lazy_series,
    should_warn_string_column,
)

if TYPE_CHECKING:
    from typing import Any

    import pandas as pd

    from anndata import AnnData

    from .registry import FormattedOutput

# Category color dot - needs inline for dynamic background color
STYLE_CAT_DOT = "width:8px;height:8px;border-radius:50%;display:inline-block;"


# -----------------------------------------------------------------------------
# Lazy import for generate_repr_html (only case requiring it)
# Used by _render_nested_anndata_entry for recursive AnnData rendering.
# All other shared functions are now imported from core.py.
# -----------------------------------------------------------------------------


def _get_generate_repr_html():
    """Lazy import of generate_repr_html to avoid circular imports."""
    from .html import generate_repr_html

    return generate_repr_html


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
    doc_url = f"{DOCS_BASE_URL}generated/anndata.AnnData.{section}.html"
    tooltip = "Observation annotations" if section == "obs" else "Variable annotations"

    if n_cols == 0:
        return render_empty_section(section, doc_url, tooltip)

    # Set section for section-specific formatters (e.g., LazyColumnFormatter)
    from dataclasses import replace

    section_context = replace(context, section=section)

    # Render entries (with truncation)
    rows = []
    for i, col_name in enumerate(df.columns):
        if i >= context.max_items:
            rows.append(render_truncation_indicator(n_cols - context.max_items))
            break
        col = df[col_name]
        rows.append(
            _render_dataframe_entry(adata, section, col_name, col, section_context)
        )

    return render_section(
        section,
        "\n".join(rows),
        n_items=n_cols,
        doc_url=doc_url,
        tooltip=tooltip,
        should_collapse=n_cols > context.fold_threshold,
        count_str=f"({n_cols} columns)",
    )


def _render_category_list(
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
    """
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

    # Calculate total hidden: from max_cats truncation + lazy truncation
    hidden_from_max_cats = max(0, len(categories) - max_cats)
    total_hidden = hidden_from_max_cats + n_hidden

    if total_hidden > 0:
        parts.append(f'<span class="adata-text-muted">...+{total_hidden}</span>')
    parts.append("</span>")
    return "".join(parts)


def _get_nunique(col: pd.Series, context: FormatterContext) -> int | None:
    """
    Get nunique count for a column, respecting unique_limit.

    Returns None if:
    - Column is lazy (to avoid triggering data loading)
    - Column length exceeds unique_limit (performance optimization)
    - nunique() fails (unhashable types)
    """
    if is_lazy_series(col):
        return None

    if context.unique_limit > 0 and len(col) > context.unique_limit:
        return None

    try:
        return col.nunique()
    except TypeError:
        # nunique() fails on unhashable types (e.g., object columns with lists/dicts)
        return None


def _render_unique_count(n_unique: int | None, *, is_lazy: bool) -> str:
    """Render unique count display for a non-categorical column."""
    if is_lazy:
        return '<span class="adata-text-muted">(lazy)</span>'
    if n_unique is not None:
        return f'<span class="adata-text-muted">({n_unique} unique)</span>'
    return ""


def _is_categorical_column(col: Any) -> bool:
    """
    Check if a column is categorical.

    Works for both pandas Series (.cat accessor) and xarray DataArray
    (CategoricalDtype). For lazy categoricals, checks for CategoricalArray
    without triggering data loading.
    """
    import pandas as pd

    # Pandas Series with categorical dtype
    if hasattr(col, "cat"):
        return True

    # Check for lazy categorical (CategoricalArray) without accessing dtype
    try:
        from anndata.experimental.backed._lazy_arrays import CategoricalArray

        if hasattr(col, "variable") and hasattr(col.variable, "_data"):
            lazy_indexed = col.variable._data
            if hasattr(lazy_indexed, "array") and isinstance(
                lazy_indexed.array, CategoricalArray
            ):
                return True
    except ImportError:
        pass

    # Fallback: xarray DataArray or other objects with CategoricalDtype
    # Note: This may trigger loading for lazy categoricals
    return hasattr(col, "dtype") and isinstance(col.dtype, pd.CategoricalDtype)


def _get_categories_from_column(col: Any) -> list:
    """
    Get categories from a categorical column.

    Works for both pandas Series (.cat.categories) and xarray DataArray
    (dtype.categories).
    """
    # Pandas Series
    if hasattr(col, "cat"):
        return list(col.cat.categories)

    # xarray DataArray or other objects with CategoricalDtype
    if hasattr(col, "dtype") and hasattr(col.dtype, "categories"):
        return list(col.dtype.categories)

    return []


def _get_categorical_array(col: Any):
    """
    Get the underlying CategoricalArray from a lazy xarray DataArray.

    Returns None if not a lazy categorical column.
    """
    try:
        from anndata.experimental.backed._lazy_arrays import CategoricalArray

        # Navigate through xarray structure to find CategoricalArray
        # DataArray -> Variable -> LazilyIndexedArray -> CategoricalArray
        if hasattr(col, "variable") and hasattr(col.variable, "_data"):
            lazy_indexed = col.variable._data
            if hasattr(lazy_indexed, "array"):
                arr = lazy_indexed.array
                if isinstance(arr, CategoricalArray):
                    return arr
    except ImportError:
        pass
    return None


def _get_lazy_category_count(col: Any) -> int | None:
    """
    Get the number of categories for a lazy categorical without loading them.

    For lazy categoricals, we access the underlying CategoricalArray directly
    and read the category count from the zarr/h5 storage metadata, avoiding
    any data loading.

    Returns None if the count cannot be determined.
    """
    # Try to get category count from CategoricalArray without loading
    cat_arr = _get_categorical_array(col)
    if cat_arr is not None:
        try:
            # Access the raw _categories group/array shape
            # For zarr: _categories is a Group with 'values' array
            # For h5: similar structure
            cats = cat_arr._categories
            if hasattr(cats, "keys"):  # It's a group
                values = cats["values"]
                return values.shape[0]
            elif hasattr(cats, "shape"):  # It's an array directly
                return cats.shape[0]
        except Exception:  # noqa: BLE001
            pass
    return None


def _get_lazy_categories(
    col: Any, context: FormatterContext
) -> tuple[list, bool, int | None]:
    """
    Get categories for a lazy categorical column, respecting limits.

    For lazy AnnData (from read_lazy()), this accesses the underlying
    CategoricalArray directly and reads only the needed categories from
    storage, avoiding loading the full categorical data.

    Parameters
    ----------
    col
        Column (potentially lazy) to get categories from
    context
        FormatterContext with max_lazy_categories limit

    Returns
    -------
    tuple of (categories_list, was_skipped, n_categories)
        categories_list: List of category values (empty if skipped)
        was_skipped: True if categories were skipped due to limit
        n_categories: Number of categories (if known), for display when skipped
    """
    # Try to get category count without loading
    n_cats = _get_lazy_category_count(col)

    # If max_lazy_categories is 0, skip loading entirely (metadata-only mode)
    if context.max_lazy_categories == 0:
        return [], True, n_cats

    # Determine if we need to truncate
    should_truncate = n_cats is not None and n_cats > context.max_lazy_categories
    n_to_read = context.max_lazy_categories if should_truncate else n_cats

    # Try to read categories directly from CategoricalArray storage
    # Use read_elem/read_elem_partial for proper decoding of string arrays
    cat_arr = _get_categorical_array(col)
    if cat_arr is not None:
        try:
            from anndata._io.specs.registry import read_elem, read_elem_partial

            cats = cat_arr._categories
            # Get values array: zarr uses group with "values" key, h5 uses array directly
            values = cats["values"] if hasattr(cats, "keys") else cats

            # Use read_elem_partial for sliced reads, read_elem for full reads
            # This ensures proper decoding of string arrays
            if n_to_read is not None and n_to_read < (n_cats or float("inf")):
                categories = list(
                    read_elem_partial(values, indices=slice(0, n_to_read))
                )
            else:
                categories = list(read_elem(values))
            return categories, should_truncate, n_cats
        except Exception:  # noqa: BLE001
            pass

    # Fallback to unified accessor (will trigger loading)
    try:
        return _get_categories_from_column(col), False, n_cats
    except Exception:  # noqa: BLE001
        return [], True, n_cats


def _get_categories_for_column(
    col: Any,
    context: FormatterContext,
    *,
    is_categorical: bool,
    is_lazy: bool,
) -> tuple[list, bool, int | None]:
    """
    Get categories for a column, handling lazy loading appropriately.

    Returns
    -------
    tuple of (categories_list, was_skipped, n_categories)
        categories_list: List of category values
        was_skipped: True if categories were skipped for lazy columns
        n_categories: Number of categories (if known)
    """
    if not is_categorical:
        return [], False, None

    if is_lazy:
        return _get_lazy_categories(col, context)

    # Non-lazy categorical - use unified accessor
    categories = _get_categories_from_column(col)
    return categories, False, len(categories) if categories else None


def _render_dataframe_entry(  # noqa: PLR0912, PLR0915
    adata: AnnData,
    section: str,
    col_name: str,
    col: pd.Series,
    context: FormatterContext,
) -> str:
    """Render a single DataFrame column entry."""
    # Format the column
    output = formatter_registry.format_value(col, context)

    # Check if categorical (determines what we show in meta cell)
    # Use _is_categorical_column to handle both pandas Series and xarray DataArray
    is_categorical = _is_categorical_column(col)
    is_lazy = is_lazy_series(col)

    # Compute nunique once, reuse for warning check and display
    # Only compute for non-categorical columns (categorical already shows categories)
    n_unique = None if is_categorical else _get_nunique(col, context)

    # Check for string->category warning (uses pre-computed n_unique)
    entry_warnings = list(output.warnings)
    should_warn, warn_msg = should_warn_string_column(col, n_unique)
    if should_warn:
        entry_warnings.append(warn_msg)

    # Check column name validity (issue #321)
    name_valid, name_reason, name_hard_error = check_column_name(col_name)
    name_error = False
    if not name_valid:
        entry_warnings.append(name_reason)
        name_error = name_hard_error

    # Get category info first - this determines if we load colors
    # For lazy categoricals, respect max_lazy_categories to avoid loading too much
    # categories_truncated is True if we only loaded first N categories (not all)
    categories, categories_truncated, n_total_cats = _get_categories_for_column(
        col, context, is_categorical=is_categorical, is_lazy=is_lazy
    )
    max_cats = context.max_categories
    n_cats_to_show = min(len(categories), max_cats) if is_categorical else 0

    # Only load colors if we're actually going to display categories
    # This avoids disk I/O in metadata-only mode (max_lazy_categories=0)
    colors = None
    if is_categorical and len(categories) > 0:
        # Check for color mismatch (loads colors from uns, may trigger disk I/O)
        color_warning = check_color_category_mismatch(adata, col_name)
        if color_warning:
            entry_warnings.append(color_warning)

        # Get colors if they match
        # For lazy categoricals with truncated categories, only load the colors
        # we need (up to n_cats_to_show) to avoid loading all colors from disk
        color_limit = n_cats_to_show if (is_lazy and categories_truncated) else None
        colors = get_matching_column_colors(adata, col_name, limit=color_limit)

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
    parts.append(render_name_cell(col_name))

    # Type cell
    parts.append('<td class="adata-entry-type">')
    parts.append(
        f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>'
    )
    is_error = not output.is_serializable or name_error
    parts.append(render_warning_icon(entry_warnings, is_error=is_error))

    # Add wrap button for categories in the type column
    if is_categorical and n_cats_to_show > 0:
        parts.append(
            '<button class="adata-cats-wrap-btn" title="Toggle multi-line view">⋯</button>'
        )
    parts.append("</td>")

    # Meta cell - show category values with colors or unique count
    parts.append('<td class="adata-entry-meta">')
    if is_categorical:
        if len(categories) == 0:
            # Metadata-only mode: no categories loaded, show just count
            if n_total_cats is not None:
                parts.append(
                    f'<span class="adata-text-muted">({n_total_cats} categories)</span>'
                )
            else:
                parts.append('<span class="adata-text-muted">(categories)</span>')
        elif categories_truncated:
            # Partial categories loaded (lazy truncation): show them with truncation
            # Calculate how many additional categories are hidden
            n_hidden = (n_total_cats - len(categories)) if n_total_cats else 0
            parts.append(
                _render_category_list(categories, colors, max_cats, n_hidden=n_hidden)
            )
        else:
            # All categories loaded: show with standard max_categories truncation
            parts.append(_render_category_list(categories, colors, max_cats))
    elif not is_lazy:
        # Non-lazy non-categorical: show unique count
        parts.append(_render_unique_count(n_unique, is_lazy=False))
    # Lazy non-categorical: nothing in meta (would require loading data)
    parts.append("</td>")

    parts.append("</tr>")
    return "\n".join(parts)


# -----------------------------------------------------------------------------
# Mapping Section (obsm, varm, layers, obsp, varp)
# -----------------------------------------------------------------------------


def _render_mapping_section(
    adata: AnnData,
    section: str,
    context: FormatterContext,
) -> str:
    """Render obsm, varm, layers, obsp, varp sections."""
    from dataclasses import replace

    mapping = getattr(adata, section, None)
    if mapping is None:
        return ""

    keys = list(mapping.keys())
    n_items = len(keys)

    # Doc URL and tooltip for this section
    doc_url = f"{DOCS_BASE_URL}generated/anndata.AnnData.{section}.html"
    tooltip = get_section_tooltip(section)

    if n_items == 0:
        return render_empty_section(section, doc_url, tooltip)

    # Set section for section-specific formatters (e.g., DaskArrayFormatter)
    section_context = replace(context, section=section)

    # Render entries (with truncation)
    rows = []
    for i, key in enumerate(keys):
        if i >= context.max_items:
            rows.append(render_truncation_indicator(n_items - context.max_items))
            break
        value = mapping[key]
        rows.append(_render_mapping_entry(key, value, section_context, section))

    return render_section(
        section,
        "\n".join(rows),
        n_items=n_items,
        doc_url=doc_url,
        tooltip=tooltip,
        should_collapse=n_items > context.fold_threshold,
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
    parts.append(render_warning_icon(all_warnings, is_error=has_error))

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
            f'<span class="adata-text-muted" title="{escape_html(full_preview)}">{escape_html(output.details["meta_preview"])}</span>'
        )
    elif "shape" in output.details and section in ("obsm", "varm"):
        shape = output.details["shape"]
        if len(shape) >= 2:
            parts.append(
                f'<span class="adata-text-muted">({format_number(shape[1])} cols)</span>'
            )

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
    parts.append(render_name_cell(key))

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


# -----------------------------------------------------------------------------
# Uns Section (unstructured annotations)
# -----------------------------------------------------------------------------


def _render_uns_section(
    adata: AnnData,
    context: FormatterContext,
) -> str:
    """Render the uns section with special handling."""
    uns = adata.uns
    keys = list(uns.keys())
    n_items = len(keys)

    # Doc URL and tooltip
    doc_url = f"{DOCS_BASE_URL}generated/anndata.AnnData.uns.html"
    tooltip = "Unstructured annotation"

    if n_items == 0:
        return render_empty_section("uns", doc_url, tooltip)

    # Render entries (with truncation)
    rows = []
    for i, key in enumerate(keys):
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
        return _render_nested_anndata_entry(key, value, context)

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

    # Generate preview based on type
    preview = generate_value_preview(value, context.max_string_length)
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
    parts.append(render_name_cell(key))

    # Type
    parts.append('<td class="adata-entry-type">')
    parts.append(
        f'<span class="{output.css_class}">{escape_html(output.type_name)}</span>'
    )
    is_error = not output.is_serializable or key_hard_error
    parts.append(render_warning_icon(all_warnings, is_error=is_error))
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
    parts.append(render_name_cell(key))

    # Type
    parts.append('<td class="adata-entry-type">')
    parts.append(f'<span class="{output.css_class}">{escape_html(type_label)}</span>')
    is_error = not output.is_serializable or key_hard_error
    parts.append(render_warning_icon(all_warnings, is_error=is_error))
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
    parts.append(render_name_cell(key))

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
) -> str:
    """Render a nested AnnData entry."""
    n_obs = getattr(value, "n_obs", "?")
    n_vars = getattr(value, "n_vars", "?")

    can_expand = context.depth < context.max_depth - 1

    # Inline styles for layout - colors via CSS
    # Expand button hidden by default - JS shows it

    parts = [
        f'<tr class="adata-entry" data-key="{escape_html(key)}" data-dtype="AnnData">'
    ]

    # Name
    parts.append(render_name_cell(key))

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
        # Recursive call (lazy import to avoid circular dependency)
        generate_repr_html = _get_generate_repr_html()
        nested_html = generate_repr_html(
            value,
            depth=context.depth + 1,
            max_depth=context.max_depth,
            show_header=True,
            show_search=False,
        )
        parts.append(nested_html)
        parts.append("</div>")
        parts.append("</td>")
        parts.append("</tr>")

    return "\n".join(parts)


# -----------------------------------------------------------------------------
# Unknown Sections (extension attributes)
# -----------------------------------------------------------------------------


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
        parts.append(render_name_cell(attr_name))
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


# -----------------------------------------------------------------------------
# Raw Section
# -----------------------------------------------------------------------------


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
    parts.append('<tr class="adata-entry" data-key="raw" data-dtype="Raw">')
    parts.append(render_name_cell("raw"))
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

        nested_html = _generate_raw_repr_html(raw, context.child("raw"))
        parts.append(nested_html)

        parts.append("</div>")
        parts.append("</td>")
        parts.append("</tr>")

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
    parts.append('<div class="anndata-hdr">')
    parts.append('<span class="adata-type">Raw</span>')
    shape_str = f"{format_number(n_obs)} obs × {format_number(n_vars)} var"
    parts.append(f'<span class="adata-shape">{shape_str}</span>')
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
            var_context = FormatterContext(
                depth=context.depth,
                max_depth=context.max_depth,
                fold_threshold=context.fold_threshold,
                max_items=context.max_items,
                adata_ref=None,
                section="var",
            )
            parts.append(_render_dataframe_section(raw, "var", var_context))
    except Exception as e:  # noqa: BLE001
        parts.append(_render_error_entry("var", str(e)))

    # varm section (like AnnData's varm)
    try:
        if hasattr(raw, "varm") and raw.varm is not None and len(raw.varm) > 0:
            varm_context = FormatterContext(
                depth=context.depth,
                max_depth=context.max_depth,
                fold_threshold=context.fold_threshold,
                max_items=context.max_items,
                adata_ref=None,
                section="varm",
            )
            parts.append(_render_mapping_section(raw, "varm", varm_context))
    except Exception as e:  # noqa: BLE001
        parts.append(_render_error_entry("varm", str(e)))

    parts.append("</div>")

    return "\n".join(parts)
