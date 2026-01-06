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

from . import (
    DEFAULT_FOLD_THRESHOLD,
    DEFAULT_MAX_CATEGORIES,
    DEFAULT_MAX_DEPTH,
    DEFAULT_MAX_FIELD_WIDTH,
    DEFAULT_MAX_ITEMS,
    DEFAULT_MAX_LAZY_CATEGORIES,
    DEFAULT_MAX_STRING_LENGTH,
    DEFAULT_PREVIEW_ITEMS,
    DEFAULT_TYPE_WIDTH,
    DEFAULT_UNIQUE_LIMIT,
    SECTION_ORDER,
)
from .._repr_constants import TOOLTIP_TRUNCATE_LENGTH
from .components import (
    render_badge,
    render_entry_preview_cell,
    render_entry_row_open,
    render_entry_type_cell,
    render_name_cell,
    render_nested_content_cell,
    render_search_box,
)
from .core import (
    render_section,
    render_truncation_indicator,
    render_x_entry,
)
from .css import get_css
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
    get_setting,
    is_backed,
    is_lazy,
    is_view,
)

if TYPE_CHECKING:
    from typing import Any

    import pandas as pd

    from anndata import AnnData

    from .registry import SectionFormatter

    from .registry import FormattedEntry, FormattedOutput

# Import formatters to register them (side-effect import)
from . import formatters as _formatters  # noqa: F401

from .._repr_constants import (
    CHAR_WIDTH_PX,
    COPY_BUTTON_PADDING_PX,
    DEFAULT_FIELD_WIDTH_PX,
    MIN_FIELD_WIDTH_PX,
)


def _calculate_field_name_width(adata: AnnData, max_width: int) -> int:
    """
    Calculate the optimal field name column width based on longest field name.

    Collects field names from obs, var, uns, obsm, varm, layers, obsp, varp
    and returns a pixel width that fits the longest name (up to max_width).

    Uses constants from _repr_constants.py tuned for the default 13px monospace font.
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
        return DEFAULT_FIELD_WIDTH_PX

    # Find longest name and convert to pixels
    max_len = max(len(name) for name in all_names)
    width_px = (max_len * CHAR_WIDTH_PX) + COPY_BUTTON_PADDING_PX

    # Clamp to reasonable range
    return max(MIN_FIELD_WIDTH_PX, min(width_px, max_width))


def _resolve_setting(override: int | None, setting_name: str, default: int) -> int:
    """Resolve a setting value with priority: explicit override > anndata.settings > default.

    Parameters
    ----------
    override
        Explicit value passed to generate_repr_html (highest priority)
    setting_name
        Name of the anndata.settings attribute to check
    default
        Fallback default value (lowest priority)
    """
    if override is not None:
        return override
    return get_setting(setting_name, default=default)


def _create_formatter_context(
    adata: AnnData,
    *,
    depth: int = 0,
    max_depth: int | None = None,
    fold_threshold: int | None = None,
    max_items: int | None = None,
    max_lazy_categories: int | None = None,
) -> FormatterContext:
    """Create a FormatterContext with settings resolution.

    Resolves None values from anndata.settings, then from defaults.
    Uses _resolve_setting() for the common pattern: override > settings > default.
    """
    return FormatterContext(
        depth=depth,
        max_depth=_resolve_setting(max_depth, "repr_html_max_depth", DEFAULT_MAX_DEPTH),
        fold_threshold=_resolve_setting(
            fold_threshold, "repr_html_fold_threshold", DEFAULT_FOLD_THRESHOLD
        ),
        max_items=_resolve_setting(max_items, "repr_html_max_items", DEFAULT_MAX_ITEMS),
        max_categories=get_setting(
            "repr_html_max_categories", default=DEFAULT_MAX_CATEGORIES
        ),
        max_lazy_categories=_resolve_setting(
            max_lazy_categories, "repr_html_max_lazy_categories", DEFAULT_MAX_LAZY_CATEGORIES
        ),
        max_string_length=get_setting(
            "repr_html_max_string_length", default=DEFAULT_MAX_STRING_LENGTH
        ),
        unique_limit=get_setting("repr_html_unique_limit", default=DEFAULT_UNIQUE_LIMIT),
        adata_ref=adata,
    )


def generate_repr_html(  # noqa: PLR0913
    adata: AnnData,
    *,
    depth: int = 0,
    max_depth: int | None = None,
    fold_threshold: int | None = None,
    max_items: int | None = None,
    max_lazy_categories: int | None = None,
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
    max_lazy_categories
        Maximum categories to load for lazy categoricals. Set to 0 to disable
        loading categories entirely (metadata-only mode). Uses settings/default if None.
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
    # Check if HTML repr is enabled
    if not get_setting("repr_html_enabled", default=True):
        return f"<pre>{escape_html(repr(adata))}</pre>"

    # Create formatter context (resolves settings)
    context = _create_formatter_context(
        adata,
        depth=depth,
        max_depth=max_depth,
        fold_threshold=fold_threshold,
        max_items=max_items,
        max_lazy_categories=max_lazy_categories,
    )

    # Check max depth
    if depth >= context.max_depth:
        return _render_max_depth_indicator(adata)

    # Generate unique container ID
    container_id = _container_id or f"anndata-repr-{uuid.uuid4().hex[:8]}"

    # Build HTML parts
    parts = []

    # CSS and JS only at top level
    if depth == 0:
        parts.append(get_css())

    # Calculate field name column width based on content
    max_field_width = get_setting(
        "repr_html_max_field_width", default=DEFAULT_MAX_FIELD_WIDTH
    )
    field_width = _calculate_field_name_width(adata, max_field_width)

    # Get type column width from settings
    type_width = get_setting("repr_html_type_width", default=DEFAULT_TYPE_WIDTH)

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
    parts.append(render_x_entry(adata, context))
    parts.extend(_render_all_sections(adata, context))
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
) -> list[str]:
    """Render all standard and custom sections."""
    parts = []
    custom_sections_after = _get_custom_sections_by_position(adata)

    for section in SECTION_ORDER:
        if section == "X":
            # X is already rendered, but check for custom sections after X
            if "X" in custom_sections_after:
                parts.extend(
                    _render_custom_section(adata, section_formatter, context)
                    for section_formatter in custom_sections_after["X"]
                )
            continue
        parts.append(_render_section(adata, section, context))

        # Render custom sections after this section
        if section in custom_sections_after:
            parts.extend(
                _render_custom_section(adata, section_formatter, context)
                for section_formatter in custom_sections_after[section]
            )

    # Custom sections at end (no specific position)
    if None in custom_sections_after:
        parts.extend(
            _render_custom_section(adata, section_formatter, context)
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
) -> str:
    """Render a single standard section."""
    try:
        if section == "X":
            return render_x_entry(adata, context)
        if section == "raw":
            return _render_raw_section(adata, context)
        if section in ("obs", "var"):
            return _render_dataframe_section(adata, section, context)
        if section == "uns":
            return _render_uns_section(adata, context)
        return _render_mapping_section(adata, section, context)
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
    adata: AnnData,
    formatter: SectionFormatter,
    context: FormatterContext,
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
        if i >= context.max_items:
            rows.append(render_truncation_indicator(n_items - context.max_items))
            break
        rows.append(render_formatted_entry(entry, section_name))

    # Use render_section for consistent structure
    return render_section(
        getattr(formatter, "display_name", section_name),
        "\n".join(rows),
        n_items=n_items,
        doc_url=getattr(formatter, "doc_url", None),
        tooltip=getattr(formatter, "tooltip", ""),
        should_collapse=n_items > context.fold_threshold,
        section_id=section_name,
    )


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
    """
    output = entry.output
    extra_warnings = extra_warnings or []

    # Compute entry CSS classes
    all_warnings = extra_warnings + list(output.warnings)
    has_error = not output.is_serializable or is_hard_error

    has_expandable_content = output.expanded_html is not None
    has_columns_list = output.details.get("has_columns_list", False)

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

    # Type cell (using consolidated helper)
    parts.append(
        render_entry_type_cell(
            output.type_name,
            output.css_class,
            type_html=output.type_html,
            tooltip=output.tooltip,
            warnings=all_warnings,
            is_error=has_error,
            has_expandable_content=has_expandable_content,
            has_columns_list=has_columns_list,
        )
    )

    # Preview cell (using consolidated helper)
    parts.append(
        render_entry_preview_cell(
            preview_html=output.preview_html,
            preview_text=output.preview,
            columns=output.details.get("columns"),
            has_columns_list=has_columns_list,
            tooltip_preview=output.details.get("meta_preview"),
            tooltip_full=output.details.get("meta_preview_full"),
            shape=output.details.get("shape"),
            section=section,
        )
    )

    parts.append("</tr>")

    # Expandable content row
    if has_expandable_content:
        parts.append(render_nested_content_cell(output.expanded_html))

    return "\n".join(parts)


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

    if is_lazy(adata):
        parts.append(render_badge("Lazy", "adata-badge-lazy"))

    # Check for extension type (not standard AnnData)
    if type_name != "AnnData":
        parts.append(render_badge(type_name, "adata-badge-extension"))

    # README icon if uns["README"] exists with a string
    readme_content = adata.uns.get("README") if hasattr(adata, "uns") else None
    if isinstance(readme_content, str) and readme_content.strip():
        escaped_readme = escape_html(readme_content)
        # Truncate for no-JS tooltip (first 500 chars)
        tooltip_text = readme_content[:TOOLTIP_TRUNCATE_LENGTH]
        if len(readme_content) > TOOLTIP_TRUNCATE_LENGTH:
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
    """Format a preview of an index.

    Handles bytes index values (from older h5ad files) by decoding them.
    """
    n = len(index)
    if n == 0:
        return "<em>empty</em>"

    def format_value(x: Any) -> str:
        """Format a single index value, decoding bytes if needed."""
        if isinstance(x, bytes):
            try:
                return x.decode("utf-8")
            except UnicodeDecodeError:
                return x.decode("latin-1")
        return str(x)

    preview_n = DEFAULT_PREVIEW_ITEMS
    if n <= preview_n * 2:
        # Show all
        items = [escape_html(format_value(x)) for x in index]
    else:
        # Show first and last
        first = [escape_html(format_value(x)) for x in index[:preview_n]]
        last = [escape_html(format_value(x)) for x in index[-preview_n:]]
        items = [*first, "...", *last]

    return ", ".join(items)


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


def _render_max_depth_indicator(adata: AnnData) -> str:
    """Render indicator when max depth is reached."""
    n_obs = getattr(adata, "n_obs", "?")
    n_vars = getattr(adata, "n_vars", "?")
    return f'<div class="adata-max-depth">AnnData ({format_number(n_obs)} × {format_number(n_vars)}) - max depth reached</div>'
