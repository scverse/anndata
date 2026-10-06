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
    CSS_BADGE_EXTENSION,
    TOOLTIP_TRUNCATE_LENGTH,
)
from .._settings import settings
from .._types import AnnDataElem
from ..utils import get_literal_members
from .components import (
    render_badge,
    render_header_badges,
    render_search_box,
)
from .core import (
    render_error_section,
    render_formatted_entry,
    render_index_preview,
    render_section,
    render_truncation_indicator,
    render_x_entry,
)
from .css import get_css
from .javascript import get_javascript
from .lazy import get_lazy_backing_info, is_lazy_adata
from .registry import (
    FormatterContext,
    formatter_registry,
)
from .sections import (
    _detect_unknown_sections,
    _render_dataframe_section,
    _render_mapping_section,
    _render_raw_section,
    _render_unknown_sections,
    _render_uns_section,
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
    from anndata import AnnData

    from .registry import SectionFormatter

# Import formatters to register them (side-effect import)
from .._repr_constants import (
    CHAR_WIDTH_PX,
    COPY_BUTTON_PADDING_PX,
    DEFAULT_FIELD_WIDTH_PX,
    MIN_FIELD_WIDTH_PX,
)
from . import formatters as _formatters  # noqa: F401

# Display order of the standard sections: the main matrix first, then the
# annotations in the order of :data:`AnnDataElem`.
_ELEMS: tuple[AnnDataElem, ...] = tuple(get_literal_members(AnnDataElem))
_SECTION_ORDER: tuple[AnnDataElem, ...] = ("X", *(e for e in _ELEMS if e != "X"))


def _collect_all_field_names(adata: AnnData) -> list[str]:
    """
    Collect all field names from standard and custom sections.

    Returns field names from obs/var columns and keys from mapping sections
    (uns, obsm, varm, layers, obsp, varp) plus any registered custom sections.
    """
    all_names: list[str] = []
    standard_sections: set[str] = set(_SECTION_ORDER)

    for section in _SECTION_ORDER:
        if section in {"X", "raw"}:
            continue
        try:
            attr = getattr(adata, section)
            if attr is None:
                continue
            if section in {"obs", "var"}:
                if hasattr(attr, "columns"):
                    all_names.extend(attr.columns.tolist())
            elif hasattr(attr, "keys"):
                # skip `layers[None]`, which is `.X`
                all_names.extend(k for k in attr if k is not None)
        except Exception:  # noqa: BLE001
            # Broken section — skip for width calculation, error placeholder is
            # rendered separately by _render_section.
            pass

    # Registered custom sections (e.g., TreeData's obst/vart)
    for section_name in formatter_registry.get_registered_sections():
        if section_name in standard_sections:
            continue
        try:
            attr = getattr(adata, section_name, None)
            if attr is not None and hasattr(attr, "keys"):
                all_names.extend(attr.keys())
        except Exception:  # noqa: BLE001
            pass

    return all_names


def _calculate_field_name_width(adata: AnnData, max_width: int) -> int:
    """
    Calculate the optimal field name column width based on longest field name.

    Uses _collect_all_field_names() to gather names from all sections,
    then converts the longest name to a pixel width (up to max_width).

    Uses constants from _repr_constants.py tuned for the default 13px monospace font.
    """
    all_names = _collect_all_field_names(adata)

    if not all_names:
        return DEFAULT_FIELD_WIDTH_PX

    # Find longest name and convert to pixels
    max_len = max(len(str(name)) for name in all_names)
    width_px = (max_len * CHAR_WIDTH_PX) + COPY_BUTTON_PADDING_PX

    # Clamp to reasonable range (max_width from user setting always wins)
    return min(max(MIN_FIELD_WIDTH_PX, width_px), max_width)


def _create_formatter_context(
    adata: AnnData,
    *,
    depth: int = 0,
    max_depth: int | None = None,
    fold_threshold: int | None = None,
    max_items: int | None = None,
    max_lazy_categories: int | None = None,
) -> FormatterContext:
    """Create a FormatterContext, using explicit overrides where given and settings otherwise."""

    def resolve(override: int | None, setting: int) -> int:
        return setting if override is None else override

    return FormatterContext(
        depth=depth,
        max_depth=resolve(max_depth, settings.repr_html_max_depth),
        fold_threshold=resolve(fold_threshold, settings.repr_html_fold_threshold),
        max_items=resolve(max_items, settings.repr_html_max_items),
        max_lazy_categories=resolve(
            max_lazy_categories, settings.repr_html_max_lazy_categories
        ),
        max_categories=settings.repr_html_max_categories,
        max_string_length=settings.repr_html_max_string_length,
        unique_limit=settings.repr_html_unique_limit,
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
    if not settings.repr_html_enabled:
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
    field_width = _calculate_field_name_width(adata, settings.repr_html_max_field_width)
    type_width = settings.repr_html_type_width

    # Container with computed column widths as CSS variables.
    # Inline font-family:monospace provides readable fallback when CSS is stripped
    # (GitHub, untrusted notebooks). CSS overrides with its own font stack.
    # Inline min-width on cells + CSS custom properties give column alignment
    # even without a stylesheet.
    style = f"font-family: monospace; --anndata-name-col-width: {field_width}px; --anndata-type-col-width: {type_width}px;"
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
        parts.append(render_index_preview(adata))

    # Sections container
    parts.append('<div class="anndata-repr__sections">')
    parts.extend(_render_all_sections(adata, context))
    parts.append("</div>")  # anndata-repr__sections

    # Footer with metadata (only at top level)
    if depth == 0:
        parts.append(_render_footer(adata))
        # Degradation hints: visible only when CSS or JS is missing.
        # No-CSS hint: visible by default, hidden by CSS.
        parts.append(
            '<div class="anndata-repr__hint-nocss">'
            "<em>Styled representation available in Jupyter and trusted notebooks "
            "(colors, search, type highlighting).</em>"
            "</div>"
        )
        # No-JS hint: hidden by default (no-CSS case already has its own hint),
        # shown by CSS (for static HTML with styles but no JS),
        # hidden again by JS on init.
        parts.append(
            '<div class="anndata-repr__hint-nojs" style="display:none">'
            "<em>Interactive features (search, copy, category wrapping) "
            "require JavaScript. Trust this notebook to enable them.</em>"
            "</div>"
        )

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
    parts: list[str] = []
    custom_sections_after = _get_custom_sections_by_position(adata)

    for section in _SECTION_ORDER:
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

    # Detect and show unknown sections (attributes not in AnnDataElem)
    unknown_sections = _detect_unknown_sections(adata)
    if unknown_sections:
        parts.append(_render_unknown_sections(unknown_sections))

    return parts


def _render_section(
    adata: AnnData,
    section: str,
    context: FormatterContext,
) -> str:
    """Render a single standard section.

    Attribute access happens inside the try/except so a broken section (one
    whose ``getattr`` raises — e.g. a corrupt aligned mapping or a subclass
    with a crashing property) renders as an error placeholder instead of
    aborting the whole repr. This is why we iterate section names directly
    via ``_SECTION_ORDER`` rather than delegating to
    ``iter_outer``, which propagates the first exception it hits.
    """
    try:
        if section == "X":
            return render_x_entry(adata, context)
        elem = getattr(adata, section)
        if section == "raw":
            return _render_raw_section(elem, context)
        if section in ("obs", "var"):
            return _render_dataframe_section(section, elem, context)
        if section == "uns":
            return _render_uns_section(elem, context)
        return _render_mapping_section(section, elem, context)
    except Exception as e:  # noqa: BLE001
        # Show error instead of hiding the section
        return render_error_section(section, f"{type(e).__name__}: {e}")


def _get_custom_sections_by_position(
    adata: object,
) -> dict[str | None, list[SectionFormatter]]:
    """
    Get registered custom section formatters grouped by their position.

    Returns a dict mapping after_section -> list of formatters.
    None key contains formatters that should appear at the end.
    """
    from collections import defaultdict

    result = defaultdict(list)
    standard_section_names: set[str] = set(_SECTION_ORDER)

    for section_name in formatter_registry.get_registered_sections():
        formatter = formatter_registry.get_section_formatter(section_name)
        if formatter is None:
            continue

        # Skip standard sections (they're handled separately)
        if section_name in standard_section_names:
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
    """Render a custom section using its registered formatter.

    If the formatter defines ``render_html(obj, context)``, it is tried
    first and the result is used as-is (no ``<details>`` wrapping).
    If ``render_html`` fails, falls back to the standard ``get_entries``
    path so formatters can provide both an enhanced and a safe representation.
    """
    # Allow formatters to produce raw HTML (e.g., compact inline rows)
    if hasattr(formatter, "render_html"):
        try:
            return formatter.render_html(adata, context)
        except Exception as e:  # noqa: BLE001
            from .._warnings import warn

            warn(
                f"Custom section formatter '{formatter.section_name}' render_html failed, "
                f"falling back to get_entries: {e}",
                UserWarning,
            )
            # Fall through to get_entries below

    try:
        entries = formatter.get_entries(adata, context)
    except Exception as e:  # noqa: BLE001
        # Intentional broad catch: custom formatters shouldn't crash the entire
        # repr. Show the failure like for built-in sections, and warn for debugging.
        from .._warnings import warn

        warn(
            f"Custom section formatter '{formatter.section_name}' failed: {e}",
            UserWarning,
        )
        return render_error_section(formatter.section_name, f"{type(e).__name__}: {e}")

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


def _render_header(
    adata: AnnData, *, show_search: bool = False, container_id: str = ""
) -> str:
    """Render the header with type, shape, badges, and optional search box."""
    parts = ['<div class="anndata-header">']

    # Type name - allow for extension types
    type_name = type(adata).__name__
    parts.append(f'<span class="anndata-header__type">{escape_html(type_name)}</span>')

    # Shape
    shape_str = f"{format_number(adata.n_obs)} obs × {format_number(adata.n_vars)} vars"
    parts.append(f'<span class="anndata-header__shape">{shape_str}</span>')

    # View / backed / lazy badges and backing file path
    backed = is_backed(adata)
    lazy = is_lazy_adata(adata)
    backing_path = backing_format = None
    is_open = None
    if backed:
        backing = get_backing_info(adata)
        backing_path = str(backing.get("filename") or "")
        backing_format = str(backing.get("format") or "")
        is_open = bool(backing.get("is_open"))
    elif lazy:
        lazy_info = get_lazy_backing_info(adata)
        backing_path = lazy_info.get("filename", "")
        backing_format = lazy_info.get("format", "")
    parts.append(
        render_header_badges(
            is_view=is_view(adata),
            is_backed=backed,
            is_lazy=lazy,
            backing_path=backing_path,
            backing_format=backing_format,
            is_open=is_open,
        )
    )

    # Mark subclasses (the type name itself is already shown above)
    if type_name != "AnnData":
        cls = type(adata)
        parts.append(
            render_badge(
                "AnnData subclass",
                CSS_BADGE_EXTENSION,
                f"{cls.__module__}.{cls.__qualname__}",
            )
        )

    # README icon if uns["README"] exists with a string
    readme_content = adata.uns.get("README") if hasattr(adata, "uns") else None
    if isinstance(readme_content, str) and readme_content.strip():
        # Check max README size setting (0 means no limit)
        max_readme_size = settings.repr_html_max_readme_size
        original_len = len(readme_content)
        if max_readme_size > 0 and original_len > max_readme_size:
            # Truncate and add note
            readme_content = readme_content[:max_readme_size]
            truncation_note = (
                f"\n\n---\n*README truncated: showing {max_readme_size:,} of "
                f"{original_len:,} characters*"
            )
            readme_content += truncation_note

        escaped_readme = escape_html(readme_content)
        # Truncate for no-JS tooltip (first 500 chars)
        tooltip_text = readme_content[:TOOLTIP_TRUNCATE_LENGTH]
        if len(readme_content) > TOOLTIP_TRUNCATE_LENGTH:
            tooltip_text += "..."
        escaped_tooltip = escape_html(tooltip_text)

        parts.append(
            f'<span class="anndata-readme__icon" '
            f'data-readme="{escaped_readme}" '
            f'title="{escaped_tooltip}" '
            f'role="button" tabindex="0" aria-label="View README">'
            f"ⓘ"
            f"</span>"
        )

    # Search box on the right (spacer pushes it right) - use render_search_box() helper
    if show_search:
        parts.append('<span class="anndata-spacer"></span>')
        parts.append(render_search_box(container_id))

    parts.append("</div>")
    return "\n".join(parts)


def _render_footer(adata: AnnData) -> str:
    """Render the footer with version and memory info."""
    parts = ['<div class="anndata-footer">']

    # Version
    version = get_anndata_version()
    parts.append(f"<span>anndata v{version}</span>")

    # Memory usage. Omitted for lazy AnnData, where everything stays on disk and
    # __sizeof__ would only count a few in-memory wrappers.
    if not is_lazy_adata(adata):
        try:
            mem_str = format_memory_size(adata.__sizeof__())
            title = (
                "Estimated in-memory size (data on disk not included)"
                if is_backed(adata)
                else "Estimated memory usage"
            )
            parts.append(f'<span title="{title}">~{mem_str}</span>')
        except Exception:  # noqa: BLE001
            # Broad catch: __sizeof__ recursively calls into user data which could raise anything
            pass

    parts.append("</div>")
    return "\n".join(parts)


def _render_max_depth_indicator(adata: AnnData) -> str:
    """Render indicator when max depth is reached."""
    n_obs = getattr(adata, "n_obs", "?")
    n_vars = getattr(adata, "n_vars", "?")
    return f'<div class="anndata-depth-limit">AnnData ({format_number(n_obs)} × {format_number(n_vars)}) - max depth reached</div>'
