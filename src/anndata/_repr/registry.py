"""
Registry pattern for extensible HTML formatting.

This module provides a registry system that allows:
1. New data types (TreeData, MuData, SpatialData) to register custom formatters
2. Graceful fallback for unknown types
3. Override of default formatters
4. Section-level and type-level customization

Usage for extending to new types:

    from anndata._repr import register_formatter, TypeFormatter, FormattedOutput

    # Format by Python type (e.g., custom array in obsm)
    @register_formatter
    class MyArrayFormatter(TypeFormatter):
        def can_format(self, obj, context):
            return isinstance(obj, MyArrayType)

        def format(self, obj, context):
            return FormattedOutput(
                type_name=f"MyArray {obj.shape}",
                css_class="anndata-dtype--myarray",
                # preview_html for rightmost column (data preview, counts, etc.)
                preview_html=f'<span class="anndata-text--muted">({obj.n_items} items)</span>',
            )

    # Format by embedded type hint (e.g., tagged data in uns)
    from anndata._repr import extract_uns_type_hint

    @register_formatter
    class MyConfigFormatter(TypeFormatter):
        priority = 100  # Higher priority to check before fallback

        def can_format(self, obj, context):
            hint, _ = extract_uns_type_hint(obj)
            return hint == "mypackage.config"

        def format(self, obj, context):
            hint, data = extract_uns_type_hint(obj)
            return FormattedOutput(
                type_name="config",
                preview_html='<span>Custom config preview</span>',
            )
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING

from .._repr_constants import (
    CSS_DTYPE_EXTENSION,
    CSS_DTYPE_UNKNOWN,
    CSS_TEXT_ERROR,
    CSS_TEXT_WARNING,
    DEFAULT_FOLD_THRESHOLD,
    DEFAULT_MAX_CATEGORIES,
    DEFAULT_MAX_DEPTH,
    DEFAULT_MAX_ITEMS,
    DEFAULT_MAX_LAZY_CATEGORIES,
    DEFAULT_MAX_STRING_LENGTH,
    DEFAULT_UNIQUE_LIMIT,
)
from .utils import escape_html, validate_key

if TYPE_CHECKING:
    from typing import Any


@dataclass
class FormattedOutput:
    """Output from a formatter.

    Visual structure of an entry row::

        ┌─────────────┬────────────────────────────┬─────────────────┐
        │ Name        │ Type                       │ Preview         │
        ├─────────────┼────────────────────────────┼─────────────────┤
        │ (from key)  │ type_html or type_name     │ preview_html or │
        │             │ + warnings + [Expand ▼]    │ preview (text)  │
        └─────────────┴────────────────────────────┴─────────────────┘
                               │ (if expanded_html provided and clicked)
                               ▼
                  ┌─────────────────────────────────────────────────┐
                  │ expanded_html content (collapsible row)         │
                  └─────────────────────────────────────────────────┘

    Field precedence rules
    ----------------------
    Some fields have precedence relationships when multiple are provided:

    **Type column** (``type_name`` vs ``type_html``):
        - ``type_name`` is always required and used for search/filter (data-dtype)
        - If ``type_html`` is provided, it replaces the visual display
        - ``type_name`` is still used for the data-dtype attribute regardless

    **Preview column** (``preview`` vs ``preview_html``):
        - If ``preview_html`` is provided, it is used (raw HTML)
        - Otherwise, ``preview`` is used as plain text (auto-escaped)
        - A warning is logged if both are provided

    Field naming convention
    -----------------------
    - ``*_html`` fields contain raw HTML (caller responsible for escaping)
    - Other string fields are plain text (auto-escaped when rendered)

    Available CSS classes
    ---------------------
    Built-in dtype classes for ``css_class`` (BEM modifiers on ``anndata-dtype``):
        - ``anndata-dtype--category``: Categorical data (purple)
        - ``anndata-dtype--int``, ``anndata-dtype--float``: Numeric types (blue)
        - ``anndata-dtype--bool``: Boolean (red)
        - ``anndata-dtype--string``: String data (dark blue)
        - ``anndata-dtype--sparse``: Sparse matrices (green)
        - ``anndata-dtype--array``, ``anndata-dtype--ndarray``: Arrays (blue)
        - ``anndata-dtype--dataframe``: DataFrames (purple)
        - ``anndata-dtype--anndata``: Nested AnnData (red, bold)
        - ``anndata-dtype--dask``: Dask arrays (orange)
        - ``anndata-dtype--gpu``: GPU arrays (lime)
        - ``anndata-dtype--awkward``: Awkward arrays (orange-red)
        - ``anndata-dtype--unknown``: Unknown types (gray, italic)
        - ``anndata-dtype--extension``: Extension types (purple)
    """

    type_name: str = "unknown"
    """Text for the type column (e.g., 'ndarray (100, 50) float32').
    Always used for data-dtype attribute (search/filter). Auto-escaped.
    Defaults to 'unknown' for resilience when type extraction fails."""

    type_html: str | None = None
    """Optional. Raw HTML to render in type column instead of type_name.
    If provided, replaces the visual rendering but type_name still used for data-dtype."""

    css_class: str = CSS_DTYPE_UNKNOWN
    """CSS class for styling the type column."""

    tooltip: str = ""
    """Tooltip text on hover."""

    warnings: list[str] = field(default_factory=list)
    """Warning messages to display with warning icon."""

    preview: str | None = None
    """Optional. Plain text for preview column (rightmost). Auto-escaped.
    Mutually exclusive with preview_html."""

    preview_html: str | None = None
    """Optional. Raw HTML for preview column (e.g., category pills with colors).
    Takes precedence over preview if both provided (with warning)."""

    expanded_html: str | None = None
    """Optional. Raw HTML for expandable content shown in collapsible row below.
    If provided, an 'Expand ▼' button is added to the type column."""

    is_serializable: bool = True
    """Whether this type can be serialized to H5AD/Zarr."""

    error: str | None = None
    """Hard error message. If set, row is highlighted red and error shown in preview.

    **Precedence**: If ``error`` is set, it takes precedence over ``preview`` and
    ``preview_html`` - the error message is displayed instead of any preview content.

    Used for: formatter failures, key validation errors, property access failures.
    Ecosystem packages can set this explicitly or just raise (caught by registry)."""


@dataclass
class FormattedEntry:
    """A single entry in a section (e.g., one column in obs)."""

    key: str
    """The key/name of this entry"""

    output: FormattedOutput
    """Formatted output for this entry"""


@dataclass
class FormatterContext:
    """Context passed to formatters for stateful formatting."""

    depth: int = 0
    """Current recursion depth"""

    max_depth: int = DEFAULT_MAX_DEPTH
    """Maximum recursion depth"""

    fold_threshold: int = DEFAULT_FOLD_THRESHOLD
    """Auto-fold sections with more than this many entries"""

    max_items: int = DEFAULT_MAX_ITEMS
    """Maximum items to show per section"""

    max_categories: int = DEFAULT_MAX_CATEGORIES
    """Maximum category values to display inline"""

    max_lazy_categories: int = DEFAULT_MAX_LAZY_CATEGORIES
    """Maximum categories to load for lazy categoricals.

    For lazy AnnData (from read_lazy()), loading categories requires reading
    from disk. This limit prevents loading categories for columns with many
    unique values (which would be slow and produce cluttered output).
    Set to 0 to skip loading categories entirely for lazy columns.
    """

    max_string_length: int = DEFAULT_MAX_STRING_LENGTH
    """Truncate strings longer than this in previews"""

    unique_limit: int = DEFAULT_UNIQUE_LIMIT
    """Max rows to compute unique counts (0 to disable)"""

    parent_keys: tuple[str, ...] = ()
    """Keys of parent objects (for building access paths)"""

    adata_ref: Any = None
    """Reference to the root AnnData object (for color lookups etc.)"""

    section: str = ""
    """Current section being formatted (obs, var, uns, etc.)"""

    key: str | None = None
    """Key/name of the current entry being formatted (column name, uns key, etc.)"""

    def child(self, key: str) -> FormatterContext:
        """Create a child context for nested formatting."""
        return replace(
            self,
            depth=self.depth + 1,
            parent_keys=(*self.parent_keys, key),
        )

    @property
    def access_path(self) -> str:
        """Build Python access path string."""
        if not self.parent_keys:
            return ""
        parts = []
        for key in self.parent_keys:
            if key.isidentifier():
                parts.append(f".{key}")
            else:
                parts.append(f"[{key!r}]")
        return "".join(parts)


class TypeFormatter(ABC):
    """
    Base class for type-specific formatters.

    Subclass this to add support for new types. The formatter will be
    called when can_format() returns True for an object.

    Attributes
    ----------
    priority : int
        Determines order of checking (higher = checked first). Default: 0.
    sections : tuple[str, ...] | None
        If set, this formatter only applies to the specified sections.
        Use standard section names: "obs", "var", "uns", "obsm", "varm",
        "layers", "obsp", "varp", "raw", "X".
        If None (default), applies to all sections.

    Examples
    --------
    Formatter that only applies to uns section::

        @register_formatter
        class MyUnsFormatter(TypeFormatter):
            sections = ("uns",)

            def can_format(self, obj, context):
                return isinstance(obj, MySpecialType)

            def format(self, obj, context):
                return FormattedOutput(type_name="MySpecial")

    Formatter that applies to obsm and varm::

        @register_formatter
        class MyMatrixFormatter(TypeFormatter):
            sections = ("obsm", "varm")

            def can_format(self, obj, context):
                return isinstance(obj, MyMatrixType)

            def format(self, obj, context):
                return FormattedOutput(type_name=f"MyMatrix {obj.shape}")

    Formatter that uses context for metadata lookup::

        @register_formatter
        class AnnotatedCategoricalFormatter(TypeFormatter):
            priority = 115  # Higher than default CategoricalFormatter
            sections = ("obs", "var")

            def can_format(self, obj, context):
                # Check if column has metadata in uns
                if not (isinstance(obj, pd.Series) and hasattr(obj, "cat")):
                    return False
                if context.adata_ref is None or context.key is None:
                    return False
                annotations = context.adata_ref.uns.get("__annotations__", {})
                return context.key in annotations.get(context.section, {})

            def format(self, obj, context):
                # Use metadata to render enhanced output
                return FormattedOutput(type_name="category[annotated]")
    """

    priority: int = 0
    sections: tuple[str, ...] | None = None

    @abstractmethod
    def can_format(self, obj: Any, context: FormatterContext) -> bool:
        """Return True if this formatter can handle the given object.

        Parameters
        ----------
        obj
            The object to check.
        context
            Formatter context (see :class:`FormatterContext`). Key attributes:

            - ``section``: Current section ("obs", "var", "uns", etc.)
            - ``key``: Current entry key (column name for obs/var, dict key for uns, etc.)
            - ``adata_ref``: Reference to root AnnData (for uns lookups)

            Use these for context-aware decisions, e.g., looking up metadata
            in ``context.adata_ref.uns`` based on ``context.key``.
            See the "Formatter that uses context" example in the class docstring.
        """
        ...

    @abstractmethod
    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        """Format the object and return FormattedOutput."""
        ...


class SectionFormatter(ABC):
    """
    Base class for section-specific formatters.

    Subclass this to customize how entire sections (obs, var, uns, etc.)
    are formatted. This allows packages like TreeData, MuData, SpatialData
    to add custom sections (e.g., obst, vart, mod, spatial).

    A single SectionFormatter can handle multiple sections by setting
    ``section_names`` (tuple). Use ``should_show()`` returning False to
    suppress sections entirely (they won't appear in "other" either).

    Example usage::

        from anndata._repr import (
            register_formatter,
            SectionFormatter,
            FormattedEntry,
            FormattedOutput,
        )


        @register_formatter
        class ObstSectionFormatter(SectionFormatter):
            section_name = "obst"
            after_section = "obsm"  # Position after obsm

            def should_show(self, obj):
                return hasattr(obj, "obst") and len(obj.obst) > 0

            def get_entries(self, obj, context):
                entries = []
                for key, value in obj.obst.items():
                    output = FormattedOutput(
                        type_name=f"Tree ({value.n_nodes} nodes)",
                        css_class="anndata-dtype--tree",
                    )
                    entries.append(FormattedEntry(key=key, output=output))
                return entries

    Example - suppress multiple sections::

        @register_formatter
        class SuppressInternalSections(SectionFormatter):
            section_names = ("obsmap", "varmap", "axis")

            @property
            def section_name(self) -> str:
                return self.section_names[0]

            def should_show(self, obj) -> bool:
                return False  # Never show

            def get_entries(self, obj, context):
                return []
    """

    @property
    @abstractmethod
    def section_name(self) -> str:
        """Primary name of the section this formatter handles."""
        ...

    @property
    def section_names(self) -> tuple[str, ...]:
        """
        All section names this formatter handles.

        Override this to handle multiple sections with one formatter.
        Defaults to a tuple containing just section_name.
        """
        return (self.section_name,)

    @property
    def display_name(self) -> str:
        """Display name (defaults to section_name)."""
        return self.section_name

    @property
    def after_section(self) -> str | None:
        """
        Section after which this custom section should appear.

        If None, appears at the end. Valid values are standard sections:
        "X", "obs", "var", "uns", "obsm", "varm", "layers", "obsp", "varp", "raw"
        """
        return None

    @property
    def doc_url(self) -> str | None:
        """URL to documentation for this section."""
        return None

    @property
    def tooltip(self) -> str:
        """Tooltip text for section header."""
        return ""

    @abstractmethod
    def get_entries(self, obj: Any, context: FormatterContext) -> list[FormattedEntry]:
        """Get all entries for this section."""
        ...

    def should_show(self, obj: Any) -> bool:
        """Return True if this section should be displayed."""
        return True


class FallbackFormatter(TypeFormatter):
    """
    Fallback formatter for unknown types.

    This is the last line of defense - it must NEVER raise an exception.
    Every single attribute access is wrapped defensively because objects may have
    malicious __getattr__, broken properties, or custom metaclasses that fail.
    """

    priority: int = -1000  # Lowest priority, always checked last

    def can_format(self, obj: Any) -> bool:
        return True  # Can format anything

    def format(  # noqa: PLR0912, PLR0915
        self,
        obj: Any,
        context: FormatterContext,
        *,
        outer_error: str | None = None,
    ) -> FormattedOutput:
        """Format any object defensively, never raising exceptions.

        Parameters
        ----------
        obj
            Object to format
        context
            Formatter context
        outer_error
            Error message from a failed formatter (passed by registry)
        """
        # === Type name (with fallback) ===
        type_name = "unknown"
        try:  # noqa: SIM105
            type_name = type(obj).__name__
        except Exception:  # noqa: BLE001
            pass

        # === Module (with fallback) ===
        module = None
        try:  # noqa: SIM105
            module = type(obj).__module__
        except Exception:  # noqa: BLE001
            pass

        # === Build full name safely ===
        full_name = type_name
        try:
            if module and module != "builtins":
                full_name = f"{module}.{type_name}"
        except Exception:  # noqa: BLE001
            pass

        # === Gather info defensively ===
        tooltip_parts: list[str] = []
        access_errors: list[str] = []

        # Type info for tooltip
        try:  # noqa: SIM105
            tooltip_parts.append(f"Type: {full_name}")
        except Exception:  # noqa: BLE001
            pass

        # Shape
        try:
            if hasattr(obj, "shape"):
                shape = obj.shape
                tooltip_parts.append(f"Shape: {shape}")
        except Exception as e:  # noqa: BLE001
            try:
                access_errors.append(f".shape raised {type(e).__name__}")
            except Exception:  # noqa: BLE001
                access_errors.append(".shape failed")

        # dtype
        try:
            if hasattr(obj, "dtype"):
                dtype = obj.dtype
                tooltip_parts.append(f"Dtype: {dtype}")
        except Exception as e:  # noqa: BLE001
            try:
                access_errors.append(f".dtype raised {type(e).__name__}")
            except Exception:  # noqa: BLE001
                access_errors.append(".dtype failed")

        # len
        try:
            if hasattr(obj, "__len__"):
                length = len(obj)
                tooltip_parts.append(f"Length: {length}")
                if length > 1_000_000_000:
                    access_errors.append(f"len() = {length:,} (suspicious)")
        except Exception as e:  # noqa: BLE001
            try:
                access_errors.append(f"len() raised {type(e).__name__}")
            except Exception:  # noqa: BLE001
                access_errors.append("len() failed")

        # repr (for tooltip only)
        repr_str = None
        try:
            repr_str = repr(obj)
            if repr_str and not repr_str.startswith("<"):
                tooltip_parts.append(f"Repr: {repr_str[:100]}")
        except Exception as e:  # noqa: BLE001
            try:
                access_errors.append(f"repr() raised {type(e).__name__}")
            except Exception:  # noqa: BLE001
                access_errors.append("repr() failed")

        # str
        try:
            str_val = str(obj)
            # Just checking it doesn't fail, not using the result
            _ = str_val
        except Exception as e:  # noqa: BLE001
            try:
                access_errors.append(f"str() raised {type(e).__name__}")
            except Exception:  # noqa: BLE001
                access_errors.append("str() failed")

        # === Combine all errors ===
        all_errors: list[str] = []
        if outer_error:
            all_errors.append(outer_error)
        all_errors.extend(access_errors)

        error = "; ".join(all_errors) if all_errors else None

        # === Check serializability ===
        is_serial = True
        serial_reason = ""
        try:
            from .utils import is_serializable

            is_serial, serial_reason = is_serializable(obj)
        except Exception:  # noqa: BLE001
            pass

        # === Build preview_html for errors ===
        # SECURITY: All text must be HTML-escaped to prevent XSS
        preview_html = None
        warnings: list[str] = []

        # Add serialization reason to warnings if not serializable
        if not is_serial and serial_reason:
            warnings.append(serial_reason)

        if all_errors:
            try:
                error_text = escape_html(", ".join(all_errors))
                preview_html = f'<span class="{CSS_TEXT_ERROR}">{error_text}</span>'
            except Exception:  # noqa: BLE001
                preview_html = f'<span class="{CSS_TEXT_ERROR}">Error</span>'
        else:
            # No errors - check if unknown type warning needed
            try:
                is_extension = module and not module.startswith((
                    "anndata",
                    "numpy",
                    "pandas",
                    "scipy",
                ))
                if not is_extension:
                    warnings.append(f"Unknown type: {full_name}")
                    warning_text = escape_html(f"Unknown type: {full_name}")
                    preview_html = (
                        f'<span class="{CSS_TEXT_WARNING}">{warning_text}</span>'
                    )
            except Exception:  # noqa: BLE001
                pass

        # === Build tooltip safely ===
        tooltip = ""
        try:  # noqa: SIM105
            tooltip = "\n".join(tooltip_parts)
        except Exception:  # noqa: BLE001
            pass

        # === Determine CSS class ===
        css_class = CSS_DTYPE_UNKNOWN
        try:
            is_extension = module and not module.startswith((
                "anndata",
                "numpy",
                "pandas",
                "scipy",
            ))
            if is_extension:
                css_class = CSS_DTYPE_EXTENSION
        except Exception:  # noqa: BLE001
            pass

        return FormattedOutput(
            type_name=type_name,
            css_class=css_class,
            tooltip=tooltip,
            warnings=warnings,
            preview_html=preview_html,
            is_serializable=is_serial,
            error=error,
        )


def _apply_key_validation(
    result: FormattedOutput, context: FormatterContext
) -> FormattedOutput:
    """Apply key validation to a formatted result.

    Checks if context.key is valid for HDF5/Zarr serialization and adds
    warnings to the result if not. This centralizes key validation so
    all formatters benefit automatically.

    Parameters
    ----------
    result
        The FormattedOutput from a formatter
    context
        The formatter context (uses context.key)

    Returns
    -------
    FormattedOutput, possibly with added warnings and updated is_serializable
    """
    if context.key is None:
        return result

    key_valid, key_reason, is_hard_error = validate_key(context.key)
    if key_valid:
        return result

    # Key has issues - add warning and possibly mark as not serializable
    new_warnings = [*result.warnings, key_reason]
    new_serializable = result.is_serializable and not is_hard_error

    return replace(result, warnings=new_warnings, is_serializable=new_serializable)


class FormatterRegistry:
    """
    Registry for type and section formatters.

    This is the central registry that manages all formatters. It supports:
    - Registering new formatters at runtime
    - Priority-based formatter selection
    - Graceful fallback for unknown types
    - Thread-safe operation
    """

    def __init__(self) -> None:
        self._type_formatters: list[TypeFormatter] = []
        self._section_formatters: dict[str, SectionFormatter] = {}
        self._fallback = FallbackFormatter()

    def register_type_formatter(self, formatter: TypeFormatter) -> None:
        """
        Register a type formatter.

        Formatters are checked in priority order (highest first).
        """
        self._type_formatters.append(formatter)
        # Keep sorted by priority (highest first)
        self._type_formatters.sort(key=lambda f: -f.priority)

    def register_section_formatter(self, formatter: SectionFormatter) -> None:
        """Register a section formatter for all its section_names."""
        for name in formatter.section_names:
            self._section_formatters[name] = formatter

    def unregister_type_formatter(self, formatter: TypeFormatter) -> bool:
        """Unregister a type formatter. Returns True if found and removed."""
        try:
            self._type_formatters.remove(formatter)
            return True
        except ValueError:
            return False

    def format_value(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        """
        Format a value using the appropriate formatter.

        Tries each registered formatter in priority order, falling back
        to the fallback formatter if none match. Formatters with a `sections`
        property are only checked if the current section matches.

        If a formatter raises an exception, we continue to try other formatters.
        Failed formatters are accumulated and:
        - If a later formatter succeeds: warn about the failures
        - If all fail: pass accumulated errors to fallback
        """
        from .._warnings import warn

        current_section = context.section
        # Track failures: (full_msg_for_warning, short_msg_for_html)
        failed_formatters: list[tuple[str, str]] = []

        for formatter in self._type_formatters:
            # Check if formatter is restricted to specific sections
            if (
                formatter.sections is not None
                and current_section not in formatter.sections
            ):
                continue

            try:
                if formatter.can_format(obj, context):
                    result = formatter.format(obj, context)
                    # Success! Warn about any earlier failures
                    if failed_formatters:
                        warn_msgs = [f[0] for f in failed_formatters]
                        try:
                            success_name = type(formatter).__name__
                        except Exception:  # noqa: BLE001
                            success_name = "Formatter"
                        warn(
                            f"Formatters failed before {success_name} succeeded: "
                            f"{'; '.join(warn_msgs)}",
                            UserWarning,
                        )
                    return _apply_key_validation(result, context)
            except Exception as e:  # noqa: BLE001
                # Formatter failed - record and continue to next
                # Build both messages defensively
                try:
                    formatter_name = type(formatter).__name__
                except Exception:  # noqa: BLE001
                    formatter_name = "Formatter"

                # Full message for warning (debugging)
                try:
                    full_msg = f"{formatter_name}: {e}"
                except Exception:  # noqa: BLE001
                    full_msg = f"{formatter_name} failed"

                # Short message for HTML (exception type only)
                try:
                    error_type = type(e).__name__
                    short_msg = f"{formatter_name} raised {error_type}"
                except Exception:  # noqa: BLE001
                    short_msg = f"{formatter_name} failed"

                failed_formatters.append((full_msg, short_msg))
                # Warn immediately for debugging
                warn(f"Formatter {full_msg}", UserWarning)

        # No formatter succeeded - pass accumulated errors to fallback
        if failed_formatters:
            short_msgs = [f[1] for f in failed_formatters]
            outer_error = "; ".join(short_msgs)
        else:
            outer_error = None

        result = self._fallback.format(obj, context, outer_error=outer_error)
        return _apply_key_validation(result, context)

    def get_section_formatter(self, section: str) -> SectionFormatter | None:
        """Get the formatter for a section, or None if not registered."""
        return self._section_formatters.get(section)

    def get_registered_sections(self) -> list[str]:
        """Get list of registered section names."""
        return list(self._section_formatters.keys())

    def get_formatter_for(
        self, obj: Any, context: FormatterContext | None = None
    ) -> tuple[TypeFormatter | None, str]:
        """
        Debug helper: find which formatter would handle an object.

        This is useful for understanding formatter priority and debugging
        why a specific formatter is or isn't being selected.

        Parameters
        ----------
        obj
            The object to find a formatter for
        context
            Optional FormatterContext. If None, creates a minimal context.

        Returns
        -------
        tuple of (formatter, reason)
            formatter: The TypeFormatter that would handle this object, or None
            reason: String explaining why this formatter was selected or why none matched

        Examples
        --------
        >>> import numpy as np
        >>> from anndata._repr import formatter_registry, FormatterContext
        >>> arr = np.array([1, 2, 3])
        >>> formatter, reason = formatter_registry.get_formatter_for(arr)
        >>> print(f"{type(formatter).__name__}: {reason}")
        NumpyArrayFormatter: can_format returned True (priority=100)
        """
        if context is None:
            context = FormatterContext()

        current_section = context.section
        for formatter in self._type_formatters:
            # Check section restriction
            if (
                formatter.sections is not None
                and current_section not in formatter.sections
            ):
                continue

            try:
                if formatter.can_format(obj, context):
                    sections_str = (
                        f", sections={formatter.sections}" if formatter.sections else ""
                    )
                    return (
                        formatter,
                        f"can_format returned True (priority={formatter.priority}{sections_str})",
                    )
            except Exception as e:  # noqa: BLE001
                return (
                    None,
                    f"{type(formatter).__name__} raised {type(e).__name__}: {e}",
                )

        return (self._fallback, "No formatter matched, using fallback")

    def list_formatters(self) -> list[dict[str, Any]]:
        """
        List all registered type formatters with their properties.

        Returns a list of dicts with formatter info, sorted by priority (highest first).
        Useful for debugging formatter registration and priority ordering.

        Returns
        -------
        List of dicts with keys: name, priority, sections, module

        Examples
        --------
        >>> from anndata._repr import formatter_registry
        >>> for f in formatter_registry.list_formatters()[:5]:
        ...     print(f"{f['priority']:4d} {f['name']}")
         150 AnnDataFormatter
         120 DaskArrayFormatter
         120 CuPyArrayFormatter
         120 AwkwardArrayFormatter
         110 NumpyMaskedArrayFormatter
        """
        return [
            {
                "name": type(formatter).__name__,
                "priority": formatter.priority,
                "sections": formatter.sections,
                "module": type(formatter).__module__,
            }
            for formatter in self._type_formatters
        ]


# Global registry instance
formatter_registry = FormatterRegistry()


# Type hint key used in uns dicts to indicate custom rendering
UNS_TYPE_HINT_KEY = "__anndata_repr__"


def extract_uns_type_hint(value: Any) -> tuple[str | None, Any]:
    """
    Extract type hint from data if present.

    This is a utility function for TypeFormatter implementations that need
    to handle tagged data. Data can be tagged with a type hint to indicate
    which package should render it, without requiring that package to be
    installed or imported.

    Supports two formats:

    1. Dict with __anndata_repr__ key::

        {"__anndata_repr__": "mypackage.mytype", "data": {"key": "value"}}
        # Returns ("mypackage.mytype", {"data": {"key": "value"}})

    2. String with prefix::

        "__anndata_repr__:mypackage.mytype::actual content here"
        # Returns ("mypackage.mytype", "actual content here")

    If no type hint found, returns (None, original_value).

    How to make a formatter available
    ---------------------------------
    When data contains a type hint but no formatter is registered for it,
    anndata shows a fallback message: "[mypackage.mytype] (import mypackage)".
    This tells users they need to import the package to see the custom rendering.

    To register a formatter that handles tagged data:

    1. In your package (e.g., mypackage/__init__.py), register a TypeFormatter::

        from anndata._repr import (
            register_formatter,
            TypeFormatter,
            FormattedOutput,
            extract_uns_type_hint,
        )


        @register_formatter
        class MyTypeFormatter(TypeFormatter):
            priority = 100  # Check before fallback formatters
            sections = ("uns",)  # Only apply to uns section

            def can_format(self, obj, context):
                hint, _ = extract_uns_type_hint(obj)
                return hint == "mypackage.mytype"

            def format(self, obj, context):
                hint, data = extract_uns_type_hint(obj)
                # Render your custom visualization
                return FormattedOutput(
                    type_name="mytype",
                    preview_html="<span>Custom rendering</span>",
                )

    2. When the user imports your package, the formatter is registered
       and will be used automatically for any data tagged with your hint.

    Parameters
    ----------
    value
        The value to check for a type hint

    Returns
    -------
    tuple of (type_hint or None, cleaned_value)
        If a type hint is found, returns (hint_string, data_without_hint).
        Otherwise returns (None, original_value).
    """
    # Check dict format
    if isinstance(value, dict) and UNS_TYPE_HINT_KEY in value:
        hint = value.get(UNS_TYPE_HINT_KEY)
        if isinstance(hint, str):
            # Return value without the type hint key
            cleaned = {k: v for k, v in value.items() if k != UNS_TYPE_HINT_KEY}
            return hint, cleaned

    # Check string prefix format
    if isinstance(value, str) and value.startswith(f"{UNS_TYPE_HINT_KEY}:"):
        # Format: "__anndata_repr__:type.hint::content"
        rest = value[len(UNS_TYPE_HINT_KEY) + 1 :]  # After "__anndata_repr__:"
        if "::" in rest:
            hint, content = rest.split("::", 1)
            return hint, content

    return None, value


def register_formatter(
    formatter: TypeFormatter | SectionFormatter,
) -> TypeFormatter | SectionFormatter:
    """
    Register a formatter with the global registry.

    Can be used as a decorator:

        @register_formatter
        class MyFormatter(TypeFormatter):
            ...

    Or called directly:

        register_formatter(MyFormatter())
    """
    if isinstance(formatter, type):
        # Called with class, instantiate it
        formatter = formatter()

    if isinstance(formatter, TypeFormatter):
        formatter_registry.register_type_formatter(formatter)
    elif isinstance(formatter, SectionFormatter):
        formatter_registry.register_section_formatter(formatter)
    else:
        msg = f"Expected TypeFormatter or SectionFormatter, got {type(formatter)}"
        raise TypeError(msg)

    return formatter
