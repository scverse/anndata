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
        def can_format(self, obj):
            return isinstance(obj, MyArrayType)

        def format(self, obj, context):
            return FormattedOutput(
                type_name=f"MyArray {obj.shape}",
                css_class="dtype-myarray",
                # preview_html for rightmost column (data preview, counts, etc.)
                preview_html=f'<span class="adata-text-muted">({obj.n_items} items)</span>',
            )

    # Format by embedded type hint (e.g., tagged data in uns)
    from anndata._repr import extract_uns_type_hint

    @register_formatter
    class MyConfigFormatter(TypeFormatter):
        priority = 100  # Higher priority to check before fallback

        def can_format(self, obj):
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

import contextlib
from abc import ABC, abstractmethod
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING

from .._repr_constants import (
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
    Built-in dtype classes for ``css_class``:
        - ``dtype-category``: Categorical data (purple)
        - ``dtype-int``, ``dtype-float``: Numeric types (blue)
        - ``dtype-bool``: Boolean (red)
        - ``dtype-string``: String data (dark blue)
        - ``dtype-sparse``: Sparse matrices (green)
        - ``dtype-array``, ``dtype-ndarray``: Arrays (blue)
        - ``dtype-dataframe``: DataFrames (purple)
        - ``dtype-anndata``: Nested AnnData (red, bold)
        - ``dtype-dask``: Dask arrays (orange)
        - ``dtype-gpu``: GPU arrays (lime)
        - ``dtype-awkward``: Awkward arrays (orange-red)
        - ``dtype-unknown``: Unknown types (gray, italic)
        - ``dtype-extension``: Extension types (purple)
    """

    type_name: str
    """Required. Text for the type column (e.g., 'ndarray (100, 50) float32').
    Always used for data-dtype attribute (search/filter). Auto-escaped."""

    type_html: str | None = None
    """Optional. Raw HTML to render in type column instead of type_name.
    If provided, replaces the visual rendering but type_name still used for data-dtype."""

    css_class: str = "dtype-unknown"
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

    column_name: str | None = None
    """Column name for obs/var entries (for color lookups in uns)"""

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

            def can_format(self, obj):
                return isinstance(obj, MySpecialType)

            def format(self, obj, context):
                return FormattedOutput(type_name="MySpecial")

    Formatter that applies to obsm and varm::

        @register_formatter
        class MyMatrixFormatter(TypeFormatter):
            sections = ("obsm", "varm")

            def can_format(self, obj):
                return isinstance(obj, MyMatrixType)

            def format(self, obj, context):
                return FormattedOutput(type_name=f"MyMatrix {obj.shape}")
    """

    priority: int = 0
    sections: tuple[str, ...] | None = None

    @abstractmethod
    def can_format(self, obj: Any) -> bool:
        """Return True if this formatter can handle the given object."""
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
                        css_class="dtype-tree",
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

    This ensures the repr never breaks, even for completely unknown types.
    """

    priority: int = -1000  # Lowest priority, always checked last

    def can_format(self, obj: Any) -> bool:
        return True  # Can format anything

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        type_name = type(obj).__name__
        module = type(obj).__module__

        # Build a useful type description
        if module and module != "builtins":
            full_name = f"{module}.{type_name}"
        else:
            full_name = type_name

        # Try to get useful info for tooltip
        # All attribute access is wrapped in try/except because objects may have
        # properties that raise exceptions, hang, or behave unexpectedly
        tooltip_parts = [f"Type: {full_name}"]
        access_errors = []

        try:
            if hasattr(obj, "shape"):
                tooltip_parts.append(f"Shape: {obj.shape}")
        except Exception as e:  # noqa: BLE001
            access_errors.append(f".shape raised {type(e).__name__}")

        try:
            if hasattr(obj, "dtype"):
                tooltip_parts.append(f"Dtype: {obj.dtype}")
        except Exception as e:  # noqa: BLE001
            access_errors.append(f".dtype raised {type(e).__name__}")

        try:
            if hasattr(obj, "__len__"):
                length = len(obj)
                tooltip_parts.append(f"Length: {length}")
                # Warn about suspiciously large lengths
                if length > 1_000_000_000:  # > 1 billion
                    access_errors.append(f"len() = {length:,} (suspicious)")
        except Exception as e:  # noqa: BLE001
            access_errors.append(f"len() raised {type(e).__name__}")

        # Try to get repr (might fail for broken objects)
        try:
            repr_str = repr(obj)
            # Only add if it's useful (not just the default object repr)
            if repr_str and not repr_str.startswith("<"):
                tooltip_parts.append(f"Repr: {repr_str[:100]}")
        except Exception as e:  # noqa: BLE001
            access_errors.append(f"repr() raised {type(e).__name__}")

        # Try to get str (might fail for broken objects)
        try:
            str_val = str(obj)
            # Only note if different from repr
            if str_val and str_val != repr_str:
                pass  # str() worked, nothing to report
        except Exception as e:  # noqa: BLE001
            access_errors.append(f"str() raised {type(e).__name__}")

        # Check if this might be from an extension package
        is_extension = module and not module.startswith((
            "anndata",
            "numpy",
            "pandas",
            "scipy",
        ))

        # Build warnings list
        warnings = []
        if not is_extension:
            warnings.append(f"Unknown type: {full_name}")
        if access_errors:
            warnings.append(f"Errors accessing attributes: {', '.join(access_errors)}")

        # Show errors/warnings visibly in preview
        preview = None
        preview_html = None
        if access_errors:
            # Make error info visible in RED
            error_text = ", ".join(access_errors)
            preview_html = f'<span class="{CSS_TEXT_ERROR}">{error_text}</span>'
        elif warnings:
            # Show warnings in ORANGE/YELLOW when there are no errors
            warning_text = "; ".join(warnings)
            preview_html = f'<span class="{CSS_TEXT_WARNING}">{warning_text}</span>'

        return FormattedOutput(
            type_name=type_name,
            css_class="dtype-unknown" if not is_extension else "dtype-extension",
            tooltip="\n".join(tooltip_parts),
            warnings=warnings,
            preview=preview,
            preview_html=preview_html,
            is_serializable=False,  # Assume unknown types aren't serializable
        )


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
        """
        current_section = context.section
        for formatter in self._type_formatters:
            # Check if formatter is restricted to specific sections
            if (
                formatter.sections is not None
                and current_section not in formatter.sections
            ):
                continue

            try:
                # Check with context if formatter supports it, else use basic check
                if hasattr(formatter, "can_format_with_context"):
                    can_fmt = formatter.can_format_with_context(obj, context)
                else:
                    can_fmt = formatter.can_format(obj)

                if can_fmt:
                    return formatter.format(obj, context)
            except Exception as e:  # noqa: BLE001
                # Intentional broad catch: formatters shouldn't crash the entire repr
                from .._warnings import warn

                warn(
                    f"Formatter {type(formatter).__name__} failed for "
                    f"{type(obj).__name__}: {e}",
                    UserWarning,
                )
                continue

        # Use fallback
        return self._fallback.format(obj, context)

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
                if hasattr(formatter, "can_format_with_context"):
                    can_fmt = formatter.can_format_with_context(obj, context)
                else:
                    can_fmt = formatter.can_format(obj)

                if can_fmt:
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

            def can_format(self, obj):
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
