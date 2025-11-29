"""
Registry pattern for extensible HTML formatting.

This module provides a registry system that allows:
1. New data types (TreeData, MuData, SpatialData) to register custom formatters
2. Graceful fallback for unknown types
3. Override of default formatters
4. Section-level and type-level customization

Usage for extending to new types:

    from anndata._repr import register_formatter, TypeFormatter

    class MyDataFormatter(TypeFormatter):
        def can_format(self, obj: Any) -> bool:
            return isinstance(obj, MyDataType)

        def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
            return FormattedOutput(
                type_name="MyDataType",
                css_class="dtype-mydata",
                details={"custom": "info"},
            )

    register_formatter(MyDataFormatter())
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Callable
import warnings

if TYPE_CHECKING:
    from collections.abc import Sequence


@dataclass
class FormattedOutput:
    """Output from a formatter."""

    type_name: str
    """Display name for the type (e.g., 'ndarray (100, 50) float32')"""

    css_class: str = "dtype-unknown"
    """CSS class for styling"""

    tooltip: str = ""
    """Tooltip text on hover"""

    details: dict[str, Any] = field(default_factory=dict)
    """Additional details (shape, dtype, sparsity, etc.)"""

    warnings: list[str] = field(default_factory=list)
    """Warning messages to display"""

    children: list[FormattedEntry] | None = None
    """Child entries for expandable types"""

    html_content: str | None = None
    """Custom HTML content (for special visualizations like colors)"""

    is_expandable: bool = False
    """Whether this entry can be expanded"""

    is_serializable: bool = True
    """Whether this type can be serialized to H5AD/Zarr"""


@dataclass
class FormattedEntry:
    """A single entry in a section (e.g., one column in obs)."""

    key: str
    """The key/name of this entry"""

    output: FormattedOutput
    """Formatted output for this entry"""

    copy_text: str | None = None
    """Text to copy to clipboard (defaults to key)"""


@dataclass
class FormatterContext:
    """Context passed to formatters for stateful formatting."""

    depth: int = 0
    """Current recursion depth"""

    max_depth: int = 3
    """Maximum recursion depth"""

    parent_keys: tuple[str, ...] = ()
    """Keys of parent objects (for building access paths)"""

    adata_ref: Any = None
    """Reference to the root AnnData object (for color lookups etc.)"""

    section: str = ""
    """Current section being formatted (obs, var, uns, etc.)"""

    def child(self, key: str) -> FormatterContext:
        """Create a child context for nested formatting."""
        return FormatterContext(
            depth=self.depth + 1,
            max_depth=self.max_depth,
            parent_keys=(*self.parent_keys, key),
            adata_ref=self.adata_ref,
            section=self.section,
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

    Priority determines order of checking (higher = checked first).
    """

    priority: int = 0

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
    are formatted. This allows packages like MuData to add custom sections.
    """

    @property
    @abstractmethod
    def section_name(self) -> str:
        """Name of the section this formatter handles."""
        ...

    @property
    def display_name(self) -> str:
        """Display name (defaults to section_name)."""
        return self.section_name

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

        # Try to get useful info
        details = {}
        tooltip_parts = [f"Type: {full_name}"]

        if hasattr(obj, "shape"):
            details["shape"] = obj.shape
            tooltip_parts.append(f"Shape: {obj.shape}")

        if hasattr(obj, "dtype"):
            details["dtype"] = str(obj.dtype)
            tooltip_parts.append(f"Dtype: {obj.dtype}")

        if hasattr(obj, "__len__"):
            try:
                details["length"] = len(obj)
                tooltip_parts.append(f"Length: {len(obj)}")
            except (TypeError, RuntimeError):
                pass

        # Check if this might be from an extension package
        is_extension = module and not module.startswith(("anndata", "numpy", "pandas", "scipy"))

        return FormattedOutput(
            type_name=type_name,
            css_class="dtype-unknown" if not is_extension else "dtype-extension",
            tooltip="\n".join(tooltip_parts),
            details=details,
            warnings=[f"Unknown type: {full_name}"] if not is_extension else [],
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
        """Register a section formatter."""
        self._section_formatters[formatter.section_name] = formatter

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
        to the fallback formatter if none match.
        """
        for formatter in self._type_formatters:
            try:
                if formatter.can_format(obj):
                    return formatter.format(obj, context)
            except Exception as e:
                # Log but don't fail - try next formatter
                warnings.warn(
                    f"Formatter {type(formatter).__name__} failed for "
                    f"{type(obj).__name__}: {e}",
                    stacklevel=2,
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


# Global registry instance
formatter_registry = FormatterRegistry()


# =============================================================================
# Uns Value Renderer Registry (for custom serialized data visualization)
# =============================================================================

# Type hint key used in uns dicts to indicate custom rendering
UNS_TYPE_HINT_KEY = "__anndata_repr__"


@dataclass
class UnsRendererOutput:
    """Output from an uns renderer."""

    html: str
    """HTML content to display (will be sanitized)"""

    collapsed: bool = True
    """Whether to show collapsed by default"""

    type_label: str = ""
    """Optional type label to show (e.g., 'analysis config')"""


class UnsRendererRegistry:
    """
    Registry for custom uns value renderers.

    This registry allows packages to register custom HTML renderers for
    serialized data stored in uns. The data must contain a type hint
    (key: '__anndata_repr__') that maps to a registered renderer.

    Security model:
    - Data in uns NEVER triggers code execution or imports
    - Packages must be imported by the user and register their renderers
    - Unrecognized type hints fall back to safe JSON/string preview
    - All HTML output is sanitized

    Example usage by a package:

        # In mypackage/__init__.py
        try:
            from anndata._repr import register_uns_renderer, UnsRendererOutput

            def render_analysis_config(value, context):
                '''Render analysis configuration stored as JSON string.'''
                import json
                data = json.loads(value) if isinstance(value, str) else value
                # Return safe HTML
                return UnsRendererOutput(
                    html='<div class="analysis-config">...</div>',
                    type_label="analysis config",
                )

            register_uns_renderer("mypackage.config", render_analysis_config)
        except ImportError:
            pass  # anndata not installed or old version

    Example data structure in uns:

        adata.uns["analysis_config"] = {
            "__anndata_repr__": "mypackage.config",
            "data": '{"params": {...}, "version": "1.0"}'
        }

    Or for simple string values with embedded hint:

        adata.uns["config"] = "__anndata_repr__:mypackage.config::{...json...}"
    """

    def __init__(self) -> None:
        self._renderers: dict[str, Callable[[Any, FormatterContext], UnsRendererOutput]] = {}

    def register(
        self,
        type_hint: str,
        renderer: Callable[[Any, FormatterContext], UnsRendererOutput],
    ) -> None:
        """
        Register a renderer for a type hint.

        Parameters
        ----------
        type_hint
            The type hint string (e.g., "kompot.history", "scanpy.settings")
        renderer
            A callable that takes (value, context) and returns UnsRendererOutput.
            The value is the uns entry (with type hint removed if it was a dict).
        """
        self._renderers[type_hint] = renderer

    def unregister(self, type_hint: str) -> bool:
        """Unregister a renderer. Returns True if found and removed."""
        return self._renderers.pop(type_hint, None) is not None

    def get_renderer(
        self, type_hint: str
    ) -> Callable[[Any, FormatterContext], UnsRendererOutput] | None:
        """Get renderer for a type hint, or None if not registered."""
        return self._renderers.get(type_hint)

    def is_registered(self, type_hint: str) -> bool:
        """Check if a type hint has a registered renderer."""
        return type_hint in self._renderers

    def get_registered_hints(self) -> list[str]:
        """Get list of registered type hints."""
        return list(self._renderers.keys())


# Global uns renderer registry
uns_renderer_registry = UnsRendererRegistry()


def register_uns_renderer(
    type_hint: str,
    renderer: Callable[[Any, FormatterContext], UnsRendererOutput],
) -> None:
    """
    Register a custom renderer for uns values with a specific type hint.

    This is the public API for packages to register custom renderers.
    See UnsRendererRegistry docstring for full documentation and examples.

    Parameters
    ----------
    type_hint
        The type hint string (e.g., "kompot.history")
    renderer
        A callable that takes (value, context) and returns UnsRendererOutput

    Examples
    --------
    >>> from anndata._repr import register_uns_renderer, UnsRendererOutput
    >>>
    >>> def render_my_data(value, context):
    ...     return UnsRendererOutput(
    ...         html='<span class="my-data">Custom view</span>',
    ...         type_label="my custom type",
    ...     )
    >>>
    >>> register_uns_renderer("mypackage.mytype", render_my_data)
    """
    uns_renderer_registry.register(type_hint, renderer)


def extract_uns_type_hint(value: Any) -> tuple[str | None, Any]:
    """
    Extract type hint from an uns value if present.

    Supports two formats:
    1. Dict with __anndata_repr__ key:
       {"__anndata_repr__": "pkg.type", "data": ...}
       Returns ("pkg.type", {"data": ...})

    2. String with prefix:
       "__anndata_repr__:pkg.type::actual content"
       Returns ("pkg.type", "actual content")

    If no type hint found, returns (None, original_value).

    Parameters
    ----------
    value
        The uns value to check

    Returns
    -------
    tuple of (type_hint or None, cleaned_value)
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
        rest = value[len(UNS_TYPE_HINT_KEY) + 1:]  # After "__anndata_repr__:"
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
