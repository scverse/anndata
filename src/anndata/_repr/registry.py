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
