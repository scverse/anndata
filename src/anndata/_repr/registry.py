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
                html_content='<span>Custom config preview</span>',
            )
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any


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
    """Custom HTML content (for expandable nested content or type column)"""

    meta_content: str | None = None
    """HTML content for the meta column (data previews, dimensions, etc.)"""

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
        is_extension = module and not module.startswith((
            "anndata",
            "numpy",
            "pandas",
            "scipy",
        ))

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
                if formatter.can_format(obj):
                    return formatter.format(obj, context)
            except Exception as e:  # noqa: BLE001
                # Intentional broad catch: formatters shouldn't crash the entire repr
                from anndata._warnings import warn

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


# Global registry instance
formatter_registry = FormatterRegistry()


# =============================================================================
# Type hint extraction for tagged data in uns
# =============================================================================

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
                    html_content="<span>Custom rendering</span>",
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
