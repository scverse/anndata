"""
Registry pattern tests for the _repr module.

Tests for FormatterRegistry, TypeFormatter, SectionFormatter,
custom formatter registration, and uns type hints.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

from anndata import AnnData

if TYPE_CHECKING:
    from typing import Any


class TestFormatterRegistry:
    """Test the formatter registry pattern for extensibility."""

    def test_registry_has_formatters(self):
        """Test registry contains registered formatters."""
        from anndata._repr.registry import formatter_registry

        assert len(formatter_registry._type_formatters) > 0

    def test_custom_formatter_registration(self):
        """Test registering a custom formatter."""
        from anndata._repr.registry import (
            FormattedOutput,
            FormatterContext,
            TypeFormatter,
            formatter_registry,
        )

        class CustomType:
            pass

        class CustomTypeFormatter(TypeFormatter):
            priority = 500

            def can_format(self, obj: Any) -> bool:
                return isinstance(obj, CustomType)

            def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
                return FormattedOutput(
                    type_name="CustomType",
                    css_class="dtype-custom",
                    is_serializable=False,
                )

        formatter = CustomTypeFormatter()
        formatter_registry.register_type_formatter(formatter)

        try:
            obj = CustomType()
            context = FormatterContext()
            result = formatter_registry.format_value(obj, context)
            assert result.type_name == "CustomType"
            assert result.css_class == "dtype-custom"
            assert result.is_serializable is False
        finally:
            formatter_registry.unregister_type_formatter(formatter)

    def test_fallback_formatter_for_unknown_types(self):
        """Test fallback formatter handles unknown types gracefully."""
        from anndata._repr.registry import FormatterContext, formatter_registry

        class UnknownType:
            pass

        obj = UnknownType()
        context = FormatterContext()
        result = formatter_registry.format_value(obj, context)

        assert result is not None
        assert "UnknownType" in result.type_name
        assert "unknown" in result.css_class or "extension" in result.css_class

    def test_formatter_priority_order(self):
        """Test formatters are checked in priority order."""
        from anndata._repr.registry import formatter_registry

        priorities = [f.priority for f in formatter_registry._type_formatters]
        assert priorities == sorted(priorities, reverse=True)

    def test_formatter_sections_filtering(self):
        """Test formatters are only applied to specified sections."""
        from anndata._repr.registry import (
            FormattedOutput,
            FormatterContext,
            TypeFormatter,
            formatter_registry,
        )

        class SectionSpecificType:
            pass

        class UnsOnlyFormatter(TypeFormatter):
            priority = 600
            sections = ("uns",)

            def can_format(self, obj: Any) -> bool:
                return isinstance(obj, SectionSpecificType)

            def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
                return FormattedOutput(
                    type_name="UnsSpecificType",
                    css_class="dtype-uns-specific",
                )

        formatter = UnsOnlyFormatter()
        formatter_registry.register_type_formatter(formatter)

        try:
            obj = SectionSpecificType()

            context_uns = FormatterContext(section="uns")
            result_uns = formatter_registry.format_value(obj, context_uns)
            assert result_uns.type_name == "UnsSpecificType"

            context_obsm = FormatterContext(section="obsm")
            result_obsm = formatter_registry.format_value(obj, context_obsm)
            assert result_obsm.type_name != "UnsSpecificType"
            assert "SectionSpecificType" in result_obsm.type_name
        finally:
            formatter_registry.unregister_type_formatter(formatter)

    def test_formatter_sections_none_applies_everywhere(self):
        """Test formatters with sections=None apply to all sections."""
        from anndata._repr.registry import (
            FormattedOutput,
            FormatterContext,
            TypeFormatter,
            formatter_registry,
        )

        class UniversalType:
            pass

        class UniversalFormatter(TypeFormatter):
            priority = 600
            sections = None

            def can_format(self, obj: Any) -> bool:
                return isinstance(obj, UniversalType)

            def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
                return FormattedOutput(type_name="UniversalType")

        formatter = UniversalFormatter()
        formatter_registry.register_type_formatter(formatter)

        try:
            obj = UniversalType()
            for section in ["uns", "obsm", "varm", "layers", "obs", "var"]:
                context = FormatterContext(section=section)
                result = formatter_registry.format_value(obj, context)
                assert result.type_name == "UniversalType"
        finally:
            formatter_registry.unregister_type_formatter(formatter)

    def test_extension_type_graceful_handling(self):
        """Test extension types are handled gracefully."""
        from anndata._repr.registry import FormatterContext, formatter_registry

        class ExtensionData:
            def __init__(self):
                self.n_obs = 100
                self.n_vars = 50
                self.shape = (100, 50)
                self.dtype = np.float32

        ExtensionData.__module__ = "treedata.core"

        obj = ExtensionData()
        context = FormatterContext()
        result = formatter_registry.format_value(obj, context)

        assert result is not None
        assert "ExtensionData" in result.type_name
        assert "Shape:" in result.tooltip  # Shape info in tooltip for fallback

    def test_anndata_in_uns_detected(self):
        """Test nested AnnData in .uns is properly detected."""
        inner = AnnData(np.zeros((5, 3)))
        outer = AnnData(np.zeros((10, 5)))
        outer.uns["inner_adata"] = inner

        html = outer._repr_html_()

        assert "inner_adata" in html
        assert "AnnData" in html

    def test_registry_formatter_exception_continues(self):
        """Test registry continues to next formatter on exception."""
        from anndata._repr.registry import (
            FormattedOutput,
            FormatterContext,
            FormatterRegistry,
            TypeFormatter,
        )

        class FailingFormatter(TypeFormatter):
            priority = 1000

            def can_format(self, obj):
                return True

            def format(self, obj, context):
                msg = "Intentional failure"
                raise RuntimeError(msg)

        class BackupFormatter(TypeFormatter):
            priority = 500

            def can_format(self, obj):
                return True

            def format(self, obj, context):
                return FormattedOutput(type_name="Backup", css_class="backup")

        registry = FormatterRegistry()
        failing = FailingFormatter()
        backup = BackupFormatter()
        registry.register_type_formatter(failing)
        registry.register_type_formatter(backup)

        try:
            with pytest.warns(UserWarning, match="Formatter.*failed"):
                result = registry.format_value("test", FormatterContext())
            assert result.type_name == "Backup"
        finally:
            registry.unregister_type_formatter(failing)
            registry.unregister_type_formatter(backup)

    def test_register_formatter_decorator_with_class(self):
        """Test register_formatter works as decorator with class."""
        from anndata._repr.registry import (
            FormattedOutput,
            FormatterContext,
            TypeFormatter,
            formatter_registry,
            register_formatter,
        )

        class DecoratorTestFormatter(TypeFormatter):
            priority = 999

            def can_format(self, obj):
                return (
                    isinstance(obj, tuple)
                    and len(obj) == 3
                    and obj[0] == "decorator_test"
                )

            def format(self, obj, context):
                return FormattedOutput(type_name="DecoratorTest", css_class="test")

        formatter_instance = DecoratorTestFormatter()
        register_formatter(formatter_instance)

        try:
            result = formatter_registry.format_value(
                ("decorator_test", 1, 2), FormatterContext()
            )
            assert result.type_name == "DecoratorTest"
        finally:
            formatter_registry.unregister_type_formatter(formatter_instance)


class TestFormatterContext:
    """Tests for FormatterContext."""

    def test_context_child_creates_nested_context(self):
        """Test FormatterContext.child() creates proper nested context."""
        from anndata._repr.registry import FormatterContext

        parent = FormatterContext(
            depth=0,
            max_depth=5,
            parent_keys=(),
            adata_ref=None,
            section="uns",
        )

        child = parent.child("nested_key")

        assert child.depth == 1
        assert child.max_depth == 5
        assert child.parent_keys == ("nested_key",)
        assert child.section == "uns"

        grandchild = child.child("deeper")
        assert grandchild.depth == 2
        assert grandchild.parent_keys == ("nested_key", "deeper")

    def test_context_access_path_empty(self):
        """Test access_path returns empty string for no parent keys."""
        from anndata._repr.registry import FormatterContext

        context = FormatterContext(parent_keys=())
        assert context.access_path == ""

    def test_context_access_path_identifier_keys(self):
        """Test access_path with valid Python identifiers."""
        from anndata._repr.registry import FormatterContext

        context = FormatterContext(parent_keys=("uns", "neighbors", "params"))
        assert context.access_path == ".uns.neighbors.params"

    def test_context_access_path_non_identifier_keys(self):
        """Test access_path with non-identifier keys."""
        from anndata._repr.registry import FormatterContext

        context = FormatterContext(parent_keys=("uns", "key with spaces", "123numeric"))
        path = context.access_path
        assert ".uns" in path
        assert "['key with spaces']" in path
        assert "['123numeric']" in path

    def test_context_access_path_mixed_keys(self):
        """Test access_path with mixed identifier and non-identifier keys."""
        from anndata._repr.registry import FormatterContext

        context = FormatterContext(parent_keys=("valid", "has-hyphen", "also_valid"))
        path = context.access_path
        assert ".valid" in path
        assert "['has-hyphen']" in path
        assert ".also_valid" in path


class TestSectionFormatter:
    """Tests for SectionFormatter abstract class."""

    def test_section_formatter_default_methods(self):
        """Test SectionFormatter default method implementations."""
        from anndata._repr.registry import SectionFormatter

        class TestSectionFormatter(SectionFormatter):
            @property
            def section_name(self) -> str:
                return "test_section"

            def get_entries(self, obj, context):
                return []

        formatter = TestSectionFormatter()

        assert formatter.display_name == "test_section"
        assert formatter.doc_url is None
        assert formatter.tooltip == ""
        assert formatter.should_show(None) is True


class TestFallbackFormatter:
    """Tests for FallbackFormatter edge cases."""

    def test_fallback_extension_type_no_warning(self):
        """Test fallback formatter for extension types doesn't add warning."""
        from anndata._repr.registry import FallbackFormatter, FormatterContext

        class ExtensionType:
            pass

        ExtensionType.__module__ = "treedata.core"

        formatter = FallbackFormatter()
        obj = ExtensionType()
        context = FormatterContext()

        result = formatter.format(obj, context)

        assert result.type_name == "ExtensionType"
        assert result.css_class == "dtype-extension"
        assert len(result.warnings) == 0

    def test_fallback_with_shape_and_dtype(self):
        """Test fallback formatter extracts shape and dtype."""
        from anndata._repr.registry import FallbackFormatter, FormatterContext

        class ShapedType:
            shape = (10, 5)
            dtype = "float32"

        formatter = FallbackFormatter()
        obj = ShapedType()
        context = FormatterContext()

        result = formatter.format(obj, context)

        assert "Shape: (10, 5)" in result.tooltip
        assert "Dtype: float32" in result.tooltip

    def test_fallback_with_len(self):
        """Test fallback formatter extracts length."""
        from anndata._repr.registry import FallbackFormatter, FormatterContext

        class LengthType:
            def __len__(self):
                return 42

        formatter = FallbackFormatter()
        obj = LengthType()
        context = FormatterContext()

        result = formatter.format(obj, context)

        assert "Length: 42" in result.tooltip

    def test_fallback_len_raises_error(self):
        """Test fallback formatter handles __len__ errors gracefully."""
        from anndata._repr.registry import FallbackFormatter, FormatterContext

        class BrokenLenType:
            def __len__(self):
                msg = "Cannot get length"
                raise TypeError(msg)

        formatter = FallbackFormatter()
        obj = BrokenLenType()
        context = FormatterContext()

        result = formatter.format(obj, context)
        assert "Length:" not in result.tooltip


class TestRegistryAbstractMethods:
    """Tests for registry abstract method checks."""

    def test_type_formatter_is_abstract(self):
        """Test TypeFormatter requires abstract methods."""
        from anndata._repr.registry import TypeFormatter

        # Should not be instantiable directly
        assert hasattr(TypeFormatter, "can_format")
        assert hasattr(TypeFormatter, "format")

    def test_section_formatter_is_abstract(self):
        """Test SectionFormatter requires abstract methods."""
        from anndata._repr.registry import SectionFormatter

        assert hasattr(SectionFormatter, "section_name")
        assert hasattr(SectionFormatter, "get_entries")


class TestUnsRendererRegistry:
    """Test the uns renderer registry for custom serialized data visualization."""

    def test_extract_type_hint_dict_format(self):
        """Test extracting type hint from dict format."""
        from anndata._repr.registry import UNS_TYPE_HINT_KEY, extract_uns_type_hint

        value = {
            UNS_TYPE_HINT_KEY: "mypackage.config",
            "data": {"setting": "value"},
            "version": "1.0",
        }
        hint, cleaned = extract_uns_type_hint(value)

        assert hint == "mypackage.config"
        assert UNS_TYPE_HINT_KEY not in cleaned
        assert cleaned["data"] == {"setting": "value"}
        assert cleaned["version"] == "1.0"

    def test_extract_type_hint_string_format(self):
        """Test extracting type hint from string prefix format."""
        from anndata._repr.registry import UNS_TYPE_HINT_KEY, extract_uns_type_hint

        value = f'{UNS_TYPE_HINT_KEY}:mypackage.config::{{"setting": "value"}}'
        hint, cleaned = extract_uns_type_hint(value)

        assert hint == "mypackage.config"
        assert cleaned == '{"setting": "value"}'

    def test_extract_type_hint_no_hint_returns_none(self):
        """Test that values without type hints return (None, original_value)."""
        from anndata._repr.registry import extract_uns_type_hint

        value = {"data": "value"}
        hint, cleaned = extract_uns_type_hint(value)
        assert hint is None
        assert cleaned == value

        value = "just a string"
        hint, cleaned = extract_uns_type_hint(value)
        assert hint is None
        assert cleaned == value

        value = 42
        hint, cleaned = extract_uns_type_hint(value)
        assert hint is None
        assert cleaned == value

    def test_extract_type_hint_invalid_dict_hint_type(self):
        """Test that non-string type hints in dict are ignored."""
        from anndata._repr.registry import UNS_TYPE_HINT_KEY, extract_uns_type_hint

        value = {UNS_TYPE_HINT_KEY: 123, "data": "value"}
        hint, cleaned = extract_uns_type_hint(value)
        assert hint is None
        assert cleaned == value

    def test_extract_type_hint_malformed_string_format(self):
        """Test that malformed string format returns no hint."""
        from anndata._repr.registry import UNS_TYPE_HINT_KEY, extract_uns_type_hint

        value = f"{UNS_TYPE_HINT_KEY}:mypackage.config:data"
        hint, cleaned = extract_uns_type_hint(value)
        assert hint is None
        assert cleaned == value

    def test_type_formatter_for_tagged_uns_data(self):
        """Test using TypeFormatter to handle tagged data in uns."""
        from anndata._repr import (
            FormattedOutput,
            TypeFormatter,
            extract_uns_type_hint,
            formatter_registry,
            register_formatter,
        )

        class TestConfigFormatter(TypeFormatter):
            priority = 100

            def can_format(self, obj):
                hint, _ = extract_uns_type_hint(obj)
                return hint == "test.config_format"

            def format(self, obj, context):
                _hint, data = extract_uns_type_hint(obj)
                items = data.get("data", {})
                return FormattedOutput(
                    type_name="test config",
                    preview_html=f'<span class="test-custom">Items: {len(items)}</span>',
                )

        formatter = TestConfigFormatter()
        register_formatter(formatter)

        try:
            adata = AnnData(np.zeros((5, 3)))
            adata.uns["my_config"] = {
                "__anndata_repr__": "test.config_format",
                "data": {"a": 1, "b": 2, "c": 3},
            }

            html = adata._repr_html_()

            assert "Items: 3" in html
            assert "test config" in html
        finally:
            formatter_registry.unregister_type_formatter(formatter)

    def test_unregistered_type_hint_shows_import_message(self):
        """Test that unregistered type hints show helpful import message."""
        adata = AnnData(np.zeros((5, 3)))
        adata.uns["external_data"] = {
            "__anndata_repr__": "externalpackage.customtype",
            "data": {"key": "value"},
        }

        html = adata._repr_html_()

        assert "externalpackage.customtype" in html
        assert "import externalpackage" in html

    def test_formatter_error_handled_gracefully(self):
        """Test that TypeFormatter errors don't crash the repr."""
        from anndata._repr import (
            TypeFormatter,
            extract_uns_type_hint,
            formatter_registry,
            register_formatter,
        )

        class FailingFormatter(TypeFormatter):
            priority = 100

            def can_format(self, obj):
                hint, _ = extract_uns_type_hint(obj)
                return hint == "test.failing_format"

            def format(self, obj, context):
                msg = "Intentional test error"
                raise ValueError(msg)

        formatter = FailingFormatter()
        register_formatter(formatter)

        try:
            adata = AnnData(np.zeros((5, 3)))
            adata.uns["will_fail"] = {
                "__anndata_repr__": "test.failing_format",
                "data": "test",
            }

            with pytest.warns(UserWarning, match="Formatter.*failed"):
                html = adata._repr_html_()

            assert html is not None
            assert "will_fail" in html
        finally:
            formatter_registry.unregister_type_formatter(formatter)

    def test_string_format_type_hint_in_html(self):
        """Test string format type hints work in HTML output."""
        adata = AnnData(np.zeros((5, 3)))
        adata.uns["string_hint"] = (
            "__anndata_repr__:somepackage.config::actual content here"
        )

        html = adata._repr_html_()

        assert "somepackage.config" in html
        assert "import somepackage" in html

    def test_type_hint_key_constant_exported(self):
        """Test that UNS_TYPE_HINT_KEY constant is properly exported."""
        from anndata._repr import UNS_TYPE_HINT_KEY

        assert UNS_TYPE_HINT_KEY == "__anndata_repr__"

    def test_security_data_never_triggers_import(self):
        """Test that data in uns NEVER triggers imports or code execution."""
        import sys

        fake_package = "definitely_not_a_real_package_12345"
        assert fake_package not in sys.modules

        adata = AnnData(np.zeros((5, 3)))
        adata.uns["malicious"] = {
            "__anndata_repr__": f"{fake_package}.evil",
            "data": "some data",
        }

        html = adata._repr_html_()

        assert fake_package not in sys.modules
        assert html is not None
        assert "malicious" in html


class TestCustomSectionFormatterCodePaths:
    """Tests for custom section formatter code paths to improve coverage."""

    def test_custom_section_with_entries(self):
        """Test custom section formatter with entries."""
        from anndata._repr.registry import (
            FormattedEntry,
            FormattedOutput,
            FormatterContext,  # noqa: TC001
            SectionFormatter,
            formatter_registry,
        )

        class TestSectionFormatter(SectionFormatter):
            @property
            def section_name(self) -> str:
                return "test_custom_section_entries"

            @property
            def display_name(self) -> str:
                return "Test Custom Entries"

            @property
            def doc_url(self) -> str:
                return "https://example.com/docs"

            @property
            def tooltip(self) -> str:
                return "A test section"

            def should_show(self, obj) -> bool:
                return hasattr(obj, "_test_marker")

            def get_entries(self, obj, context: FormatterContext):
                return [
                    FormattedEntry(
                        key="entry1",
                        output=FormattedOutput(type_name="TestType", css_class="test"),
                    ),
                    FormattedEntry(
                        key="entry2",
                        output=FormattedOutput(
                            type_name="TestType2",
                            css_class="test",
                            warnings=["Test warning"],
                        ),
                    ),
                ]

        formatter = TestSectionFormatter()
        formatter_registry.register_section_formatter(formatter)

        try:
            adata = AnnData(np.zeros((10, 5)))
            adata._test_marker = True

            html = adata._repr_html_()
            assert (
                "Test Custom Entries" in html or "test_custom_section_entries" in html
            )
            assert "entry1" in html
            assert "entry2" in html
        finally:
            formatter_registry._section_formatters.pop(
                "test_custom_section_entries", None
            )

    def test_custom_section_exception_handling(self):
        """Test custom section formatter handles exceptions gracefully."""
        from anndata._repr.registry import (
            FormatterContext,  # noqa: TC001
            SectionFormatter,
            formatter_registry,
        )

        class FailingSectionFormatter(SectionFormatter):
            @property
            def section_name(self) -> str:
                return "failing_section"

            def should_show(self, obj) -> bool:
                return True

            def get_entries(self, obj, context: FormatterContext):
                msg = "Intentional failure"
                raise ValueError(msg)

        formatter = FailingSectionFormatter()
        formatter_registry.register_section_formatter(formatter)

        try:
            adata = AnnData(np.zeros((10, 5)))
            # Should not crash, just skip the failing section
            with pytest.warns(UserWarning, match="Custom section.*failed"):
                html = adata._repr_html_()
            assert html is not None
        finally:
            formatter_registry._section_formatters.pop("failing_section", None)

    def test_custom_section_should_show_exception(self):
        """Test custom section handles should_show exception."""
        from anndata._repr.registry import (
            FormatterContext,  # noqa: TC001
            SectionFormatter,
            formatter_registry,
        )

        class ShouldShowFailingFormatter(SectionFormatter):
            @property
            def section_name(self) -> str:
                return "should_show_failing"

            def should_show(self, obj) -> bool:
                msg = "Intentional failure in should_show"
                raise RuntimeError(msg)

            def get_entries(self, obj, context: FormatterContext):
                return []

        formatter = ShouldShowFailingFormatter()
        formatter_registry.register_section_formatter(formatter)

        try:
            adata = AnnData(np.zeros((10, 5)))
            # Should not crash
            html = adata._repr_html_()
            assert html is not None
        finally:
            formatter_registry._section_formatters.pop("should_show_failing", None)


class TestFormattedEntryRendering:
    """Tests for FormattedEntry rendering in custom sections."""

    def test_formatted_entry_with_expandable_html(self):
        """Test formatted entry with expandable HTML content."""
        from anndata._repr.registry import (
            FormattedEntry,
            FormattedOutput,
            FormatterContext,  # noqa: TC001
            SectionFormatter,
            formatter_registry,
        )

        class ExpandableEntryFormatter(SectionFormatter):
            @property
            def section_name(self) -> str:
                return "expandable_section"

            def should_show(self, obj) -> bool:
                return hasattr(obj, "_expandable_marker")

            def get_entries(self, obj, context: FormatterContext):
                return [
                    FormattedEntry(
                        key="expandable_entry",
                        output=FormattedOutput(
                            type_name="Expandable",
                            css_class="test",
                            expanded_html="<div>Expanded content here</div>",
                        ),
                    ),
                ]

        formatter = ExpandableEntryFormatter()
        formatter_registry.register_section_formatter(formatter)

        try:
            adata = AnnData(np.zeros((10, 5)))
            adata._expandable_marker = True

            html = adata._repr_html_()
            assert "expandable_entry" in html
            assert "Expand" in html or "expand" in html.lower()
        finally:
            formatter_registry._section_formatters.pop("expandable_section", None)

    def test_formatted_entry_with_inline_html(self):
        """Test formatted entry with inline (non-expandable) HTML."""
        from anndata._repr.registry import (
            FormattedEntry,
            FormattedOutput,
            FormatterContext,  # noqa: TC001
            SectionFormatter,
            formatter_registry,
        )

        class InlineEntryFormatter(SectionFormatter):
            @property
            def section_name(self) -> str:
                return "inline_section"

            def should_show(self, obj) -> bool:
                return hasattr(obj, "_inline_marker")

            def get_entries(self, obj, context: FormatterContext):
                return [
                    FormattedEntry(
                        key="inline_entry",
                        output=FormattedOutput(
                            type_name="Inline",
                            css_class="test",
                            preview_html="<span class='preview'>Preview content</span>",
                        ),
                    ),
                ]

        formatter = InlineEntryFormatter()
        formatter_registry.register_section_formatter(formatter)

        try:
            adata = AnnData(np.zeros((10, 5)))
            adata._inline_marker = True

            html = adata._repr_html_()
            assert "inline_entry" in html
            assert "Preview content" in html
        finally:
            formatter_registry._section_formatters.pop("inline_section", None)

    def test_formatted_entry_not_serializable(self):
        """Test formatted entry with not serializable marker."""
        from anndata._repr.registry import (
            FormattedEntry,
            FormattedOutput,
            FormatterContext,  # noqa: TC001
            SectionFormatter,
            formatter_registry,
        )

        class NotSerializableFormatter(SectionFormatter):
            @property
            def section_name(self) -> str:
                return "not_serializable_section"

            def should_show(self, obj) -> bool:
                return hasattr(obj, "_not_serializable_marker")

            def get_entries(self, obj, context: FormatterContext):
                return [
                    FormattedEntry(
                        key="bad_entry",
                        output=FormattedOutput(
                            type_name="BadType",
                            css_class="test",
                            is_serializable=False,
                            warnings=["Cannot be saved"],
                        ),
                    ),
                ]

        formatter = NotSerializableFormatter()
        formatter_registry.register_section_formatter(formatter)

        try:
            adata = AnnData(np.zeros((10, 5)))
            adata._not_serializable_marker = True

            html = adata._repr_html_()
            assert "bad_entry" in html
            assert "⚠" in html or "warning" in html.lower()
        finally:
            formatter_registry._section_formatters.pop("not_serializable_section", None)
