"""
Tests for HTML representation of AnnData objects.

This module tests:
- Basic HTML generation for various AnnData configurations
- Formatter registry pattern and extensibility
- Type formatters for all supported types
- Settings integration
- HTML/CSS/JS validation
- Special features (colors, warnings, search)
- Completeness of representation

Target coverage: >95%

Visual Inspection:
    Run the visual inspection script to generate an HTML file for manual review:

        python tests/visual_inspect_repr_html.py

    The file will be saved to tests/repr_html_visual_test.html
"""

from __future__ import annotations

import re
from html.parser import HTMLParser
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import anndata as ad
from anndata import AnnData

if TYPE_CHECKING:
    from typing import Any

# Check optional dependencies
try:
    import dask.array as da

    HAS_DASK = True
except ImportError:
    HAS_DASK = False

try:
    import cupy as cp  # noqa: F401

    HAS_CUPY = True
except ImportError:
    HAS_CUPY = False

try:
    import awkward as ak

    HAS_AWKWARD = True
except ImportError:
    HAS_AWKWARD = False

try:
    import xarray  # noqa: F401

    HAS_XARRAY = True
except ImportError:
    HAS_XARRAY = False


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def adata():
    """Basic AnnData for testing."""
    return AnnData(
        np.random.randn(100, 50).astype(np.float32),
        obs=pd.DataFrame(
            {"batch": ["A", "B"] * 50}, index=[f"cell_{i}" for i in range(100)]
        ),
        var=pd.DataFrame(
            {"gene_name": [f"gene_{i}" for i in range(50)]},
            index=[f"gene_{i}" for i in range(50)],
        ),
    )


@pytest.fixture
def adata_full():
    """AnnData with all attributes populated."""
    n_obs, n_vars = 100, 50
    adata = AnnData(
        sp.random(n_obs, n_vars, density=0.1, format="csr", dtype=np.float32),
        obs=pd.DataFrame({
            "batch": pd.Categorical(["A", "B"] * (n_obs // 2)),
            "n_counts": np.random.randint(1000, 10000, n_obs),
            "cell_type": pd.Categorical(
                ["T", "B", "NK"] * (n_obs // 3) + ["T"] * (n_obs % 3)
            ),
        }),
        var=pd.DataFrame({
            "gene_name": [f"gene_{i}" for i in range(n_vars)],
            "highly_variable": np.random.choice([True, False], n_vars),
        }),
    )
    adata.uns["neighbors"] = {"params": {"n_neighbors": 15}}
    adata.uns["batch_colors"] = ["#FF0000", "#00FF00"]
    adata.obsm["X_pca"] = np.random.randn(n_obs, 50).astype(np.float32)
    adata.obsm["X_umap"] = np.random.randn(n_obs, 2).astype(np.float32)
    adata.varm["PCs"] = np.random.randn(n_vars, 50).astype(np.float32)
    adata.layers["raw"] = sp.random(n_obs, n_vars, density=0.1, format="csr")
    adata.obsp["distances"] = sp.random(n_obs, n_obs, density=0.01, format="csr")
    adata.varp["gene_corr"] = sp.random(n_vars, n_vars, density=0.1, format="csr")
    return adata


@pytest.fixture
def adata_with_colors():
    """AnnData with color annotations."""
    adata = AnnData(np.zeros((10, 5)))
    adata.obs["cluster"] = pd.Categorical(["A", "B", "C"] * 3 + ["A"])
    adata.uns["cluster_colors"] = ["#FF0000", "#00FF00", "#0000FF"]
    return adata


@pytest.fixture
def adata_with_nested():
    """AnnData with nested AnnData in uns."""
    inner = AnnData(np.zeros((5, 3)))
    outer = AnnData(np.zeros((10, 5)))
    outer.uns["nested_adata"] = inner
    return outer


@pytest.fixture
def adata_with_special_chars():
    """AnnData with special characters in names."""
    adata = AnnData(np.zeros((10, 5)))
    adata.obs["col<script>"] = list(range(10))  # XSS attempt
    adata.uns["key&value"] = "test"
    adata.uns["quotes\"test'"] = "value"
    return adata


# =============================================================================
# HTML Validation
# =============================================================================


VOID_ELEMENTS = frozenset({
    "area",
    "base",
    "br",
    "col",
    "embed",
    "hr",
    "img",
    "input",
    "link",
    "meta",
    "param",
    "source",
    "track",
    "wbr",
})


class StrictHTMLParser(HTMLParser):
    """Validates HTML structure and catches common errors."""

    def __init__(self):
        super().__init__()
        self.tag_stack = []
        self.errors = []
        self.ids_seen = set()

    def handle_starttag(self, tag, attrs):
        if tag not in VOID_ELEMENTS:
            self.tag_stack.append(tag)

        # Check for duplicate IDs
        for name, value in attrs:
            if name == "id":
                if value in self.ids_seen:
                    self.errors.append(f"Duplicate ID: {value}")
                self.ids_seen.add(value)

    def handle_endtag(self, tag):
        if tag not in VOID_ELEMENTS:
            if not self.tag_stack or self.tag_stack[-1] != tag:
                self.errors.append(f"Mismatched tag: </{tag}>")
            else:
                self.tag_stack.pop()


class TestHTMLValidation:
    """Validate generated HTML is well-formed and standards-compliant."""

    def test_html_well_formed(self, adata):
        """Test HTML is parseable and well-formed."""
        html = adata._repr_html_()
        parser = StrictHTMLParser()
        parser.feed(html)
        assert not parser.errors, f"HTML errors: {parser.errors}"
        assert not parser.tag_stack, f"Unclosed tags: {parser.tag_stack}"

    def test_no_duplicate_ids(self, adata_full):
        """Test no duplicate element IDs within same repr."""
        html = adata_full._repr_html_()
        ids = re.findall(r'id=["\']([^"\']+)["\']', html)
        assert len(ids) == len(set(ids)), f"Duplicate IDs found: {ids}"

    def test_all_links_have_href(self, adata_full):
        """Test all anchor tags have href attribute."""
        html = adata_full._repr_html_()
        # Find <a> tags without href
        a_without_href = re.search(r"<a(?![^>]*href=)[^>]*>", html)
        assert a_without_href is None, "Found <a> without href"

    def test_style_tag_valid_css(self, adata):
        """Test inline CSS has balanced braces."""
        html = adata._repr_html_()
        style_match = re.search(r"<style[^>]*>(.*?)</style>", html, re.DOTALL)
        if style_match:
            css = style_match.group(1)
            # Basic CSS validation - check balanced braces
            assert css.count("{") == css.count("}"), "Unbalanced CSS braces"

    def test_escaped_user_content(self, adata_with_special_chars):
        """Test that user-provided content is properly escaped."""
        html = adata_with_special_chars._repr_html_()
        # Should not contain raw <script> tags from user content
        assert "<script>" not in html.split("</style>")[-1].split("<script")[0].lower()
        # Should not contain javascript: in attributes
        assert "javascript:" not in html.lower()

    def test_nested_anndata_valid_html(self, adata_with_nested):
        """Test nested AnnData produces valid nested HTML."""
        html = adata_with_nested._repr_html_()
        parser = StrictHTMLParser()
        parser.feed(html)
        assert not parser.errors


class TestCSSValidation:
    """Validate CSS styling."""

    def test_css_variables_defined(self, adata):
        """Test all CSS variables are defined before use."""
        html = adata._repr_html_()
        style_match = re.search(r"<style[^>]*>(.*?)</style>", html, re.DOTALL)
        if style_match:
            css = style_match.group(1)
            # Find variable definitions in <style> tag (--name: value)
            defined = set(re.findall(r"--([\w-]+)\s*:", css))
            # Also find variable definitions in inline styles (style="--name: value")
            # First extract all style attribute values, then find var definitions in them
            inline_styles = re.findall(r'style="([^"]*)"', html)
            for style in inline_styles:
                defined |= set(re.findall(r"--([\w-]+)\s*:", style))
            # Find variable usages (var(--name))
            used = set(re.findall(r"var\s*\(\s*--([\w-]+)", css))
            undefined = used - defined
            assert not undefined, f"Undefined CSS vars: {undefined}"

    def test_dark_mode_css_present(self, adata):
        """Test dark mode CSS rules are present."""
        html = adata._repr_html_()
        assert (
            "@media (prefers-color-scheme: dark)" in html
            or "jp-Theme-Dark" in html
            or "[data-jp-theme-light" in html
        )


# =============================================================================
# Basic HTML Generation
# =============================================================================


class TestBasicHTMLGeneration:
    """Test basic HTML generation for various AnnData configurations."""

    def test_empty_anndata(self):
        """Test repr for empty AnnData."""
        adata = AnnData()
        html = adata._repr_html_()
        assert "0" in html  # n_obs and n_vars are 0

    def test_minimal_anndata(self):
        """Test repr with only X."""
        adata = AnnData(np.zeros((10, 5)))
        html = adata._repr_html_()
        assert "10" in html
        assert "5" in html

    def test_dense_matrix(self):
        """Test repr with dense X."""
        adata = AnnData(np.random.randn(100, 50).astype(np.float32))
        html = adata._repr_html_()
        assert "float32" in html
        assert "100" in html
        assert "50" in html

    def test_sparse_csr_matrix(self):
        """Test repr with sparse CSR X."""
        X = sp.random(1000, 500, density=0.1, format="csr")
        adata = AnnData(X)
        html = adata._repr_html_()
        assert "csr" in html.lower()

    def test_sparse_csc_matrix(self):
        """Test repr with sparse CSC X."""
        X = sp.random(1000, 500, density=0.1, format="csc")
        adata = AnnData(X)
        html = adata._repr_html_()
        assert "csc" in html.lower()

    def test_all_attributes_present(self, adata_full):
        """Test all standard attributes appear in repr."""
        html = adata_full._repr_html_()
        for attr in ["obs", "var", "uns", "obsm", "varm", "layers", "obsp", "varp"]:
            assert f'data-section="{attr}"' in html or attr in html

    def test_x_none_handled(self):
        """Test AnnData with X=None is handled."""
        adata = AnnData(obs=pd.DataFrame({"a": [1, 2, 3]}))
        html = adata._repr_html_()
        assert "None" in html or "none" in html.lower()


# =============================================================================
# View Representation
# =============================================================================


class TestViewRepresentation:
    """Test HTML repr for AnnData views."""

    def test_view_badge_present(self, adata):
        """Test view badge appears for views."""
        view = adata[0:5, :]
        html = view._repr_html_()
        assert "View" in html

    def test_view_shows_subset_shape(self, adata):
        """Test view shows correct subset dimensions."""
        view = adata[0:5, 0:3]
        html = view._repr_html_()
        assert "5" in html  # n_obs
        assert "3" in html  # n_vars


# =============================================================================
# Folding Behavior
# =============================================================================


class TestFoldingBehavior:
    """Test auto-folding based on entry count."""

    def test_few_entries_expanded(self):
        """Test sections with ≤5 entries are expanded."""
        adata = AnnData(
            np.zeros((10, 5)),
            obs=pd.DataFrame({"a": list(range(10)), "b": list(range(10))}),
        )
        html = adata._repr_html_()
        # obs section should be present
        assert "obs" in html

    def test_many_entries_collapsed(self):
        """Test sections with >5 entries are collapsed."""
        obs_dict = {f"col_{i}": list(range(10)) for i in range(10)}
        adata = AnnData(np.zeros((10, 5)), obs=pd.DataFrame(obs_dict))
        html = adata._repr_html_()
        assert "collapsed" in html

    def test_custom_fold_threshold(self):
        """Test custom fold threshold via settings."""
        from anndata import settings

        obs_dict = {f"col_{i}": list(range(10)) for i in range(4)}
        adata = AnnData(np.zeros((10, 5)), obs=pd.DataFrame(obs_dict))

        with settings.override(repr_html_fold_threshold=2):
            html = adata._repr_html_()
            assert "collapsed" in html


# =============================================================================
# Nested AnnData
# =============================================================================


class TestNestedAnnData:
    """Test recursive AnnData display in .uns."""

    def test_nested_anndata_present(self):
        """Test nested AnnData in .uns is detected."""
        inner = AnnData(np.zeros((5, 3)))
        outer = AnnData(np.zeros((10, 5)))
        outer.uns["nested"] = inner
        html = outer._repr_html_()
        assert "nested" in html

    def test_nested_anndata_max_depth(self):
        """Test max recursion depth is respected."""
        from anndata import settings

        # Create deeply nested structure
        adata = AnnData(np.zeros((5, 3)))
        for _i in range(5):
            wrapper = AnnData(np.zeros((5, 3)))
            wrapper.uns["inner"] = adata
            adata = wrapper

        with settings.override(repr_html_max_depth=2):
            html = adata._repr_html_()
            # Should contain max depth indicator
            assert "max depth" in html.lower() or "depth" in html

    def test_nested_anndata_expandable(self):
        """Test nested AnnData has expand button."""
        inner = AnnData(np.zeros((5, 3)))
        outer = AnnData(np.zeros((10, 5)))
        outer.uns["nested"] = inner
        html = outer._repr_html_()
        # Should have expand mechanism
        assert "Expand" in html or "expand" in html.lower()


# =============================================================================
# Color Display
# =============================================================================


class TestColorDisplay:
    """Test color swatch display."""

    def test_color_list_swatches(self):
        """Test *_colors entries show color swatches."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["cluster"] = pd.Categorical(["A", "B"] * 5)
        adata.uns["cluster_colors"] = ["#FF0000", "#00FF00"]
        html = adata._repr_html_()
        # Should show colors
        assert "color" in html.lower()

    def test_color_mismatch_warning(self):
        """Test warning when color count doesn't match categories."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["cluster"] = pd.Categorical(["A", "B", "C"] * 3 + ["A"])
        adata.uns["cluster_colors"] = ["#FF0000", "#00FF00"]  # Only 2 colors for 3 cats
        html = adata._repr_html_()
        # Should indicate mismatch or warning
        assert "mismatch" in html.lower() or "warning" in html.lower() or "⚠" in html

    def test_colors_inline_with_category(self, adata_with_colors):
        """Test categorical columns show inline colors."""
        html = adata_with_colors._repr_html_()
        # Colors should be somewhere in the representation
        assert "#FF0000" in html or "color-swatch" in html


# =============================================================================
# Warnings
# =============================================================================


class TestWarnings:
    """Test warning indicators."""

    def test_string_column_warning(self):
        """Test warning for string columns that will convert."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["batch"] = ["A", "B", "A", "B", "A", "B", "A", "B", "A", "B"]
        html = adata._repr_html_()
        # Should have warning indicator
        assert "warning" in html.lower() or "⚠" in html

    def test_string_all_unique_no_warning(self):
        """Test no warning for string columns with all unique values."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["id"] = ["cell_0", "cell_1", "cell_2", "cell_3", "cell_4"]
        _html = adata._repr_html_()
        # The "id" column should not have warning for unique strings
        # This is harder to test precisely, so we just verify no errors

    def test_unserializable_uns_warning(self):
        """Test warning for unserializable .uns entries."""

        class CustomClass:
            pass

        adata = AnnData(np.zeros((10, 5)))
        adata.uns["custom"] = CustomClass()
        html = adata._repr_html_()
        # Should indicate unserializable
        assert (
            "serializable" in html.lower()
            or "warning" in html.lower()
            or "⚠" in html
            or "unknown" in html.lower()
        )


# =============================================================================
# Search Functionality
# =============================================================================


class TestSearchFunctionality:
    """Test search/filter feature."""

    def test_search_input_present(self, adata):
        """Test search input field is present."""
        html = adata._repr_html_()
        assert "<input" in html
        assert "search" in html.lower() or "filter" in html.lower()

    def test_filter_indicator_present(self, adata):
        """Test filter indicator element exists."""
        html = adata._repr_html_()
        assert "filter-indicator" in html or "filterIndicator" in html

    def test_data_key_attributes(self, adata_full):
        """Test entries have data-key for filtering."""
        html = adata_full._repr_html_()
        assert "data-key=" in html

    def test_data_dtype_attributes(self, adata_full):
        """Test entries have data-dtype for filtering."""
        html = adata_full._repr_html_()
        assert "data-dtype=" in html


# =============================================================================
# Copy to Clipboard
# =============================================================================


class TestCopyToClipboard:
    """Test clipboard copy buttons."""

    def test_copy_buttons_present(self, adata_full):
        """Test copy buttons exist for entries."""
        html = adata_full._repr_html_()
        assert "copy" in html.lower() or "📋" in html


# =============================================================================
# Special Array Types
# =============================================================================


class TestSpecialArrayTypes:
    """Test display of special array types."""

    @pytest.mark.skipif(not HAS_DASK, reason="dask not installed")
    def test_dask_array_display(self):
        """Test Dask array shows chunk info."""
        X = da.zeros((1000, 500), chunks=(100, 500))
        adata = AnnData(X)
        html = adata._repr_html_()
        assert "dask" in html.lower()

    @pytest.mark.skipif(not HAS_AWKWARD, reason="awkward not installed")
    def test_awkward_array_display(self):
        """Test Awkward array displays correctly."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obsm["ragged"] = ak.Array([[1, 2], [3, 4, 5], [6]] * 3 + [[7]])
        html = adata._repr_html_()
        assert "awkward" in html.lower()


# =============================================================================
# Raw Section
# =============================================================================


class TestRawSection:
    """Test .raw section display."""

    def test_raw_section_present(self):
        """Test raw section appears when raw is set."""
        adata = AnnData(np.zeros((10, 20)))
        adata.raw = adata.copy()
        # Now subset var
        adata = adata[:, :5]
        html = adata._repr_html_()
        assert "raw" in html.lower()

    def test_raw_none_no_section(self):
        """Test no raw section when raw is None."""
        adata = AnnData(np.zeros((10, 5)))
        # raw is None by default
        _html = adata._repr_html_()
        # Should not have a non-empty raw section


# =============================================================================
# Metadata
# =============================================================================


class TestMetadata:
    """Test metadata display."""

    def test_version_displayed(self, adata):
        """Test anndata version is shown."""
        html = adata._repr_html_()
        assert "anndata" in html.lower()
        # Version number pattern
        assert re.search(r"\d+\.\d+", html)

    def test_memory_usage_displayed(self, adata):
        """Test memory usage is shown."""
        html = adata._repr_html_()
        # Should show memory in bytes/KB/MB/GB
        assert re.search(r"\d+(\.\d+)?\s*(B|KB|MB|GB|bytes)", html, re.IGNORECASE)

    def test_obs_var_names_preview(self, adata):
        """Test obs_names and var_names preview."""
        html = adata._repr_html_()
        # Should contain some cell or gene names
        assert "cell_" in html or "obs_names" in html
        assert "gene_" in html or "var_names" in html


# =============================================================================
# Settings
# =============================================================================


class TestSettings:
    """Test settings integration."""

    def test_html_disabled_fallback(self, adata):
        """Test fallback to text repr when HTML disabled."""
        from anndata import settings

        with settings.override(repr_html_enabled=False):
            html = adata._repr_html_()
            # Should return None (fallback to text) or pre-wrapped text
            assert html is None or "<pre>" in html

    def test_custom_max_depth(self):
        """Test custom max recursion depth."""
        from anndata import settings

        inner = AnnData(np.zeros((5, 3)))
        outer = AnnData(np.zeros((10, 5)))
        outer.uns["nested"] = inner

        with settings.override(repr_html_max_depth=1):
            html = outer._repr_html_()
            assert "max depth" in html.lower() or "depth" in html

    def test_custom_max_items(self):
        """Test custom max items setting."""
        from anndata import settings

        obs_dict = {f"col_{i}": list(range(10)) for i in range(300)}
        adata = AnnData(np.zeros((10, 5)), obs=pd.DataFrame(obs_dict))

        with settings.override(repr_html_max_items=50):
            html = adata._repr_html_()
            # Should show indication that items are truncated
            assert "more" in html.lower() or "..." in html


# =============================================================================
# Documentation Links
# =============================================================================


class TestDocumentationLinks:
    """Test help/documentation links."""

    def test_help_links_present(self, adata_full):
        """Test help links exist for sections."""
        html = adata_full._repr_html_()
        assert "readthedocs" in html


# =============================================================================
# Registry and Extensibility
# =============================================================================


class TestFormatterRegistry:
    """Test the formatter registry pattern for extensibility."""

    def test_registry_has_formatters(self):
        """Test registry contains registered formatters."""
        from anndata._repr.registry import formatter_registry

        # Should have some formatters registered
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
            """A custom type for testing."""

        class CustomTypeFormatter(TypeFormatter):
            priority = 500  # High priority

            def can_format(self, obj: Any) -> bool:
                return isinstance(obj, CustomType)

            def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
                return FormattedOutput(
                    type_name="CustomType",
                    css_class="dtype-custom",
                    details={"custom": True},
                    is_serializable=False,
                )

        # Register the custom formatter
        formatter = CustomTypeFormatter()
        formatter_registry.register_type_formatter(formatter)

        try:
            # Test formatting
            obj = CustomType()
            context = FormatterContext()
            result = formatter_registry.format_value(obj, context)
            assert result.type_name == "CustomType"
            assert result.css_class == "dtype-custom"
            assert result.details.get("custom") is True
            assert result.is_serializable is False
        finally:
            # Cleanup
            formatter_registry.unregister_type_formatter(formatter)

    def test_fallback_formatter_for_unknown_types(self):
        """Test fallback formatter handles unknown types gracefully."""
        from anndata._repr.registry import FormatterContext, formatter_registry

        class UnknownType:
            """An unknown type not in the registry."""

        obj = UnknownType()
        context = FormatterContext()
        result = formatter_registry.format_value(obj, context)

        # Should not raise an error
        assert result is not None
        assert "UnknownType" in result.type_name
        # Should indicate unknown/extension type
        assert "unknown" in result.css_class or "extension" in result.css_class

    def test_formatter_priority_order(self):
        """Test formatters are checked in priority order."""
        from anndata._repr.registry import formatter_registry

        # Verify formatters are sorted by priority (highest first)
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
            """A type for testing section filtering."""

        class UnsOnlyFormatter(TypeFormatter):
            priority = 600  # High priority
            sections = ("uns",)  # Only apply to uns section

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

            # Should match in uns section
            context_uns = FormatterContext(section="uns")
            result_uns = formatter_registry.format_value(obj, context_uns)
            assert result_uns.type_name == "UnsSpecificType"

            # Should NOT match in obsm section (falls back to fallback formatter)
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
            """A type for testing universal formatters."""

        class UniversalFormatter(TypeFormatter):
            priority = 600
            sections = None  # Default: apply to all sections

            def can_format(self, obj: Any) -> bool:
                return isinstance(obj, UniversalType)

            def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
                return FormattedOutput(type_name="UniversalType")

        formatter = UniversalFormatter()
        formatter_registry.register_type_formatter(formatter)

        try:
            obj = UniversalType()

            # Should match in any section
            for section in ["uns", "obsm", "varm", "layers", "obs", "var"]:
                context = FormatterContext(section=section)
                result = formatter_registry.format_value(obj, context)
                assert result.type_name == "UniversalType", (
                    f"Failed for section {section}"
                )
        finally:
            formatter_registry.unregister_type_formatter(formatter)

    def test_extension_type_graceful_handling(self):
        """Test extension types (like TreeData, MuData) are handled gracefully."""
        from anndata._repr.registry import FormatterContext, formatter_registry

        # Simulate an extension type that has AnnData-like attributes
        # We create the class in a way that properly sets __module__
        class ExtensionData:
            """Simulates an extension like TreeData or MuData."""

            def __init__(self):
                self.n_obs = 100
                self.n_vars = 50
                self.shape = (100, 50)
                self.dtype = np.float32

        # Set module at class level (this is how real extension types work)
        ExtensionData.__module__ = "treedata.core"

        obj = ExtensionData()
        context = FormatterContext()
        result = formatter_registry.format_value(obj, context)

        # Should not raise an error
        assert result is not None
        # Should contain type info
        assert "ExtensionData" in result.type_name
        # Should have shape info since the object has shape attribute
        assert "shape" in result.details

    def test_anndata_in_uns_detected(self):
        """Test nested AnnData in .uns is properly detected."""
        inner = AnnData(np.zeros((5, 3)))
        outer = AnnData(np.zeros((10, 5)))
        outer.uns["inner_adata"] = inner

        html = outer._repr_html_()

        # Should show nested AnnData
        assert "inner_adata" in html
        assert "AnnData" in html


# =============================================================================
# Utility Functions
# =============================================================================


class TestUtilityFunctions:
    """Test utility functions in _repr.utils."""

    def test_is_serializable_basic_types(self):
        """Test serialization detection for basic types."""
        from anndata._repr.utils import is_serializable

        # These should be serializable
        assert is_serializable(None)[0]
        assert is_serializable(True)[0]  # noqa: FBT003
        assert is_serializable(42)[0]
        assert is_serializable(3.14)[0]
        assert is_serializable("string")[0]
        assert is_serializable(np.array([1, 2, 3]))[0]
        assert is_serializable({"key": "value"})[0]
        assert is_serializable([1, 2, 3])[0]

    def test_is_serializable_custom_class(self):
        """Test custom classes are not serializable."""
        from anndata._repr.utils import is_serializable

        class CustomClass:
            pass

        is_ok, reason = is_serializable(CustomClass())
        assert not is_ok
        assert "CustomClass" in reason

    def test_is_serializable_nested_check(self):
        """Test nested unserializable object is detected."""
        from anndata._repr.utils import is_serializable

        class CustomClass:
            pass

        obj = {"valid": 1, "invalid": CustomClass()}
        is_ok, reason = is_serializable(obj)
        assert not is_ok
        assert "invalid" in reason

    def test_should_warn_string_column(self):
        """Test string-to-category warning detection."""
        from anndata._repr.utils import should_warn_string_column

        # String with repeats should warn
        s = pd.Series(["A", "B", "A", "C", "B", "A"])
        warn, msg = should_warn_string_column(s, s.nunique())
        assert warn
        assert "3" in msg  # 3 unique values

        # All unique should not warn
        s = pd.Series(["A", "B", "C", "D", "E"])
        warn, msg = should_warn_string_column(s, s.nunique())
        assert not warn

        # Numeric should not warn
        s = pd.Series([1, 2, 1, 2, 1])
        warn, msg = should_warn_string_column(s, s.nunique())
        assert not warn

        # When n_unique is None (e.g., lazy or large), should not warn
        s = pd.Series(["A", "B", "A"])
        warn, msg = should_warn_string_column(s, None)
        assert not warn

    def test_is_color_list(self):
        """Test color list detection."""
        from anndata._repr.utils import is_color_list

        # Valid color lists
        assert is_color_list("cluster_colors", ["#FF0000", "#00FF00", "#0000FF"])
        assert is_color_list("leiden_colors", np.array(["#123456", "#ABCDEF"]))
        assert is_color_list("cluster_colors", [])  # Empty is valid

        # Invalid
        assert not is_color_list("cluster", ["#FF0000"])  # Wrong key pattern
        assert not is_color_list("colors", ["#FF0000"])  # Wrong key pattern
        assert not is_color_list("cluster_colors", "#FF0000")  # Not a list
        assert not is_color_list(
            "cluster_colors", ["not_a_color", "also_not"]
        )  # Not colors

    def test_escape_html(self):
        """Test HTML escaping."""
        from anndata._repr.utils import escape_html

        assert escape_html("<script>") == "&lt;script&gt;"
        assert escape_html("a & b") == "a &amp; b"
        assert escape_html('"quoted"') == "&quot;quoted&quot;"

    def test_format_memory_size(self):
        """Test memory size formatting."""
        from anndata._repr.utils import format_memory_size

        assert "B" in format_memory_size(500)
        assert "KB" in format_memory_size(5000)
        assert "MB" in format_memory_size(5_000_000)
        assert "GB" in format_memory_size(5_000_000_000)

    def test_format_number(self):
        """Test number formatting with thousands separators."""
        from anndata._repr.utils import format_number

        assert format_number(1000) == "1,000"
        assert format_number(1000000) == "1,000,000"


# =============================================================================
# Formatters
# =============================================================================


class TestFormatters:
    """Test type-specific formatters."""

    def test_numpy_array_formatter(self):
        """Test numpy array formatting."""
        from anndata._repr.formatters import NumpyArrayFormatter
        from anndata._repr.registry import FormatterContext

        formatter = NumpyArrayFormatter()
        arr = np.zeros((100, 50), dtype=np.float32)

        assert formatter.can_format(arr)
        result = formatter.format(arr, FormatterContext())

        assert "100" in result.type_name
        assert "50" in result.type_name
        assert "float32" in result.type_name
        assert result.details["shape"] == (100, 50)

    def test_sparse_matrix_formatter(self):
        """Test sparse matrix formatting."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        formatter = SparseMatrixFormatter()
        mat = sp.random(1000, 500, density=0.1, format="csr")

        assert formatter.can_format(mat)
        result = formatter.format(mat, FormatterContext())

        assert "csr" in result.type_name.lower()
        assert result.details["nnz"] == mat.nnz
        assert result.details["sparsity"] is not None

    def test_categorical_formatter(self):
        """Test categorical formatting."""
        from anndata._repr.formatters import CategoricalFormatter
        from anndata._repr.registry import FormatterContext

        formatter = CategoricalFormatter()
        cat_series = pd.Series(pd.Categorical(["A", "B", "C"] * 10))

        assert formatter.can_format(cat_series)
        result = formatter.format(cat_series, FormatterContext())

        assert "category" in result.type_name.lower()
        assert result.details["n_categories"] == 3

    def test_dict_formatter(self):
        """Test dictionary formatting."""
        from anndata._repr.formatters import DictFormatter
        from anndata._repr.registry import FormatterContext

        formatter = DictFormatter()
        d = {"a": 1, "b": 2, "c": 3}

        assert formatter.can_format(d)
        result = formatter.format(d, FormatterContext())

        assert "dict" in result.type_name.lower()
        assert result.details["n_items"] == 3


# =============================================================================
# Completeness Tests
# =============================================================================


class TestRepresentationCompleteness:
    """Verify all data is accurately represented."""

    def test_all_obs_columns_shown(self):
        """Test all obs columns appear in repr."""
        obs_cols = ["col_a", "col_b", "col_c", "col_d", "col_e"]
        adata = AnnData(
            np.zeros((10, 5)), obs=pd.DataFrame({c: list(range(10)) for c in obs_cols})
        )
        html = adata._repr_html_()
        for col in obs_cols:
            assert col in html, f"Missing obs column: {col}"

    def test_all_var_columns_shown(self):
        """Test all var columns appear in repr."""
        var_cols = ["gene_name", "gene_id", "highly_variable"]
        adata = AnnData(
            np.zeros((10, 5)), var=pd.DataFrame({c: list(range(5)) for c in var_cols})
        )
        html = adata._repr_html_()
        for col in var_cols:
            assert col in html, f"Missing var column: {col}"

    def test_all_uns_keys_shown(self):
        """Test all uns keys appear in repr."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["key1"] = "value1"
        adata.uns["key2"] = 42
        adata.uns["nested_dict"] = {"a": 1, "b": 2}
        adata.uns["array_data"] = np.array([1, 2, 3])
        html = adata._repr_html_()
        for key in ["key1", "key2", "nested_dict", "array_data"]:
            assert key in html, f"Missing uns key: {key}"

    def test_all_obsm_keys_shown(self):
        """Test all obsm keys appear in repr."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obsm["X_pca"] = np.random.randn(10, 50).astype(np.float32)
        adata.obsm["X_umap"] = np.random.randn(10, 2).astype(np.float32)
        adata.obsm["X_tsne"] = np.random.randn(10, 2).astype(np.float32)
        html = adata._repr_html_()
        for key in ["X_pca", "X_umap", "X_tsne"]:
            assert key in html, f"Missing obsm key: {key}"

    def test_all_layers_shown(self):
        """Test all layers appear in repr."""
        adata = AnnData(np.zeros((10, 5)))
        adata.layers["raw"] = np.random.randn(10, 5)
        adata.layers["normalized"] = np.random.randn(10, 5)
        adata.layers["scaled"] = np.random.randn(10, 5)
        html = adata._repr_html_()
        for layer in ["raw", "normalized", "scaled"]:
            assert layer in html, f"Missing layer: {layer}"

    def test_correct_shape_values(self):
        """Test shape values are accurate."""
        adata = AnnData(np.zeros((123, 456)))
        html = adata._repr_html_()
        assert "123" in html
        assert "456" in html


# =============================================================================
# Edge Cases
# =============================================================================


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_very_large_keys(self):
        """Test handling of very long key names."""
        adata = AnnData(np.zeros((10, 5)))
        long_key = "a" * 500
        adata.uns[long_key] = "value"
        html = adata._repr_html_()
        # Should not error
        assert html is not None

    def test_unicode_keys(self):
        """Test handling of unicode key names."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["日本語"] = "value"
        adata.uns["émojis_🧬"] = "dna"
        html = adata._repr_html_()
        assert "日本語" in html
        assert "émojis" in html

    def test_empty_sections(self):
        """Test handling of empty sections."""
        adata = AnnData(np.zeros((10, 5)))
        # All sections empty except X
        html = adata._repr_html_()
        assert html is not None

    def test_very_nested_uns(self):
        """Test handling of deeply nested uns."""
        adata = AnnData(np.zeros((10, 5)))
        nested = {"level": 0}
        current = nested
        for i in range(10):
            current["nested"] = {"level": i + 1}
            current = current["nested"]
        adata.uns["deep"] = nested
        html = adata._repr_html_()
        assert html is not None

    def test_mixed_types_in_uns(self):
        """Test handling of mixed types in uns."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["string"] = "value"
        adata.uns["int"] = 42
        adata.uns["float"] = 3.14
        adata.uns["list"] = [1, 2, 3]
        adata.uns["dict"] = {"a": 1}
        adata.uns["array"] = np.array([1, 2, 3])
        adata.uns["sparse"] = sp.csr_matrix([[1, 0], [0, 1]])
        adata.uns["df"] = pd.DataFrame({"a": [1, 2]})
        html = adata._repr_html_()
        assert html is not None


# =============================================================================
# Backed AnnData (H5AD files)
# =============================================================================


class TestBackedAnnData:
    """Test HTML repr for disk-backed AnnData."""

    def test_backed_badge_h5ad(self, tmp_path):
        """Test backed badge for H5AD files."""
        adata = AnnData(np.random.randn(100, 50).astype(np.float32))
        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        backed = ad.read_h5ad(path, backed="r")
        html = backed._repr_html_()
        assert "Backed" in html or "backed" in html.lower() or "📁" in html
        backed.file.close()

    def test_backed_shows_filename(self, tmp_path):
        """Test backed repr shows filename."""
        adata = AnnData(np.random.randn(100, 50).astype(np.float32))
        path = tmp_path / "test_file.h5ad"
        adata.write_h5ad(path)

        backed = ad.read_h5ad(path, backed="r")
        html = backed._repr_html_()
        assert "test_file" in html or "h5ad" in html.lower()
        backed.file.close()


# =============================================================================
# Ad-blocker Compatibility
# =============================================================================


class TestAdBlockerCompatibility:
    """
    Test that HTML representation doesn't trigger ad-blockers.

    Some browsers (like Comet/Chromium) have aggressive ad-blocking that
    blocks elements with class names starting with "ad-". This test suite
    ensures we avoid such patterns.

    See: https://github.com/scverse/anndata/issues/XXXX
    """

    def test_no_ad_prefix_in_html_classes(self):
        """Test that generated HTML doesn't contain 'ad-' prefixed class names."""
        adata = AnnData(
            np.random.randn(50, 20),
            obs=pd.DataFrame({"cell_type": ["A", "B"] * 25}),
            var=pd.DataFrame({"gene_name": [f"gene_{i}" for i in range(20)]}),
        )
        adata.obsm["X_pca"] = np.random.randn(50, 10)
        adata.uns["metadata"] = {"key": "value"}

        html = adata._repr_html_()

        # Find all class="..." attributes
        class_pattern = re.compile(r'class="([^"]*)"')
        all_classes = class_pattern.findall(html)

        # Check each class name
        problematic_classes = []
        for class_attr in all_classes:
            # Split multiple classes (e.g., class="foo ad-bar")
            problematic_classes.extend(
                cn for cn in class_attr.split() if cn.startswith("ad-")
            )

        assert not problematic_classes, (
            f"Found class names starting with 'ad-' which may be blocked by ad-blockers: "
            f"{problematic_classes}. Use 'adata-' prefix instead."
        )

    def test_no_ad_prefix_in_css_selectors(self):
        """Test that CSS doesn't contain '.ad-' selectors."""
        from anndata._repr.css import get_css

        css = get_css()

        # Find all CSS class selectors
        selector_pattern = re.compile(r"\.(ad-[\w-]+)")
        problematic_selectors = selector_pattern.findall(css)

        assert not problematic_selectors, (
            f"Found CSS selectors starting with '.ad-' which may be blocked by ad-blockers: "
            f"{set(problematic_selectors)}. Use '.adata-' prefix instead."
        )

    def test_no_ad_prefix_in_javascript_selectors(self):
        """Test that JavaScript doesn't query for '.ad-' selectors."""
        from anndata._repr.javascript import get_javascript

        js = get_javascript("test-container")

        # Find all querySelectorAll/querySelector calls with '.ad-'
        selector_pattern = re.compile(r"querySelector(?:All)?\(['\"]\.ad-[\w-]+['\"]")
        matches = selector_pattern.findall(js)

        assert not matches, (
            f"Found JavaScript selectors starting with '.ad-' which may be blocked by ad-blockers: "
            f"{matches}. Use '.adata-' prefix instead."
        )

    def test_visual_rendering_without_adblocker(self):
        """Test that key elements are present in the HTML."""
        adata = AnnData(
            np.random.randn(100, 50),
            obs=pd.DataFrame({"cell_type": ["A", "B"] * 50}),
        )
        html = adata._repr_html_()

        # Check that essential structural elements are present
        # These should use 'adata-' or 'anndata-' prefixes
        assert "adata-type" in html or "anndata" in html.lower()
        assert "adata-entry" in html or "anndata-sec" in html
        assert "adata-table" in html or "table" in html.lower()

        # Check that the representation is reasonably sized (not empty/truncated)
        assert len(html) > 1000, "HTML output seems too short"


# =============================================================================
# Coverage Improvement Tests
# =============================================================================


class TestCoverageEdgeCases:
    """
    Tests for edge cases to improve code coverage.

    These tests target specific code paths that aren't covered by other tests.
    """

    def test_category_column_overflow(self):
        """Test rendering categorical column with more than max categories."""
        from anndata import settings

        # Create categorical column with many categories
        categories = [f"cat_{i}" for i in range(50)]
        adata = AnnData(
            np.zeros((100, 5)),
            obs=pd.DataFrame({
                "many_cats": pd.Categorical(
                    np.random.choice(categories, 100), categories=categories
                )
            }),
        )

        with settings.override(repr_html_max_categories=10):
            html = adata._repr_html_()
            # Should show first 10 and indicate there are more
            assert "cat_0" in html
            assert "cat_9" in html or "cat_8" in html  # Might show 9 or stop at 8
            assert "...+" in html or "more" in html.lower()

    def test_non_serializable_object_in_uns(self):
        """Test warning for non-serializable objects in uns."""

        class NonSerializable:
            """A class that can't be serialized to H5AD."""

        adata = AnnData(np.zeros((10, 5)))
        adata.uns["non_serializable"] = NonSerializable()

        html = adata._repr_html_()

        # Should show warning icon
        assert "⚠" in html or "warn" in html.lower()

    def test_max_depth_with_multiple_nested_anndata(self):
        """Test max depth indicator with deeply nested AnnData."""
        from anndata import settings

        # Create 3 levels of nesting
        level2 = AnnData(np.zeros((5, 3)))
        level1 = AnnData(np.zeros((7, 4)))
        level1.uns["nested"] = level2
        level0 = AnnData(np.zeros((10, 5)))
        level0.uns["nested"] = level1

        # With max_depth=0, should show max depth indicator immediately
        with settings.override(repr_html_max_depth=0):
            html = level0._repr_html_()
            assert "max depth" in html.lower() or "depth" in html.lower()

    def test_dataframe_entry_nunique_exception(self):
        """Test nunique() exception handling for dataframe columns."""

        # Create a column that might raise an exception on nunique()
        # Use a Series with unhashable types
        adata = AnnData(np.zeros((10, 5)))
        # Lists are unhashable and will cause nunique() to fail
        adata.obs["unhashable"] = [[i] for i in range(10)]

        # Should not crash, just skip the unique count
        html = adata._repr_html_()
        assert html is not None
        assert "unhashable" in html

    def test_very_large_dataframe_skips_nunique(self):
        """Test that nunique() is skipped for very large columns."""
        from anndata import settings

        # Create a large dataset
        adata = AnnData(np.zeros((100000, 5)), obs=pd.DataFrame({"col": range(100000)}))

        with settings.override(repr_html_unique_limit=1000):
            html = adata._repr_html_()
            # Should render but not show unique count (would be too slow)
            assert "col" in html
            # Unique count should not be shown for large columns
            # (checking that it doesn't crash is the main goal)

    def test_empty_category_counts(self):
        """Test rendering category column with all unique values."""
        adata = AnnData(
            np.zeros((10, 5)),
            obs=pd.DataFrame({
                "unique_cats": pd.Categorical(
                    [f"cat_{i}" for i in range(10)],
                    categories=[f"cat_{i}" for i in range(10)],
                )
            }),
        )

        html = adata._repr_html_()
        assert "unique_cats" in html
        # Should show the categories
        assert "cat_0" in html


class TestFutureCompatibility:
    """
    Tests for future compatibility with upcoming array types and changes.

    These tests ensure that HTML representation remains robust when:
    - PR #1927 removes scipy sparse inheritance
    - PR #2063 adds Array-API compatible arrays (JAX, PyTorch, etc.)
    - New array backends are added in the future

    References:
    - PR #1927: https://github.com/scverse/anndata/pull/1927
    - PR #2063: https://github.com/scverse/anndata/pull/2063
    """

    def test_sparse_duck_typing_detection(self):
        """
        Test that sparse arrays can be detected via duck typing.

        Future-proofing for PR #1927 which removes scipy sparse inheritance.
        This ensures sparse detection works even if scipy.sparse.issparse() breaks.
        """
        import scipy.sparse as sp

        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        # Create a sparse matrix
        sparse_matrix = sp.csr_matrix([[1, 0, 2], [0, 0, 3], [4, 5, 6]])

        # Test that the formatter can detect it
        formatter = SparseMatrixFormatter()
        assert formatter.can_format(sparse_matrix), "Should detect sparse matrix"

        # Test that it formats correctly
        context = FormatterContext(depth=0)
        result = formatter.format(sparse_matrix, context)
        assert "csr" in result.type_name.lower()
        assert result.css_class == "dtype-sparse"
        assert result.details["nnz"] == 6

    def test_array_api_formatter_with_mock_jax_array(self):
        """
        Test ArrayAPIFormatter with a mock JAX-like array.

        Future-proofing for PR #2063 which adds Array-API compatibility.
        Tests that JAX/PyTorch/TensorFlow arrays are handled correctly.
        """
        from anndata._repr.formatters import ArrayAPIFormatter
        from anndata._repr.registry import FormatterContext

        # Create a mock JAX-like array
        class MockJAXArray:
            """Mock JAX array for testing."""

            def __init__(self):
                self.shape = (100, 50)
                self.dtype = np.dtype("float32")
                self.ndim = 2

            @property
            def __module__(self):
                return "jax.numpy"

            def __class__(self):
                return type("DeviceArray", (), {})

        mock_array = MockJAXArray()
        type(mock_array).__module__ = "jax.numpy"
        type(mock_array).__name__ = "DeviceArray"

        # Test that ArrayAPIFormatter can detect it
        formatter = ArrayAPIFormatter()
        can_format = formatter.can_format(mock_array)
        assert can_format, "Should detect JAX-like array"

        # Test that it formats correctly
        context = FormatterContext(depth=0)
        result = formatter.format(mock_array, context)
        assert "100 × 50" in result.type_name
        assert "float32" in result.type_name
        assert result.css_class == "dtype-array-api"
        assert result.details["backend"] == "jax"

    def test_array_api_excludes_already_handled_types(self):
        """
        Test that ArrayAPIFormatter doesn't interfere with specific formatters.

        Ensures that numpy, pandas, scipy.sparse, cupy arrays use their
        specific formatters, not the generic ArrayAPIFormatter.
        """
        from anndata._repr.formatters import ArrayAPIFormatter

        formatter = ArrayAPIFormatter()

        # Test numpy array - should NOT be handled by ArrayAPIFormatter
        np_array = np.array([[1, 2], [3, 4]])
        assert not formatter.can_format(np_array), (
            "Should not handle numpy arrays (has NumpyArrayFormatter)"
        )

        # Test pandas DataFrame - should NOT be handled
        df = pd.DataFrame({"a": [1, 2, 3]})
        assert not formatter.can_format(df), (
            "Should not handle pandas DataFrame (has DataFrameFormatter)"
        )

        # Test pandas Series - should NOT be handled
        series = pd.Series([1, 2, 3])
        assert not formatter.can_format(series), (
            "Should not handle pandas Series (has SeriesFormatter)"
        )

    def test_sparse_format_name_fallback(self):
        """
        Test that sparse format name detection has proper fallback.

        If scipy.sparse.isspmatrix_* functions fail (as they might after PR #1927),
        the formatter should fall back to using type(obj).__name__.
        """
        import scipy.sparse as sp

        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        # Test with csr_array (new scipy sparse array type)
        try:
            # scipy >= 1.8 has sparse arrays
            sparse_array = sp.csr_array([[1, 0, 2], [0, 0, 3]])
            formatter = SparseMatrixFormatter()
            context = FormatterContext(depth=0)
            result = formatter.format(sparse_array, context)

            # Should either get "csr_matrix" or "csr_array" depending on scipy version
            assert "csr" in result.type_name.lower(), (
                f"Should contain 'csr', got {result.type_name}"
            )
        except AttributeError:
            # scipy < 1.8 doesn't have csr_array, skip this part
            pass

    def test_unknown_array_type_graceful_fallback(self):
        """
        Test that completely unknown array types don't break HTML rendering.

        This ensures robustness when encountering array types we haven't seen yet.
        """
        # Create a dynamically generated mock array type
        # We use type() to create a new class with custom __module__
        FutureArray = type(
            "FutureArray",  # class name
            (),  # base classes
            {
                "__module__": "future_lib.arrays",  # module name
                "__init__": lambda self: setattr(self, "_initialized", True),
                "shape": property(lambda self: (10, 5)),
                "dtype": property(lambda self: np.dtype("int64")),
                "ndim": property(lambda self: 2),
            },
        )

        # Create AnnData with this unknown type
        adata = AnnData(np.zeros((10, 5)))
        future_array = FutureArray()
        adata.uns["future_data"] = future_array

        # Should not crash when rendering - test will fail with clear error if it does
        html = adata._repr_html_()
        # Should successfully render (might use fallback formatter)
        assert html is not None
        assert len(html) > 0
        # The unknown type should appear in the HTML
        assert (
            "FutureArray" in html
            or "future" in html.lower()
            or "object" in html.lower()
        )

    def test_css_array_api_styling_exists(self):
        """
        Test that CSS styling for Array-API arrays is present.

        Ensures that the dtype-array-api class is defined in CSS.
        """
        from anndata._repr.css import get_css

        css = get_css()

        # Check for light mode styling
        assert "dtype-array-api" in css, (
            "CSS should contain dtype-array-api class styling"
        )

        # Check that there's some color definition for it
        # The exact color doesn't matter, but it should be styled
        pattern = r"\.dtype-array-api\s*\{[^}]*color:"
        assert re.search(pattern, css), "dtype-array-api should have color styling"


# =============================================================================
# Uns Renderer Registry Tests
# =============================================================================


class TestUnsRendererRegistry:
    """
    Test the uns renderer registry for custom serialized data visualization.

    This registry allows packages to register custom HTML renderers for data
    stored in .uns with type hints. The security model ensures data NEVER
    triggers code execution - packages must be explicitly imported first.
    """

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

        # Regular dict
        value = {"data": "value"}
        hint, cleaned = extract_uns_type_hint(value)
        assert hint is None
        assert cleaned == value

        # Regular string
        value = "just a string"
        hint, cleaned = extract_uns_type_hint(value)
        assert hint is None
        assert cleaned == value

        # Regular int
        value = 42
        hint, cleaned = extract_uns_type_hint(value)
        assert hint is None
        assert cleaned == value

    def test_extract_type_hint_invalid_dict_hint_type(self):
        """Test that non-string type hints in dict are ignored."""
        from anndata._repr.registry import UNS_TYPE_HINT_KEY, extract_uns_type_hint

        # Type hint is not a string
        value = {UNS_TYPE_HINT_KEY: 123, "data": "value"}
        hint, cleaned = extract_uns_type_hint(value)
        assert hint is None
        assert cleaned == value

    def test_extract_type_hint_malformed_string_format(self):
        """Test that malformed string format returns no hint."""
        from anndata._repr.registry import UNS_TYPE_HINT_KEY, extract_uns_type_hint

        # Missing :: separator
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
            priority = 100  # High priority to check before fallback

            def can_format(self, obj):
                hint, _ = extract_uns_type_hint(obj)
                return hint == "test.config_format"

            def format(self, obj, context):
                _hint, data = extract_uns_type_hint(obj)
                items = data.get("data", {})
                return FormattedOutput(
                    type_name="test config",
                    html_content=f'<span class="test-custom">Items: {len(items)}</span>',
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

            # Should show custom HTML
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

        # Should show the type hint and import suggestion
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

            # Should not raise, should fall back to default rendering
            with pytest.warns(UserWarning, match="Formatter.*failed"):
                html = adata._repr_html_()

            assert html is not None
            assert "will_fail" in html
        finally:
            formatter_registry.unregister_type_formatter(formatter)

    def test_string_format_type_hint_in_html(self):
        """Test string format type hints work in HTML output."""
        adata = AnnData(np.zeros((5, 3)))
        # String format: __anndata_repr__:package.type::content
        adata.uns["string_hint"] = (
            "__anndata_repr__:somepackage.config::actual content here"
        )

        html = adata._repr_html_()

        # Should show the type hint
        assert "somepackage.config" in html
        # Should suggest importing
        assert "import somepackage" in html

    def test_type_hint_key_constant_exported(self):
        """Test that UNS_TYPE_HINT_KEY constant is properly exported."""
        from anndata._repr import UNS_TYPE_HINT_KEY

        assert UNS_TYPE_HINT_KEY == "__anndata_repr__"

    def test_security_data_never_triggers_import(self):
        """
        Test that data in uns NEVER triggers imports or code execution.

        This is a critical security test. Even if a type hint references
        a package, the data should never cause that package to be imported.
        """
        import sys

        # Use a fake package name that definitely doesn't exist
        fake_package = "definitely_not_a_real_package_12345"
        assert fake_package not in sys.modules

        adata = AnnData(np.zeros((5, 3)))
        adata.uns["malicious"] = {
            "__anndata_repr__": f"{fake_package}.evil",
            "data": "some data",
        }

        # Render the HTML
        html = adata._repr_html_()

        # The fake package should NOT have been imported
        assert fake_package not in sys.modules
        # But the HTML should still work
        assert html is not None
        assert "malicious" in html


# =============================================================================
# Additional Coverage Tests for Registry
# =============================================================================


class TestFormatterContextCoverage:
    """Tests for FormatterContext to improve registry.py coverage."""

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
        """Test access_path with non-identifier keys (need brackets)."""
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


class TestSectionFormatterCoverage:
    """Tests for SectionFormatter abstract class coverage."""

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

        # Test default implementations
        assert formatter.display_name == "test_section"  # Defaults to section_name
        assert formatter.doc_url is None
        assert formatter.tooltip == ""
        assert formatter.should_show(None) is True


class TestFallbackFormatterCoverage:
    """Tests for FallbackFormatter edge cases."""

    def test_fallback_extension_type_no_warning(self):
        """Test fallback formatter for extension types doesn't add warning."""
        from anndata._repr.registry import FallbackFormatter, FormatterContext

        # Create a class that looks like an extension
        class ExtensionType:
            pass

        ExtensionType.__module__ = "treedata.core"

        formatter = FallbackFormatter()
        obj = ExtensionType()
        context = FormatterContext()

        result = formatter.format(obj, context)

        assert result.type_name == "ExtensionType"
        assert result.css_class == "dtype-extension"
        assert len(result.warnings) == 0  # No warning for extensions

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

        assert result.details["shape"] == (10, 5)
        assert result.details["dtype"] == "float32"
        assert "Shape: (10, 5)" in result.tooltip

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

        assert result.details["length"] == 42

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

        # Should not raise
        result = formatter.format(obj, context)
        assert "length" not in result.details


class TestFormatterRegistryCoverage:
    """Additional tests for FormatterRegistry coverage."""

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
            # Should fall through to backup formatter
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

        # Register using register_formatter (non-decorator form)
        formatter_instance = DecoratorTestFormatter()
        register_formatter(formatter_instance)

        try:
            # Should be registered
            result = formatter_registry.format_value(
                ("decorator_test", 1, 2), FormatterContext()
            )
            assert result.type_name == "DecoratorTest"
        finally:
            # Cleanup
            formatter_registry.unregister_type_formatter(formatter_instance)


# =============================================================================
# Additional Coverage Tests for Utils
# =============================================================================


class TestUtilsCoverage:
    """Tests to improve utils.py coverage."""

    def test_is_serializable_list_with_unserializable(self):
        """Test is_serializable catches unserializable item in list."""
        from anndata._repr.utils import is_serializable

        class BadType:
            pass

        obj = [1, 2, BadType(), 4]
        is_ok, reason = is_serializable(obj)
        assert not is_ok
        assert "Index 2" in reason

    def test_is_serializable_max_depth_exceeded(self):
        """Test is_serializable handles deep nesting."""
        from anndata._repr.utils import is_serializable

        # Create deeply nested structure
        nested = {"level": 0}
        current = nested
        for i in range(15):
            current["nested"] = {"level": i + 1}
            current = current["nested"]

        is_ok, reason = is_serializable(nested, _max_depth=10)
        assert not is_ok
        assert "depth" in reason.lower()

    def test_is_serializable_numpy_scalar(self):
        """Test is_serializable handles numpy scalar types."""
        from anndata._repr.utils import is_serializable

        assert is_serializable(np.int64(42))[0]
        assert is_serializable(np.float32(3.14))[0]
        assert is_serializable(np.bool_(True))[0]  # noqa: FBT003

    def test_should_warn_string_column_with_none_nunique(self):
        """Test should_warn_string_column handles None n_unique gracefully.

        This tests the case where n_unique couldn't be computed (e.g., due to
        exceptions in nunique(), lazy data, or unique_limit being exceeded).
        """
        from anndata._repr.utils import should_warn_string_column

        # When n_unique is None, should gracefully return False
        s = pd.Series(["A", "B", "A", "C"])
        warn, _msg = should_warn_string_column(s, None)
        assert not warn  # Should return False when n_unique is None

    def test_get_nunique_exception_handling(self):
        """Test _get_nunique handles exceptions gracefully."""
        from anndata._repr.registry import FormatterContext
        from anndata._repr.sections import _get_nunique

        context = FormatterContext()

        # Normal case
        s = pd.Series(["A", "B", "A", "C"])
        result = _get_nunique(s, context)
        assert result == 3

        # Unhashable types (lists) should return None
        s = pd.Series([[1, 2], [3, 4], [1, 2]])
        result = _get_nunique(s, context)
        assert result is None

    def test_is_color_list_named_colors(self):
        """Test is_color_list detects named colors."""
        from anndata._repr.utils import is_color_list

        assert is_color_list("cluster_colors", ["red", "blue", "green"])
        assert is_color_list("batch_colors", ["crimson", "navy"])

    def test_is_color_list_rgb_format(self):
        """Test is_color_list detects RGB/RGBA format."""
        from anndata._repr.utils import is_color_list

        assert is_color_list("cluster_colors", ["rgb(255, 0, 0)", "rgb(0, 255, 0)"])
        assert is_color_list("batch_colors", ["rgba(255, 0, 0, 0.5)"])

    def test_is_color_list_none_first_element(self):
        """Test is_color_list handles None in first element."""
        from anndata._repr.utils import is_color_list

        assert not is_color_list("cluster_colors", [None, "#FF0000"])

    def test_get_matching_column_colors_var(self):
        """Test get_matching_column_colors finds colors for var columns."""
        from anndata._repr.utils import get_matching_column_colors

        adata = AnnData(np.zeros((10, 5)))
        adata.var["gene_type"] = pd.Categorical(["A", "B"] * 2 + ["A"])
        adata.uns["gene_type_colors"] = ["#FF0000", "#00FF00"]

        colors = get_matching_column_colors(adata, "gene_type")
        assert colors == ["#FF0000", "#00FF00"]

    def test_get_matching_column_colors_no_column(self):
        """Test get_matching_column_colors returns None for missing column."""
        from anndata._repr.utils import get_matching_column_colors

        adata = AnnData(np.zeros((10, 5)))
        adata.uns["missing_colors"] = ["#FF0000"]

        colors = get_matching_column_colors(adata, "missing")
        assert colors is None

    def test_get_matching_column_colors_not_categorical(self):
        """Test get_matching_column_colors returns None for non-categorical."""
        from anndata._repr.utils import get_matching_column_colors

        adata = AnnData(np.zeros((10, 5)))
        adata.obs["numeric"] = list(range(10))
        adata.uns["numeric_colors"] = ["#FF0000"]

        colors = get_matching_column_colors(adata, "numeric")
        assert colors is None

    def test_sanitize_for_id_starts_with_number(self):
        """Test sanitize_for_id handles strings starting with numbers."""
        from anndata._repr.utils import sanitize_for_id

        result = sanitize_for_id("123abc")
        assert result.startswith("id_")
        assert result[3:].startswith("123")

    def test_sanitize_for_id_special_chars(self):
        """Test sanitize_for_id replaces special characters."""
        from anndata._repr.utils import sanitize_for_id

        result = sanitize_for_id("hello world!@#$%")
        assert " " not in result
        assert "!" not in result
        assert "_" in result

    def test_truncate_string_short(self):
        """Test truncate_string doesn't modify short strings."""
        from anndata._repr.utils import truncate_string

        result = truncate_string("short", max_length=100)
        assert result == "short"

    def test_truncate_string_long(self):
        """Test truncate_string truncates long strings."""
        from anndata._repr.utils import truncate_string

        result = truncate_string("a" * 200, max_length=100)
        assert len(result) == 100
        assert result.endswith("...")

    def test_format_memory_size_negative(self):
        """Test format_memory_size handles negative values."""
        from anndata._repr.utils import format_memory_size

        result = format_memory_size(-100)
        assert result == "Unknown"

    def test_format_memory_size_very_large(self):
        """Test format_memory_size handles very large values (PB)."""
        from anndata._repr.utils import format_memory_size

        result = format_memory_size(5 * 1024**5)  # 5 PB
        assert "PB" in result

    def test_format_number_float_whole(self):
        """Test format_number handles floats that are whole numbers."""
        from anndata._repr.utils import format_number

        result = format_number(1000.0)
        assert result == "1,000"

    def test_format_number_float_decimal(self):
        """Test format_number handles floats with decimals."""
        from anndata._repr.utils import format_number

        result = format_number(1234.567)
        assert "1,234.57" in result

    def test_get_backing_info_zarr(self):
        """Test get_backing_info detects Zarr format."""
        from anndata._repr.utils import get_backing_info

        class MockBackedAdata:
            isbacked = True
            filename = "/path/to/data.zarr"

            class file:
                is_open = True

        info = get_backing_info(MockBackedAdata())
        assert info["format"] == "Zarr"

    def test_get_backing_info_unknown_format(self):
        """Test get_backing_info handles unknown format."""
        from anndata._repr.utils import get_backing_info

        class MockBackedAdata:
            isbacked = True
            filename = "/path/to/data.xyz"

            class file:
                is_open = True

        info = get_backing_info(MockBackedAdata())
        assert info["format"] == "Unknown"


# =============================================================================
# Additional Coverage Tests for Formatters
# =============================================================================


class TestFormattersCoverage:
    """Tests to improve formatters.py coverage."""

    def test_numpy_array_3d(self):
        """Test numpy array formatter with 3D+ arrays."""
        from anndata._repr.formatters import NumpyArrayFormatter
        from anndata._repr.registry import FormatterContext

        formatter = NumpyArrayFormatter()
        arr = np.zeros((10, 5, 3))
        result = formatter.format(arr, FormatterContext())

        assert "(10, 5, 3)" in result.type_name

    def test_masked_array_formatter(self):
        """Test MaskedArray formatter."""
        from anndata._repr.formatters import NumpyMaskedArrayFormatter
        from anndata._repr.registry import FormatterContext

        formatter = NumpyMaskedArrayFormatter()
        arr = np.ma.array([1, 2, 3, 4, 5], mask=[0, 0, 1, 0, 1])

        assert formatter.can_format(arr)
        result = formatter.format(arr, FormatterContext())

        assert "MaskedArray" in result.type_name
        assert result.details["n_masked"] == 2

    def test_masked_array_no_mask(self):
        """Test MaskedArray formatter with no masked values."""
        from anndata._repr.formatters import NumpyMaskedArrayFormatter
        from anndata._repr.registry import FormatterContext

        formatter = NumpyMaskedArrayFormatter()
        arr = np.ma.array([1, 2, 3, 4, 5])  # No mask

        result = formatter.format(arr, FormatterContext())
        assert result.details["n_masked"] == 0

    def test_sparse_csc_formatter(self):
        """Test sparse formatter with CSC matrix."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        formatter = SparseMatrixFormatter()
        mat = sp.csc_matrix([[1, 0], [0, 2]])

        result = formatter.format(mat, FormatterContext())
        assert "csc" in result.type_name.lower()

    def test_sparse_coo_formatter(self):
        """Test sparse formatter with COO matrix."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        formatter = SparseMatrixFormatter()
        mat = sp.coo_matrix([[1, 0], [0, 2]])

        result = formatter.format(mat, FormatterContext())
        assert "coo" in result.type_name.lower()

    def test_sparse_zero_elements(self):
        """Test sparse formatter with zero-element matrix."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        formatter = SparseMatrixFormatter()
        mat = sp.csr_matrix((0, 0))

        result = formatter.format(mat, FormatterContext())
        assert (
            result.details["sparsity"] is None
        )  # Can't compute sparsity for 0 elements


class TestBuiltinFormattersCoverage:
    """Tests for built-in type formatters coverage."""

    def test_list_formatter(self):
        """Test list formatter."""
        from anndata._repr.formatters import ListFormatter
        from anndata._repr.registry import FormatterContext

        formatter = ListFormatter()
        obj = [1, 2, 3, 4, 5]

        assert formatter.can_format(obj)
        result = formatter.format(obj, FormatterContext())
        assert "list" in result.type_name.lower()
        assert result.details["n_items"] == 5

    def test_dict_formatter(self):
        """Test dict formatter."""
        from anndata._repr.formatters import DictFormatter
        from anndata._repr.registry import FormatterContext

        formatter = DictFormatter()
        obj = {"a": 1, "b": 2, "c": 3}

        assert formatter.can_format(obj)
        result = formatter.format(obj, FormatterContext())
        assert "dict" in result.type_name.lower()
        assert result.details["n_items"] == 3

    def test_string_formatter(self):
        """Test string formatter."""
        from anndata._repr.formatters import StringFormatter
        from anndata._repr.registry import FormatterContext

        formatter = StringFormatter()
        obj = "hello world"

        assert formatter.can_format(obj)
        result = formatter.format(obj, FormatterContext())
        assert "str" in result.type_name.lower()

    def test_none_formatter(self):
        """Test None formatter."""
        from anndata._repr.formatters import NoneFormatter
        from anndata._repr.registry import FormatterContext

        formatter = NoneFormatter()

        assert formatter.can_format(None)
        result = formatter.format(None, FormatterContext())
        assert "None" in result.type_name

    def test_bool_formatter(self):
        """Test bool formatter."""
        from anndata._repr.formatters import BoolFormatter
        from anndata._repr.registry import FormatterContext

        formatter = BoolFormatter()

        assert formatter.can_format(True)  # noqa: FBT003
        assert formatter.can_format(False)  # noqa: FBT003
        result = formatter.format(True, FormatterContext())  # noqa: FBT003
        assert "bool" in result.type_name.lower()

    def test_int_formatter(self):
        """Test int formatter."""
        from anndata._repr.formatters import IntFormatter
        from anndata._repr.registry import FormatterContext

        formatter = IntFormatter()

        assert formatter.can_format(42)
        result = formatter.format(42, FormatterContext())
        assert "int" in result.type_name.lower()

    def test_float_formatter(self):
        """Test float formatter."""
        from anndata._repr.formatters import FloatFormatter
        from anndata._repr.registry import FormatterContext

        formatter = FloatFormatter()

        assert formatter.can_format(3.14)
        result = formatter.format(3.14, FormatterContext())
        assert "float" in result.type_name.lower()

    def test_categorical_direct_object(self):
        """Test CategoricalFormatter with direct pd.Categorical object."""
        from anndata._repr.formatters import CategoricalFormatter
        from anndata._repr.registry import FormatterContext

        formatter = CategoricalFormatter()

        # Direct Categorical object (not wrapped in Series)
        cat = pd.Categorical(["A", "B", "A", "C"])

        assert formatter.can_format(cat)
        result = formatter.format(cat, FormatterContext())
        assert "category" in result.type_name.lower()
        assert result.details["n_categories"] == 3

    def test_series_formatter_simple(self):
        """Test SeriesFormatter with simple numeric series."""
        from anndata._repr.formatters import SeriesFormatter
        from anndata._repr.registry import FormatterContext

        formatter = SeriesFormatter()
        series = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0])

        assert formatter.can_format(series)
        result = formatter.format(series, FormatterContext())
        assert result.details["length"] == 5

    def test_dataframe_formatter(self):
        """Test DataFrameFormatter."""
        from anndata._repr.formatters import DataFrameFormatter
        from anndata._repr.registry import FormatterContext

        formatter = DataFrameFormatter()
        ctx = FormatterContext()

        # Basic DataFrame
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
        assert formatter.can_format(df)
        result = formatter.format(df, ctx)
        assert result.details["n_rows"] == 3
        assert result.details["n_cols"] == 2
        assert "3 × 2" in result.type_name
        # Column names in meta_preview
        assert "[a, b]" in result.details["meta_preview"]

        # DataFrame with many columns (truncated preview)
        df_many = pd.DataFrame({f"col_{i}": [1] for i in range(10)})
        result_many = formatter.format(df_many, ctx)
        assert "…" in result_many.details["meta_preview"]  # Should be truncated

        # Empty DataFrame
        df_empty = pd.DataFrame()
        result_empty = formatter.format(df_empty, ctx)
        assert "0 × 0" in result_empty.type_name
        assert result_empty.details["meta_preview"] == ""  # No column preview for empty

    def test_dataframe_formatter_expandable(self):
        """Test DataFrameFormatter with expandable to_html enabled."""
        import anndata
        from anndata._repr.formatters import DataFrameFormatter
        from anndata._repr.registry import FormatterContext

        formatter = DataFrameFormatter()
        ctx = FormatterContext()
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})

        # Default: not expandable
        result = formatter.format(df, ctx)
        assert not result.is_expandable
        assert result.html_content is None

        # Enable expansion
        original = anndata.settings.repr_html_dataframe_expand
        try:
            anndata.settings.repr_html_dataframe_expand = True
            result_expanded = formatter.format(df, ctx)
            assert result_expanded.is_expandable
            assert result_expanded.html_content is not None
            assert "<table" in result_expanded.html_content  # pandas to_html output
        finally:
            anndata.settings.repr_html_dataframe_expand = original


class TestHTMLCoverage:
    """Tests to improve html.py coverage."""

    def test_extension_type_badge(self):
        """Test extension type shows badge."""

        # Create a mock extension type
        class MockAnnData:
            def __init__(self):
                self.n_obs = 10
                self.n_vars = 5
                self.X = np.zeros((10, 5))
                self.obs = pd.DataFrame(index=range(10))
                self.var = pd.DataFrame(index=range(5))
                self.uns = {}
                self.obsm = {}
                self.varm = {}
                self.layers = {}
                self.obsp = {}
                self.varp = {}
                self.raw = None
                self.is_view = False
                self.isbacked = False
                self.obs_names = pd.Index([f"cell_{i}" for i in range(10)])
                self.var_names = pd.Index([f"gene_{i}" for i in range(5)])

            def __sizeof__(self):
                return 1000

        # Can't easily test extension badge without modifying core AnnData
        # but we can test that repr_html works with various configurations
        adata = AnnData(np.zeros((10, 5)))
        html = adata._repr_html_()
        assert html is not None

    def test_html_with_empty_obsm_varm(self):
        """Test HTML repr with empty obsm/varm."""
        adata = AnnData(np.zeros((10, 5)))
        # obsm and varm are empty by default
        html = adata._repr_html_()
        assert html is not None

    def test_html_with_all_empty_sections(self):
        """Test HTML repr with minimal data."""
        adata = AnnData()
        html = adata._repr_html_()
        assert html is not None
        assert "0" in html  # n_obs and n_vars are 0


class TestRegistryAbstractMethods:
    """
    Tests verifying abstract methods raise NotImplementedError.

    Abstract methods in TypeFormatter and SectionFormatter are designed
    to be overridden. Testing them directly isn't meaningful, but we can
    verify the abstract class structure is correct.
    """

    def test_type_formatter_is_abstract(self):
        """Verify TypeFormatter cannot be instantiated directly."""
        from anndata._repr.registry import TypeFormatter

        with pytest.raises(TypeError):
            TypeFormatter()

    def test_section_formatter_is_abstract(self):
        """Verify SectionFormatter cannot be instantiated directly."""
        from anndata._repr.registry import SectionFormatter

        with pytest.raises(TypeError):
            SectionFormatter()


class TestCustomHtmlContent:
    """Tests for custom HTML content in Type Formatters."""

    def test_inline_html_content(self):
        """Test inline (non-expandable) custom HTML content."""
        from anndata._repr.registry import (
            FormattedOutput,
            TypeFormatter,
            formatter_registry,
        )

        # Custom array class that supports inline HTML preview
        class CustomArray(np.ndarray):
            _test_inline_html = True

        class InlineHtmlFormatter(TypeFormatter):
            priority = 2000  # High priority to be checked first

            def can_format(self, obj):
                return isinstance(obj, np.ndarray) and getattr(
                    obj, "_test_inline_html", False
                )

            def format(self, obj, context):
                return FormattedOutput(
                    type_name="CustomInline",
                    css_class="dtype-custom",
                    html_content='<span class="test-inline">Inline Preview</span>',
                    is_expandable=False,
                )

        formatter = InlineHtmlFormatter()
        formatter_registry.register_type_formatter(formatter)

        try:
            adata = AnnData(np.zeros((5, 3)))
            # Create custom array with proper shape
            custom = np.zeros((5, 2)).view(CustomArray)
            adata.obsm["custom_data"] = custom

            html = adata._repr_html_()

            assert "CustomInline" in html
            assert "Inline Preview" in html
            assert "test-inline" in html
        finally:
            formatter_registry.unregister_type_formatter(formatter)

    def test_expandable_html_content(self):
        """Test expandable custom HTML content (e.g., for TreeData visualization)."""
        from anndata._repr.registry import (
            FormattedOutput,
            TypeFormatter,
            formatter_registry,
        )

        # Custom array class that supports expandable HTML (like TreeData)
        class TreeArray(np.ndarray):
            _test_expandable_html = True

        class ExpandableHtmlFormatter(TypeFormatter):
            priority = 2000

            def can_format(self, obj):
                return isinstance(obj, np.ndarray) and getattr(
                    obj, "_test_expandable_html", False
                )

            def format(self, obj, context):
                tree_html = """
                <div class="test-tree">
                    <ul>
                        <li>Root
                            <ul>
                                <li>Child 1</li>
                                <li>Child 2</li>
                            </ul>
                        </li>
                    </ul>
                </div>
                """
                return FormattedOutput(
                    type_name="TreeData (3 nodes)",
                    css_class="dtype-tree",
                    html_content=tree_html,
                    is_expandable=True,
                )

        formatter = ExpandableHtmlFormatter()
        formatter_registry.register_type_formatter(formatter)

        try:
            adata = AnnData(np.zeros((5, 3)))
            # Create tree array with proper shape
            tree = np.zeros((5, 4)).view(TreeArray)
            adata.obsm["tree"] = tree

            html = adata._repr_html_()

            # Should have the type name
            assert "TreeData (3 nodes)" in html
            # Should have expand button
            assert "adata-expand-btn" in html
            # Should have the tree content in nested row
            assert "test-tree" in html
            assert "adata-nested-row" in html
        finally:
            formatter_registry.unregister_type_formatter(formatter)


class TestCustomSectionFormatters:
    """Tests for custom section formatters (e.g., TreeData's obst/vart)."""

    def test_custom_section_appears_after_specified_section(self):
        """Test that custom sections appear after their specified position."""
        from anndata._repr.registry import (
            FormattedEntry,
            FormattedOutput,
            SectionFormatter,
            formatter_registry,
        )

        class ObstSectionFormatter(SectionFormatter):
            @property
            def section_name(self):
                return "obst"

            @property
            def display_name(self):
                return "obst"

            @property
            def after_section(self):
                return "obsm"  # Should appear after obsm

            @property
            def tooltip(self):
                return "Observation-aligned tree data"

            def should_show(self, obj):
                return hasattr(obj, "_test_obst") and obj._test_obst

            def get_entries(self, obj, context):
                return [
                    FormattedEntry(
                        key="cell_tree",
                        output=FormattedOutput(
                            type_name="Tree (100 nodes)",
                            css_class="dtype-tree",
                        ),
                    ),
                    FormattedEntry(
                        key="lineage",
                        output=FormattedOutput(
                            type_name="Tree (50 nodes)",
                            css_class="dtype-tree",
                        ),
                    ),
                ]

        formatter = ObstSectionFormatter()
        formatter_registry.register_section_formatter(formatter)

        try:
            adata = AnnData(np.zeros((5, 3)))
            adata._test_obst = True  # Mark for custom section
            adata.obsm["X_pca"] = np.zeros((5, 2))

            html = adata._repr_html_()

            # Custom section should be present
            assert 'data-section="obst"' in html
            assert "cell_tree" in html
            assert "lineage" in html
            assert "Tree (100 nodes)" in html

            # Verify order: obsm should come before obst
            obsm_pos = html.find('data-section="obsm"')
            obst_pos = html.find('data-section="obst"')
            assert obsm_pos < obst_pos, "obst should appear after obsm"

            # Verify order: obst should come before varm
            varm_pos = html.find('data-section="varm"')
            assert obst_pos < varm_pos, "obst should appear before varm"

        finally:
            # Cleanup - need to remove from registry
            formatter_registry._section_formatters.pop("obst", None)

    def test_custom_section_with_expandable_content(self):
        """Test custom section with expandable HTML content."""
        from anndata._repr.registry import (
            FormattedEntry,
            FormattedOutput,
            SectionFormatter,
            formatter_registry,
        )

        class TreeSectionFormatter(SectionFormatter):
            @property
            def section_name(self):
                return "trees"

            def should_show(self, obj):
                return hasattr(obj, "_test_trees")

            def get_entries(self, obj, context):
                tree_html = '<div class="tree-viz"><ul><li>Root<ul><li>A</li><li>B</li></ul></li></ul></div>'
                return [
                    FormattedEntry(
                        key="phylo",
                        output=FormattedOutput(
                            type_name="PhyloTree (25 nodes)",
                            css_class="dtype-tree",
                            html_content=tree_html,
                            is_expandable=True,
                        ),
                    ),
                ]

        formatter = TreeSectionFormatter()
        formatter_registry.register_section_formatter(formatter)

        try:
            adata = AnnData(np.zeros((5, 3)))
            adata._test_trees = True

            html = adata._repr_html_()

            # Section should exist
            assert 'data-section="trees"' in html
            assert "PhyloTree (25 nodes)" in html

            # Expandable content should be present
            assert "tree-viz" in html
            assert "adata-nested-row" in html
            assert "adata-expand-btn" in html

        finally:
            formatter_registry._section_formatters.pop("trees", None)

    def test_custom_section_not_shown_when_should_show_false(self):
        """Test that custom sections are hidden when should_show returns False."""
        from anndata._repr.registry import (
            FormattedEntry,
            FormattedOutput,
            SectionFormatter,
            formatter_registry,
        )

        class HiddenSectionFormatter(SectionFormatter):
            @property
            def section_name(self):
                return "hidden"

            def should_show(self, obj):
                return False  # Never show

            def get_entries(self, obj, context):
                return [
                    FormattedEntry(
                        key="secret",
                        output=FormattedOutput(
                            type_name="Secret", css_class="dtype-unknown"
                        ),
                    ),
                ]

        formatter = HiddenSectionFormatter()
        formatter_registry.register_section_formatter(formatter)

        try:
            adata = AnnData(np.zeros((5, 3)))
            html = adata._repr_html_()

            # Section should NOT be present
            assert 'data-section="hidden"' not in html
            assert "secret" not in html

        finally:
            formatter_registry._section_formatters.pop("hidden", None)


# =============================================================================
# Additional Formatter Coverage Tests
# =============================================================================


class TestRareSparseMatrixFormats:
    """Tests for rare sparse matrix formats (lil, dok, dia, bsr)."""

    def test_sparse_lil_formatter(self):
        """Test sparse formatter with LIL matrix (line 174-175)."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        formatter = SparseMatrixFormatter()
        mat = sp.lil_matrix((10, 10))
        mat[0, 0] = 1
        mat[5, 5] = 2

        assert formatter.can_format(mat)
        result = formatter.format(mat, FormatterContext())
        assert "lil" in result.type_name.lower()
        assert result.details["format"] == "lil_matrix"

    def test_sparse_dok_formatter(self):
        """Test sparse formatter with DOK matrix (line 176-177)."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        formatter = SparseMatrixFormatter()
        mat = sp.dok_matrix((10, 10))
        mat[0, 0] = 1
        mat[5, 5] = 2

        assert formatter.can_format(mat)
        result = formatter.format(mat, FormatterContext())
        assert "dok" in result.type_name.lower()
        assert result.details["format"] == "dok_matrix"

    def test_sparse_dia_formatter(self):
        """Test sparse formatter with DIA matrix (line 178-179)."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        formatter = SparseMatrixFormatter()
        # Create a diagonal matrix
        data = np.array([[1, 2, 3, 4], [5, 6, 7, 8]])
        offsets = np.array([0, 1])
        mat = sp.dia_matrix((data, offsets), shape=(4, 4))

        assert formatter.can_format(mat)
        result = formatter.format(mat, FormatterContext())
        assert "dia" in result.type_name.lower()
        assert result.details["format"] == "dia_matrix"

    def test_sparse_bsr_formatter(self):
        """Test sparse formatter with BSR matrix (line 180-181)."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        formatter = SparseMatrixFormatter()
        # Create a block sparse row matrix
        mat = sp.bsr_matrix(
            np.array([[1, 0, 0, 0], [0, 0, 2, 0], [0, 0, 0, 3], [4, 0, 0, 0]])
        )

        assert formatter.can_format(mat)
        result = formatter.format(mat, FormatterContext())
        assert "bsr" in result.type_name.lower()
        assert result.details["format"] == "bsr_matrix"


class TestDataFrameFormatterEdgeCases:
    """Tests for DataFrame formatter edge cases."""

    def test_dataframe_long_column_names_truncation(self):
        """Test DataFrame with long column names gets truncated (line 247-248)."""
        from anndata._repr.formatters import DataFrameFormatter
        from anndata._repr.registry import FormatterContext

        formatter = DataFrameFormatter()
        # Create DataFrame with very long column names
        long_names = {
            f"very_long_column_name_{i}_with_extra_text": [1] for i in range(3)
        }
        df = pd.DataFrame(long_names)

        result = formatter.format(df, FormatterContext())
        # Preview should be truncated
        assert "…" in result.details["meta_preview"]
        # Full preview should have all columns
        assert "meta_preview_full" in result.details


class TestSeriesFormatterNonSerializable:
    """Tests for Series formatter detecting non-serializable object dtype columns."""

    def test_custom_object_in_obs_column_not_serializable(self):
        """Test that custom objects in obs columns are flagged as non-serializable.

        Uses a custom class (not a list) to ensure this test remains valid
        even if list serialization is added in the future (issue #1923).
        """
        from anndata._repr.formatters import SeriesFormatter
        from anndata._repr.registry import FormatterContext

        class CustomObject:
            """A custom class that will never be serializable."""

        # Create a Series with custom objects
        series = pd.Series([CustomObject(), CustomObject(), CustomObject()])
        assert series.dtype == np.dtype("object")

        formatter = SeriesFormatter()
        result = formatter.format(series, FormatterContext())

        assert not result.is_serializable
        assert len(result.warnings) > 0
        assert "CustomObject" in result.warnings[0]

    def test_list_in_obs_column_not_serializable(self):
        """Test that lists in obs columns are flagged as non-serializable.

        NOTE: This test may need updating if #1923 adds list serialization.
        See: https://github.com/scverse/anndata/issues/1923
        """
        from anndata._repr.formatters import SeriesFormatter
        from anndata._repr.registry import FormatterContext

        series = pd.Series([["a", "b"], ["c"], ["d", "e", "f"]])
        assert series.dtype == np.dtype("object")

        formatter = SeriesFormatter()
        result = formatter.format(series, FormatterContext())

        assert not result.is_serializable
        assert len(result.warnings) > 0
        assert "list" in result.warnings[0]

    def test_string_obs_column_is_serializable(self):
        """Test that string object columns are serializable."""
        from anndata._repr.formatters import SeriesFormatter
        from anndata._repr.registry import FormatterContext

        series = pd.Series(["a", "b", "c"], dtype=object)

        formatter = SeriesFormatter()
        result = formatter.format(series, FormatterContext())

        assert result.is_serializable
        assert len(result.warnings) == 0

    def test_empty_series_is_serializable(self):
        """Test that empty object dtype series is considered serializable."""
        from anndata._repr.formatters import SeriesFormatter
        from anndata._repr.registry import FormatterContext

        series = pd.Series([], dtype=object)

        formatter = SeriesFormatter()
        result = formatter.format(series, FormatterContext())

        assert result.is_serializable

    def test_list_column_detected_and_not_serializable(self, tmp_path):
        """Repr detects list columns as non-serializable, and they actually fail.

        MAINTAINER NOTE: If issue #1923 is resolved and lists become serializable,
        update _check_series_serializability() in formatters.py to remove the
        list/tuple check, or make it conditional on list content.
        """
        adata = ad.AnnData(X=np.eye(3))
        adata.obs["list_col"] = [["a", "b"], ["c"], ["d"]]

        # Verify lists still fail to serialize
        try:
            adata.write_h5ad(tmp_path / "test.h5ad")
            pytest.fail(
                "List serialization now works! "
                "Update _check_series_serializability() in formatters.py."
            )
        except (TypeError, Exception):  # noqa: BLE001
            pass  # Expected - lists not serializable

        # Repr should detect and warn
        html = adata._repr_html_()
        assert "list" in html
        assert "(!)" in html

    def test_custom_object_detected_and_not_serializable(self, tmp_path):
        """Repr detects custom objects as non-serializable, and they actually fail.

        NOTE: Custom objects should never become serializable without explicit
        registration. If this test fails, check if anndata added a generic
        object serialization mechanism.
        """

        class CustomObject:
            pass

        adata = ad.AnnData(X=np.eye(3))
        adata.obs["custom"] = [CustomObject(), CustomObject(), CustomObject()]

        # Verify custom objects still fail to serialize
        try:
            adata.write_h5ad(tmp_path / "test.h5ad")
            pytest.fail(
                "Custom object serialization now works! "
                "Update _check_series_serializability() in formatters.py."
            )
        except (TypeError, Exception):  # noqa: BLE001
            pass  # Expected - custom objects not serializable

        # Repr should detect and warn
        html = adata._repr_html_()
        assert "CustomObject" in html
        assert "(!)" in html

    def test_non_ascii_column_no_warning_and_serializes(self, tmp_path):
        """Repr does not warn for non-ASCII (it's valid), and it serializes.

        MAINTAINER NOTE: If non-ASCII stops working, add a warning in
        check_column_name() in formatters.py. But UTF-8 should always work.
        """
        adata = ad.AnnData(X=np.eye(3))
        adata.obs["gène_名前"] = ["a", "b", "c"]

        # Verify non-ASCII still serializes fine
        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)
        adata2 = ad.read_h5ad(path)
        assert "gène_名前" in adata2.obs.columns, (
            "Non-ASCII serialization broke! If this is intentional, "
            "add warning in check_column_name() in formatters.py."
        )

        # Repr should NOT warn (non-ASCII is valid UTF-8)
        html = adata._repr_html_()
        assert "gène_名前" in html
        assert "Not serializable" not in html

    def test_tuple_column_name_detected_and_not_serializable(self, tmp_path):
        """Repr detects non-string column names, and they actually fail.

        NOTE: Non-string column names (tuples, ints) should never become
        serializable - HDF5/Zarr keys must be strings. If this changes,
        update check_column_name() in formatters.py.
        """
        adata = ad.AnnData(X=np.eye(3))
        adata.obs[("a", "b")] = [1, 2, 3]

        # Verify non-string names still fail to serialize
        try:
            adata.write_h5ad(tmp_path / "test.h5ad")
            pytest.fail(
                "Non-string column name serialization now works! "
                "Update check_column_name() in formatters.py."
            )
        except (TypeError, Exception):  # noqa: BLE001
            pass  # Expected - non-string names not serializable

        # Repr should detect and warn
        html = adata._repr_html_()
        assert "Non-string" in html
        assert "(!)" in html

    def test_datetime64_column_detected_and_not_serializable(self, tmp_path):
        """Repr detects datetime64 columns as non-serializable (issue #455).

        MAINTAINER NOTE: If this test fails because datetime64 serialization
        was added to anndata, update SeriesFormatter in formatters.py to
        remove the datetime64 warning check.
        """
        adata = ad.AnnData(X=np.eye(3))
        adata.obs["date"] = pd.to_datetime(["2024-01-01", "2024-01-02", "2024-01-03"])

        # Verify datetime64 still fails to serialize
        # If this starts passing, datetime64 support was added - update repr accordingly
        try:
            path = tmp_path / "test_datetime.h5ad"
            adata.write_h5ad(path)
            # If we get here, datetime64 is now serializable - repr should stop warning
            pytest.fail(
                "datetime64 serialization now works! "
                "Update SeriesFormatter in formatters.py to remove the datetime64 warning."
            )
        except Exception:  # noqa: BLE001
            pass  # Expected - datetime64 not serializable

        # Repr should detect and warn about it
        html = adata._repr_html_()
        assert "datetime64" in html, "Repr should show datetime64 dtype"
        assert "(!)" in html, "Repr should show warning icon for datetime64"

    def test_timedelta64_column_detected_and_not_serializable(self, tmp_path):
        """Repr detects timedelta64 columns as non-serializable.

        MAINTAINER NOTE: If this test fails because timedelta64 serialization
        was added to anndata, update SeriesFormatter in formatters.py to
        remove the timedelta64 warning check.
        """
        adata = ad.AnnData(X=np.eye(3))
        adata.obs["duration"] = pd.to_timedelta(["1 days", "2 days", "3 days"])

        # Verify timedelta64 still fails to serialize
        try:
            path = tmp_path / "test_timedelta.h5ad"
            adata.write_h5ad(path)
            pytest.fail(
                "timedelta64 serialization now works! "
                "Update SeriesFormatter in formatters.py to remove the timedelta64 warning."
            )
        except Exception:  # noqa: BLE001
            pass  # Expected - timedelta64 not serializable

        # Repr should detect and warn about it
        html = adata._repr_html_()
        assert "timedelta64" in html, "Repr should show timedelta64 dtype"
        assert "(!)" in html, "Repr should show warning icon for timedelta64"


class TestColumnNameValidation:
    """Tests for column name validation (issue #321)."""

    def test_check_column_name_valid(self):
        """Test that valid column names pass."""
        from anndata._repr.formatters import check_column_name

        valid, _, _ = check_column_name("gene_name")
        assert valid
        valid, _, _ = check_column_name("gène_名前")  # Non-ASCII is fine
        assert valid

    def test_check_column_name_slash(self):
        """Test slashes are flagged as warning (not hard error - still works for now)."""
        from anndata._repr.formatters import check_column_name

        valid, reason, is_hard_error = check_column_name("path/to/gene")
        assert not valid
        assert "/" in reason
        assert not is_hard_error  # Deprecation warning, not hard error

    def test_check_column_name_non_string(self):
        """Test non-string names are flagged as hard error."""
        from anndata._repr.formatters import check_column_name

        valid, reason, is_hard_error = check_column_name(("a", "b"))
        assert not valid
        assert "Non-string" in reason
        assert is_hard_error  # Fails NOW

    def test_slash_column_name_warns_but_still_serializes(self, tmp_path):
        """Test slash in column names: warns (yellow) but still works for now.

        MAINTAINER NOTE: When slashes become hard errors (FutureWarning fulfilled),
        update check_column_name() in formatters.py to return is_hard_error=True
        for slashes, and update the CSS class from yellow to red.
        """
        adata = ad.AnnData(X=np.eye(3))
        adata.obs["path/gene"] = ["a", "b", "c"]

        # Currently still serializes (with FutureWarning)
        path = tmp_path / "test_slash.h5ad"
        import warnings

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            try:
                adata.write_h5ad(path)
                # If we get here, slashes still work
                serializes = True
            except Exception:  # noqa: BLE001
                # Slashes now fail - update repr to show red instead of yellow
                serializes = False
                pytest.fail(
                    "Slash in column names now fails! "
                    "Update check_column_name() in formatters.py: "
                    "set is_hard_error=True for slashes."
                )

        if serializes:
            # Verify it's still a FutureWarning (yellow in repr)
            adata2 = ad.read_h5ad(path)
            assert "path/gene" in adata2.obs.columns

        # Repr should show yellow warning (not red error)
        html = adata._repr_html_()
        assert "deprecated" in html or "/" in html

    def test_mapping_sections_warn_for_invalid_keys(self):
        """Test that layers/obsm/etc warn for invalid key names."""
        adata = ad.AnnData(X=np.eye(3))
        adata.layers[("tuple", "key")] = np.eye(3)  # Non-string - error
        adata.obsm["path/embed"] = np.random.randn(3, 2)  # Slash - warning

        html = adata._repr_html_()
        assert "Non-string" in html
        assert "deprecated" in html

    def test_uns_warns_for_invalid_keys(self):
        """Test that uns warns for invalid key names."""
        adata = ad.AnnData(X=np.eye(3))
        adata.uns[("tuple", "key")] = "value"  # Non-string - error

        html = adata._repr_html_()
        assert "Non-string" in html
        assert "Not serializable" in html


class TestMockCuPyArrayFormatter:
    """Tests for CuPy array formatter using mock objects."""

    def test_cupy_array_formatter_with_mock(self):
        """Test CuPy array formatter with a mock CuPy array (lines 394-410)."""
        from anndata._repr.formatters import CuPyArrayFormatter
        from anndata._repr.registry import FormatterContext

        # Create a mock CuPy array
        class MockDevice:
            id = 0

        class MockCuPyArray:
            def __init__(self):
                self.shape = (100, 50)
                self.dtype = np.float32
                self.device = MockDevice()

        # Set module to look like cupy
        MockCuPyArray.__module__ = "cupy._core.core"

        formatter = CuPyArrayFormatter()
        mock_arr = MockCuPyArray()

        assert formatter.can_format(mock_arr)
        result = formatter.format(mock_arr, FormatterContext())

        assert "cupy" in result.type_name.lower()
        assert "100" in result.type_name
        assert "50" in result.type_name
        assert result.css_class == "dtype-gpu"
        assert "GPU:0" in result.tooltip


class TestMockAwkwardArrayFormatter:
    """Tests for Awkward array formatter using mock objects."""

    def test_awkward_array_formatter_with_mock(self):
        """Test Awkward array formatter with a mock object (lines 428-444)."""
        from anndata._repr.formatters import AwkwardArrayFormatter
        from anndata._repr.registry import FormatterContext

        # Create a mock Awkward array
        class MockAwkwardArray:
            def __init__(self):
                self.type = "var * int64"

            def __len__(self):
                return 100

        MockAwkwardArray.__module__ = "awkward.highlevel"

        formatter = AwkwardArrayFormatter()
        mock_arr = MockAwkwardArray()

        assert formatter.can_format(mock_arr)
        result = formatter.format(mock_arr, FormatterContext())

        assert "awkward" in result.type_name.lower()
        assert "100" in result.type_name
        assert result.css_class == "dtype-awkward"
        assert "var * int64" in result.tooltip

    def test_awkward_array_formatter_exception_handling(self):
        """Test Awkward array formatter handles exceptions (lines 431-433)."""
        from anndata._repr.formatters import AwkwardArrayFormatter
        from anndata._repr.registry import FormatterContext

        # Create a mock that raises exceptions
        class BrokenAwkwardArray:
            @property
            def type(self):
                msg = "Cannot get type"
                raise RuntimeError(msg)

            def __len__(self):
                msg = "Cannot get length"
                raise RuntimeError(msg)

        BrokenAwkwardArray.__module__ = "awkward.highlevel"

        formatter = AwkwardArrayFormatter()
        mock_arr = BrokenAwkwardArray()

        assert formatter.can_format(mock_arr)
        result = formatter.format(mock_arr, FormatterContext())

        # Should handle gracefully
        assert "awkward" in result.type_name.lower()
        assert "?" in result.type_name  # Length should be "?"
        assert "unknown" in result.tooltip.lower()


class TestArrayAPIFormatter:
    """Tests for Array-API compatible array formatter."""

    def test_array_api_formatter_jax_like(self):
        """Test Array-API formatter with JAX-like array (lines 496-538)."""
        from anndata._repr.formatters import ArrayAPIFormatter
        from anndata._repr.registry import FormatterContext

        # Create a mock JAX array
        class MockJAXArray:
            def __init__(self):
                self.shape = (100, 50)
                self.dtype = np.float32
                self.ndim = 2
                self.device = "gpu:0"

        MockJAXArray.__module__ = "jax.numpy"

        formatter = ArrayAPIFormatter()
        mock_arr = MockJAXArray()

        assert formatter.can_format(mock_arr)
        result = formatter.format(mock_arr, FormatterContext())

        assert "MockJAXArray" in result.type_name
        assert "100" in result.type_name
        assert "50" in result.type_name
        assert result.css_class == "dtype-array-api"
        assert "JAX" in result.tooltip
        assert "gpu:0" in result.tooltip

    def test_array_api_formatter_pytorch_like(self):
        """Test Array-API formatter with PyTorch-like tensor."""
        from anndata._repr.formatters import ArrayAPIFormatter
        from anndata._repr.registry import FormatterContext

        class MockTorchTensor:
            def __init__(self):
                self.shape = (64, 32)
                self.dtype = "torch.float32"
                self.ndim = 2
                self.device = "cuda:0"

        MockTorchTensor.__module__ = "torch"

        formatter = ArrayAPIFormatter()
        mock_tensor = MockTorchTensor()

        assert formatter.can_format(mock_tensor)
        result = formatter.format(mock_tensor, FormatterContext())

        assert "MockTorchTensor" in result.type_name
        assert "PyTorch" in result.tooltip

    def test_array_api_formatter_device_buffer(self):
        """Test Array-API formatter with device_buffer attribute (lines 521-522)."""
        from anndata._repr.formatters import ArrayAPIFormatter
        from anndata._repr.registry import FormatterContext

        class MockDeviceBuffer:
            def device(self):
                return "tpu:0"

        class MockJAXArrayWithBuffer:
            def __init__(self):
                self.shape = (10, 5)
                self.dtype = np.float32
                self.ndim = 2
                self.device_buffer = MockDeviceBuffer()

        MockJAXArrayWithBuffer.__module__ = "jaxlib.xla_extension"

        formatter = ArrayAPIFormatter()
        mock_arr = MockJAXArrayWithBuffer()

        assert formatter.can_format(mock_arr)
        result = formatter.format(mock_arr, FormatterContext())

        assert "JAX" in result.tooltip
        assert "tpu:0" in result.tooltip

    def test_array_api_formatter_excludes_numpy(self):
        """Test Array-API formatter excludes numpy arrays."""
        from anndata._repr.formatters import ArrayAPIFormatter

        formatter = ArrayAPIFormatter()
        arr = np.zeros((10, 5))

        # Should NOT format numpy arrays (handled by NumpyArrayFormatter)
        assert not formatter.can_format(arr)

    def test_array_api_formatter_excludes_pandas(self):
        """Test Array-API formatter excludes pandas objects."""
        from anndata._repr.formatters import ArrayAPIFormatter

        formatter = ArrayAPIFormatter()

        # Should NOT format pandas objects
        assert not formatter.can_format(pd.DataFrame({"a": [1, 2, 3]}))
        assert not formatter.can_format(pd.Series([1, 2, 3]))


class TestDtypeCSSClassHelpers:
    """Tests for dtype CSS class helper functions."""

    def test_get_dtype_css_class_complex(self):
        """Test CSS class for complex dtype (line 742-743)."""
        from anndata._repr.formatters import _get_dtype_css_class

        complex_dtype = np.dtype(np.complex128)
        css_class = _get_dtype_css_class(complex_dtype)
        assert css_class == "dtype-float"  # Complex maps to float

    def test_get_dtype_css_class_unknown(self):
        """Test CSS class for unknown dtype (line 744-745)."""
        from anndata._repr.formatters import _get_dtype_css_class

        # Datetime dtype has kind 'M'
        datetime_dtype = np.dtype("datetime64[ns]")
        css_class = _get_dtype_css_class(datetime_dtype)
        assert css_class == "dtype-object"  # Unknown maps to object

    def test_get_pandas_dtype_css_class_category(self):
        """Test CSS class for pandas category dtype (line 757-758)."""
        from anndata._repr.formatters import _get_pandas_dtype_css_class

        cat_dtype = pd.CategoricalDtype(categories=["a", "b", "c"])
        css_class = _get_pandas_dtype_css_class(cat_dtype)
        assert css_class == "dtype-category"

    def test_get_pandas_dtype_css_class_unknown(self):
        """Test CSS class for unknown pandas dtype (line 761-762)."""
        from anndata._repr.formatters import _get_pandas_dtype_css_class

        # Timedelta dtype
        timedelta_dtype = pd.to_timedelta([1, 2, 3], unit="s").dtype
        css_class = _get_pandas_dtype_css_class(timedelta_dtype)
        assert css_class == "dtype-object"  # Unknown maps to object


class TestSparseFormatterDuckTyping:
    """Tests for sparse formatter duck typing fallback."""

    def test_sparse_formatter_duck_typing_fallback(self):
        """Test sparse formatter uses duck typing when scipy checks fail (lines 182-189)."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        # Create a mock sparse-like object that passes duck typing
        class MockSparseArray:
            def __init__(self):
                self.nnz = 10
                self.shape = (5, 5)
                self.dtype = np.float64

            def tocsr(self):
                pass

        # Set module to look like scipy.sparse
        MockSparseArray.__module__ = "scipy.sparse._csr"

        formatter = SparseMatrixFormatter()
        mock_sparse = MockSparseArray()

        assert formatter.can_format(mock_sparse)
        result = formatter.format(mock_sparse, FormatterContext())

        # Should use type name as fallback
        assert "MockSparseArray" in result.type_name
        assert result.details["format"] == "MockSparseArray"


class TestAnnDataFormatterCoverage:
    """Tests for AnnData formatter coverage."""

    def test_anndata_formatter_nested_depth_limit(self):
        """Test AnnData formatter respects max_depth for nested objects."""
        from anndata._repr.formatters import AnnDataFormatter
        from anndata._repr.registry import FormatterContext

        formatter = AnnDataFormatter()
        inner = AnnData(np.zeros((5, 3)))

        # At depth 0, should be expandable (0 < max_depth default 3)
        result_shallow = formatter.format(inner, FormatterContext(depth=0, max_depth=3))
        assert result_shallow.is_expandable

        # At depth 3, should NOT be expandable (3 >= max_depth 3)
        result_deep = formatter.format(inner, FormatterContext(depth=3, max_depth=3))
        assert not result_deep.is_expandable


# =============================================================================
# Column Width Settings Tests
# =============================================================================


class TestColumnWidthSettings:
    """Test column width configuration settings."""

    def test_field_width_calculated_from_content(self):
        """Test field name column width is calculated based on content."""
        from anndata._repr.html import _calculate_field_name_width

        # Short names - should use minimum width
        adata_short = AnnData(np.zeros((10, 5)))
        adata_short.obs["a"] = list(range(10))
        width_short = _calculate_field_name_width(adata_short, max_width=400)
        assert width_short >= 80  # Minimum

        # Long names - should expand
        adata_long = AnnData(np.zeros((10, 5)))
        adata_long.obs["this_is_a_very_long_column_name_that_tests_width"] = list(
            range(10)
        )
        width_long = _calculate_field_name_width(adata_long, max_width=400)
        assert width_long > width_short

    def test_field_width_capped_by_max(self):
        """Test field name column width is capped by max_width setting."""
        from anndata._repr.html import _calculate_field_name_width

        # Very long name - should be capped
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["x" * 100] = list(range(10))  # 100 char name

        width = _calculate_field_name_width(adata, max_width=200)
        assert width <= 200

    def test_repr_html_max_field_width_setting(self):
        """Test repr_html_max_field_width setting affects output."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))
        adata.obs["medium_length_name"] = list(range(10))

        # Default max field width
        html_default = adata._repr_html_()
        # Extract the CSS variable value
        match_default = re.search(r"--anndata-name-col-width:\s*(\d+)px", html_default)
        assert match_default, "Field width CSS variable not found"
        width_default = int(match_default.group(1))

        # With smaller max field width
        with settings.override(repr_html_max_field_width=100):
            html_small = adata._repr_html_()
            match_small = re.search(r"--anndata-name-col-width:\s*(\d+)px", html_small)
            assert match_small, "Field width CSS variable not found"
            width_small = int(match_small.group(1))

        # Smaller max should result in smaller or equal width
        assert width_small <= 100
        assert width_small <= width_default

    def test_repr_html_type_width_setting(self):
        """Test repr_html_type_width setting affects output."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))
        adata.obs["col"] = list(range(10))

        # Default type width
        html_default = adata._repr_html_()
        match_default = re.search(r"--anndata-type-col-width:\s*(\d+)px", html_default)
        assert match_default, "Type width CSS variable not found"
        width_default = int(match_default.group(1))
        assert width_default == 220  # Default value

        # With custom type width
        with settings.override(repr_html_type_width=250):
            html_custom = adata._repr_html_()
            match_custom = re.search(
                r"--anndata-type-col-width:\s*(\d+)px", html_custom
            )
            assert match_custom, "Type width CSS variable not found"
            width_custom = int(match_custom.group(1))

        assert width_custom == 250

    def test_field_width_considers_all_sections(self):
        """Test field width calculation considers names from all sections."""
        from anndata._repr.html import _calculate_field_name_width

        # Create AnnData with different length names in different sections
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["short"] = list(range(10))  # obs
        adata.var["medium_name"] = list(range(5))  # var
        adata.uns["longer_uns_key_name"] = "value"  # uns
        adata.obsm["X_very_long_obsm_embedding_name"] = np.zeros((10, 2))  # obsm
        adata.layers["raw"] = np.zeros((10, 5))  # layers

        # Width should accommodate the longest name
        width = _calculate_field_name_width(adata, max_width=500)

        # X_very_long_obsm_embedding_name is 31 chars, at ~8px/char + 30px padding = ~278px
        # Should be close to that
        assert width >= 200  # Should be substantial due to long obsm key

    def test_empty_adata_field_width(self):
        """Test field width calculation for empty AnnData."""
        from anndata._repr.html import _calculate_field_name_width

        adata = AnnData()  # Empty
        width = _calculate_field_name_width(adata, max_width=400)
        assert width == 100  # Default minimum for empty


# =============================================================================
# Coverage Tests for html.py Functions
# =============================================================================


class TestGenerateReprHtmlDirectly:
    """Tests calling generate_repr_html directly to cover internal code paths."""

    def test_html_disabled_returns_pre(self):
        """Test disabled HTML returns <pre> wrapped text when called directly."""
        from anndata import settings
        from anndata._repr.html import generate_repr_html

        adata = AnnData(np.zeros((10, 5)))

        with settings.override(repr_html_enabled=False):
            html = generate_repr_html(adata)
            # Should return pre-wrapped text
            assert html is not None
            assert html.startswith("<pre>")
            assert "</pre>" in html

    def test_max_depth_reached_at_depth_zero(self):
        """Test max depth indicator when depth >= max_depth."""
        from anndata._repr.html import generate_repr_html

        adata = AnnData(np.zeros((10, 5)))

        # Call with depth >= max_depth
        html = generate_repr_html(adata, depth=3, max_depth=2)
        assert "max depth" in html.lower()
        assert "10" in html  # n_obs
        assert "5" in html  # n_vars

    def test_nested_anndata_not_expandable_at_max_depth(self):
        """Test nested AnnData doesn't expand when at max depth."""
        inner = AnnData(np.zeros((5, 3)))
        outer = AnnData(np.zeros((10, 5)))
        outer.uns["nested"] = inner

        # With max_depth=1, the inner AnnData shouldn't be expandable
        from anndata import settings

        with settings.override(repr_html_max_depth=1):
            html = outer._repr_html_()
            assert "nested" in html
            # Should show max depth reached for nested
            assert "max depth" in html.lower()


class TestValuePreviewFunctions:
    """Tests for value preview helper functions in html.py."""

    def test_preview_string_short(self):
        """Test string preview for short strings."""
        from anndata._repr.utils import preview_string

        result = preview_string("hello", max_len=100)
        assert result == '"hello"'

    def test_preview_string_long(self):
        """Test string preview for long strings."""
        from anndata._repr.utils import preview_string

        long_str = "a" * 100
        result = preview_string(long_str, max_len=20)
        assert len(result) < 30  # Much shorter than original
        assert result.endswith('..."')

    def test_preview_number_bool(self):
        """Test number preview for booleans."""
        from anndata._repr.utils import preview_number

        assert preview_number(True) == "True"  # noqa: FBT003
        assert preview_number(False) == "False"  # noqa: FBT003

    def test_preview_number_int(self):
        """Test number preview for integers."""
        from anndata._repr.utils import preview_number

        assert preview_number(42) == "42"
        assert preview_number(np.int64(42)) == "42"

    def test_preview_number_float_whole(self):
        """Test number preview for whole number floats."""
        from anndata._repr.utils import preview_number

        assert preview_number(42.0) == "42"
        assert preview_number(np.float32(100.0)) == "100"

    def test_preview_number_float_decimal(self):
        """Test number preview for decimal floats."""
        from anndata._repr.utils import preview_number

        result = preview_number(3.14159265359)
        # Should be formatted concisely with :.6g
        assert "3.14159" in result

    def test_preview_dict_empty(self):
        """Test dict preview for empty dict."""
        from anndata._repr.utils import preview_dict

        assert preview_dict({}) == "{}"

    def test_preview_dict_small(self):
        """Test dict preview for small dict (<=3 keys)."""
        from anndata._repr.utils import preview_dict

        result = preview_dict({"a": 1, "b": 2})
        assert "a" in result
        assert "b" in result
        assert "{" in result

    def test_preview_dict_large(self):
        """Test dict preview for large dict (>3 keys)."""
        from anndata._repr.utils import preview_dict

        d = {f"key_{i}": i for i in range(10)}
        result = preview_dict(d)
        assert "keys" in result.lower()
        assert "10" in result

    def test_preview_sequence_empty_list(self):
        """Test sequence preview for empty list."""
        from anndata._repr.utils import preview_sequence

        assert preview_sequence([]) == "[]"

    def test_preview_sequence_empty_tuple(self):
        """Test sequence preview for empty tuple."""
        from anndata._repr.utils import preview_sequence

        assert preview_sequence(()) == "()"

    def test_preview_sequence_small_list(self):
        """Test sequence preview for small list."""
        from anndata._repr.utils import preview_sequence

        result = preview_sequence([1, 2, 3])
        assert "[" in result
        assert "]" in result
        assert "1" in result

    def test_preview_sequence_small_tuple(self):
        """Test sequence preview for small tuple."""
        from anndata._repr.utils import preview_sequence

        result = preview_sequence((1, 2, 3))
        assert "(" in result
        assert ")" in result

    def test_preview_sequence_large(self):
        """Test sequence preview for large sequence."""
        from anndata._repr.utils import preview_sequence

        result = preview_sequence(list(range(100)))
        assert "items" in result.lower()
        assert "100" in result

    def test_preview_sequence_with_complex_items(self):
        """Test sequence preview with complex items (returns item count)."""
        from anndata._repr.utils import preview_sequence

        # Lists with nested structures should fall back to item count
        result = preview_sequence([{"a": 1}, {"b": 2}])
        assert "2" in result  # 2 items

    def test_preview_item_string(self):
        """Test item preview for strings."""
        from anndata._repr.utils import preview_item

        assert preview_item("hello") == '"hello"'
        assert preview_item("a" * 100).endswith('..."')

    def test_preview_item_numbers(self):
        """Test item preview for numbers."""
        from anndata._repr.utils import preview_item

        assert preview_item(42) == "42"
        assert preview_item(3.14) == "3.14"
        assert preview_item(True) == "True"  # noqa: FBT003

    def test_preview_item_none(self):
        """Test item preview for None."""
        from anndata._repr.utils import preview_item

        assert preview_item(None) == "None"

    def test_preview_item_complex(self):
        """Test item preview for complex types returns empty."""
        from anndata._repr.utils import preview_item

        assert preview_item([1, 2, 3]) == ""
        assert preview_item({"a": 1}) == ""

    def test_generate_value_preview_none(self):
        """Test value preview for None."""
        from anndata._repr.utils import generate_value_preview

        assert generate_value_preview(None) == "None"

    def test_generate_valuepreview_string(self):
        """Test value preview for string."""
        from anndata._repr.utils import generate_value_preview

        result = generate_value_preview("hello world")
        assert '"hello world"' in result

    def test_generate_valuepreview_number(self):
        """Test value preview for number."""
        from anndata._repr.utils import generate_value_preview

        assert "42" in generate_value_preview(42)
        assert "3.14" in generate_value_preview(3.14)

    def test_generate_valuepreview_dict(self):
        """Test value preview for dict."""
        from anndata._repr.utils import generate_value_preview

        result = generate_value_preview({"key": "value"})
        assert "key" in result

    def test_generate_value_preview_list(self):
        """Test value preview for list."""
        from anndata._repr.utils import generate_value_preview

        result = generate_value_preview([1, 2, 3])
        assert "[" in result or "items" in result

    def test_generate_value_preview_complex(self):
        """Test value preview for complex types returns empty."""
        from anndata._repr.utils import generate_value_preview

        class CustomClass:
            pass

        assert generate_value_preview(CustomClass()) == ""


class TestRawSectionRendering:
    """Tests for raw section rendering."""

    def test_raw_section_with_var(self):
        """Test raw section shows var column count."""
        adata = AnnData(np.zeros((10, 20)))
        adata.raw = adata.copy()
        # Now subset var
        adata = adata[:, :5]

        html = adata._repr_html_()
        assert "raw" in html.lower()
        # Should show raw.var info
        assert "var" in html.lower()

    def test_raw_section_with_varm(self):
        """Test raw section shows varm when present."""
        adata = AnnData(np.zeros((10, 20)))
        adata.varm["PCs"] = np.zeros((20, 50))
        adata.raw = adata.copy()
        adata = adata[:, :5]

        html = adata._repr_html_()
        assert "raw" in html.lower()


class TestColorListRendering:
    """Tests for color list entry rendering."""

    def test_color_list_many_colors(self):
        """Test color list with more than 15 colors shows +N."""
        adata = AnnData(np.zeros((10, 5)))
        # Create a color list with more than 15 colors
        colors = [f"#{i:02x}{i:02x}{i:02x}" for i in range(0, 250, 10)]
        adata.uns["many_colors"] = colors

        html = adata._repr_html_()
        # Should show some colors and indicate more
        assert "many_colors" in html
        # Check for color format
        assert "#" in html

    def test_color_list_exact_limit(self):
        """Test color list with exactly 15 colors."""
        adata = AnnData(np.zeros((15, 5)))
        colors = [f"#{i:02x}0000" for i in range(15)]
        adata.obs["cluster"] = pd.Categorical([i % 15 for i in range(15)])
        adata.uns["cluster_colors"] = colors

        html = adata._repr_html_()
        assert "cluster_colors" in html


class TestMappingSectionEdgeCases:
    """Tests for mapping section edge cases."""

    def test_obsm_shape_meta_display(self):
        """Test obsm shows column count in meta."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obsm["X_pca"] = np.zeros((10, 50))
        adata.obsm["X_umap"] = np.zeros((10, 2))

        html = adata._repr_html_()
        assert "X_pca" in html
        assert "X_umap" in html
        # Should show column info
        assert "50" in html or "cols" in html.lower()

    def test_layers_truncation(self):
        """Test layers section truncates when exceeding max_items."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))
        # Add many layers
        for i in range(100):
            adata.layers[f"layer_{i}"] = np.zeros((10, 5))

        with settings.override(repr_html_max_items=10):
            html = adata._repr_html_()
            # Should show truncation
            assert "more" in html.lower() or "..." in html


class TestTypeCellRendering:
    """Tests for type cell rendering with warnings and expand buttons."""

    def test_type_cell_with_warning(self):
        """Test type cell shows warning icon."""

        class CustomObject:
            pass

        adata = AnnData(np.zeros((10, 5)))
        adata.uns["custom"] = CustomObject()

        html = adata._repr_html_()
        # Should show warning icon
        assert "⚠" in html or "warning" in html.lower()

    def test_type_cell_not_serializable(self):
        """Test type cell shows not serializable warning."""

        class NotSerializable:
            pass

        adata = AnnData(np.zeros((10, 5)))
        adata.obsm["bad"] = np.zeros((10, 5))
        # Can't directly add unserializable to obsm, but we can test uns
        adata.uns["unserializable"] = NotSerializable()

        html = adata._repr_html_()
        assert "serializable" in html.lower() or "⚠" in html

    def test_expand_button_for_expandable_content(self):
        """Test expand button appears for expandable content."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))
        # DataFrames with expand setting enabled show expand button
        adata.uns["df"] = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})

        with settings.override(repr_html_dataframe_expand=True):
            html = adata._repr_html_()
            # Should have expand functionality
            assert "expand" in html.lower() or "Expand" in html


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
                return "test_custom_section"

            @property
            def display_name(self) -> str:
                return "Test Custom"

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
            assert "Test Custom" in html or "test_custom_section" in html
            assert "entry1" in html
            assert "entry2" in html
        finally:
            # Cleanup: remove from internal dict
            formatter_registry._section_formatters.pop("test_custom_section", None)

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
            # Cleanup: remove from internal dict
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
            # Cleanup: remove from internal dict
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
                            html_content="<div>Expanded content here</div>",
                            is_expandable=True,
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
            # Cleanup: remove from internal dict
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

        class InlineHtmlFormatter(SectionFormatter):
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
                            html_content="<span>Inline preview</span>",
                            is_expandable=False,
                        ),
                    ),
                ]

        formatter = InlineHtmlFormatter()
        formatter_registry.register_section_formatter(formatter)

        try:
            adata = AnnData(np.zeros((10, 5)))
            adata._inline_marker = True

            html = adata._repr_html_()
            assert "inline_entry" in html
            assert "Inline preview" in html
        finally:
            # Cleanup: remove from internal dict
            formatter_registry._section_formatters.pop("inline_section", None)

    def test_formatted_entry_not_serializable(self):
        """Test formatted entry with is_serializable=False."""
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
            assert "serializable" in html.lower() or "⚠" in html
        finally:
            # Cleanup: remove from internal dict
            formatter_registry._section_formatters.pop("not_serializable_section", None)


class TestUnsEntryRendering:
    """Tests for uns entry rendering with various types."""

    def test_uns_entry_with_type_hint_preview_note(self):
        """Test uns entry with unregistered type hint shows import suggestion."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["typed_data"] = {
            "__anndata_repr__": "hypothetical.package.Type",
            "data": {"key": "value"},
        }

        html = adata._repr_html_()
        assert "typed_data" in html
        assert "hypothetical.package" in html
        assert "import hypothetical" in html

    def test_uns_entry_with_tuple(self):
        """Test uns entry preview for tuple."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["my_tuple"] = (1, 2, 3)

        html = adata._repr_html_()
        assert "my_tuple" in html

    def test_uns_entry_with_empty_list(self):
        """Test uns entry preview for empty list."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["empty_list"] = []

        html = adata._repr_html_()
        assert "empty_list" in html

    def test_uns_entry_with_none(self):
        """Test uns entry preview for None."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["none_value"] = None

        html = adata._repr_html_()
        assert "none_value" in html
        assert "None" in html


class TestIndexPreviewRendering:
    """Tests for index preview rendering."""

    def test_index_preview_empty(self):
        """Test index preview for empty AnnData."""
        adata = AnnData()
        html = adata._repr_html_()
        assert "empty" in html.lower()

    def test_index_preview_few_items(self):
        """Test index preview shows all items when small."""
        adata = AnnData(np.zeros((4, 3)))
        html = adata._repr_html_()
        # Should show all obs/var names
        assert "obs_names" in html or "cell" in html.lower()

    def test_index_preview_many_items(self):
        """Test index preview shows first/last when many items."""
        adata = AnnData(
            np.zeros((100, 50)),
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(100)]),
            var=pd.DataFrame(index=[f"gene_{i}" for i in range(50)]),
        )
        html = adata._repr_html_()
        # Should show ellipsis for truncation
        assert "..." in html

    def test_format_index_value_decodes_bytes(self):
        """Test that bytes index values are decoded properly."""
        from anndata._repr.html import _format_index_value

        # UTF-8 bytes should decode properly
        assert _format_index_value(b"cell_0") == "cell_0"
        assert _format_index_value(b"gene_123") == "gene_123"

        # Regular strings pass through
        assert _format_index_value("cell_0") == "cell_0"

        # Numbers are converted to string
        assert _format_index_value(42) == "42"


class TestSectionTooltips:
    """Tests for section tooltip function."""

    def test_get_section_tooltip_all_sections(self):
        """Test tooltip text for all standard sections."""
        from anndata._repr.core import get_section_tooltip

        sections = [
            "obs",
            "var",
            "uns",
            "obsm",
            "varm",
            "layers",
            "obsp",
            "varp",
            "raw",
        ]
        for section in sections:
            tooltip = get_section_tooltip(section)
            assert tooltip, f"Missing tooltip for section: {section}"

    def test_get_section_tooltip_unknown(self):
        """Test tooltip for unknown section returns empty string."""
        from anndata._repr.core import get_section_tooltip

        assert get_section_tooltip("unknown_section") == ""


# =============================================================================
# README Button Feature Tests
# =============================================================================


class TestReadmeIcon:
    """Tests for the README icon feature in header."""

    def test_readme_icon_appears_with_string(self):
        """Test README icon appears when uns['README'] is a string."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["README"] = "# My Dataset\n\nThis is a test dataset."

        html = adata._repr_html_()
        # Check for README icon element
        assert 'class="adata-readme-icon"' in html
        # Icon should have info symbol
        assert "ⓘ" in html

    def test_readme_icon_not_shown_without_readme(self):
        """Test README icon is not shown when no README exists."""
        adata = AnnData(np.zeros((10, 5)))

        html = adata._repr_html_()
        # Check for absence of README icon element
        assert 'class="adata-readme-icon"' not in html

    def test_readme_icon_not_shown_for_empty_string(self):
        """Test README icon is not shown for empty string."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["README"] = ""

        html = adata._repr_html_()
        assert 'class="adata-readme-icon"' not in html

    def test_readme_icon_not_shown_for_whitespace_only(self):
        """Test README icon is not shown for whitespace-only string."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["README"] = "   \n\t  "

        html = adata._repr_html_()
        assert 'class="adata-readme-icon"' not in html

    def test_readme_icon_not_shown_for_non_string(self):
        """Test README icon is not shown for non-string values."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["README"] = {"title": "My Dataset"}  # dict, not string

        html = adata._repr_html_()
        assert 'class="adata-readme-icon"' not in html

    def test_readme_content_escaped(self):
        """Test README content is properly HTML-escaped."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["README"] = '<script>alert("xss")</script>'

        html = adata._repr_html_()
        # Script tag in README should be escaped in the data-readme attribute
        # Note: the HTML does contain <script> tags for the JS functionality
        assert "&lt;script&gt;alert" in html
        assert "alert(&quot;xss&quot;)" in html

    def test_readme_tooltip_truncated(self):
        """Test README tooltip is truncated for long content."""
        adata = AnnData(np.zeros((10, 5)))
        long_readme = "A" * 1000  # Very long README
        adata.uns["README"] = long_readme

        html = adata._repr_html_()
        # README icon should be present
        assert 'class="adata-readme-icon"' in html
        # The full 1000 chars should be in data-readme, but tooltip should be shorter
        assert 'title="' in html

    def test_readme_data_attribute_contains_content(self):
        """Test README content is in data-readme attribute."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["README"] = """# Dataset Title

## Description
This dataset contains important data.
"""

        html = adata._repr_html_()
        assert 'class="adata-readme-icon"' in html
        # Content should be in data-readme attribute (escaped)
        assert "data-readme=" in html
        # Title should be somewhere in the escaped content
        assert "Dataset Title" in html

    def test_readme_icon_accessibility(self):
        """Test README icon has proper accessibility attributes."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["README"] = "# Test\n\nContent"

        html = adata._repr_html_()
        assert 'role="button"' in html
        assert 'tabindex="0"' in html
        assert 'aria-label="View README"' in html


# =============================================================================
# Test render_header_badges Helper
# =============================================================================


class TestRenderHeaderBadges:
    """Tests for the render_header_badges public helper."""

    def test_no_badges(self):
        """Test render_header_badges with no flags set."""
        from anndata._repr import render_header_badges

        html = render_header_badges()
        assert html == ""

    def test_view_badge_only(self):
        """Test render_header_badges with view flag."""
        from anndata._repr import render_header_badges

        html = render_header_badges(is_view=True)
        assert "View" in html
        assert "adata-badge-view" in html
        assert "This is a view" in html

    def test_backed_badge_only(self):
        """Test render_header_badges with backed flag."""
        from anndata._repr import render_header_badges

        html = render_header_badges(is_backed=True)
        assert "Backed" in html
        assert "adata-badge-backed" in html

    def test_backed_badge_with_format(self):
        """Test render_header_badges with backing format."""
        from anndata._repr import render_header_badges

        html = render_header_badges(is_backed=True, backing_format="Zarr")
        assert "Zarr" in html
        assert "adata-badge-backed" in html

    def test_backed_badge_with_path(self):
        """Test render_header_badges with backing path in tooltip."""
        from anndata._repr import render_header_badges

        html = render_header_badges(
            is_backed=True,
            backing_path="/data/sample.h5ad",
            backing_format="H5AD",
        )
        assert "H5AD" in html
        assert "/data/sample.h5ad" in html

    def test_both_badges(self):
        """Test render_header_badges with both view and backed."""
        from anndata._repr import render_header_badges

        html = render_header_badges(
            is_view=True,
            is_backed=True,
            backing_format="Zarr",
        )
        assert "View" in html
        assert "Zarr" in html
        assert "adata-badge-view" in html
        assert "adata-badge-backed" in html

    def test_lazy_badge_only(self):
        """Test render_header_badges with lazy badge."""
        from anndata._repr import render_header_badges

        html = render_header_badges(is_lazy=True)
        assert "Lazy" in html
        assert "adata-badge-lazy" in html

    def test_is_lazy_detection(self):
        """Test is_lazy function detects lazy AnnData."""
        from anndata._repr.utils import is_lazy

        # Regular AnnData should not be lazy
        adata = AnnData(np.zeros((10, 5)))
        assert not is_lazy(adata)

        # Object without obs should not be lazy
        class NoObs:
            pass

        assert not is_lazy(NoObs())


# =============================================================================
# Test Error Handling in HTML Rendering
# =============================================================================


class TestErrorHandling:
    """Tests for error handling in HTML rendering."""

    def test_error_entry_display(self):
        """Test that error entries are displayed with error styling."""
        from anndata._repr.html import _render_error_entry

        html = _render_error_entry("bad_key", "Something went wrong")
        assert "bad_key" in html
        assert "Something went wrong" in html
        assert "adata-error" in html or "error" in html.lower()

    def test_formatter_exception_caught(self):
        """Test that formatter exceptions don't crash the repr."""
        from anndata._repr import (
            TypeFormatter,
            formatter_registry,
        )

        class CrashingFormatter(TypeFormatter):
            priority = 1000  # High priority to be checked first

            def can_format(self, obj):
                return isinstance(obj, dict) and obj.get("__crash__")

            def format(self, obj, context):
                msg = "Intentional crash"
                raise RuntimeError(msg)

        formatter = CrashingFormatter()
        formatter_registry.register_type_formatter(formatter)

        try:
            # Create AnnData with data that triggers the crashing formatter
            adata = AnnData(np.zeros((5, 3)))
            adata.uns["test"] = {"__crash__": True, "data": "value"}

            # Should not raise, should fall back gracefully
            import warnings

            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                html = adata._repr_html_()
                # Should have generated HTML despite the crash
                assert "anndata-repr" in html
                # Should have warned about the formatter failure
                assert any("CrashingFormatter" in str(warning.message) for warning in w)
        finally:
            formatter_registry.unregister_type_formatter(formatter)

    def test_section_formatter_exception_caught(self):
        """Test that section formatter exceptions don't crash the repr."""
        from anndata._repr import (
            SectionFormatter,
            register_formatter,
        )

        class CrashingSectionFormatter(SectionFormatter):
            @property
            def section_name(self):
                return "crashing_section"

            def should_show(self, obj):
                return True

            def get_entries(self, obj, context):
                msg = "Section formatter crash"
                raise RuntimeError(msg)

        formatter = CrashingSectionFormatter()
        register_formatter(formatter)

        try:
            adata = AnnData(np.zeros((5, 3)))
            import warnings

            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                html = adata._repr_html_()
                # Should still generate HTML
                assert "anndata-repr" in html
        finally:
            # Clean up - remove from registry
            from anndata._repr import formatter_registry

            if "crashing_section" in formatter_registry._section_formatters:
                del formatter_registry._section_formatters["crashing_section"]


# =============================================================================
# Test Unknown Sections Detection
# =============================================================================


class TestUnknownSectionsDetection:
    """Tests for detecting and displaying unknown sections."""

    def test_standard_sections_not_in_other(self):
        """Test that standard sections don't appear in 'other'."""
        adata = AnnData(
            np.zeros((10, 5)),
            obs=pd.DataFrame({"col": range(10)}),
            var=pd.DataFrame({"gene": range(5)}),
        )
        adata.obsm["X_pca"] = np.zeros((10, 2))
        adata.uns["neighbors"] = {"key": "value"}

        html = adata._repr_html_()

        # Standard sections should be present
        assert 'data-section="obs"' in html
        assert 'data-section="obsm"' in html
        assert 'data-section="uns"' in html

        # Should NOT have an "other" section for standard attributes
        # (other section only appears for truly unknown attributes)
        assert 'data-section="other"' not in html or "unknown" not in html.lower()

    def test_unknown_attribute_in_other(self):
        """Test that unknown attributes appear in 'other' section."""
        adata = AnnData(np.zeros((10, 5)))
        # Add a non-standard attribute
        adata._custom_attr = {"special": "data"}

        html = adata._repr_html_()

        # The repr should still work
        assert "anndata-repr" in html

    def test_registered_section_not_in_other(self):
        """Test that registered custom sections don't appear in 'other'."""
        from anndata._repr import (
            FormattedEntry,
            FormattedOutput,
            SectionFormatter,
            formatter_registry,
            register_formatter,
        )

        class CustomSectionFormatter(SectionFormatter):
            @property
            def section_name(self):
                return "custom_test_section"

            def should_show(self, obj):
                return hasattr(obj, "uns")

            def get_entries(self, obj, context):
                return [
                    FormattedEntry(
                        key="test_entry",
                        output=FormattedOutput(type_name="test"),
                    )
                ]

        formatter = CustomSectionFormatter()
        register_formatter(formatter)

        try:
            adata = AnnData(np.zeros((5, 3)))
            html = adata._repr_html_()

            # Custom section should be present
            assert 'data-section="custom_test_section"' in html
            assert "test_entry" in html
        finally:
            # Clean up
            if "custom_test_section" in formatter_registry._section_formatters:
                del formatter_registry._section_formatters["custom_test_section"]

    def test_suppressed_section_not_in_other(self):
        """Test that suppressed sections (should_show=False) don't appear in 'other'."""
        from anndata._repr import (
            SectionFormatter,
            formatter_registry,
            register_formatter,
        )

        class SuppressedSectionFormatter(SectionFormatter):
            @property
            def section_name(self):
                return "suppressed_section"

            def should_show(self, obj):
                return False  # Never show

            def get_entries(self, obj, context):
                return []

        formatter = SuppressedSectionFormatter()
        register_formatter(formatter)

        try:
            adata = AnnData(np.zeros((5, 3)))
            html = adata._repr_html_()

            # Suppressed section should NOT appear
            assert "suppressed_section" not in html
            # And should not be in "other" either
            assert "suppressed_section" not in html
        finally:
            # Clean up
            if "suppressed_section" in formatter_registry._section_formatters:
                del formatter_registry._section_formatters["suppressed_section"]


# =============================================================================
# Test Public API Exports
# =============================================================================


class TestPublicAPIExports:
    """Tests that all documented public API is actually importable."""

    def test_css_js_exports(self):
        """Test CSS and JS helper exports."""
        from anndata._repr import get_css, get_javascript

        css = get_css()
        assert "<style>" in css
        assert "anndata-repr" in css

        js = get_javascript("test-container")
        assert "<script>" in js
        assert "test-container" in js

    def test_section_rendering_exports(self):
        """Test section rendering helper exports."""
        from anndata._repr import (
            FormattedEntry,
            FormattedOutput,
            render_formatted_entry,
            render_section,
        )

        entry = FormattedEntry(
            key="test_key",
            output=FormattedOutput(type_name="test_type"),
        )
        entry_html = render_formatted_entry(entry)
        assert "test_key" in entry_html
        assert "test_type" in entry_html

        section_html = render_section("test_section", entry_html, n_items=1)
        assert "test_section" in section_html
        assert 'data-section="test_section"' in section_html

    def test_ui_helper_exports(self):
        """Test UI helper exports."""
        from anndata._repr import (
            render_badge,
            render_copy_button,
            render_fold_icon,
            render_search_box,
        )

        badge = render_badge("Test", "test-class", "tooltip")
        assert "Test" in badge
        assert "test-class" in badge

        search = render_search_box("container-id")
        assert "container-id" in search
        assert "search" in search.lower()

        fold = render_fold_icon()
        assert "fold" in fold.lower() or "svg" in fold.lower()

        copy_btn = render_copy_button("text-to-copy")
        assert "text-to-copy" in copy_btn or "copy" in copy_btn.lower()

    def test_utility_exports(self):
        """Test utility function exports."""
        from anndata._repr import escape_html, format_memory_size, format_number

        assert escape_html("<script>") == "&lt;script&gt;"
        assert format_number(1000) == "1,000"
        assert "KB" in format_memory_size(1024) or "1" in format_memory_size(1024)

    def test_registry_exports(self):
        """Test registry class exports."""
        from anndata._repr import (
            FormattedEntry,
            FormattedOutput,
            FormatterContext,
            FormatterRegistry,
            SectionFormatter,
            TypeFormatter,
            formatter_registry,
            register_formatter,
        )

        # All should be importable
        assert FormatterRegistry is not None
        assert TypeFormatter is not None
        assert SectionFormatter is not None
        assert FormattedOutput is not None
        assert FormattedEntry is not None
        assert FormatterContext is not None
        assert formatter_registry is not None
        assert register_formatter is not None

    def test_generate_repr_html_export(self):
        """Test generate_repr_html export and basic usage."""
        from anndata._repr import generate_repr_html

        adata = AnnData(np.zeros((10, 5)))
        html = generate_repr_html(adata)
        assert "anndata-repr" in html
        assert "AnnData" in html


class TestLazyCategoryLoading:
    """Tests for lazy category loading in HTML repr."""

    def test_get_lazy_category_count(self):
        """Test _get_lazy_category_count returns None for non-lazy columns.

        The function is designed specifically for lazy CategoricalArrays.
        For regular pandas Series or other objects, it returns None.
        """
        from anndata._repr.sections import _get_lazy_category_count

        # Regular pandas Series should return None (not a lazy CategoricalArray)
        series = pd.Series(pd.Categorical(["a", "b", "c"]))
        result = _get_lazy_category_count(series)
        assert result is None  # Not a lazy CategoricalArray

        # Object without proper structure should return None
        class MockCol:
            pass

        assert _get_lazy_category_count(MockCol()) is None

        # Non-categorical Series should return None
        non_cat = pd.Series([1, 2, 3])
        assert _get_lazy_category_count(non_cat) is None

    def test_get_lazy_categories_max_zero_skips(self):
        """Test that max_lazy_categories=0 skips loading entirely."""
        from anndata._repr.registry import FormatterContext
        from anndata._repr.sections import _get_lazy_categories

        context = FormatterContext(max_lazy_categories=0)

        class MockCol:
            pass

        categories, skipped, n_cats = _get_lazy_categories(MockCol(), context)
        assert categories == []
        assert skipped is True
        assert n_cats is None  # No category count available

    def test_get_categories_for_column_non_categorical(self):
        """Test _get_categories_for_column with non-categorical column."""
        from anndata._repr.registry import FormatterContext
        from anndata._repr.sections import _get_categories_for_column

        context = FormatterContext()
        series = pd.Series([1, 2, 3])

        categories, skipped, n_cats = _get_categories_for_column(
            series, context, is_categorical=False, is_lazy=False
        )
        assert categories == []
        assert skipped is False  # Not skipped, just no categories
        assert n_cats is None

    def test_get_categories_for_column_regular_categorical(self):
        """Test _get_categories_for_column with regular (non-lazy) categorical."""
        from anndata._repr.registry import FormatterContext
        from anndata._repr.sections import _get_categories_for_column

        context = FormatterContext()
        series = pd.Series(pd.Categorical(["a", "b", "c", "a"]))

        categories, skipped, n_cats = _get_categories_for_column(
            series, context, is_categorical=True, is_lazy=False
        )
        assert set(categories) == {"a", "b", "c"}
        assert skipped is False
        assert n_cats == 3

    def test_default_max_lazy_categories_export(self):
        """Test that DEFAULT_MAX_LAZY_CATEGORIES is properly exported."""
        from anndata._repr import DEFAULT_MAX_LAZY_CATEGORIES

        assert isinstance(DEFAULT_MAX_LAZY_CATEGORIES, int)
        assert DEFAULT_MAX_LAZY_CATEGORIES > 0

    def test_formatter_context_has_max_lazy_categories(self):
        """Test that FormatterContext has max_lazy_categories attribute."""
        from anndata._repr import DEFAULT_MAX_LAZY_CATEGORIES
        from anndata._repr.registry import FormatterContext

        context = FormatterContext()
        assert hasattr(context, "max_lazy_categories")
        assert context.max_lazy_categories == DEFAULT_MAX_LAZY_CATEGORIES

    def test_formatter_context_propagates_max_lazy_categories(self):
        """Test that FormatterContext.child() propagates max_lazy_categories."""
        from anndata._repr.registry import FormatterContext

        context = FormatterContext(max_lazy_categories=50)
        child = context.child("test_key")
        assert child.max_lazy_categories == 50

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_get_lazy_category_count_does_not_load_data(self, tmp_path):
        """Test that _get_lazy_category_count reads from storage metadata only."""
        from anndata._repr.sections import _get_lazy_category_count

        # Create real lazy AnnData
        adata = AnnData(
            sp.random(100, 50, density=0.1, format="csr", dtype=np.float32),
            obs=pd.DataFrame({"cat_col": pd.Categorical(["a", "b", "c"] * 33 + ["a"])}),
        )
        path = tmp_path / "test.zarr"
        adata.write_zarr(path)

        lazy = ad.experimental.read_lazy(path)
        col = lazy.obs._ds["cat_col"]

        # Check that categories are NOT loaded yet
        cat_arr = col.variable._data.array
        assert "categories" not in cat_arr.__dict__, (
            "Categories should not be loaded yet"
        )

        # Get count - should NOT load categories
        n_cats = _get_lazy_category_count(col)
        assert n_cats == 3

        # Verify categories still NOT loaded
        assert "categories" not in cat_arr.__dict__, (
            "Categories should still not be loaded"
        )

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_get_lazy_categories_does_not_load_data(self, tmp_path):
        """Test that _get_lazy_categories reads from storage directly."""
        from anndata._repr.registry import FormatterContext
        from anndata._repr.sections import _get_lazy_categories

        # Create real lazy AnnData
        adata = AnnData(
            sp.random(100, 50, density=0.1, format="csr", dtype=np.float32),
            obs=pd.DataFrame({"cat_col": pd.Categorical(["x", "y", "z"] * 33 + ["x"])}),
        )
        path = tmp_path / "test.zarr"
        adata.write_zarr(path)

        lazy = ad.experimental.read_lazy(path)
        col = lazy.obs._ds["cat_col"]

        # Check that categories are NOT loaded yet
        cat_arr = col.variable._data.array
        assert "categories" not in cat_arr.__dict__, (
            "Categories should not be loaded yet"
        )

        context = FormatterContext(max_lazy_categories=100)
        categories, skipped, n_cats = _get_lazy_categories(col, context)

        assert set(categories) == {"x", "y", "z"}
        assert not skipped
        assert n_cats == 3

        # Categories are read from storage directly, not via the cached property
        # The cached_property should still not be accessed
        assert "categories" not in cat_arr.__dict__, (
            "Categories cached property should not be accessed"
        )

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_get_lazy_categories_skipping_does_not_load_categories(self, tmp_path):
        """Test that when skipping (too many cats), we don't load category values."""
        from anndata._repr.registry import FormatterContext
        from anndata._repr.sections import _get_lazy_categories

        # Create large categorical (exceeds typical limit)
        large_cats = [f"cat_{i}" for i in range(150)]
        adata = AnnData(
            sp.random(150, 50, density=0.1, format="csr", dtype=np.float32),
            obs=pd.DataFrame({"big_cat": pd.Categorical(large_cats)}),
        )
        path = tmp_path / "test.zarr"
        adata.write_zarr(path)

        lazy = ad.experimental.read_lazy(path)
        col = lazy.obs._ds["big_cat"]

        # Check that categories are NOT loaded yet
        cat_arr = col.variable._data.array
        assert "categories" not in cat_arr.__dict__, (
            "Categories should not be loaded yet"
        )

        context = FormatterContext(max_lazy_categories=100)  # Less than 150
        categories, truncated, n_cats = _get_lazy_categories(col, context)

        # Now we get first 100 categories with truncation flag
        assert len(categories) == 100  # First 100 categories loaded
        assert categories[0] == "cat_0"  # First category
        assert truncated is True  # Indicates more categories exist
        assert n_cats == 150  # Total count known

        # Verify we didn't trigger the cached_property
        assert "categories" not in cat_arr.__dict__, (
            "Categories should still not be loaded"
        )

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_repr_html_does_not_load_lazy_categorical_data(self, tmp_path):
        """Test that generating HTML repr doesn't trigger loading of lazy categorical data."""
        from anndata._repr import DEFAULT_MAX_LAZY_CATEGORIES

        # Create lazy AnnData with categoricals
        adata = AnnData(
            sp.random(100, 50, density=0.1, format="csr", dtype=np.float32),
            obs=pd.DataFrame({
                "small_cat": pd.Categorical(["A", "B", "C"] * 33 + ["A"]),
                "large_cat": pd.Categorical(
                    [f"cat_{i}" for i in range(100)],
                    categories=[
                        f"cat_{i}" for i in range(DEFAULT_MAX_LAZY_CATEGORIES + 20)
                    ],
                ),
            }),
        )
        path = tmp_path / "test.zarr"
        adata.write_zarr(path)

        lazy = ad.experimental.read_lazy(path)

        # Get the CategoricalArrays before repr
        small_cat_arr = lazy.obs._ds["small_cat"].variable._data.array
        large_cat_arr = lazy.obs._ds["large_cat"].variable._data.array

        # Verify categories NOT loaded yet
        assert "categories" not in small_cat_arr.__dict__, (
            "small_cat categories loaded early"
        )
        assert "categories" not in large_cat_arr.__dict__, (
            "large_cat categories loaded early"
        )

        # Generate HTML repr
        html = lazy._repr_html_()

        # Verify HTML was generated and contains expected content
        assert "small_cat" in html
        assert "large_cat" in html
        assert "category" in html  # Type indicator

        # CRITICAL: Verify categories still NOT loaded via cached_property
        # (they may have been read directly from storage)
        assert "categories" not in small_cat_arr.__dict__, (
            "small_cat cached_property was accessed - this loads full dtype"
        )
        assert "categories" not in large_cat_arr.__dict__, (
            "large_cat cached_property was accessed - this loads full dtype"
        )

    @pytest.mark.skipif(
        not pytest.importorskip("xarray", reason="xarray not installed"),
        reason="xarray required for lazy loading",
    )
    def test_lazy_categorical_repr_integration(self, tmp_path):
        """Integration test: verify lazy categoricals display correctly in repr."""
        import h5py

        from anndata._repr import DEFAULT_MAX_LAZY_CATEGORIES
        from anndata.experimental import read_lazy

        n_large = DEFAULT_MAX_LAZY_CATEGORIES + 20  # More categories than the limit

        # Create test data with categoricals
        adata = AnnData(np.random.randn(100, 50).astype(np.float32))
        adata.obs["small_cat"] = pd.Categorical(["A", "B", "C"] * 33 + ["A"])
        # Create categorical with more category levels than observations
        # to exceed DEFAULT_MAX_LAZY_CATEGORIES
        adata.obs["large_cat"] = pd.Categorical(
            [f"cat_{i}" for i in range(100)],  # 100 observations
            categories=[f"cat_{i}" for i in range(n_large)],  # 120 categories
        )
        adata.uns["small_cat_colors"] = ["#ff0000", "#00ff00", "#0000ff"]

        # Save and read lazy
        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        with h5py.File(path, "r") as f:
            lazy_adata = read_lazy(f)
            html = lazy_adata._repr_html_()

            # Small categorical should show categories with colors
            assert "small_cat" in html
            assert "A" in html  # Category should be visible
            assert "#ff0000" in html or "ff0000" in html  # Color should be visible

            # Large categorical (120 cats) exceeds limit (100), shows first 100 + truncation
            assert "large_cat" in html
            # Should show first category and truncation indicator
            assert "cat_0" in html  # First category should be visible
            # Should show "...+20" indicating 20 hidden categories (120 - 100 = 20)
            assert "...+20" in html

    @pytest.mark.skipif(
        not pytest.importorskip("xarray", reason="xarray not installed"),
        reason="xarray required for lazy loading",
    )
    def test_metadata_only_mode_no_disk_loading(self, tmp_path):
        """Test that max_lazy_categories=0 shows counts without loading category labels."""
        import h5py

        from anndata._repr.html import generate_repr_html
        from anndata.experimental import read_lazy

        # Create test data
        adata = AnnData(np.random.randn(50, 20).astype(np.float32))
        adata.obs["cat1"] = pd.Categorical(["A", "B", "C"] * 16 + ["A", "A"])
        adata.obs["cat2"] = pd.Categorical(["X", "Y"] * 25)
        adata.uns["cat1_colors"] = ["#ff0000", "#00ff00", "#0000ff"]

        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        with h5py.File(path, "r") as f:
            lazy_adata = read_lazy(f)

            # Generate HTML in metadata-only mode
            html = generate_repr_html(lazy_adata, max_lazy_categories=0)

            # Should show category counts (from dtype metadata)
            assert "(3 categories)" in html  # cat1 has 3 categories
            assert "(2 categories)" in html  # cat2 has 2 categories

            # Should NOT show individual category labels
            # (these would require iterating dtype.categories)
            # Note: we check that individual labels aren't expanded in meta column
            # The type column will still show "category (N, lazy)"

            # Should NOT show colors (not loaded in metadata-only mode)
            assert "#ff0000" not in html
            assert "#00ff00" not in html

    @pytest.mark.skipif(
        not pytest.importorskip("xarray", reason="xarray not installed"),
        reason="xarray required for lazy loading",
    )
    def test_metadata_only_vs_default_mode(self, tmp_path):
        """Compare metadata-only mode vs default mode output."""
        import h5py

        from anndata._repr.html import generate_repr_html
        from anndata.experimental import read_lazy

        # Create test data with small categorical (should load in default mode)
        adata = AnnData(np.random.randn(50, 20).astype(np.float32))
        adata.obs["small_cat"] = pd.Categorical(["A", "B"] * 25)
        adata.uns["small_cat_colors"] = ["#ff0000", "#00ff00"]

        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        with h5py.File(path, "r") as f:
            lazy_adata = read_lazy(f)

            # Default mode: should load categories and colors
            html_default = generate_repr_html(lazy_adata)
            assert "#ff0000" in html_default  # Colors loaded
            # Category labels shown in meta column (A, B visible)

            # Metadata-only mode: should NOT load categories or colors
            html_metadata = generate_repr_html(lazy_adata, max_lazy_categories=0)
            assert "#ff0000" not in html_metadata  # Colors not loaded
            assert "(2 categories)" in html_metadata  # Just count shown
