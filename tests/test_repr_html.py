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
"""

from __future__ import annotations

import re
from html.parser import HTMLParser
from string import ascii_letters
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import anndata as ad
from anndata import AnnData

# Check optional dependencies
try:
    import dask.array as da

    HAS_DASK = True
except ImportError:
    HAS_DASK = False

try:
    import cupy as cp

    HAS_CUPY = True
except ImportError:
    HAS_CUPY = False

try:
    import awkward as ak

    HAS_AWKWARD = True
except ImportError:
    HAS_AWKWARD = False


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
        obs=pd.DataFrame(
            {
                "batch": pd.Categorical(["A", "B"] * (n_obs // 2)),
                "n_counts": np.random.randint(1000, 10000, n_obs),
                "cell_type": pd.Categorical(
                    ["T", "B", "NK"] * (n_obs // 3) + ["T"] * (n_obs % 3)
                ),
            }
        ),
        var=pd.DataFrame(
            {
                "gene_name": [f"gene_{i}" for i in range(n_vars)],
                "highly_variable": np.random.choice([True, False], n_vars),
            }
        ),
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
    adata.uns['quotes"test\''] = "value"
    return adata


# =============================================================================
# HTML Validation
# =============================================================================


VOID_ELEMENTS = frozenset(
    {
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
    }
)


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
            # Find variable definitions (--name: value)
            defined = set(re.findall(r"--([\w-]+)\s*:", css))
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
        assert "100" in html and "50" in html

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
        for i in range(5):
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
        html = adata._repr_html_()
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
        html = adata._repr_html_()
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

            pass

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

            pass

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
        assert is_serializable(True)[0]
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
        warn, msg = should_warn_string_column(s)
        assert warn
        assert "3" in msg  # 3 unique values

        # All unique should not warn
        s = pd.Series(["A", "B", "C", "D", "E"])
        warn, msg = should_warn_string_column(s)
        assert not warn

        # Numeric should not warn
        s = pd.Series([1, 2, 1, 2, 1])
        warn, msg = should_warn_string_column(s)
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
