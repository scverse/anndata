"""
UI feature tests for the _repr module.

Tests for folding behavior, color display, search functionality,
clipboard, ad-blocker compatibility, and other UI features.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
import pytest

from anndata import AnnData


class TestFoldingBehavior:
    """Test auto-folding based on entry count."""

    def test_few_entries_expanded(self, validate_html):
        """Test sections with few entries are NOT collapsed."""
        adata = AnnData(
            np.zeros((10, 5)),
            obs=pd.DataFrame({"a": list(range(10)), "b": list(range(10))}),
        )
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("obs")
        v.assert_section_contains_entry("obs", "a")
        v.assert_section_contains_entry("obs", "b")
        # With 2 entries (below threshold of 5), section should NOT be collapsed
        v.assert_section_not_initially_collapsed("obs")

    def test_at_threshold_not_collapsed(self, validate_html):
        """Test section at exactly the threshold is NOT collapsed."""
        # Default threshold is 5, so 5 entries should NOT be collapsed
        obs_dict = {f"col_{i}": list(range(10)) for i in range(5)}
        adata = AnnData(np.zeros((10, 5)), obs=pd.DataFrame(obs_dict))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("obs")
        v.assert_section_not_initially_collapsed("obs")

    def test_above_threshold_collapsed(self, validate_html):
        """Test section above threshold IS collapsed."""
        # Default threshold is 5, so 6 entries SHOULD be collapsed
        obs_dict = {f"col_{i}": list(range(10)) for i in range(6)}
        adata = AnnData(np.zeros((10, 5)), obs=pd.DataFrame(obs_dict))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("obs")
        v.assert_collapse_functionality_present()
        v.assert_section_initially_collapsed("obs")

    def test_many_entries_collapsed(self, validate_html):
        """Test sections with many entries have collapse functionality."""
        obs_dict = {f"col_{i}": list(range(10)) for i in range(10)}
        adata = AnnData(np.zeros((10, 5)), obs=pd.DataFrame(obs_dict))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("obs")
        v.assert_collapse_functionality_present()
        v.assert_section_initially_collapsed("obs")

    def test_custom_fold_threshold(self, validate_html):
        """Test custom fold threshold via settings."""
        from anndata import settings

        obs_dict = {f"col_{i}": list(range(10)) for i in range(4)}
        adata = AnnData(np.zeros((10, 5)), obs=pd.DataFrame(obs_dict))

        with settings.override(repr_html_fold_threshold=2):
            html = adata._repr_html_()
            v = validate_html(html)
            v.assert_section_exists("obs")
            v.assert_collapse_functionality_present()
            v.assert_section_initially_collapsed("obs")


class TestNestedAnnData:
    """Test recursive AnnData display in .uns."""

    def test_nested_anndata_present(self, validate_html):
        """Test nested AnnData in .uns is detected."""
        inner = AnnData(np.zeros((5, 3)))
        outer = AnnData(np.zeros((10, 5)))
        outer.uns["nested"] = inner
        html = outer._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("uns")
        v.assert_section_contains_entry("uns", "nested")

    def test_nested_anndata_max_depth(self, validate_html):
        """Test max recursion depth is respected."""
        from anndata import settings

        adata = AnnData(np.zeros((5, 3)))
        for _i in range(5):
            wrapper = AnnData(np.zeros((5, 3)))
            wrapper.uns["inner"] = adata
            adata = wrapper

        with settings.override(repr_html_max_depth=2):
            html = adata._repr_html_()
            v = validate_html(html)
            v.assert_element_exists(".anndata-repr")
            assert "max depth" in html.lower() or "depth" in html

    def test_nested_anndata_expandable(self, validate_html):
        """Test nested AnnData has expand button."""
        inner = AnnData(np.zeros((5, 3)))
        outer = AnnData(np.zeros((10, 5)))
        outer.uns["nested"] = inner
        html = outer._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("uns", "nested")
        assert "Expand" in html or "expand" in html.lower()


class TestColorDisplay:
    """Test color swatch display."""

    def test_color_list_swatches(self, validate_html):
        """Test *_colors entries show color swatches."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["cluster"] = pd.Categorical(["A", "B"] * 5)
        adata.uns["cluster_colors"] = ["#FF0000", "#00FF00"]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "cluster")
        v.assert_color_swatch("#FF0000")

    def test_color_mismatch_warning(self, validate_html):
        """Test warning when color count doesn't match categories."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["cluster"] = pd.Categorical(["A", "B", "C"] * 3 + ["A"])
        adata.uns["cluster_colors"] = ["#FF0000", "#00FF00"]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "cluster")
        # Should show mismatch warning
        assert "mismatch" in html.lower() or "warning" in html.lower() or "⚠" in html

    def test_colors_inline_with_category(self, validate_html):
        """Test categorical columns show inline colors."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["cluster"] = pd.Categorical(["A", "B", "C"] * 3 + ["A"])
        adata.uns["cluster_colors"] = ["#FF0000", "#00FF00", "#0000FF"]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "cluster")
        v.assert_color_swatch("#FF0000")

    def test_color_list_many_colors(self, validate_html):
        """Test color list with many colors."""
        adata = AnnData(np.zeros((20, 5)))
        cats = [f"cat_{i}" for i in range(20)]
        adata.obs["cluster"] = pd.Categorical(cats)
        adata.uns["cluster_colors"] = [f"#{i:02x}0000" for i in range(20)]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "cluster")

    def test_color_list_exact_limit(self, validate_html):
        """Test color list at exact limit."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["cluster"] = pd.Categorical([f"cat_{i}" for i in range(10)])
        adata.uns["cluster_colors"] = [f"#{i:02x}0000" for i in range(10)]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")


class TestSearchFunctionality:
    """Test search/filter feature."""

    def test_search_input_present(self, adata, validate_html):
        """Test search input field is present."""
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists("input")
        assert "search" in html.lower() or "filter" in html.lower()

    def test_filter_indicator_present(self, adata, validate_html):
        """Test filter indicator element exists."""
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        assert "filter-indicator" in html or "filterIndicator" in html

    def test_data_key_attributes(self, adata_full, validate_html):
        """Test entries have data-key for filtering."""
        html = adata_full._repr_html_()
        v = validate_html(html)
        v.assert_element_exists("[data-key]")

    def test_data_dtype_attributes(self, adata_full, validate_html):
        """Test entries have data-dtype for filtering."""
        html = adata_full._repr_html_()
        v = validate_html(html)
        v.assert_element_exists("[data-dtype]")


class TestCopyToClipboard:
    """Test clipboard copy buttons."""

    def test_copy_buttons_present(self, adata_full, validate_html):
        """Test copy buttons exist for entries."""
        html = adata_full._repr_html_()
        v = validate_html(html)
        # Validate copy button elements exist with proper class
        v.assert_element_exists(".adata-copy-btn")
        # Validate copy buttons have data-copy attribute for the key name
        v.assert_element_exists("[data-copy]")

    def test_copy_button_has_accessibility_attributes(self, adata_full, validate_html):
        """Test copy buttons have proper accessibility attributes."""
        html = adata_full._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".adata-copy-btn")
        # Should have title and aria-label for accessibility
        assert 'aria-label="Copy name"' in html or 'title="Copy' in html


class TestDocumentationLinks:
    """Test help/documentation links."""

    def test_help_links_present(self, adata_full, validate_html):
        """Test help links exist for sections."""
        html = adata_full._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        # Check for readthedocs links in href attributes
        assert "readthedocs" in html
        # Should have links to documentation for standard sections
        assert "anndata.AnnData.obs" in html or "anndata.readthedocs.io" in html


class TestAdBlockerCompatibility:
    """Test that HTML representation doesn't trigger ad-blockers."""

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

        class_pattern = re.compile(r'class="([^"]*)"')
        all_classes = class_pattern.findall(html)

        problematic_classes = []
        for class_attr in all_classes:
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

        assert "adata-type" in html or "anndata" in html.lower()
        assert "adata-entry" in html or "anndata-sec" in html
        assert "adata-table" in html or "table" in html.lower()
        assert len(html) > 1000


class TestSectionTooltips:
    """Tests for section tooltips."""

    def test_get_section_tooltip_all_sections(self):
        """Test tooltips exist for all standard sections."""
        from anndata._repr.core import get_section_tooltip

        sections = [
            "obs",
            "var",
            "uns",
            "obsm",
            "varm",
            "obsp",
            "varp",
            "layers",
            "raw",
        ]
        for section in sections:
            tooltip = get_section_tooltip(section)
            assert isinstance(tooltip, str)

    def test_get_section_tooltip_unknown(self):
        """Test tooltip for unknown section."""
        from anndata._repr.core import get_section_tooltip

        tooltip = get_section_tooltip("unknown_section")
        assert tooltip == ""


class TestColorSwatchesAndCategories:
    """Tests for color swatch display and category handling."""

    def test_matching_colors_show_swatches(self, validate_html):
        """Test matching *_colors show color swatches."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["cluster"] = pd.Categorical(["A", "B"] * 5)
        adata.uns["cluster_colors"] = ["#FF0000", "#00FF00"]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "cluster")
        v.assert_color_swatch("#FF0000")

    def test_color_mismatch_shows_warning(self, validate_html):
        """Test color count mismatch shows warning."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["cluster"] = pd.Categorical(["A", "B", "C"] * 3 + ["A"])
        # Only 2 colors for 3 categories
        adata.uns["cluster_colors"] = ["#FF0000", "#00FF00"]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "cluster")
        # Should indicate mismatch warning
        assert "mismatch" in html.lower() or "⚠" in html or "warning" in html.lower()

    def test_categories_shown_inline(self, validate_html):
        """Test category values are shown inline."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["cluster"] = pd.Categorical(
            ["TypeA", "TypeB", "TypeC"] * 3 + ["TypeA"]
        )
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "cluster")
        v.assert_text_visible("TypeA")
        v.assert_text_visible("TypeB")
        v.assert_text_visible("TypeC")

    def test_many_categories_truncated_with_count(self, validate_html):
        """Test many categories are truncated with remaining count."""
        from anndata import settings

        categories = [f"Type_{i}" for i in range(30)]
        adata = AnnData(
            np.zeros((30, 5)),
            obs=pd.DataFrame({"cat": pd.Categorical(categories)}),
        )

        with settings.override(repr_html_max_categories=10):
            html = adata._repr_html_()
            v = validate_html(html)
            v.assert_section_contains_entry("obs", "cat")
            v.assert_text_visible("Type_0")
            assert "...+" in html or "20" in html or "more" in html.lower()


class TestNestedStructuresAndDepth:
    """Tests for nested structures and depth handling."""

    def test_nested_anndata_shows_expand_button(self):
        """Test nested AnnData has expand functionality."""
        inner = AnnData(np.zeros((5, 3)))
        outer = AnnData(np.zeros((10, 5)))
        outer.uns["nested"] = inner
        html = outer._repr_html_()
        # Should have expand functionality
        assert "expand" in html.lower() or "Expand" in html

    def test_max_depth_shows_indicator(self):
        """Test max depth shows indicator."""
        from anndata import settings

        level2 = AnnData(np.zeros((5, 3)))
        level1 = AnnData(np.zeros((7, 4)))
        level1.uns["inner"] = level2
        level0 = AnnData(np.zeros((10, 5)))
        level0.uns["inner"] = level1

        with settings.override(repr_html_max_depth=1):
            html = level0._repr_html_()
            # Should indicate max depth reached
            assert "max depth" in html.lower() or "depth" in html

    def test_deeply_nested_dict_handled(self):
        """Test deeply nested dicts are handled."""
        adata = AnnData(np.zeros((5, 3)))
        nested = {"level": 0}
        current = nested
        for i in range(20):
            current["nested"] = {"level": i + 1}
            current = current["nested"]
        adata.uns["deep"] = nested
        html = adata._repr_html_()
        assert html is not None
        assert "deep" in html


# Fixtures needed by tests above


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
    import scipy.sparse as sp

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
