"""
Tests for HTMLValidator and examples of proper HTML testing patterns.

These tests demonstrate the recommended approach for testing HTML repr output:
- Use structured assertions instead of string matching
- Validate element presence and attributes
- Check text appears in correct elements
- Verify sections and entries are properly rendered
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from anndata import AnnData

from .conftest import HTMLValidator


class TestHTMLValidatorBasics:
    """Tests for HTMLValidator functionality."""

    def test_assert_element_exists_by_class(self):
        """Test finding elements by class selector."""
        html = '<div class="my-class">content</div>'
        v = HTMLValidator(html)
        v.assert_element_exists(".my-class")

    def test_assert_element_exists_by_id(self):
        """Test finding elements by ID selector."""
        html = '<div id="my-id">content</div>'
        v = HTMLValidator(html)
        v.assert_element_exists("#my-id")

    def test_assert_element_exists_by_tag(self):
        """Test finding elements by tag selector."""
        html = "<span>content</span>"
        v = HTMLValidator(html)
        v.assert_element_exists("span")

    def test_assert_element_exists_by_attribute(self):
        """Test finding elements by attribute selector."""
        html = '<div data-section="obs">content</div>'
        v = HTMLValidator(html)
        v.assert_element_exists("[data-section=obs]")

    def test_assert_element_not_exists(self):
        """Test asserting element doesn't exist."""
        html = '<div class="other">content</div>'
        v = HTMLValidator(html)
        v.assert_element_not_exists(".missing-class")

    def test_assert_element_exists_fails(self):
        """Test assertion fails when element missing."""
        html = "<div>content</div>"
        v = HTMLValidator(html)
        with pytest.raises(AssertionError, match="not found"):
            v.assert_element_exists(".missing")

    def test_assert_text_visible(self):
        """Test asserting text is visible (not in style/script)."""
        html = "<style>.foo{}</style><div>visible text</div><script>hidden</script>"
        v = HTMLValidator(html)
        v.assert_text_visible("visible text")
        with pytest.raises(AssertionError):
            v.assert_text_visible("hidden")

    def test_assert_text_in_element(self):
        """Test asserting text appears in specific element."""
        html = '<div class="target">expected text</div><div class="other">other</div>'
        v = HTMLValidator(html)
        v.assert_text_in_element(".target", "expected text")

    def test_assert_section_exists(self):
        """Test asserting data section exists."""
        html = '<div data-section="obs">content</div>'
        v = HTMLValidator(html)
        v.assert_section_exists("obs")

    def test_assert_section_contains_entry(self):
        """Test asserting section contains entry."""
        html = '<div data-section="obs">batch cell_type</div>'
        v = HTMLValidator(html)
        v.assert_section_contains_entry("obs", "batch")

    def test_assert_badge_shown(self):
        """Test asserting badge is displayed."""
        html = '<span class="adata-badge-view">View</span>'
        v = HTMLValidator(html)
        v.assert_badge_shown("view")

    def test_assert_badge_not_shown(self):
        """Test asserting badge is not displayed."""
        html = '<span class="other">content</span>'
        v = HTMLValidator(html)
        v.assert_badge_not_shown("view")

    def test_assert_shape_displayed(self):
        """Test asserting shape values are displayed."""
        html = "<div>100 obs × 50 var</div>"
        v = HTMLValidator(html)
        v.assert_shape_displayed(100, 50)

    def test_assert_dtype_displayed(self):
        """Test asserting dtype is displayed."""
        html = '<span class="dtype">float32</span>'
        v = HTMLValidator(html)
        v.assert_dtype_displayed("float32")

    def test_count_elements(self):
        """Test counting matching elements."""
        html = '<div class="item">1</div><div class="item">2</div><div class="other">3</div>'
        v = HTMLValidator(html)
        assert v.count_elements(".item") == 2

    def test_chaining(self):
        """Test method chaining works."""
        html = '<div class="a" data-section="obs">text</div>'
        v = HTMLValidator(html)
        (
            v
            .assert_element_exists(".a")
            .assert_section_exists("obs")
            .assert_text_visible("text")
        )


class TestHTMLValidatorWithAnnData:
    """Tests demonstrating HTMLValidator with actual AnnData repr."""

    def test_basic_structure(self, validate_html):
        """Test basic AnnData repr structure."""
        adata = AnnData(np.zeros((100, 50)))
        html = adata._repr_html_()
        v = validate_html(html)

        # Validate structure, not just string presence
        v.assert_element_exists(".anndata-repr")
        v.assert_shape_displayed(100, 50)

    def test_sections_present(self, validate_html):
        """Test all populated sections are present."""
        adata = AnnData(
            np.zeros((10, 5)),
            obs=pd.DataFrame({"batch": ["A", "B"] * 5}),
            var=pd.DataFrame({"gene": range(5)}),
        )
        adata.uns["key"] = "value"
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_section_exists("obs")
        v.assert_section_exists("var")
        v.assert_section_exists("uns")

    def test_obs_entries_in_correct_section(self, validate_html):
        """Test obs column names appear in obs section."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["cell_type"] = pd.Categorical(["A", "B"] * 5)
        adata.obs["n_counts"] = range(10)
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_section_contains_entry("obs", "cell_type")
        v.assert_section_contains_entry("obs", "n_counts")

    def test_view_badge_displayed_correctly(self, validate_html):
        """Test View badge is shown for views."""
        adata = AnnData(np.zeros((10, 5)))
        view = adata[0:5, :]
        html = view._repr_html_()
        v = validate_html(html)

        v.assert_badge_shown("view")
        v.assert_shape_displayed(5, 5)  # View shape, not original

    def test_view_badge_not_shown_for_non_view(self, validate_html):
        """Test View badge is NOT shown for non-views."""
        adata = AnnData(np.zeros((10, 5)))
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_badge_not_shown("view")

    def test_sparse_matrix_dtype_displayed(self, validate_html):
        """Test sparse matrix shows dtype."""
        X = sp.random(100, 50, density=0.1, format="csr", dtype=np.float32)
        adata = AnnData(X)
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_dtype_displayed("float32")
        v.assert_text_visible("csr")

    def test_categorical_with_colors(self, validate_html):
        """Test categorical with colors shows color values."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["cluster"] = pd.Categorical(["A", "B"] * 5)
        adata.uns["cluster_colors"] = ["#FF0000", "#00FF00"]
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_section_contains_entry("obs", "cluster")
        v.assert_color_swatch("#FF0000")
        v.assert_color_swatch("#00FF00")

    def test_warning_for_unserializable(self, validate_html):
        """Test warning indicator for unserializable objects."""

        class CustomClass:
            pass

        adata = AnnData(np.zeros((5, 3)))
        adata.uns["custom"] = CustomClass()
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_section_contains_entry("uns", "custom")
        v.assert_warning_indicator()

    def test_backed_badge(self, validate_html, tmp_path):
        """Test backed badge is shown for backed AnnData."""
        import anndata as ad

        adata = AnnData(np.zeros((10, 5)))
        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        backed = ad.read_h5ad(path, backed="r")
        html = backed._repr_html_()
        v = validate_html(html)

        v.assert_badge_shown("backed")
        backed.file.close()

    def test_raw_section_visible(self, validate_html):
        """Test raw section is visible when raw is set."""
        adata = AnnData(np.zeros((10, 20)))
        adata.raw = adata.copy()
        adata = adata[:, :5]
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_section_exists("raw")
        # Raw should show original var count
        v.assert_text_visible("20")

    def test_layers_section(self, validate_html):
        """Test layers section shows layer names."""
        adata = AnnData(np.zeros((10, 5)))
        adata.layers["counts"] = np.ones((10, 5))
        adata.layers["normalized"] = np.zeros((10, 5))
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_section_exists("layers")
        v.assert_section_contains_entry("layers", "counts")
        v.assert_section_contains_entry("layers", "normalized")


class TestMigrationExamples:
    """Examples showing how to migrate old-style tests to HTMLValidator.

    OLD STYLE (string matching):
        html = adata._repr_html_()
        assert "batch" in html
        assert "View" in html

    NEW STYLE (structured validation):
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "batch")
        v.assert_badge_shown("view")
    """

    def test_old_vs_new_section_check(self, validate_html):
        """Compare old vs new section checking."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["batch"] = ["A", "B"] * 5
        html = adata._repr_html_()

        # OLD: String in HTML (could match CSS class name, comment, etc.)
        assert "batch" in html  # Weak assertion

        # NEW: Entry in correct section
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "batch")  # Strong assertion

    def test_old_vs_new_badge_check(self, validate_html):
        """Compare old vs new badge checking."""
        adata = AnnData(np.zeros((10, 5)))
        view = adata[0:5, :]
        html = view._repr_html_()

        # OLD: String matching (could match "View" in documentation, comments, etc.)
        assert "View" in html  # Weak assertion

        # NEW: Check for actual badge element
        v = validate_html(html)
        v.assert_badge_shown("view")  # Strong assertion - checks CSS class

    def test_old_vs_new_shape_check(self, validate_html):
        """Compare old vs new shape checking."""
        adata = AnnData(np.zeros((123, 456)))
        html = adata._repr_html_()

        # OLD: Numbers anywhere in HTML
        assert "123" in html  # Could match anything
        assert "456" in html

        # NEW: Explicit shape check
        v = validate_html(html)
        v.assert_shape_displayed(123, 456)
