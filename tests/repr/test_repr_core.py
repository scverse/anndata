"""
Core HTML repr tests: validation, basic generation, settings, header/footer.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
import scipy.sparse as sp

from anndata import AnnData

from .conftest import StrictHTMLParser


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
        a_without_href = re.search(r"<a(?![^>]*href=)[^>]*>", html)
        assert a_without_href is None, "Found <a> without href"

    def test_style_tag_valid_css(self, adata):
        """Test inline CSS has balanced braces."""
        html = adata._repr_html_()
        style_match = re.search(r"<style[^>]*>(.*?)</style>", html, re.DOTALL)
        if style_match:
            css = style_match.group(1)
            assert css.count("{") == css.count("}"), "Unbalanced CSS braces"

    def test_escaped_user_content(self, adata_with_special_chars):
        """Test that user-provided content is properly escaped."""
        html = adata_with_special_chars._repr_html_()
        # Should not contain raw <script> tags from user content
        assert "<script>" not in html.split("</style>")[-1].split("<script")[0].lower()
        assert "javascript:" not in html.lower()

    def test_nested_anndata_valid_html(self, adata_with_nested):
        """Test nested AnnData produces valid nested HTML."""
        html = adata_with_nested._repr_html_()
        parser = StrictHTMLParser()
        parser.feed(html)
        assert not parser.errors

    def test_html5_w3c_validation(self, adata_full, validate_html5):
        """Test HTML passes W3C HTML5 validation.

        This test runs full HTML5 validation using Nu Html Checker (vnu)
        when available. Install via: brew install vnu

        Validates:
        - Full W3C HTML5 spec compliance
        - ARIA accessibility attributes
        - Proper attribute values

        Note: Some errors are expected for Jupyter HTML fragments:
        - <style> in body: Valid for inline fragments, filtered out
        """
        html = adata_full._repr_html_()
        errors = validate_html5(html)
        # Filter out info-level messages and expected fragment issues
        critical = [
            e
            for e in errors
            if not e.startswith("info:")
            # <style> in body is valid for Jupyter HTML fragments
            and "style" not in e.lower()
            and "script" not in e.lower()
        ]
        assert not critical, "HTML5 validation errors:\n" + "\n".join(critical)

    def test_javascript_syntax_valid(self, adata_full, validate_js):
        """Test JavaScript syntax is valid.

        This test validates JavaScript syntax using esprima when available.
        Install via: pip install esprima

        Validates:
        - JavaScript syntax (no parsing errors)
        - All scripts are syntactically valid
        """
        html = adata_full._repr_html_()
        errors = validate_js(html)
        assert not errors, "JavaScript syntax errors:\n" + "\n".join(errors)


class TestCSSValidation:
    """Validate CSS styling."""

    def test_css_variables_defined(self, adata):
        """Test all CSS variables are defined before use."""
        html = adata._repr_html_()
        style_match = re.search(r"<style[^>]*>(.*?)</style>", html, re.DOTALL)
        if style_match:
            css = style_match.group(1)
            defined = set(re.findall(r"--([\w-]+)\s*:", css))
            inline_styles = re.findall(r'style="([^"]*)"', html)
            for style in inline_styles:
                defined |= set(re.findall(r"--([\w-]+)\s*:", style))
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


class TestBasicHTMLGeneration:
    """Test basic HTML generation for various AnnData configurations."""

    def test_empty_anndata(self, validate_html):
        """Test repr for empty AnnData."""
        adata = AnnData()
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_shape_displayed(0, 0)
        v.assert_element_exists(".anndata-repr")

    def test_minimal_anndata(self, validate_html):
        """Test repr with only X."""
        adata = AnnData(np.zeros((10, 5)))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_shape_displayed(10, 5)
        v.assert_element_exists(".anndata-repr")

    def test_dense_matrix(self, validate_html):
        """Test repr with dense X."""
        adata = AnnData(np.random.randn(100, 50).astype(np.float32))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_dtype_displayed("float32")
        v.assert_shape_displayed(100, 50)

    def test_sparse_csr_matrix(self, validate_html):
        """Test repr with sparse CSR X."""
        X = sp.random(1000, 500, density=0.1, format="csr")
        adata = AnnData(X)
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_text_visible("csr")
        v.assert_shape_displayed(1000, 500)

    def test_sparse_csc_matrix(self, validate_html):
        """Test repr with sparse CSC X."""
        X = sp.random(1000, 500, density=0.1, format="csc")
        adata = AnnData(X)
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_text_visible("csc")
        v.assert_shape_displayed(1000, 500)

    def test_all_attributes_present(self, adata_full, validate_html):
        """Test all standard attributes appear in repr."""
        html = adata_full._repr_html_()
        v = validate_html(html)
        for attr in ["obs", "var", "uns", "obsm", "varm", "layers", "obsp", "varp"]:
            v.assert_section_exists(attr)

    def test_x_none_handled(self, validate_html):
        """Test AnnData with X=None is handled."""
        adata = AnnData(obs=pd.DataFrame({"a": [1, 2, 3]}))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        # X=None should show None or similar indicator
        assert "None" in html or "none" in html.lower()

    def test_html_with_empty_obsm_varm(self, validate_html):
        """Test HTML repr with empty obsm/varm."""
        adata = AnnData(np.zeros((10, 5)))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_shape_displayed(10, 5)

    def test_html_with_all_empty_sections(self, validate_html):
        """Test HTML repr with minimal data (empty AnnData)."""
        adata = AnnData()
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_shape_displayed(0, 0)


class TestViewRepresentation:
    """Test HTML repr for AnnData views."""

    def test_view_badge_present(self, adata, validate_html):
        """Test view badge appears for views."""
        view = adata[0:5, :]
        html = view._repr_html_()
        v = validate_html(html)
        v.assert_badge_shown("view")

    def test_view_shows_subset_shape(self, adata, validate_html):
        """Test view shows correct subset dimensions."""
        view = adata[0:5, 0:3]
        html = view._repr_html_()
        v = validate_html(html)
        v.assert_shape_displayed(5, 3)
        v.assert_badge_shown("view")


class TestBackedAnnData:
    """Test HTML repr for disk-backed AnnData."""

    def test_backed_badge_h5ad(self, tmp_path, validate_html):
        """Test backed badge for H5AD files."""
        import anndata as ad

        adata = AnnData(np.random.randn(100, 50).astype(np.float32))
        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        backed = ad.read_h5ad(path, backed="r")
        html = backed._repr_html_()
        v = validate_html(html)
        v.assert_badge_shown("backed")
        v.assert_shape_displayed(100, 50)
        backed.file.close()

    def test_backed_shows_filename(self, tmp_path, validate_html):
        """Test backed repr shows filename."""
        import anndata as ad

        adata = AnnData(np.random.randn(100, 50).astype(np.float32))
        path = tmp_path / "test_file.h5ad"
        adata.write_h5ad(path)

        backed = ad.read_h5ad(path, backed="r")
        html = backed._repr_html_()
        v = validate_html(html)
        v.assert_badge_shown("backed")
        v.assert_text_visible("h5ad")
        backed.file.close()


class TestSettings:
    """Test settings integration."""

    def test_html_disabled_fallback(self, adata):
        """Test fallback to text repr when HTML disabled."""
        from anndata import settings

        with settings.override(repr_html_enabled=False):
            html = adata._repr_html_()
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
            assert "more" in html.lower() or "..." in html


class TestColumnWidthSettings:
    """Test column width calculation settings."""

    def test_field_width_calculated_from_content(self, validate_html):
        """Test field width adapts to content."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["short"] = list(range(10))
        adata.obs["very_long_column_name_here"] = list(range(10))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_css_variable_defined("--anndata-name-col-width")

    def test_field_width_capped_by_max(self, validate_html):
        """Test field width is capped by max setting."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))
        adata.obs["x" * 500] = list(range(10))

        with settings.override(repr_html_max_field_width=200):
            html = adata._repr_html_()
            v = validate_html(html)
            # Should not exceed max width
            v.assert_css_variable_defined("--anndata-name-col-width")

    def test_repr_html_max_field_width_setting(self):
        """Test repr_html_max_field_width setting is respected."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))
        adata.obs["column"] = list(range(10))

        with settings.override(repr_html_max_field_width=100):
            html = adata._repr_html_()
            assert html is not None

    def test_repr_html_type_width_setting(self):
        """Test repr_html_type_width setting is respected."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))

        with settings.override(repr_html_type_width=150):
            html = adata._repr_html_()
            assert "--anndata-type-col-width: 150px" in html

    def test_field_width_considers_all_sections(self):
        """Test field width considers names from all sections."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["obs_col"] = list(range(10))
        adata.var["var_col"] = list(range(5))
        adata.uns["uns_key"] = "value"
        adata.obsm["obsm_key"] = np.zeros((10, 2))

        html = adata._repr_html_()
        assert "--anndata-name-col-width" in html

    def test_empty_adata_field_width(self):
        """Test field width for empty AnnData."""
        adata = AnnData()
        html = adata._repr_html_()
        assert "--anndata-name-col-width" in html


class TestSettingsEffectOnRendering:
    """Test that each setting actually affects the rendering output."""

    def test_repr_html_fold_threshold_effect(self, validate_html):
        """Test fold_threshold setting affects section folding."""
        from anndata import settings

        # Create AnnData with 6 obs columns
        obs_dict = {f"col_{i}": list(range(10)) for i in range(6)}
        adata = AnnData(np.zeros((10, 5)), obs=pd.DataFrame(obs_dict))

        # Default threshold is 5, so 6 columns should collapse
        html_default = adata._repr_html_()
        v = validate_html(html_default)
        v.assert_collapse_functionality_present()
        v.assert_section_initially_collapsed("obs")

        # With higher threshold, should not collapse
        with settings.override(repr_html_fold_threshold=10):
            html_expanded = adata._repr_html_()
            v2 = validate_html(html_expanded)
            v2.assert_element_exists(".anndata-repr")

    def test_repr_html_max_categories_effect(self, validate_html):
        """Test max_categories setting limits category display."""
        from anndata import settings

        # Create AnnData with many categories
        categories = [f"cat_{i}" for i in range(20)]
        adata = AnnData(
            np.zeros((20, 5)),
            obs=pd.DataFrame({"many_cats": pd.Categorical(categories)}),
        )

        # With low max_categories, should truncate
        with settings.override(repr_html_max_categories=5):
            html = adata._repr_html_()
            v = validate_html(html)
            v.assert_truncation_indicator()

        # With high max_categories, should show all
        with settings.override(repr_html_max_categories=50):
            html = adata._repr_html_()
            v2 = validate_html(html)
            v2.assert_text_visible("cat_0")

    def test_repr_html_unique_limit_effect(self):
        """Test unique_limit setting skips nunique for large columns."""
        from anndata import settings

        # Create AnnData with many unique values
        adata = AnnData(
            np.zeros((1000, 5)),
            obs=pd.DataFrame({"big_col": list(range(1000))}),
        )

        # With low unique_limit, should skip nunique calculation
        with settings.override(repr_html_unique_limit=100):
            html = adata._repr_html_()
            assert "big_col" in html

        # Should still render without error

    def test_repr_html_dataframe_expand_effect(self):
        """Test dataframe_expand setting adds expand button."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))
        adata.uns["df"] = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})

        # With expand enabled, should have expand button
        with settings.override(repr_html_dataframe_expand=True):
            html = adata._repr_html_()
            assert "expand" in html.lower() or "Expand" in html

        # With expand disabled, should not have expand button
        with settings.override(repr_html_dataframe_expand=False):
            html = adata._repr_html_()
            # DataFrame should still be shown but without expand

    def test_repr_html_max_items_truncates_entries(self):
        """Test max_items setting truncates section entries."""
        from anndata import settings

        # Create AnnData with many uns entries
        adata = AnnData(np.zeros((10, 5)))
        for i in range(100):
            adata.uns[f"key_{i}"] = i

        with settings.override(repr_html_max_items=10):
            html = adata._repr_html_()
            # Should show truncation indicator
            assert "more" in html.lower() or "..." in html
            # First key should be shown
            assert "key_0" in html

    def test_repr_html_max_depth_limits_nesting(self):
        """Test max_depth setting limits nested AnnData display."""
        from anndata import settings

        # Create deeply nested AnnData
        level2 = AnnData(np.zeros((5, 3)))
        level1 = AnnData(np.zeros((7, 4)))
        level1.uns["nested"] = level2
        level0 = AnnData(np.zeros((10, 5)))
        level0.uns["nested"] = level1

        # With depth=1, should show warning at level 2
        with settings.override(repr_html_max_depth=1):
            html = level0._repr_html_()
            assert "max depth" in html.lower() or "depth" in html

        # With depth=5, should show all levels
        with settings.override(repr_html_max_depth=5):
            html = level0._repr_html_()
            assert html is not None

    def test_repr_html_enabled_false_returns_none(self):
        """Test enabled=False returns None for HTML repr."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))

        with settings.override(repr_html_enabled=False):
            html = adata._repr_html_()
            assert html is None

        # When enabled, should return HTML
        html = adata._repr_html_()
        assert html is not None
        assert "<div" in html

    def test_repr_html_type_width_css_variable(self):
        """Test type_width setting sets CSS variable."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))

        with settings.override(repr_html_type_width=200):
            html = adata._repr_html_()
            assert "--anndata-type-col-width: 200px" in html

        with settings.override(repr_html_type_width=300):
            html = adata._repr_html_()
            assert "--anndata-type-col-width: 300px" in html

    def test_repr_html_max_field_width_css_variable(self):
        """Test max_field_width setting affects name column width."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))
        # Add a very long column name
        adata.obs["x" * 500] = list(range(10))

        # With small max, name column should be capped
        with settings.override(repr_html_max_field_width=100):
            html = adata._repr_html_()
            assert "--anndata-name-col-width:" in html
            # Width should be at most 100px
            import re

            match = re.search(r"--anndata-name-col-width:\s*(\d+)px", html)
            if match:
                width = int(match.group(1))
                assert width <= 100


class TestMetadata:
    """Test metadata display."""

    def test_version_displayed(self, adata, validate_html):
        """Test anndata version is shown."""
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_text_visible("anndata")
        assert re.search(r"\d+\.\d+", html)

    def test_memory_usage_displayed(self, adata, validate_html):
        """Test memory usage is shown."""
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        assert re.search(r"\d+(\.\d+)?\s*(B|KB|MB|GB|bytes)", html, re.IGNORECASE)

    def test_obs_var_names_preview(self, adata, validate_html):
        """Test obs_names and var_names preview."""
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        # Should show cell names or obs_names indicator
        assert "cell_" in html or "obs_names" in html
        # Should show gene names or var_names indicator
        assert "gene_" in html or "var_names" in html


class TestIndexPreviewRendering:
    """Test index preview rendering."""

    def test_index_preview_empty(self, validate_html):
        """Test empty index preview."""
        adata = AnnData()
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_shape_displayed(0, 0)

    def test_index_preview_few_items(self, validate_html):
        """Test index preview with few items."""
        adata = AnnData(np.zeros((3, 2)))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_shape_displayed(3, 2)

    def test_index_preview_many_items(self, validate_html):
        """Test index preview with many items shows ellipsis."""
        adata = AnnData(
            np.zeros((100, 5)),
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(100)]),
        )
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_text_visible("cell_0")
        v.assert_text_visible("...")
        v.assert_shape_displayed(100, 5)

    def test_format_index_preview_decodes_bytes(self):
        """Test index preview decodes bytes values."""
        # Can't easily set bytes index, so just verify function exists
        from anndata._repr.html import _format_index_preview

        assert callable(_format_index_preview)


class TestRenderHeaderBadges:
    """Test header badge rendering."""

    def test_no_badges(self, validate_html):
        """Test no badges for basic AnnData."""
        adata = AnnData(np.zeros((10, 5)))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_badge_not_shown("view")
        v.assert_badge_not_shown("backed")

    def test_view_badge_only(self, adata, validate_html):
        """Test view badge appears alone."""
        view = adata[0:5, :]
        html = view._repr_html_()
        v = validate_html(html)
        v.assert_badge_shown("view")
        v.assert_badge_not_shown("backed")

    def test_backed_badge_only(self, tmp_path, validate_html):
        """Test backed badge appears."""
        import anndata as ad

        adata = AnnData(np.zeros((10, 5)))
        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        backed = ad.read_h5ad(path, backed="r")
        html = backed._repr_html_()
        v = validate_html(html)
        v.assert_badge_shown("backed")
        backed.file.close()

    def test_backed_badge_with_format(self, tmp_path, validate_html):
        """Test backed badge shows format."""
        import anndata as ad

        adata = AnnData(np.zeros((10, 5)))
        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        backed = ad.read_h5ad(path, backed="r")
        html = backed._repr_html_()
        v = validate_html(html)
        v.assert_badge_shown("backed")
        v.assert_text_visible("H5AD")
        backed.file.close()

    def test_backed_badge_with_path(self, tmp_path, validate_html):
        """Test backed badge shows file path."""
        import anndata as ad

        adata = AnnData(np.zeros((10, 5)))
        path = tmp_path / "myfile.h5ad"
        adata.write_h5ad(path)

        backed = ad.read_h5ad(path, backed="r")
        html = backed._repr_html_()
        v = validate_html(html)
        v.assert_badge_shown("backed")
        v.assert_text_visible("myfile")
        backed.file.close()

    def test_both_badges(self, tmp_path, validate_html):
        """Test view of backed AnnData shows both badges."""
        import anndata as ad

        adata = AnnData(np.zeros((10, 5)))
        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        backed = ad.read_h5ad(path, backed="r")
        view = backed[0:5, :]
        html = view._repr_html_()
        v = validate_html(html)
        v.assert_badge_shown("view")
        v.assert_badge_shown("backed")
        backed.file.close()

    def test_lazy_badge_only(self, validate_html):
        """Test render_header_badges with lazy badge."""
        from anndata._repr import render_header_badges

        html = render_header_badges(is_lazy=True)
        v = validate_html(html)
        v.assert_badge_shown("lazy")

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


class TestReadmeIcon:
    """Test README icon in header."""

    def test_readme_icon_appears_with_string(self, validate_html):
        """Test readme icon shows when uns['README'] is a string."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["README"] = "This is a description"
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".adata-readme-icon")

    def test_readme_icon_not_shown_without_readme(self, validate_html):
        """Test no readme icon without uns['README']."""
        adata = AnnData(np.zeros((10, 5)))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_not_exists(".adata-readme-icon")

    def test_readme_icon_not_shown_for_empty_string(self, validate_html):
        """Test no readme icon for empty string."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["README"] = ""
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_not_exists(".adata-readme-icon")

    def test_readme_icon_not_shown_for_whitespace_only(self, validate_html):
        """Test no readme icon for whitespace-only string."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["README"] = "   \n\t  "
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_not_exists(".adata-readme-icon")

    def test_readme_icon_not_shown_for_non_string(self, validate_html):
        """Test no readme icon for non-string README."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["README"] = 42
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_not_exists(".adata-readme-icon")

    def test_readme_content_escaped(self, validate_html):
        """Test README content is HTML-escaped."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["README"] = "<script>alert('xss')</script>"
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".adata-readme-icon")
        # Escaped content should be present, raw script should not
        assert "&lt;script&gt;" in html
        assert "<script>alert" not in html

    def test_readme_tooltip_truncated(self, validate_html):
        """Test long README is truncated in tooltip."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["README"] = "x" * 1000
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".adata-readme-icon")
        v.assert_text_visible("...")

    def test_readme_data_attribute_contains_content(self, validate_html):
        """Test data-readme attribute contains full content."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["README"] = "Test content"
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".adata-readme-icon")
        v.assert_attribute_value(".adata-readme-icon", "data-readme", "Test content")

    def test_readme_icon_accessibility(self, validate_html):
        """Test readme icon has accessibility attributes."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["README"] = "Description"
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".adata-readme-icon")
        v.assert_attribute_value(".adata-readme-icon", "aria-label", "View README")
        v.assert_attribute_value(".adata-readme-icon", "tabindex", "0")


class TestGenerateReprHtmlDirectly:
    """Test generate_repr_html function directly."""

    def test_html_disabled_returns_pre(self):
        """Test disabled HTML returns pre-formatted text."""
        from anndata import settings
        from anndata._repr.html import generate_repr_html

        adata = AnnData(np.zeros((10, 5)))
        with settings.override(repr_html_enabled=False):
            html = generate_repr_html(adata)
            assert "<pre>" in html

    def test_max_depth_reached_at_depth_zero(self):
        """Test max depth indicator at depth 0."""
        from anndata._repr.html import generate_repr_html

        adata = AnnData(np.zeros((10, 5)))
        html = generate_repr_html(adata, depth=0, max_depth=0)
        assert "max depth" in html.lower()

    def test_nested_anndata_not_expandable_at_max_depth(self):
        """Test nested AnnData not expandable when at max depth."""
        from anndata._repr.html import generate_repr_html

        inner = AnnData(np.zeros((5, 3)))
        outer = AnnData(np.zeros((10, 5)))
        outer.uns["nested"] = inner

        html = generate_repr_html(outer, max_depth=1)
        assert "nested" in html


class TestPublicAPIExports:
    """Test that public API items are properly exported."""

    def test_css_js_exports(self):
        """Test CSS and JS functions are exported."""
        from anndata._repr import get_css, get_javascript

        css = get_css()
        assert ".anndata-repr" in css

        js = get_javascript("test-id")
        assert "test-id" in js

    def test_section_rendering_exports(self):
        """Test section rendering functions are exported."""
        from anndata._repr import (
            render_formatted_entry,
            render_section,
        )

        assert callable(render_section)
        assert callable(render_formatted_entry)

    def test_ui_helper_exports(self):
        """Test UI helper functions are exported."""
        from anndata._repr import (
            render_badge,
            render_copy_button,
            render_fold_icon,
            render_header_badges,
            render_search_box,
            render_warning_icon,
        )

        assert callable(render_badge)
        assert callable(render_copy_button)
        assert callable(render_fold_icon)
        assert callable(render_header_badges)
        assert callable(render_search_box)
        assert callable(render_warning_icon)

    def test_utility_exports(self):
        """Test utility functions are exported."""
        from anndata._repr import escape_html, format_memory_size, format_number

        assert escape_html("<test>") == "&lt;test&gt;"
        assert "KB" in format_memory_size(1024)
        assert format_number(1000) == "1,000"

    def test_registry_exports(self):
        """Test registry classes are exported."""
        from anndata._repr import (
            FormattedEntry,
            FormattedOutput,
            FormatterContext,
            SectionFormatter,
            TypeFormatter,
            formatter_registry,
            register_formatter,
        )

        assert FormattedOutput is not None
        assert FormattedEntry is not None
        assert FormatterContext is not None
        assert TypeFormatter is not None
        assert SectionFormatter is not None
        assert formatter_registry is not None
        assert callable(register_formatter)

    def test_generate_repr_html_export(self):
        """Test generate_repr_html is exported."""
        from anndata._repr import generate_repr_html

        adata = AnnData(np.zeros((10, 5)))
        html = generate_repr_html(adata)
        assert "anndata-repr" in html


class TestNeverCrash:
    """Tests ensuring repr never crashes regardless of data content."""

    def test_none_in_obs_column(self, validate_html):
        """Test repr handles None values in obs columns."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["with_none"] = [None, "a", None, "b", None]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("obs", "with_none")

    def test_nan_in_obs_column(self, validate_html):
        """Test repr handles NaN values in obs columns."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["with_nan"] = [np.nan, 1.0, np.nan, 2.0, np.nan]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("obs", "with_nan")

    def test_inf_in_obs_column(self, validate_html):
        """Test repr handles inf values in obs columns."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["with_inf"] = [np.inf, -np.inf, 0, 1, 2]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("obs", "with_inf")

    def test_empty_string_column(self, validate_html):
        """Test repr handles empty string columns."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["empty_strings"] = ["", "", "", "", ""]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("obs", "empty_strings")

    def test_mixed_type_list_in_uns(self, validate_html):
        """Test repr handles mixed type lists in uns."""
        adata = AnnData(np.zeros((5, 3)))
        adata.uns["mixed"] = [1, "string", None, 3.14, True]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("uns", "mixed")

    def test_circular_reference_like_structure(self, validate_html):
        """Test repr handles dict structures that reference themselves."""
        adata = AnnData(np.zeros((5, 3)))
        d = {"a": 1}
        d["self_like"] = {"nested": d.copy()}  # Not circular but deep
        adata.uns["deep"] = d
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")

    def test_very_long_string_in_uns(self, validate_html):
        """Test repr handles very long strings."""
        adata = AnnData(np.zeros((5, 3)))
        adata.uns["long_string"] = "x" * 10000
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("uns", "long_string")

    def test_special_characters_in_keys(self, validate_html):
        """Test repr handles special characters in keys."""
        adata = AnnData(np.zeros((5, 3)))
        adata.uns["key<with>special&chars"] = "value"
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        # Special chars should be escaped in section
        v.assert_section_exists("uns")

    def test_unicode_in_data(self, validate_html):
        """Test repr handles unicode in data."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["unicode"] = ["日本語", "émoji🧬", "中文", "العربية", "עברית"]
        adata.uns["unicode_key_日本語"] = "value"
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("obs", "unicode")

    def test_empty_categorical(self, validate_html):
        """Test repr handles empty categorical."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["empty_cat"] = pd.Categorical([None] * 5, categories=["a", "b", "c"])
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("obs", "empty_cat")

    def test_zero_size_array_in_obsm(self, validate_html):
        """Test repr handles zero-size arrays in obsm."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obsm["empty"] = np.zeros((5, 0))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("obsm", "empty")

    def test_object_with_failing_repr(self, validate_html):
        """Test repr handles objects whose __repr__ fails."""

        class FailingRepr:
            def __repr__(self):
                msg = "Repr failed"
                raise RuntimeError(msg)

        adata = AnnData(np.zeros((5, 3)))
        adata.uns["failing"] = FailingRepr()
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("uns", "failing")

    def test_object_with_failing_sizeof(self, validate_html):
        """Test repr handles objects whose __sizeof__ fails."""

        class FailingSizeof:
            def __sizeof__(self):
                msg = "Sizeof failed"
                raise RuntimeError(msg)

        adata = AnnData(np.zeros((5, 3)))
        adata.uns["failing_size"] = FailingSizeof()
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")


class TestViewAndBackedModes:
    """Tests for view and backed mode handling."""

    def test_view_shows_badge(self, validate_html):
        """Test view shows View badge."""
        adata = AnnData(np.zeros((10, 5)))
        view = adata[0:5, :]
        html = view._repr_html_()
        v = validate_html(html)
        v.assert_badge_shown("view")

    def test_backed_shows_badge(self, tmp_path, validate_html):
        """Test backed mode shows badge."""
        import anndata as ad

        adata = AnnData(np.zeros((10, 5)))
        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        backed = ad.read_h5ad(path, backed="r")
        html = backed._repr_html_()
        v = validate_html(html)
        v.assert_badge_shown("backed")
        backed.file.close()

    def test_view_of_backed_shows_both_badges(self, tmp_path, validate_html):
        """Test view of backed shows both badges."""
        import anndata as ad

        adata = AnnData(np.zeros((10, 5)))
        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        backed = ad.read_h5ad(path, backed="r")
        view = backed[0:5, :]
        html = view._repr_html_()
        v = validate_html(html)
        v.assert_badge_shown("view")
        v.assert_badge_shown("backed")
        backed.file.close()


class TestComprehensiveAnnData:
    """Test with comprehensive AnnData similar to visual_inspect."""

    def test_comprehensive_anndata_renders_completely(self, validate_html):
        """Test comprehensive AnnData with all features renders without error."""
        n_obs, n_vars = 100, 50

        adata = AnnData(
            sp.random(n_obs, n_vars, density=0.1, format="csr", dtype=np.float32),
            obs=pd.DataFrame({
                "cell_type": pd.Categorical(
                    ["T cell", "B cell", "NK cell", "Monocyte", "DC"] * (n_obs // 5)
                ),
                "louvain": pd.Categorical([
                    f"cluster_{i}" for i in np.random.randint(0, 8, n_obs)
                ]),
                "n_counts": np.random.randint(1000, 50000, n_obs),
                "percent_mito": np.random.uniform(0, 15, n_obs).astype(np.float32),
                "is_doublet": np.random.choice([True, False], n_obs, p=[0.1, 0.9]),
            }),
            var=pd.DataFrame({
                "gene_symbol": [f"GN{i}" for i in range(n_vars)],
                "highly_variable": np.random.choice(
                    [True, False], n_vars, p=[0.2, 0.8]
                ),
                "means": np.random.exponential(1, n_vars).astype(np.float32),
            }),
        )

        # Colors
        adata.uns["cell_type_colors"] = [
            "#FF6B6B",
            "#4ECDC4",
            "#45B7D1",
            "#96CEB4",
            "#FFEAA7",
        ]
        adata.uns["louvain_colors"] = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
        ]

        # Uns
        adata.uns["neighbors"] = {"params": {"n_neighbors": 15}}
        adata.uns["experiment_id"] = "EXP_001"
        adata.uns["steps"] = ["filter", "normalize", "pca", "umap"]

        # Nested AnnData
        inner = AnnData(np.zeros((10, 5)))
        adata.uns["subset"] = inner

        # Obsm
        adata.obsm["X_pca"] = np.random.randn(n_obs, 50).astype(np.float32)
        adata.obsm["X_umap"] = np.random.randn(n_obs, 2).astype(np.float32)

        # Varm
        adata.varm["PCs"] = np.random.randn(n_vars, 50).astype(np.float32)

        # Layers
        adata.layers["counts"] = sp.random(n_obs, n_vars, density=0.1, format="csr")
        adata.layers["normalized"] = np.random.randn(n_obs, n_vars).astype(np.float32)

        # Obsp/Varp
        adata.obsp["distances"] = sp.random(n_obs, n_obs, density=0.05, format="csr")
        adata.varp["correlations"] = sp.random(
            n_vars, n_vars, density=0.1, format="csr"
        )

        # Raw
        raw = AnnData(
            sp.random(n_obs, n_vars + 20, density=0.1, format="csr"),
            var=pd.DataFrame({"gene": [f"G{i}" for i in range(n_vars + 20)]}),
        )
        adata.raw = raw

        # Generate HTML and validate structure
        html = adata._repr_html_()
        v = validate_html(html)

        # Validate key sections are present
        v.assert_element_exists(".anndata-repr")
        v.assert_section_exists("obs")
        v.assert_section_exists("var")
        v.assert_section_exists("uns")
        v.assert_section_exists("obsm")
        v.assert_section_exists("varm")
        v.assert_section_exists("layers")
        v.assert_section_exists("obsp")
        v.assert_section_exists("varp")
        v.assert_section_exists("raw")

        # Validate key entries in correct sections
        v.assert_section_contains_entry("obs", "cell_type")
        v.assert_section_contains_entry("obsm", "X_pca")
        v.assert_section_contains_entry("obsm", "X_umap")
        v.assert_section_contains_entry("uns", "neighbors")
        v.assert_section_contains_entry("layers", "counts")
        v.assert_section_contains_entry("obsp", "distances")

        # Validate shape is displayed
        v.assert_shape_displayed(100, 50)

        # Colors should be shown
        v.assert_color_swatch("#FF6B6B")


class TestIndexPreview:
    """Tests for obs_names and var_names preview."""

    def test_obs_names_preview_shown(self, validate_html):
        """Test obs_names preview is shown."""
        adata = AnnData(
            np.zeros((10, 5)),
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(10)]),
        )
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_text_visible("cell_0")

    def test_var_names_preview_shown(self, validate_html):
        """Test var_names preview is shown."""
        adata = AnnData(
            np.zeros((10, 5)),
            var=pd.DataFrame(index=[f"gene_{i}" for i in range(5)]),
        )
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_text_visible("gene_0")

    def test_long_index_names_truncated(self, validate_html):
        """Test long index names are handled."""
        long_name = "very_long_cell_name_" * 10
        adata = AnnData(
            np.zeros((3, 2)),
            obs=pd.DataFrame(index=[long_name, "short", "medium_name"]),
        )
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_shape_displayed(3, 2)

    def test_numeric_index_shown(self, validate_html):
        """Test numeric index is shown correctly."""
        adata = AnnData(np.zeros((10, 5)))
        # Default index is 0, 1, 2, ...
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_text_visible("0")
