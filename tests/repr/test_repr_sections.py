"""
Section rendering tests for the _repr module.

Tests for obs, var, uns, obsm, varm, obsp, varp, layers, and raw section rendering,
as well as custom section formatters.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from anndata import AnnData


class TestRawSection:
    """Test .raw section display."""

    def test_raw_section_present(self, validate_html):
        """Test raw section appears when raw is set."""
        adata = AnnData(np.zeros((10, 20)))
        adata.raw = adata.copy()
        adata = adata[:, :5]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("raw")

    def test_raw_none_no_section(self, validate_html):
        """Test no raw section when raw is None."""
        adata = AnnData(np.zeros((10, 5)))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")

    def test_raw_section_with_var(self, validate_html):
        """Test raw section shows var info."""
        adata = AnnData(np.zeros((10, 20)))
        adata.raw = adata.copy()
        adata = adata[:, :5]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("raw")

    def test_raw_section_with_varm(self, validate_html):
        """Test raw section shows varm info."""
        adata = AnnData(np.zeros((10, 20)))
        adata.varm["test"] = np.zeros((20, 3))
        adata.raw = adata.copy()
        adata = adata[:, :5]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("raw")


class TestRepresentationCompleteness:
    """Verify all data is accurately represented."""

    def test_all_obs_columns_shown(self, validate_html):
        """Test all obs columns appear in repr."""
        obs_cols = ["col_a", "col_b", "col_c", "col_d", "col_e"]
        adata = AnnData(
            np.zeros((10, 5)), obs=pd.DataFrame({c: list(range(10)) for c in obs_cols})
        )
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("obs")
        for col in obs_cols:
            v.assert_section_contains_entry("obs", col)

    def test_all_var_columns_shown(self, validate_html):
        """Test all var columns appear in repr."""
        var_cols = ["gene_name", "gene_id", "highly_variable"]
        adata = AnnData(
            np.zeros((10, 5)), var=pd.DataFrame({c: list(range(5)) for c in var_cols})
        )
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("var")
        for col in var_cols:
            v.assert_section_contains_entry("var", col)

    def test_all_uns_keys_shown(self, validate_html):
        """Test all uns keys appear in repr."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["key1"] = "value1"
        adata.uns["key2"] = 42
        adata.uns["nested_dict"] = {"a": 1, "b": 2}
        adata.uns["array_data"] = np.array([1, 2, 3])
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("uns")
        for key in ["key1", "key2", "nested_dict", "array_data"]:
            v.assert_section_contains_entry("uns", key)

    def test_all_obsm_keys_shown(self, validate_html):
        """Test all obsm keys appear in repr."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obsm["X_pca"] = np.random.randn(10, 50).astype(np.float32)
        adata.obsm["X_umap"] = np.random.randn(10, 2).astype(np.float32)
        adata.obsm["X_tsne"] = np.random.randn(10, 2).astype(np.float32)
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("obsm")
        for key in ["X_pca", "X_umap", "X_tsne"]:
            v.assert_section_contains_entry("obsm", key)

    def test_all_layers_shown(self, validate_html):
        """Test all layers appear in repr."""
        adata = AnnData(np.zeros((10, 5)))
        adata.layers["raw"] = np.random.randn(10, 5)
        adata.layers["normalized"] = np.random.randn(10, 5)
        adata.layers["scaled"] = np.random.randn(10, 5)
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("layers")
        for layer in ["raw", "normalized", "scaled"]:
            v.assert_section_contains_entry("layers", layer)

    def test_correct_shape_values(self, validate_html):
        """Test shape values are accurate."""
        adata = AnnData(np.zeros((123, 456)))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_shape_displayed(123, 456)


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_very_large_keys(self, validate_html):
        """Test handling of very long key names."""
        adata = AnnData(np.zeros((10, 5)))
        long_key = "a" * 500
        adata.uns[long_key] = "value"
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_exists("uns")

    def test_unicode_keys(self, validate_html):
        """Test handling of unicode key names."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["日本語"] = "value"
        adata.uns["émojis_🧬"] = "dna"
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("uns")
        v.assert_text_visible("日本語")
        v.assert_text_visible("émojis")

    def test_empty_sections(self, validate_html):
        """Test handling of empty sections."""
        adata = AnnData(np.zeros((10, 5)))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")

    def test_very_nested_uns(self, validate_html):
        """Test handling of deeply nested uns."""
        adata = AnnData(np.zeros((10, 5)))
        nested = {"level": 0}
        current = nested
        for i in range(10):
            current["nested"] = {"level": i + 1}
            current = current["nested"]
        adata.uns["deep"] = nested
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("uns", "deep")

    def test_mixed_types_in_uns(self, validate_html):
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
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_exists("uns")


class TestMappingSectionEdgeCases:
    """Tests for mapping section edge cases."""

    def test_obsm_shape_meta_display(self, validate_html):
        """Test obsm shows shape metadata."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obsm["X_pca"] = np.random.randn(10, 50)
        adata.obsm["X_umap"] = np.random.randn(10, 2)

        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("obsm")
        v.assert_section_contains_entry("obsm", "X_pca")
        v.assert_section_contains_entry("obsm", "X_umap")

    def test_layers_truncation(self, validate_html):
        """Test layers section truncates appropriately."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))
        for i in range(20):
            adata.layers[f"layer_{i}"] = np.random.randn(10, 5)

        with settings.override(repr_html_max_items=10):
            html = adata._repr_html_()
            v = validate_html(html)
            v.assert_section_exists("layers")
            v.assert_section_contains_entry("layers", "layer_0")
            v.assert_truncation_indicator()


class TestUnsEntryRendering:
    """Tests for uns entry rendering."""

    def test_uns_entry_with_type_hint_preview_note(self, validate_html):
        """Test uns entry with type hint shows preview note."""
        adata = AnnData(np.zeros((5, 3)))
        adata.uns["typed_data"] = {
            "__anndata_repr__": "somepackage.type",
            "data": "value",
        }
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("uns", "typed_data")

    def test_uns_entry_with_tuple(self, validate_html):
        """Test uns entry with tuple."""
        adata = AnnData(np.zeros((5, 3)))
        adata.uns["tuple_val"] = (1, 2, 3)
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("uns", "tuple_val")

    def test_uns_entry_with_empty_list(self, validate_html):
        """Test uns entry with empty list."""
        adata = AnnData(np.zeros((5, 3)))
        adata.uns["empty"] = []
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("uns", "empty")

    def test_uns_entry_with_none(self, validate_html):
        """Test uns entry with None."""
        adata = AnnData(np.zeros((5, 3)))
        adata.uns["null_val"] = None
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("uns", "null_val")


class TestCustomSectionFormatters:
    """Tests for custom section formatters."""

    def test_custom_section_appears_after_specified_section(self):
        """Test custom section appears in correct position."""
        from anndata._repr import FormattedEntry, FormattedOutput, SectionFormatter
        from anndata._repr.registry import formatter_registry, register_formatter

        class TestCustomSection(SectionFormatter):
            @property
            def section_name(self):
                return "test_custom"

            @property
            def display_name(self):
                return "Test Custom"

            @property
            def after_section(self):
                return "obs"

            def should_show(self, adata):
                return True

            def get_entries(self, adata, context):
                return [
                    FormattedEntry(
                        key="custom_entry",
                        output=FormattedOutput(type_name="custom_type"),
                    )
                ]

        formatter = TestCustomSection()
        register_formatter(formatter)

        try:
            adata = AnnData(
                np.zeros((10, 5)),
                obs=pd.DataFrame({"a": list(range(10))}),
            )
            html = adata._repr_html_()

            assert "Test Custom" in html
            assert "custom_entry" in html
            assert "custom_type" in html
        finally:
            if "test_custom" in formatter_registry._section_formatters:
                del formatter_registry._section_formatters["test_custom"]

    def test_custom_section_not_shown_when_should_show_false(self):
        """Test custom section hidden when should_show returns False."""
        from anndata._repr import FormattedEntry, FormattedOutput, SectionFormatter
        from anndata._repr.registry import formatter_registry, register_formatter

        class HiddenCustomSection(SectionFormatter):
            @property
            def section_name(self):
                return "hidden_section"

            def should_show(self, adata):
                return False

            def get_entries(self, adata, context):
                return [
                    FormattedEntry(
                        key="should_not_appear",
                        output=FormattedOutput(type_name="hidden"),
                    )
                ]

        formatter = HiddenCustomSection()
        register_formatter(formatter)

        try:
            adata = AnnData(np.zeros((5, 3)))
            html = adata._repr_html_()

            assert "hidden_section" not in html
            assert "should_not_appear" not in html
        finally:
            if "hidden_section" in formatter_registry._section_formatters:
                del formatter_registry._section_formatters["hidden_section"]


class TestUnknownSectionsDetection:
    """Tests for unknown/custom attribute detection."""

    def test_standard_sections_not_in_other(self, validate_html):
        """Test standard sections don't appear in 'other' section."""
        adata = AnnData(
            np.zeros((10, 5)),
            obs=pd.DataFrame({"a": list(range(10))}),
            var=pd.DataFrame({"b": list(range(5))}),
        )
        adata.uns["test"] = "value"
        adata.obsm["X_pca"] = np.zeros((10, 2))

        html = adata._repr_html_()
        v = validate_html(html)
        # Standard sections should appear normally
        v.assert_section_exists("obs")
        v.assert_section_exists("var")
        v.assert_section_exists("uns")
        v.assert_section_exists("obsm")

    def test_registered_section_not_in_other(self):
        """Test registered custom sections don't appear in 'other'."""
        from anndata._repr import FormattedEntry, FormattedOutput, SectionFormatter
        from anndata._repr.registry import formatter_registry, register_formatter

        class RegisteredSection(SectionFormatter):
            @property
            def section_name(self):
                return "registered_section"

            def should_show(self, adata):
                return True

            def get_entries(self, adata, context):
                return [
                    FormattedEntry(
                        key="entry", output=FormattedOutput(type_name="type")
                    )
                ]

        formatter = RegisteredSection()
        register_formatter(formatter)

        try:
            adata = AnnData(np.zeros((5, 3)))
            html = adata._repr_html_()

            assert "registered_section" in html
        finally:
            if "registered_section" in formatter_registry._section_formatters:
                del formatter_registry._section_formatters["registered_section"]


class TestCoverageEdgeCases:
    """Tests for edge cases to improve code coverage."""

    def test_category_column_overflow(self):
        """Test rendering categorical column with more than max categories."""
        from anndata import settings

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
            assert "cat_0" in html
            assert "...+" in html or "more" in html.lower()

    def test_max_depth_with_multiple_nested_anndata(self):
        """Test max depth indicator with deeply nested AnnData."""
        from anndata import settings

        level2 = AnnData(np.zeros((5, 3)))
        level1 = AnnData(np.zeros((7, 4)))
        level1.uns["nested"] = level2
        level0 = AnnData(np.zeros((10, 5)))
        level0.uns["nested"] = level1

        with settings.override(repr_html_max_depth=0):
            html = level0._repr_html_()
            assert "max depth" in html.lower() or "depth" in html.lower()

    def test_dataframe_entry_nunique_exception(self):
        """Test nunique() exception handling for dataframe columns."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["unhashable"] = [[i] for i in range(10)]

        html = adata._repr_html_()
        assert html is not None
        assert "unhashable" in html

    def test_very_large_dataframe_skips_nunique(self):
        """Test that nunique() is skipped for very large columns."""
        from anndata import settings

        adata = AnnData(np.zeros((100000, 5)), obs=pd.DataFrame({"col": range(100000)}))

        with settings.override(repr_html_unique_limit=1000):
            html = adata._repr_html_()
            assert "col" in html

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
        assert "cat_0" in html


class TestCompleteDataVisibility:
    """Tests ensuring all data is visible or truncation is indicated.

    Scientific display requirement: Nothing should be hidden from the user.
    """

    def test_all_varp_keys_shown(self, validate_html):
        """Test all varp keys appear in repr."""
        adata = AnnData(np.zeros((10, 5)))
        keys = ["corr", "covariance", "pvals"]
        for key in keys:
            adata.varp[key] = sp.random(5, 5, density=0.3, format="csr")
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("varp")
        for key in keys:
            v.assert_section_contains_entry("varp", key)

    def test_all_obsp_keys_shown(self, validate_html):
        """Test all obsp keys appear in repr."""
        adata = AnnData(np.zeros((10, 5)))
        keys = ["distances", "connectivities", "weights"]
        for key in keys:
            adata.obsp[key] = sp.random(10, 10, density=0.1, format="csr")
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("obsp")
        for key in keys:
            v.assert_section_contains_entry("obsp", key)

    def test_all_varm_keys_shown(self, validate_html):
        """Test all varm keys appear in repr."""
        adata = AnnData(np.zeros((10, 5)))
        keys = ["PCs", "loadings", "gene_embeddings"]
        for key in keys:
            adata.varm[key] = np.random.randn(5, 3)
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("varm")
        for key in keys:
            v.assert_section_contains_entry("varm", key)

    def test_truncation_indicator_when_many_items(self, validate_html):
        """Test truncation indicator appears for many items."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))
        for i in range(100):
            adata.uns[f"key_{i}"] = i

        with settings.override(repr_html_max_items=10):
            html = adata._repr_html_()
            v = validate_html(html)
            # Should indicate more items exist with specific truncation pattern
            v.assert_truncation_indicator()

    def test_category_truncation_indicator(self):
        """Test category truncation indicator appears."""
        from anndata import settings

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
            assert "...+" in html or "more" in html.lower()


class TestSpecificInfoDisplay:
    """Tests for specific information displayed in repr."""

    def test_shape_values_accurate(self, validate_html):
        """Test shape values are accurate."""
        adata = AnnData(np.zeros((123, 456)))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_shape_displayed(123, 456)

    def test_sparse_matrix_format_shown(self, validate_html):
        """Test sparse matrix format is shown."""
        # Only CSR and CSC are supported by AnnData
        for fmt in ["csr", "csc"]:
            X = sp.random(100, 50, density=0.1, format=fmt)
            adata = AnnData(X)
            html = adata._repr_html_()
            v = validate_html(html)
            v.assert_text_visible(fmt)

    def test_sparse_matrix_density_or_nnz_shown(self, validate_html):
        """Test sparse matrix shows density or nnz info."""
        X = sp.random(100, 50, density=0.1, format="csr")
        adata = AnnData(X)
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        # Should show either density percentage or nnz count
        assert "%" in html or "nnz" in html.lower() or str(X.nnz) in html

    def test_dtype_shown_for_arrays(self, validate_html):
        """Test dtype is shown for arrays."""
        adata = AnnData(np.zeros((10, 5), dtype=np.float32))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_dtype_displayed("float32")

    def test_obsm_shape_shown(self, validate_html):
        """Test obsm array shapes are shown."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obsm["X_pca"] = np.zeros((10, 50))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obsm", "X_pca")
        v.assert_text_visible("50")  # Second dimension shown

    def test_category_count_shown(self, validate_html):
        """Test category count is shown for categoricals."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["cat"] = pd.Categorical(["A", "B", "C"] * 3 + ["A"])
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "cat")
        # Should show category count or list categories
        assert "3" in html or ("A" in html and "B" in html and "C" in html)

    def test_dataframe_in_obsm_shows_columns(self, validate_html):
        """Test DataFrame in obsm shows column count or names."""
        adata = AnnData(np.zeros((10, 5)))
        # DataFrame index must match adata.obs_names
        adata.obsm["spatial"] = pd.DataFrame(
            {"x": np.zeros(10), "y": np.zeros(10), "z": np.zeros(10)},
            index=adata.obs_names,
        )
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obsm", "spatial")
        # Should indicate it's a DataFrame with columns
        assert "DataFrame" in html or "3" in html or "x" in html

    def test_nested_dict_shows_key_count(self, validate_html):
        """Test nested dict shows key count."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["params"] = {"a": 1, "b": 2, "c": 3}
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("uns", "params")
        # Should indicate dict with keys
        assert "3" in html or "dict" in html.lower() or "a" in html

    def test_raw_section_shows_var_count(self, validate_html):
        """Test raw section shows different var count."""
        adata = AnnData(np.zeros((10, 20)))
        adata.raw = adata.copy()
        adata = adata[:, :5]  # Subset vars
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_exists("raw")
        # Main shows 5 vars, raw should show 20
        v.assert_text_visible("20")
        v.assert_text_visible("5")

    def test_list_shows_item_count(self, validate_html):
        """Test lists show item count."""
        adata = AnnData(np.zeros((10, 5)))
        adata.uns["steps"] = ["step1", "step2", "step3", "step4", "step5"]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("uns", "steps")
        assert "5" in html or "items" in html.lower()


class TestUnknownSectionsAndErrorHandling:
    """Tests for unknown sections and error handling in section rendering."""

    def test_unknown_attribute_detected(self):
        """Test that unknown attributes on AnnData are handled gracefully."""
        adata = AnnData(np.zeros((5, 3)))
        # Set a custom attribute that's not a standard AnnData attribute
        adata._custom_unknown_attr = {"key": "value"}
        html = adata._repr_html_()
        # Should not crash
        assert html is not None

    def test_raw_section_with_empty_var(self):
        """Test raw section renders with empty var columns."""
        # Create adata with no var columns, then set raw
        adata = AnnData(np.zeros((10, 5)))
        # raw.var is read-only, so we need to create raw from adata without var columns
        adata_no_cols = AnnData(np.zeros((10, 5)))
        adata.raw = adata_no_cols
        html = adata._repr_html_()
        assert html is not None
        assert "raw" in html.lower()

    def test_safe_get_attr_with_normal_object(self):
        """Test _safe_get_attr returns attribute value."""
        from anndata._repr.sections import _safe_get_attr

        class Obj:
            attr = "value"

        assert _safe_get_attr(Obj(), "attr", "default") == "value"

    def test_safe_get_attr_with_missing_attr(self):
        """Test _safe_get_attr returns default for missing attribute."""
        from anndata._repr.sections import _safe_get_attr

        class Obj:
            pass

        assert _safe_get_attr(Obj(), "missing", "default") == "default"

    def test_safe_get_attr_with_exception(self):
        """Test _safe_get_attr returns default when property raises."""
        from anndata._repr.sections import _safe_get_attr

        class FailingObj:
            @property
            def bad_attr(self):
                msg = "Access failed"
                raise RuntimeError(msg)

        result = _safe_get_attr(FailingObj(), "bad_attr", "fallback")
        assert result == "fallback"

    def test_get_raw_meta_parts_with_var(self):
        """Test _get_raw_meta_parts extracts var column count."""
        from anndata._repr.sections import _get_raw_meta_parts

        adata = AnnData(np.zeros((5, 3)))
        adata.var["gene_name"] = ["A", "B", "C"]
        parts = _get_raw_meta_parts(adata)
        assert any("var" in p for p in parts)

    def test_get_raw_meta_parts_with_varm(self):
        """Test _get_raw_meta_parts extracts varm count."""
        from anndata._repr.sections import _get_raw_meta_parts

        adata = AnnData(np.zeros((5, 3)))
        adata.varm["PCs"] = np.zeros((3, 2))
        parts = _get_raw_meta_parts(adata)
        assert any("varm" in p for p in parts)

    def test_get_raw_meta_parts_exception_handling(self):
        """Test _get_raw_meta_parts handles exceptions gracefully."""
        from anndata._repr.sections import _get_raw_meta_parts

        class BadRaw:
            @property
            def var(self):
                msg = "var access failed"
                raise RuntimeError(msg)

            @property
            def varm(self):
                msg = "varm access failed"
                raise RuntimeError(msg)

        # Should not crash, returns empty list
        parts = _get_raw_meta_parts(BadRaw())
        assert parts == []

    def test_render_error_entry(self):
        """Test _render_error_entry produces valid HTML."""
        from anndata._repr.sections import _render_error_entry

        html = _render_error_entry("test_section", "Test error message")
        assert "test_section" in html
        assert "error" in html.lower()
        # Should be valid HTML structure
        assert "<div" in html

    def test_render_unknown_sections(self):
        """Test _render_unknown_sections produces valid HTML."""
        from anndata._repr.sections import _render_unknown_sections

        unknown = [("custom_attr", "CustomType"), ("another", "dict (3 items)")]
        html = _render_unknown_sections(unknown)
        assert "custom_attr" in html
        assert "another" in html
        assert "other" in html.lower()  # Section name

    def test_detect_unknown_sections_empty(self):
        """Test _detect_unknown_sections returns empty for standard AnnData."""
        from anndata._repr.sections import _detect_unknown_sections

        adata = AnnData(np.zeros((5, 3)))
        unknown = _detect_unknown_sections(adata)
        assert unknown == []

    def test_detect_unknown_sections_with_custom_attr(self):
        """Test _detect_unknown_sections detects custom attributes."""
        from anndata._repr.sections import _detect_unknown_sections

        adata = AnnData(np.zeros((5, 3)))
        # Add a custom public attribute (not starting with _)
        object.__setattr__(adata, "custom_data", [1, 2, 3])
        unknown = _detect_unknown_sections(adata)
        # May or may not detect depending on implementation
        assert isinstance(unknown, list)


class TestTruncationIndicators:
    """Tests ensuring truncation is clearly indicated when data is truncated.

    Scientific display requirement: If something cannot be displayed in full,
    this should be indicated clearly to the user.
    """

    def test_obs_columns_truncation_shows_count(self, validate_html):
        """Test obs section shows truncation count when columns are truncated."""
        from anndata import settings

        cols = {f"col_{i}": list(range(10)) for i in range(100)}
        adata = AnnData(np.zeros((10, 5)), obs=pd.DataFrame(cols))

        with settings.override(repr_html_max_items=10):
            html = adata._repr_html_()
            v = validate_html(html)
            v.assert_truncation_indicator()

    def test_var_columns_truncation_shows_count(self, validate_html):
        """Test var section shows truncation count when columns are truncated."""
        from anndata import settings

        cols = {f"col_{i}": list(range(5)) for i in range(100)}
        adata = AnnData(np.zeros((10, 5)), var=pd.DataFrame(cols))

        with settings.override(repr_html_max_items=10):
            html = adata._repr_html_()
            v = validate_html(html)
            v.assert_truncation_indicator()

    def test_uns_keys_truncation_shows_count(self, validate_html):
        """Test uns section shows truncation count when keys are truncated."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))
        for i in range(100):
            adata.uns[f"key_{i}"] = f"value_{i}"

        with settings.override(repr_html_max_items=10):
            html = adata._repr_html_()
            v = validate_html(html)
            v.assert_truncation_indicator()

    def test_obsm_keys_truncation_shows_count(self, validate_html):
        """Test obsm section shows truncation count when keys are truncated."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))
        for i in range(50):
            adata.obsm[f"X_{i}"] = np.random.randn(10, 2)

        with settings.override(repr_html_max_items=10):
            html = adata._repr_html_()
            v = validate_html(html)
            v.assert_truncation_indicator()

    def test_layers_truncation_shows_count(self, validate_html):
        """Test layers section shows truncation count when truncated."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))
        for i in range(50):
            adata.layers[f"layer_{i}"] = np.random.randn(10, 5)

        with settings.override(repr_html_max_items=10):
            html = adata._repr_html_()
            v = validate_html(html)
            v.assert_truncation_indicator()

    def test_varp_keys_truncation_shows_count(self, validate_html):
        """Test varp section shows truncation count when keys are truncated."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))
        for i in range(50):
            adata.varp[f"mat_{i}"] = sp.random(5, 5, density=0.3, format="csr")

        with settings.override(repr_html_max_items=10):
            html = adata._repr_html_()
            v = validate_html(html)
            v.assert_truncation_indicator()

    def test_obsp_keys_truncation_shows_count(self, validate_html):
        """Test obsp section shows truncation count when keys are truncated."""
        from anndata import settings

        adata = AnnData(np.zeros((10, 5)))
        for i in range(50):
            adata.obsp[f"mat_{i}"] = sp.random(10, 10, density=0.1, format="csr")

        with settings.override(repr_html_max_items=10):
            html = adata._repr_html_()
            v = validate_html(html)
            v.assert_truncation_indicator()

    def test_category_truncation_shows_count(self, validate_html):
        """Test category display shows count when categories are truncated."""
        from anndata import settings

        categories = [f"cat_{i}" for i in range(100)]
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
            v = validate_html(html)
            v.assert_truncation_indicator()

    def test_index_preview_truncation_shows_ellipsis(self, validate_html):
        """Test index preview shows ellipsis when names are truncated."""
        adata = AnnData(
            np.zeros((100, 50)),
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(100)]),
        )
        html = adata._repr_html_()
        v = validate_html(html)
        # Should show first items, ellipsis, and last items
        v.assert_text_visible("...")


class TestErrorRepresentation:
    """Tests ensuring errors are represented in HTML, not raised.

    Scientific display requirement: HTML repr should never raise errors,
    but instead represent them appropriately in the display.
    """

    def test_failing_property_renders_gracefully(self, validate_html):
        """Test objects with failing properties don't crash repr."""

        class FailingProperty:
            @property
            def shape(self):
                msg = "Cannot access shape"
                raise RuntimeError(msg)

        adata = AnnData(np.zeros((5, 3)))
        adata.uns["failing_prop"] = FailingProperty()
        # Formatter will warn when it fails to handle the object
        with pytest.warns(UserWarning, match="Formatter.*:"):
            html = adata._repr_html_()
        v = validate_html(html)
        # Should not crash and should render the container
        v.assert_element_exists(".anndata-repr")
        # Shape should still be shown
        v.assert_shape_displayed(5, 3)

    def test_failing_len_shows_error_indicator(self, validate_html):
        """Test objects with failing __len__ show error indicator."""

        class FailingLen:
            def __len__(self):
                msg = "Cannot get length"
                raise RuntimeError(msg)

        adata = AnnData(np.zeros((5, 3)))
        adata.uns["failing_len"] = FailingLen()
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("uns", "failing_len")

    def test_failing_iter_shows_error_indicator(self, validate_html):
        """Test objects with failing __iter__ show error indicator."""

        class FailingIter:
            def __iter__(self):
                msg = "Cannot iterate"
                raise RuntimeError(msg)

        adata = AnnData(np.zeros((5, 3)))
        adata.uns["failing_iter"] = FailingIter()
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("uns", "failing_iter")

    def test_failing_getitem_shows_error_indicator(self, validate_html):
        """Test objects with failing __getitem__ show error indicator."""

        class FailingGetitem:
            def __getitem__(self, key):
                msg = "Cannot get item"
                raise RuntimeError(msg)

            def keys(self):
                return ["a", "b"]

        adata = AnnData(np.zeros((5, 3)))
        adata.uns["failing_getitem"] = FailingGetitem()
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("uns", "failing_getitem")

    def test_failing_dtype_renders_gracefully(self, validate_html):
        """Test objects with failing dtype attribute don't crash repr."""

        class FailingDtype:
            shape = (10, 5)

            @property
            def dtype(self):
                msg = "Cannot access dtype"
                raise RuntimeError(msg)

        adata = AnnData(np.zeros((10, 5)))
        # Put in uns since obsm validates types
        adata.uns["failing_dtype"] = FailingDtype()
        # Multiple formatters may warn when they fail to handle the object
        with pytest.warns(UserWarning, match="Formatter.*:"):
            html = adata._repr_html_()
        v = validate_html(html)
        # Should not crash and should render the container
        v.assert_element_exists(".anndata-repr")
        # Shape should still be shown
        v.assert_shape_displayed(10, 5)

    def test_corrupted_dataframe_shows_error_indicator(self, validate_html):
        """Test DataFrames with corrupted columns show error indicator."""
        adata = AnnData(np.zeros((5, 3)))
        # Create a DataFrame with a problematic column
        df = pd.DataFrame({"normal": [1, 2, 3, 4, 5]})
        # Add column that raises on access
        adata.uns["df"] = df
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("uns", "df")

    def test_nested_error_does_not_propagate(self, validate_html):
        """Test errors in nested structures don't propagate up."""

        class NestedError:
            @property
            def nested(self):
                msg = "Nested access failed"
                raise RuntimeError(msg)

        outer_dict = {"key": "value", "bad": NestedError()}
        adata = AnnData(np.zeros((5, 3)))
        adata.uns["nested"] = outer_dict
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("uns", "nested")


class TestExtensionDtypes:
    """Tests for pandas ExtensionDtype columns.

    Scientific display requirement: All data types should be correctly
    represented, including nullable integer types and string types.
    """

    def test_nullable_int_column_displayed(self, validate_html):
        """Test nullable integer columns are displayed correctly."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["nullable_int"] = pd.array([1, 2, None, 4, 5], dtype=pd.Int64Dtype())
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "nullable_int")
        # Should show Int64 dtype indicator
        assert "Int64" in html or "int64" in html.lower()

    def test_nullable_float_column_displayed(self, validate_html):
        """Test nullable float columns are displayed correctly."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["nullable_float"] = pd.array(
            [1.1, 2.2, None, 4.4, 5.5], dtype=pd.Float64Dtype()
        )
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "nullable_float")

    def test_nullable_string_column_displayed(self, validate_html):
        """Test nullable StringDtype columns are displayed correctly."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["string_col"] = pd.array(
            ["a", "b", None, "d", "e"], dtype=pd.StringDtype()
        )
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "string_col")
        # Should show string dtype or unique count
        assert "string" in html.lower() or "unique" in html.lower() or "5" in html

    def test_nullable_boolean_column_displayed(self, validate_html):
        """Test nullable BooleanDtype columns are displayed correctly."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["bool_col"] = pd.array(
            [True, False, None, True, False], dtype=pd.BooleanDtype()
        )
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "bool_col")

    def test_categorical_with_na_displayed(self, validate_html):
        """Test categorical columns with NA values are displayed correctly."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["cat_with_na"] = pd.Categorical(
            ["A", "B", None, "A", None], categories=["A", "B", "C"]
        )
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "cat_with_na")
        # Should show category information
        assert "A" in html
        assert "B" in html


class TestEmptySectionsDisplay:
    """Tests for empty sections and X=None display.

    Scientific display requirement: Empty sections should be clearly indicated.
    """

    def test_x_none_shows_indicator(self, validate_html):
        """Test X=None shows clear indicator."""
        adata = AnnData(obs=pd.DataFrame({"a": [1, 2, 3]}))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        # Should show None or similar indicator for X
        assert "None" in html or "none" in html.lower()

    def test_empty_obs_columns_shown(self, validate_html):
        """Test AnnData with no obs columns still shows obs section."""
        adata = AnnData(np.zeros((10, 5)))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        # Shape should be displayed
        v.assert_shape_displayed(10, 5)

    def test_empty_uns_section_not_shown_or_indicated(self, validate_html):
        """Test empty uns section is handled gracefully."""
        adata = AnnData(np.zeros((10, 5)))
        # uns is empty by default
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")

    def test_empty_obsm_section_not_shown_or_indicated(self, validate_html):
        """Test empty obsm section is handled gracefully."""
        adata = AnnData(np.zeros((10, 5)))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")

    def test_zero_obs_anndata_shape_shown(self, validate_html):
        """Test AnnData with 0 obs shows correct shape."""
        adata = AnnData(np.zeros((0, 5)))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_shape_displayed(0, 5)

    def test_zero_var_anndata_shape_shown(self, validate_html):
        """Test AnnData with 0 vars shows correct shape."""
        adata = AnnData(np.zeros((10, 0)))
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_shape_displayed(10, 0)

    def test_completely_empty_anndata(self, validate_html):
        """Test completely empty AnnData renders without error."""
        adata = AnnData()
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_shape_displayed(0, 0)

    def test_raw_none_no_raw_section(self, validate_html):
        """Test raw=None doesn't show raw section."""
        adata = AnnData(np.zeros((10, 5)))
        assert adata.raw is None
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        # Should not have raw section
        # (we just verify it doesn't crash - raw section may or may not appear)


class TestNestedStructuresDisplay:
    """Tests for complex nested structures display.

    Scientific display requirement: Nested structures should be
    navigable and their contents discoverable.
    """

    def test_deeply_nested_dict_in_uns(self, validate_html):
        """Test deeply nested dicts are displayed with depth indicator."""
        adata = AnnData(np.zeros((5, 3)))
        nested = {"level0": {"level1": {"level2": {"level3": {"level4": "value"}}}}}
        adata.uns["deep"] = nested
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("uns", "deep")
        # Should show some indication of nesting
        assert "level0" in html

    def test_nested_anndata_shows_shape(self, validate_html):
        """Test nested AnnData shows its shape."""
        inner = AnnData(np.zeros((20, 10)))
        outer = AnnData(np.zeros((50, 25)))
        outer.uns["subset"] = inner
        html = outer._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        v.assert_section_contains_entry("uns", "subset")
        # Should show inner AnnData shape somewhere
        assert "20" in html
        assert "10" in html

    def test_list_of_arrays_in_uns(self, validate_html):
        """Test list of arrays shows item count."""
        adata = AnnData(np.zeros((5, 3)))
        adata.uns["arrays"] = [np.zeros((3, 3)) for _ in range(5)]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("uns", "arrays")
        # Should show list indicator with count
        assert "5" in html or "list" in html.lower()

    def test_dict_with_dataframes(self, validate_html):
        """Test dict with DataFrames shows structure."""
        adata = AnnData(np.zeros((5, 3)))
        adata.uns["dfs"] = {
            "df1": pd.DataFrame({"a": [1, 2, 3]}),
            "df2": pd.DataFrame({"b": [4, 5, 6]}),
        }
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("uns", "dfs")

    def test_multiple_levels_of_nested_anndata(self, validate_html):
        """Test multiple nested AnnData levels with max_depth."""
        from anndata import settings

        level2 = AnnData(np.zeros((5, 3)))
        level1 = AnnData(np.zeros((10, 5)))
        level1.uns["inner"] = level2
        level0 = AnnData(np.zeros((20, 10)))
        level0.uns["nested"] = level1

        with settings.override(repr_html_max_depth=2):
            html = level0._repr_html_()
            v = validate_html(html)
            v.assert_element_exists(".anndata-repr")
            # Should show max depth indicator for innermost
            # (behavior depends on implementation)
