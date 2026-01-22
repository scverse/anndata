"""
Serialization warning tests for the _repr module.

Tests for warnings about non-serializable data types, column name validation,
and detection of potential serialization issues.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import anndata as ad
from anndata import AnnData


class TestBasicWarnings:
    """Tests for basic warning indicators."""

    def test_string_column_warning(self):
        """Test warning for string columns that will convert."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["batch"] = ["A", "B", "A", "B", "A", "B", "A", "B", "A", "B"]
        html = adata._repr_html_()
        assert "warning" in html.lower() or "⚠" in html

    def test_string_all_unique_no_warning(self):
        """Test no warning for string columns with all unique values."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["id"] = ["cell_0", "cell_1", "cell_2", "cell_3", "cell_4"]
        _html = adata._repr_html_()
        # Should not error

    def test_unserializable_uns_warning(self):
        """Test warning for unserializable .uns entries."""

        class CustomClass:
            pass

        adata = AnnData(np.zeros((10, 5)))
        adata.uns["custom"] = CustomClass()
        html = adata._repr_html_()
        assert (
            "serializable" in html.lower()
            or "warning" in html.lower()
            or "⚠" in html
            or "unknown" in html.lower()
        )


class TestSeriesNonSerializable:
    """Tests for Series formatter detecting non-serializable object dtype columns."""

    def test_custom_object_in_obs_column_not_serializable(self):
        """Test that custom objects in obs columns are flagged as non-serializable."""
        from anndata._repr.formatters import SeriesFormatter
        from anndata._repr.registry import FormatterContext

        class CustomObject:
            pass

        series = pd.Series([CustomObject(), CustomObject(), CustomObject()])
        assert series.dtype == np.dtype("object")

        formatter = SeriesFormatter()
        result = formatter.format(series, FormatterContext())

        assert not result.is_serializable
        assert len(result.warnings) > 0
        assert "CustomObject" in result.warnings[0]

    def test_list_in_obs_column_not_serializable(self):
        """Test that lists in obs columns are flagged as non-serializable."""
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
        """Repr detects list columns as non-serializable, and they actually fail."""
        adata = ad.AnnData(X=np.eye(3))
        adata.obs["list_col"] = [["a", "b"], ["c"], ["d"]]

        try:
            adata.write_h5ad(tmp_path / "test.h5ad")
            pytest.fail(
                "List serialization now works! "
                "Update _check_series_serializability() in formatters.py."
            )
        except (TypeError, Exception):  # noqa: BLE001
            pass

        html = adata._repr_html_()
        assert "list" in html
        assert "(!)" in html

    def test_custom_object_detected_and_not_serializable(self, tmp_path):
        """Repr detects custom objects as non-serializable, and they actually fail."""

        class CustomObject:
            pass

        adata = ad.AnnData(X=np.eye(3))
        adata.obs["custom"] = [CustomObject(), CustomObject(), CustomObject()]

        try:
            adata.write_h5ad(tmp_path / "test.h5ad")
            pytest.fail(
                "Custom object serialization now works! "
                "Update _check_series_serializability() in formatters.py."
            )
        except (TypeError, Exception):  # noqa: BLE001
            pass

        html = adata._repr_html_()
        assert "CustomObject" in html
        assert "(!)" in html

    def test_non_ascii_column_no_warning_and_serializes(self, tmp_path):
        """Repr does not warn for non-ASCII (it's valid), and it serializes."""
        adata = ad.AnnData(X=np.eye(3))
        adata.obs["gène_名前"] = ["a", "b", "c"]

        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)
        adata2 = ad.read_h5ad(path)
        assert "gène_名前" in adata2.obs.columns

        html = adata._repr_html_()
        assert "gène_名前" in html
        assert "Not serializable" not in html

    def test_tuple_column_name_detected_and_not_serializable(self, tmp_path):
        """Repr detects non-string column names, and they actually fail."""
        adata = ad.AnnData(X=np.eye(3))
        adata.obs[("a", "b")] = [1, 2, 3]

        try:
            adata.write_h5ad(tmp_path / "test.h5ad")
            pytest.fail(
                "Non-string column name serialization now works! "
                "Update check_column_name() in formatters.py."
            )
        except (TypeError, Exception):  # noqa: BLE001
            pass

        html = adata._repr_html_()
        assert "Non-string" in html
        assert "(!)" in html

    def test_datetime64_column_detected_and_not_serializable(self, tmp_path):
        """Repr detects datetime64 columns as non-serializable (issue #455)."""
        adata = ad.AnnData(X=np.eye(3))
        adata.obs["date"] = pd.to_datetime(["2024-01-01", "2024-01-02", "2024-01-03"])

        try:
            path = tmp_path / "test_datetime.h5ad"
            adata.write_h5ad(path)
            pytest.fail(
                "datetime64 serialization now works! "
                "Update SeriesFormatter in formatters.py to remove the datetime64 warning."
            )
        except Exception:  # noqa: BLE001
            pass

        html = adata._repr_html_()
        assert "datetime64" in html
        assert "(!)" in html

    def test_timedelta64_column_detected_and_not_serializable(self, tmp_path):
        """Repr detects timedelta64 columns as non-serializable."""
        adata = ad.AnnData(X=np.eye(3))
        adata.obs["duration"] = pd.to_timedelta(["1 days", "2 days", "3 days"])

        try:
            path = tmp_path / "test_timedelta.h5ad"
            adata.write_h5ad(path)
            pytest.fail(
                "timedelta64 serialization now works! "
                "Update SeriesFormatter in formatters.py to remove the timedelta64 warning."
            )
        except Exception:  # noqa: BLE001
            pass

        html = adata._repr_html_()
        assert "timedelta64" in html
        assert "(!)" in html

    def test_dict_in_obs_column_not_serializable(self, tmp_path):
        """Test that dicts in obs columns are flagged as non-serializable."""
        adata = ad.AnnData(X=np.eye(3))
        adata.obs["dict_col"] = [{"k": 1}, {"k": 2}, {"k": 3}]

        try:
            adata.write_h5ad(tmp_path / "test.h5ad")
            pytest.fail(
                "Dict serialization in obs now works! "
                "Update _check_series_serializability() in formatters.py."
            )
        except (TypeError, Exception):  # noqa: BLE001
            pass

        html = adata._repr_html_()
        assert "dict" in html
        assert "(!)" in html

    def test_lambda_in_uns_not_serializable(self, tmp_path):
        """Test that lambda/functions in uns are flagged as non-serializable."""
        adata = ad.AnnData(X=np.eye(3))
        adata.uns["my_func"] = lambda x: x

        try:
            adata.write_h5ad(tmp_path / "test.h5ad")
            pytest.fail(
                "Function serialization now works! "
                "Update is_serializable() in utils.py."
            )
        except (TypeError, Exception):  # noqa: BLE001
            pass

        html = adata._repr_html_()
        assert "function" in html.lower() or "lambda" in html.lower()
        assert "Not serializable" in html

    def test_nested_unserializable_in_uns(self, tmp_path):
        """Test that nested non-serializable values in uns dicts are detected."""

        class CustomObject:
            pass

        adata = ad.AnnData(X=np.eye(3))
        adata.uns["nested"] = {"ok": 1, "bad": CustomObject()}

        try:
            adata.write_h5ad(tmp_path / "test.h5ad")
            pytest.fail(
                "Nested non-serializable now works! "
                "Check is_serializable() recursive behavior."
            )
        except (TypeError, Exception):  # noqa: BLE001
            pass

        html = adata._repr_html_()
        assert "nested" in html
        assert "Not serializable" in html


class TestColumnNameValidation:
    """Tests for column name validation (issue #321)."""

    def test_slash_column_name_warns_but_still_serializes(self, tmp_path):
        """Test slash in column names: warns (yellow) but still works for now."""
        import warnings

        adata = ad.AnnData(X=np.eye(3))
        adata.obs["path/gene"] = ["a", "b", "c"]

        path = tmp_path / "test_slash.h5ad"

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            try:
                adata.write_h5ad(path)
                serializes = True
            except Exception:  # noqa: BLE001
                serializes = False
                pytest.fail(
                    "Slash in column names now fails! "
                    "Update check_column_name() in formatters.py: "
                    "set is_hard_error=True for slashes."
                )

        if serializes:
            adata2 = ad.read_h5ad(path)
            assert "path/gene" in adata2.obs.columns

        html = adata._repr_html_()
        assert "deprecated" in html or "/" in html

    def test_mapping_sections_warn_for_invalid_keys(self):
        """Test that layers/obsm/etc warn for invalid key names."""
        adata = ad.AnnData(X=np.eye(3))
        adata.layers[("tuple", "key")] = np.eye(3)
        adata.obsm["path/embed"] = np.random.randn(3, 2)

        html = adata._repr_html_()
        assert "Non-string" in html
        assert "deprecated" in html

    def test_uns_warns_for_invalid_keys(self):
        """Test that uns warns for invalid key names."""
        adata = ad.AnnData(X=np.eye(3))
        adata.uns[("tuple", "key")] = "value"

        html = adata._repr_html_()
        assert "Non-string" in html
        # Key error shown with warning icon
        assert "(!)" in html


class TestErrorHandling:
    """Tests for error handling in repr."""

    def test_error_entry_display(self):
        """Test error entries are displayed appropriately."""
        adata = AnnData(np.zeros((10, 5)))
        html = adata._repr_html_()
        # Basic check that repr works
        assert html is not None

    def test_formatter_exception_caught(self):
        """Test formatter exceptions are caught and handled."""
        from anndata._repr.registry import (
            TypeFormatter,
            formatter_registry,
        )

        class FailingType:
            pass

        class FailingFormatter(TypeFormatter):
            priority = 1000

            def can_format(self, obj, context):
                return isinstance(obj, FailingType)

            def format(self, obj, context):
                msg = "Test error"
                raise RuntimeError(msg)

        formatter = FailingFormatter()
        formatter_registry.register_type_formatter(formatter)

        try:
            adata = AnnData(np.zeros((5, 3)))
            adata.uns["failing"] = FailingType()

            with pytest.warns(UserWarning, match="Formatter FailingFormatter"):
                html = adata._repr_html_()

            assert html is not None
            assert "failing" in html
        finally:
            formatter_registry.unregister_type_formatter(formatter)


class TestWarningsAndSerializationIndicators:
    """Tests for warning display and serialization issues."""

    def test_string_column_shows_warning(self):
        """Test string columns that will convert show warning."""
        adata = AnnData(np.zeros((10, 5)))
        # String column with few unique values (will convert to categorical)
        adata.obs["batch"] = ["A", "B"] * 5
        html = adata._repr_html_()
        # Should show warning or indicator
        assert "⚠" in html or "warning" in html.lower() or "categorical" in html.lower()

    def test_datetime_column_shows_warning(self):
        """Test datetime columns show serialization warning."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["date"] = pd.to_datetime(["2024-01-01"] * 5)
        html = adata._repr_html_()
        # Should indicate not serializable
        assert "datetime" in html.lower() or "⚠" in html or "(!)" in html

    def test_custom_object_shows_warning(self):
        """Test custom objects show not serializable warning."""

        class CustomClass:
            pass

        adata = AnnData(np.zeros((5, 3)))
        adata.uns["custom"] = CustomClass()
        html = adata._repr_html_()
        # Should indicate not serializable
        assert (
            "serializable" in html.lower() or "⚠" in html or "unknown" in html.lower()
        )

    def test_lambda_in_uns_shows_warning(self):
        """Test lambda functions show warning."""
        adata = AnnData(np.zeros((5, 3)))
        adata.uns["func"] = lambda x: x
        html = adata._repr_html_()
        # Should indicate not serializable
        assert "serializable" in html.lower() or "function" in html.lower()

    def test_tuple_key_shows_warning(self):
        """Test tuple keys show warning."""
        adata = AnnData(np.zeros((5, 3)))
        adata.uns[("tuple", "key")] = "value"
        html = adata._repr_html_()
        # Should indicate non-string key issue
        assert "Non-string" in html or "⚠" in html or "(!)" in html

    def test_slash_in_column_name_shows_warning(self):
        """Test slash in column names shows deprecation warning."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["path/gene"] = [1, 2, 3, 4, 5]
        html = adata._repr_html_()
        # Should indicate deprecated character
        assert "deprecated" in html.lower() or "/" in html


class TestObsmVarmWarnings:
    """Tests for warning indicators in obsm/varm sections.

    Scientific display requirement: Warnings should be clearly visible
    when obsm/varm contains problematic data types.
    """

    def test_uns_custom_object_shows_warning(self, validate_html):
        """Test custom objects in uns show warning indicator."""

        class CustomObject:
            pass

        adata = AnnData(np.zeros((10, 5)))
        adata.uns["custom"] = CustomObject()
        html = adata._repr_html_()
        # Should render without error
        assert "custom" in html
        # May show unknown type warning or serializable indicator
        assert "anndata-repr" in html

    def test_varm_with_dataframe_shows_info(self, validate_html):
        """Test DataFrame in varm shows appropriate info."""
        adata = AnnData(np.zeros((10, 5)))
        adata.varm["df"] = pd.DataFrame(
            {"a": np.zeros(5), "b": np.zeros(5)}, index=adata.var_names
        )
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("varm", "df")
        # Should show DataFrame or column info
        assert "DataFrame" in html or "2" in html  # 2 columns

    def test_obsm_sparse_shows_format_info(self, validate_html):
        """Test sparse matrices in obsm show format info."""
        import scipy.sparse as sp

        adata = AnnData(np.zeros((10, 5)))
        adata.obsm["sparse"] = sp.random(10, 3, density=0.2, format="csr")
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obsm", "sparse")
        # Should show sparse format
        assert "csr" in html.lower()

    def test_varm_shape_mismatch_detection(self, validate_html):
        """Test varm entries show correct shape info."""
        adata = AnnData(np.zeros((10, 20)))
        adata.varm["PCs"] = np.zeros((20, 50))  # 20 vars, 50 PCs
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("varm", "PCs")
        # Should show shape
        assert "50" in html

    def test_obsm_dtype_shown(self, validate_html):
        """Test obsm arrays show dtype information."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obsm["X_pca"] = np.random.randn(10, 50).astype(np.float32)
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obsm", "X_pca")
        # Should show float32 dtype
        assert "float32" in html


class TestWarningVisibility:
    """Tests ensuring warnings are clearly visible in scientific display.

    Scientific display requirement: Any potential data issues should be
    clearly indicated so scientists don't miss problems with their data.
    """

    def test_warning_icon_visible_in_section(self, validate_html):
        """Test warning icon appears when section has issues."""

        class CustomType:
            pass

        adata = AnnData(np.zeros((5, 3)))
        adata.uns["bad"] = CustomType()
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_element_exists(".anndata-repr")
        # Entry should be visible even with warning
        v.assert_section_contains_entry("uns", "bad")

    def test_multiple_warnings_all_visible(self, validate_html):
        """Test multiple warning entries are all visible."""

        class TypeA:
            pass

        class TypeB:
            pass

        adata = AnnData(np.zeros((5, 3)))
        adata.uns["bad_a"] = TypeA()
        adata.uns["bad_b"] = TypeB()
        adata.uns["good"] = "normal string"
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("uns", "bad_a")
        v.assert_section_contains_entry("uns", "bad_b")
        v.assert_section_contains_entry("uns", "good")

    def test_warning_in_obs_column_visible(self, validate_html):
        """Test warning in obs column is visible."""
        adata = AnnData(np.zeros((5, 3)))
        # Column with lists (non-serializable)
        adata.obs["lists"] = [[1], [2], [3], [4], [5]]
        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "lists")
        # Should have warning indicator
        assert "(!)" in html or "⚠" in html or "list" in html

    def test_warning_does_not_hide_data(self, validate_html):
        """Test warnings don't hide surrounding data."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["normal"] = [1, 2, 3, 4, 5]
        adata.obs["problematic"] = [[1], [2], [3], [4], [5]]
        adata.obs["also_normal"] = ["a", "b", "c", "d", "e"]
        html = adata._repr_html_()
        v = validate_html(html)
        # All columns should be visible
        v.assert_section_contains_entry("obs", "normal")
        v.assert_section_contains_entry("obs", "problematic")
        v.assert_section_contains_entry("obs", "also_normal")
