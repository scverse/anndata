"""
Utility function tests for the _repr module.

Tests for serialization checks, color detection, formatting helpers,
preview functions, and other utilities.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from anndata import AnnData


class TestSerializability:
    """Tests for serialization detection utilities."""

    def test_is_serializable_basic_types(self):
        """Test serialization detection for basic types."""
        from anndata._repr.utils import is_serializable

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


class TestStringWarnings:
    """Tests for string-to-category warning detection."""

    def test_should_warn_string_column(self):
        """Test string-to-category warning detection."""
        from anndata._repr.utils import should_warn_string_column

        s = pd.Series(["A", "B", "A", "C", "B", "A"])
        warn, msg = should_warn_string_column(s, s.nunique())
        assert warn
        assert "3" in msg

        s = pd.Series(["A", "B", "C", "D", "E"])
        warn, msg = should_warn_string_column(s, s.nunique())
        assert not warn

        s = pd.Series([1, 2, 1, 2, 1])
        warn, msg = should_warn_string_column(s, s.nunique())
        assert not warn

        s = pd.Series(["A", "B", "A"])
        warn, msg = should_warn_string_column(s, None)
        assert not warn

    def test_should_warn_string_column_with_none_nunique(self):
        """Test should_warn_string_column handles None n_unique gracefully."""
        from anndata._repr.utils import should_warn_string_column

        s = pd.Series(["A", "B", "A", "C"])
        warn, _msg = should_warn_string_column(s, None)
        assert not warn


class TestColorDetection:
    """Tests for color list detection."""

    def test_is_color_list(self):
        """Test color list detection."""
        from anndata._repr.utils import is_color_list

        assert is_color_list("cluster_colors", ["#FF0000", "#00FF00", "#0000FF"])
        assert is_color_list("leiden_colors", np.array(["#123456", "#ABCDEF"]))
        assert is_color_list("cluster_colors", [])

        assert not is_color_list("cluster", ["#FF0000"])
        assert not is_color_list("colors", ["#FF0000"])
        assert not is_color_list("cluster_colors", "#FF0000")
        assert not is_color_list("cluster_colors", ["not_a_color", "also_not"])

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

    def test_get_matching_column_colors_no_uns_key(self):
        """Test get_matching_column_colors returns None when no colors in uns."""
        from anndata._repr.utils import get_matching_column_colors

        adata = AnnData(np.zeros((10, 5)))
        adata.obs["cell_type"] = pd.Categorical(["A", "B"] * 5)

        colors = get_matching_column_colors(adata, "cell_type")
        assert colors is None


class TestFormatting:
    """Tests for formatting utilities."""

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

    def test_format_memory_size_negative(self):
        """Test format_memory_size handles negative values."""
        from anndata._repr.utils import format_memory_size

        result = format_memory_size(-100)
        assert result == "Unknown"

    def test_format_memory_size_very_large(self):
        """Test format_memory_size handles very large values (PB)."""
        from anndata._repr.utils import format_memory_size

        result = format_memory_size(5 * 1024**5)
        assert "PB" in result

    def test_format_number(self):
        """Test number formatting with thousands separators."""
        from anndata._repr.utils import format_number

        assert format_number(1000) == "1,000"
        assert format_number(1000000) == "1,000,000"

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


class TestBackingInfo:
    """Tests for backing info detection."""

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


class TestColumnNameValidation:
    """Tests for column name validation."""

    def test_check_column_name_valid(self):
        """Test that valid column names pass."""
        from anndata._repr.utils import check_column_name

        valid, _, _ = check_column_name("gene_name")
        assert valid
        valid, _, _ = check_column_name("gène_名前")
        assert valid

    def test_check_column_name_slash(self):
        """Test slashes are flagged as warning."""
        from anndata._repr.utils import check_column_name

        valid, reason, is_hard_error = check_column_name("path/to/gene")
        assert not valid
        assert "/" in reason
        assert not is_hard_error

    def test_check_column_name_non_string(self):
        """Test non-string names are flagged as hard error."""
        from anndata._repr.utils import check_column_name

        valid, reason, is_hard_error = check_column_name(("a", "b"))
        assert not valid
        assert "Non-string" in reason
        assert is_hard_error


class TestValuePreviewFunctions:
    """Tests for value preview helper functions."""

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
        assert len(result) < 30
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
        assert "3.14159" in result

    def test_preview_dict_empty(self):
        """Test dict preview for empty dict."""
        from anndata._repr.utils import preview_dict

        assert preview_dict({}) == "{}"

    def test_preview_dict_small(self):
        """Test dict preview for small dict."""
        from anndata._repr.utils import preview_dict

        result = preview_dict({"a": 1, "b": 2})
        assert "a" in result
        assert "b" in result
        assert "{" in result

    def test_preview_dict_large(self):
        """Test dict preview for large dict."""
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
        assert "1" in result
        assert "2" in result
        assert "3" in result

    def test_preview_item_string(self):
        """Test item preview for strings."""
        from anndata._repr.utils import preview_item

        result = preview_item("hello")
        assert '"hello"' in result

    def test_preview_item_numbers(self):
        """Test item preview for numbers."""
        from anndata._repr.utils import preview_item

        assert "42" in preview_item(42)
        assert "3.14" in preview_item(3.14)

    def test_preview_item_none(self):
        """Test item preview for None."""
        from anndata._repr.utils import preview_item

        assert "None" in preview_item(None)

    def test_generate_value_preview_none(self):
        """Test generate_value_preview for None."""
        from anndata._repr.utils import generate_value_preview

        result = generate_value_preview(None)
        assert "None" in result

    def test_generate_value_preview_string(self):
        """Test generate_value_preview for strings."""
        from anndata._repr.utils import generate_value_preview

        result = generate_value_preview("test")
        assert "test" in result

    def test_generate_value_preview_dict(self):
        """Test generate_value_preview for dicts."""
        from anndata._repr.utils import generate_value_preview

        result = generate_value_preview({"a": 1})
        assert "a" in result or "1" in result

    def test_generate_value_preview_list(self):
        """Test generate_value_preview for lists."""
        from anndata._repr.utils import generate_value_preview

        result = generate_value_preview([1, 2, 3])
        assert "1" in result

    def test_generate_value_preview_complex(self):
        """Test value preview for complex types returns empty."""
        from anndata._repr.utils import generate_value_preview

        class CustomClass:
            pass

        assert generate_value_preview(CustomClass()) == ""

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

    def test_preview_item_complex(self):
        """Test item preview for complex types returns empty."""
        from anndata._repr.utils import preview_item

        assert preview_item([1, 2, 3]) == ""
        assert preview_item({"a": 1}) == ""

    def test_preview_sequence_small_tuple(self):
        """Test sequence preview for small tuple."""
        from anndata._repr.utils import preview_sequence

        result = preview_sequence((1, 2, 3))
        assert "(" in result
        assert ")" in result
