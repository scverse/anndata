"""
Lazy loading tests for the _repr module.

Tests for module lazy loading, lazy categorical handling,
and lazy AnnData representation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import anndata as ad
from anndata import AnnData

from .conftest import HAS_XARRAY


class TestModuleLazyLoading:
    """Tests to ensure _repr module doesn't load on import anndata."""

    def test_repr_module_not_loaded_on_import(self):
        """Verify that importing anndata doesn't load the full _repr module."""
        import subprocess
        import sys

        code = """
import sys
import anndata
repr_modules = [m for m in sys.modules if 'anndata._repr' in m and m != 'anndata._repr']
repr_modules = [m for m in repr_modules if '_repr_constants' not in m]
if repr_modules:
    print(f"FAIL: These _repr modules were loaded on import: {repr_modules}")
    sys.exit(1)
else:
    print("OK: _repr module not loaded on import")
    sys.exit(0)
"""
        result = subprocess.run(
            [sys.executable, "-c", code],
            check=False,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, (
            f"Lazy loading failed: {result.stdout}\n{result.stderr}"
        )

    def test_repr_module_loads_on_repr_html_call(self):
        """Verify that _repr module loads when _repr_html_() is called."""
        import subprocess
        import sys

        code = """
import sys
import anndata as ad
import numpy as np

adata = ad.AnnData(np.eye(3))
_ = adata._repr_html_()

repr_modules = [m for m in sys.modules if 'anndata._repr.' in m]
if not repr_modules:
    print("FAIL: _repr modules not loaded after _repr_html_()")
    sys.exit(1)
else:
    print(f"OK: _repr modules loaded: {len(repr_modules)} submodules")
    sys.exit(0)
"""
        result = subprocess.run(
            [sys.executable, "-c", code],
            check=False,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, (
            f"Module loading failed: {result.stdout}\n{result.stderr}"
        )


class TestLazyCategoryLoading:
    """Tests for lazy category loading in HTML repr."""

    def test_get_lazy_category_count(self):
        """Test get_lazy_category_count returns None for non-lazy columns."""
        from anndata._repr.lazy import get_lazy_category_count

        series = pd.Series(pd.Categorical(["a", "b", "c"]))
        result = get_lazy_category_count(series)
        assert result is None

        class MockCol:
            pass

        assert get_lazy_category_count(MockCol()) is None

        non_cat = pd.Series([1, 2, 3])
        assert get_lazy_category_count(non_cat) is None

    def test_get_lazy_categories_max_zero_skips(self):
        """Test that max_lazy_categories=0 skips loading entirely."""
        from anndata._repr.lazy import get_lazy_categories
        from anndata._repr.registry import FormatterContext

        context = FormatterContext(max_lazy_categories=0)

        class MockCol:
            pass

        categories, skipped, n_cats = get_lazy_categories(MockCol(), context)
        assert categories == []
        assert skipped is True
        assert n_cats is None

    def test_get_categories_for_display_non_lazy(self):
        """Test get_categories_for_display with regular (non-lazy) categorical."""
        from anndata._repr.registry import FormatterContext
        from anndata._repr.utils import get_categories_for_display

        context = FormatterContext()
        series = pd.Series(pd.Categorical(["a", "b", "c", "a"]))

        categories, skipped, n_cats = get_categories_for_display(
            series, context, is_lazy=False
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
        """Test that get_lazy_category_count reads from storage metadata only."""
        from anndata._repr.lazy import get_lazy_category_count

        adata = AnnData(
            sp.random(100, 50, density=0.1, format="csr", dtype=np.float32),
            obs=pd.DataFrame({"cat_col": pd.Categorical(["a", "b", "c"] * 33 + ["a"])}),
        )
        path = tmp_path / "test.zarr"
        adata.write_zarr(path)

        lazy = ad.experimental.read_lazy(path)
        col = lazy.obs._ds["cat_col"]

        cat_arr = col.variable._data.array
        assert "categories" not in cat_arr.__dict__

        n_cats = get_lazy_category_count(col)
        assert n_cats == 3

        assert "categories" not in cat_arr.__dict__

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_get_lazy_categories_does_not_load_data(self, tmp_path):
        """Test that get_lazy_categories reads from storage directly."""
        from anndata._repr.lazy import get_lazy_categories
        from anndata._repr.registry import FormatterContext

        adata = AnnData(
            sp.random(100, 50, density=0.1, format="csr", dtype=np.float32),
            obs=pd.DataFrame({"cat_col": pd.Categorical(["x", "y", "z"] * 33 + ["x"])}),
        )
        path = tmp_path / "test.zarr"
        adata.write_zarr(path)

        lazy = ad.experimental.read_lazy(path)
        col = lazy.obs._ds["cat_col"]

        cat_arr = col.variable._data.array
        assert "categories" not in cat_arr.__dict__

        context = FormatterContext(max_lazy_categories=100)
        categories, skipped, n_cats = get_lazy_categories(col, context)

        assert set(categories) == {"x", "y", "z"}
        assert not skipped
        assert n_cats == 3

        assert "categories" not in cat_arr.__dict__

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_get_lazy_categories_skipping_does_not_load_categories(self, tmp_path):
        """Test that when skipping (too many cats), we don't load category values."""
        from anndata._repr.lazy import get_lazy_categories
        from anndata._repr.registry import FormatterContext

        large_cats = [f"cat_{i}" for i in range(150)]
        adata = AnnData(
            sp.random(150, 50, density=0.1, format="csr", dtype=np.float32),
            obs=pd.DataFrame({"big_cat": pd.Categorical(large_cats)}),
        )
        path = tmp_path / "test.zarr"
        adata.write_zarr(path)

        lazy = ad.experimental.read_lazy(path)
        col = lazy.obs._ds["big_cat"]

        cat_arr = col.variable._data.array
        assert "categories" not in cat_arr.__dict__

        context = FormatterContext(max_lazy_categories=100)
        categories, truncated, n_cats = get_lazy_categories(col, context)

        assert len(categories) == 100
        assert categories[0] == "cat_0"
        assert truncated is True
        assert n_cats == 150

        assert "categories" not in cat_arr.__dict__

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_get_lazy_categories_h5ad(self, tmp_path):
        """Test get_lazy_categories works with H5AD files."""
        import h5py

        from anndata._repr.lazy import get_lazy_categories, get_lazy_category_count
        from anndata._repr.registry import FormatterContext

        adata = AnnData(np.random.randn(100, 50).astype(np.float32))
        adata.obs["cat_col"] = pd.Categorical(["x", "y", "z"] * 33 + ["x"])

        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        with h5py.File(path, "r") as f:
            lazy = ad.experimental.read_lazy(f)
            col = lazy.obs._ds["cat_col"]

            n_cats = get_lazy_category_count(col)
            assert n_cats == 3

            context = FormatterContext(max_lazy_categories=100)
            categories, skipped, n = get_lazy_categories(col, context)

            assert set(categories) == {"x", "y", "z"}
            assert not skipped
            assert n == 3

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_repr_html_does_not_load_lazy_categorical_data(self, tmp_path):
        """Test that generating HTML repr doesn't trigger loading of lazy categorical data."""
        from anndata._repr import DEFAULT_MAX_LAZY_CATEGORIES

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

        small_cat_arr = lazy.obs._ds["small_cat"].variable._data.array
        large_cat_arr = lazy.obs._ds["large_cat"].variable._data.array

        assert "categories" not in small_cat_arr.__dict__
        assert "categories" not in large_cat_arr.__dict__

        html = lazy._repr_html_()

        assert "small_cat" in html
        assert "large_cat" in html
        assert "category" in html

        assert "categories" not in small_cat_arr.__dict__
        assert "categories" not in large_cat_arr.__dict__

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_lazy_categorical_repr_integration(self, tmp_path):
        """Integration test: verify lazy categoricals display correctly in repr."""
        import h5py

        from anndata._repr import DEFAULT_MAX_LAZY_CATEGORIES
        from anndata.experimental import read_lazy

        n_large = DEFAULT_MAX_LAZY_CATEGORIES + 20

        adata = AnnData(np.random.randn(100, 50).astype(np.float32))
        adata.obs["small_cat"] = pd.Categorical(["A", "B", "C"] * 33 + ["A"])
        adata.obs["large_cat"] = pd.Categorical(
            [f"cat_{i}" for i in range(100)],
            categories=[f"cat_{i}" for i in range(n_large)],
        )
        adata.uns["small_cat_colors"] = ["#ff0000", "#00ff00", "#0000ff"]

        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        with h5py.File(path, "r") as f:
            lazy_adata = read_lazy(f)
            html = lazy_adata._repr_html_()

            assert "small_cat" in html
            assert "A" in html
            assert "#ff0000" in html or "ff0000" in html

            assert "large_cat" in html
            assert "cat_0" in html
            assert "...+20" in html

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_metadata_only_mode_no_disk_loading(self, tmp_path):
        """Test that max_lazy_categories=0 shows counts without loading category labels."""
        import h5py

        from anndata._repr.html import generate_repr_html
        from anndata.experimental import read_lazy

        adata = AnnData(np.random.randn(50, 20).astype(np.float32))
        adata.obs["cat1"] = pd.Categorical(["A", "B", "C"] * 16 + ["A", "A"])
        adata.obs["cat2"] = pd.Categorical(["X", "Y"] * 25)
        adata.uns["cat1_colors"] = ["#ff0000", "#00ff00", "#0000ff"]

        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        with h5py.File(path, "r") as f:
            lazy_adata = read_lazy(f)

            html = generate_repr_html(lazy_adata, max_lazy_categories=0)

            assert "(3 categories)" in html
            assert "(2 categories)" in html

            assert "#ff0000" not in html
            assert "#00ff00" not in html

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_metadata_only_vs_default_mode(self, tmp_path):
        """Compare metadata-only mode vs default mode output."""
        import h5py

        from anndata._repr.html import generate_repr_html
        from anndata.experimental import read_lazy

        adata = AnnData(np.random.randn(50, 20).astype(np.float32))
        adata.obs["small_cat"] = pd.Categorical(["A", "B"] * 25)
        adata.uns["small_cat_colors"] = ["#ff0000", "#00ff00"]

        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        with h5py.File(path, "r") as f:
            lazy_adata = read_lazy(f)

            html_default = generate_repr_html(lazy_adata)
            assert "#ff0000" in html_default

            html_metadata = generate_repr_html(lazy_adata, max_lazy_categories=0)
            assert "#ff0000" not in html_metadata
            assert "(2 categories)" in html_metadata


class TestIsLazyColumn:
    """Tests for is_lazy_column detection."""

    def test_is_lazy_column_regular_series(self):
        """Test that regular pandas Series is not detected as lazy."""
        from anndata._repr.lazy import is_lazy_column

        series = pd.Series([1, 2, 3])
        assert is_lazy_column(series) is False

    def test_is_lazy_column_categorical_series(self):
        """Test that categorical pandas Series is not detected as lazy."""
        from anndata._repr.lazy import is_lazy_column

        series = pd.Series(pd.Categorical(["a", "b", "c"]))
        assert is_lazy_column(series) is False

    def test_is_lazy_column_numpy_array(self):
        """Test that numpy array is not detected as lazy."""
        from anndata._repr.lazy import is_lazy_column

        arr = np.array([1, 2, 3])
        assert is_lazy_column(arr) is False

    def test_is_lazy_column_mock_xarray_like(self):
        """Test that object with xarray-like attributes is detected as lazy."""
        from anndata._repr.lazy import is_lazy_column

        class MockXarrayColumn:
            variable = "something"
            dims = ("x",)

        assert is_lazy_column(MockXarrayColumn()) is True

    def test_is_lazy_column_mock_variable_backed(self):
        """Test that object with _variable attribute is detected as lazy."""
        from anndata._repr.lazy import is_lazy_column

        class MockVariableBacked:
            _variable = "something"

        assert is_lazy_column(MockVariableBacked()) is True

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_is_lazy_column_real_xarray(self, tmp_path):
        """Test is_lazy_column with real xarray DataArray from lazy AnnData."""
        from anndata._repr.lazy import is_lazy_column

        adata = AnnData(np.random.randn(50, 20).astype(np.float32))
        adata.obs["cat"] = pd.Categorical(["A", "B"] * 25)

        path = tmp_path / "test.zarr"
        adata.write_zarr(path)

        lazy = ad.experimental.read_lazy(path)
        col = lazy.obs._ds["cat"]

        assert is_lazy_column(col) is True


class TestGetLazyCategoricalInfo:
    """Tests for get_lazy_categorical_info function."""

    def test_get_lazy_categorical_info_non_lazy(self):
        """Test that non-lazy objects return (None, False)."""
        from anndata._repr.lazy import get_lazy_categorical_info

        series = pd.Series(pd.Categorical(["a", "b", "c"]))
        n_cats, ordered = get_lazy_categorical_info(series)
        assert n_cats is None
        assert ordered is False

    def test_get_lazy_categorical_info_plain_object(self):
        """Test that plain objects return (None, False)."""
        from anndata._repr.lazy import get_lazy_categorical_info

        n_cats, ordered = get_lazy_categorical_info("not a categorical")
        assert n_cats is None
        assert ordered is False

    def test_get_lazy_categorical_info_mock_without_categorical(self):
        """Test object with xarray structure but no CategoricalArray."""
        from anndata._repr.lazy import get_lazy_categorical_info

        class MockVariable:
            _data = None

        class MockCol:
            variable = MockVariable()

        n_cats, ordered = get_lazy_categorical_info(MockCol())
        assert n_cats is None
        assert ordered is False

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_get_lazy_categorical_info_zarr(self, tmp_path):
        """Test get_lazy_categorical_info with Zarr-backed categorical."""
        from anndata._repr.lazy import get_lazy_categorical_info

        adata = AnnData(sp.random(100, 50, density=0.1, format="csr", dtype=np.float32))
        adata.obs["cat"] = pd.Categorical(["a", "b", "c", "d", "e"] * 20)

        path = tmp_path / "test.zarr"
        adata.write_zarr(path)

        lazy = ad.experimental.read_lazy(path)
        col = lazy.obs._ds["cat"]

        n_cats, ordered = get_lazy_categorical_info(col)
        assert n_cats == 5
        assert not ordered

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_get_lazy_categorical_info_h5ad(self, tmp_path):
        """Test get_lazy_categorical_info with H5AD-backed categorical."""
        import h5py

        from anndata._repr.lazy import get_lazy_categorical_info

        adata = AnnData(np.random.randn(100, 50).astype(np.float32))
        adata.obs["cat"] = pd.Categorical(["x", "y", "z"] * 33 + ["x"], ordered=False)

        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        with h5py.File(path, "r") as f:
            lazy = ad.experimental.read_lazy(f)
            col = lazy.obs._ds["cat"]

            n_cats, ordered = get_lazy_categorical_info(col)
            assert n_cats == 3
            assert not ordered


class TestGetCategoricalArrayHelper:
    """Tests for _get_categorical_array helper function."""

    def test_get_categorical_array_non_lazy(self):
        """Test that non-lazy objects return None."""
        from anndata._repr.lazy import _get_categorical_array

        series = pd.Series(pd.Categorical(["a", "b", "c"]))
        assert _get_categorical_array(series) is None

    def test_get_categorical_array_plain_object(self):
        """Test that plain objects return None."""
        from anndata._repr.lazy import _get_categorical_array

        assert _get_categorical_array("string") is None
        assert _get_categorical_array(123) is None
        assert _get_categorical_array(None) is None

    def test_get_categorical_array_mock_structure_no_categorical(self):
        """Test object with partial xarray structure returns None."""
        from anndata._repr.lazy import _get_categorical_array

        class MockLazyIndexed:
            array = "not a CategoricalArray"

        class MockVariable:
            _data = MockLazyIndexed()

        class MockCol:
            variable = MockVariable()

        assert _get_categorical_array(MockCol()) is None

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_get_categorical_array_real(self, tmp_path):
        """Test _get_categorical_array with real lazy categorical."""
        from anndata._repr.lazy import _get_categorical_array
        from anndata.experimental.backed._lazy_arrays import CategoricalArray

        adata = AnnData(sp.random(50, 20, density=0.1, format="csr", dtype=np.float32))
        adata.obs["cat"] = pd.Categorical(["A", "B", "C"] * 16 + ["A", "A"])

        path = tmp_path / "test.zarr"
        adata.write_zarr(path)

        lazy = ad.experimental.read_lazy(path)
        col = lazy.obs._ds["cat"]

        result = _get_categorical_array(col)
        assert isinstance(result, CategoricalArray)


class TestLazyAdataDetection:
    """Tests for is_lazy_adata detection with edge cases."""

    def test_is_lazy_adata_none(self):
        """Test that None returns False."""
        from anndata._repr.lazy import is_lazy_adata

        assert is_lazy_adata(None) is False

    def test_is_lazy_adata_no_obs(self):
        """Test object without obs attribute returns False."""
        from anndata._repr.lazy import is_lazy_adata

        class NoObs:
            pass

        assert is_lazy_adata(NoObs()) is False

    def test_is_lazy_adata_obs_raises(self):
        """Test object where .obs raises returns False."""
        from anndata._repr.lazy import is_lazy_adata

        class RaisingObs:
            @property
            def obs(self):
                msg = "Cannot access obs"
                raise RuntimeError(msg)

        assert is_lazy_adata(RaisingObs()) is False

    def test_is_lazy_adata_obs_none(self):
        """Test object with obs=None returns False."""
        from anndata._repr.lazy import is_lazy_adata

        class NoneObs:
            obs = None

        assert is_lazy_adata(NoneObs()) is False


class TestLazyBackingInfo:
    """Tests for lazy AnnData backing file info extraction."""

    def test_get_lazy_backing_info_non_lazy_returns_empty(self):
        """Test that non-lazy AnnData returns empty backing info."""
        from anndata._repr.lazy import get_lazy_backing_info

        adata = AnnData(np.random.randn(10, 5).astype(np.float32))
        info = get_lazy_backing_info(adata)

        assert info == {"filename": "", "format": ""}

    def test_is_lazy_adata_false_for_regular(self):
        """Test that is_lazy_adata returns False for regular AnnData."""
        from anndata._repr.lazy import is_lazy_adata

        adata = AnnData(np.random.randn(10, 5).astype(np.float32))
        assert is_lazy_adata(adata) is False

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_get_lazy_backing_info_h5ad(self, tmp_path):
        """Test backing info extraction from lazy H5AD."""
        import h5py

        from anndata._repr.lazy import get_lazy_backing_info, is_lazy_adata
        from anndata.experimental import read_lazy

        adata = AnnData(np.random.randn(50, 20).astype(np.float32))
        adata.obs["cat"] = pd.Categorical(["A", "B", "C"] * 16 + ["A", "A"])

        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        with h5py.File(path, "r") as f:
            lazy_adata = read_lazy(f)

            assert is_lazy_adata(lazy_adata) is True

            info = get_lazy_backing_info(lazy_adata)
            assert info["format"] == "H5AD"
            assert str(path) in info["filename"] or info["filename"] != ""

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_get_lazy_backing_info_zarr(self, tmp_path):
        """Test backing info extraction from lazy Zarr."""
        from anndata._repr.lazy import get_lazy_backing_info, is_lazy_adata
        from anndata.experimental import read_lazy

        adata = AnnData(sp.random(50, 20, density=0.1, format="csr", dtype=np.float32))
        adata.obs["cat"] = pd.Categorical(["A", "B", "C"] * 16 + ["A", "A"])

        path = tmp_path / "test.zarr"
        adata.write_zarr(path)

        lazy_adata = read_lazy(path)

        assert is_lazy_adata(lazy_adata) is True

        info = get_lazy_backing_info(lazy_adata)
        assert info["format"] == "Zarr"
        # Zarr path should be extracted
        assert str(path) in info["filename"] or info["filename"] != ""

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_lazy_badge_h5ad_in_html(self, tmp_path, validate_html):
        """Test that lazy H5AD shows correct badge in HTML repr."""
        import h5py

        from anndata.experimental import read_lazy

        adata = AnnData(np.random.randn(50, 20).astype(np.float32))
        adata.obs["cat"] = pd.Categorical(["A", "B"] * 25)

        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        with h5py.File(path, "r") as f:
            lazy_adata = read_lazy(f)
            html = lazy_adata._repr_html_()

            v = validate_html(html)
            # Should show "Lazy (H5AD)" badge with lazy badge styling
            v.assert_badge_shown("lazy")
            v.assert_text_visible("Lazy")
            v.assert_text_visible("H5AD")

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_lazy_badge_zarr_in_html(self, tmp_path, validate_html):
        """Test that lazy Zarr shows correct badge in HTML repr."""
        from anndata.experimental import read_lazy

        adata = AnnData(sp.random(50, 20, density=0.1, format="csr", dtype=np.float32))
        adata.obs["cat"] = pd.Categorical(["A", "B"] * 25)

        path = tmp_path / "test.zarr"
        adata.write_zarr(path)

        lazy_adata = read_lazy(path)
        html = lazy_adata._repr_html_()

        v = validate_html(html)
        # Should show "Lazy (Zarr)" badge with lazy badge styling
        v.assert_badge_shown("lazy")
        v.assert_text_visible("Lazy")
        v.assert_text_visible("Zarr")

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_lazy_file_path_in_html_h5ad(self, tmp_path, validate_html):
        """Test that lazy H5AD file path appears in HTML repr."""
        import h5py

        from anndata.experimental import read_lazy

        adata = AnnData(np.random.randn(50, 20).astype(np.float32))
        adata.obs["cat"] = pd.Categorical(["A", "B"] * 25)

        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        with h5py.File(path, "r") as f:
            lazy_adata = read_lazy(f)
            html = lazy_adata._repr_html_()

            v = validate_html(html)
            # File path should appear in visible content
            v.assert_text_visible(str(path))
            # Also verify it's in the file path element
            v.assert_element_exists(".adata-file-path")

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_lazy_file_path_in_html_zarr(self, tmp_path, validate_html):
        """Test that lazy Zarr file path appears in HTML repr."""
        from anndata.experimental import read_lazy

        adata = AnnData(sp.random(50, 20, density=0.1, format="csr", dtype=np.float32))
        adata.obs["cat"] = pd.Categorical(["A", "B"] * 25)

        path = tmp_path / "test.zarr"
        adata.write_zarr(path)

        lazy_adata = read_lazy(path)
        html = lazy_adata._repr_html_()

        v = validate_html(html)
        # File path should appear in visible content
        v.assert_text_visible(str(path))
        # Also verify it's in the file path element
        v.assert_element_exists(".adata-file-path")
