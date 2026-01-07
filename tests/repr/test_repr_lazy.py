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
        """Test _get_lazy_category_count returns None for non-lazy columns."""
        from anndata._repr.utils import _get_lazy_category_count

        series = pd.Series(pd.Categorical(["a", "b", "c"]))
        result = _get_lazy_category_count(series)
        assert result is None

        class MockCol:
            pass

        assert _get_lazy_category_count(MockCol()) is None

        non_cat = pd.Series([1, 2, 3])
        assert _get_lazy_category_count(non_cat) is None

    def test_get_lazy_categories_max_zero_skips(self):
        """Test that max_lazy_categories=0 skips loading entirely."""
        from anndata._repr.registry import FormatterContext
        from anndata._repr.utils import _get_lazy_categories

        context = FormatterContext(max_lazy_categories=0)

        class MockCol:
            pass

        categories, skipped, n_cats = _get_lazy_categories(MockCol(), context)
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
        """Test that _get_lazy_category_count reads from storage metadata only."""
        from anndata._repr.utils import _get_lazy_category_count

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

        n_cats = _get_lazy_category_count(col)
        assert n_cats == 3

        assert "categories" not in cat_arr.__dict__

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_get_lazy_categories_does_not_load_data(self, tmp_path):
        """Test that _get_lazy_categories reads from storage directly."""
        from anndata._repr.registry import FormatterContext
        from anndata._repr.utils import _get_lazy_categories

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
        categories, skipped, n_cats = _get_lazy_categories(col, context)

        assert set(categories) == {"x", "y", "z"}
        assert not skipped
        assert n_cats == 3

        assert "categories" not in cat_arr.__dict__

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_get_lazy_categories_skipping_does_not_load_categories(self, tmp_path):
        """Test that when skipping (too many cats), we don't load category values."""
        from anndata._repr.registry import FormatterContext
        from anndata._repr.utils import _get_lazy_categories

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
        categories, truncated, n_cats = _get_lazy_categories(col, context)

        assert len(categories) == 100
        assert categories[0] == "cat_0"
        assert truncated is True
        assert n_cats == 150

        assert "categories" not in cat_arr.__dict__

    @pytest.mark.skipif(not HAS_XARRAY, reason="xarray not installed")
    def test_get_lazy_categories_h5ad(self, tmp_path):
        """Test _get_lazy_categories works with H5AD files."""
        import h5py

        from anndata._repr.registry import FormatterContext
        from anndata._repr.utils import _get_lazy_categories, _get_lazy_category_count

        adata = AnnData(np.random.randn(100, 50).astype(np.float32))
        adata.obs["cat_col"] = pd.Categorical(["x", "y", "z"] * 33 + ["x"])

        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        with h5py.File(path, "r") as f:
            lazy = ad.experimental.read_lazy(f)
            col = lazy.obs._ds["cat_col"]

            n_cats = _get_lazy_category_count(col)
            assert n_cats == 3

            context = FormatterContext(max_lazy_categories=100)
            categories, skipped, n = _get_lazy_categories(col, context)

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
