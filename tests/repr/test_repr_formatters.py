"""
Type formatter tests for the _repr module.

Tests for NumpyArrayFormatter, SparseMatrixFormatter, CategoricalFormatter,
DaskArrayFormatter, CuPyArrayFormatter, AwkwardArrayFormatter, ArrayAPIFormatter,
and all built-in type formatters.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from anndata import AnnData

from .conftest import HAS_AWKWARD, HAS_DASK


class TestNumpyFormatters:
    """Tests for NumPy array formatters."""

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
        assert "2 masked values" in result.tooltip

    def test_masked_array_no_mask(self):
        """Test MaskedArray formatter with no masked values."""
        from anndata._repr.formatters import NumpyMaskedArrayFormatter
        from anndata._repr.registry import FormatterContext

        formatter = NumpyMaskedArrayFormatter()
        arr = np.ma.array([1, 2, 3, 4, 5])

        result = formatter.format(arr, FormatterContext())
        assert result.tooltip == ""  # No masked values means empty tooltip


class TestSparseFormatters:
    """Tests for sparse matrix formatters."""

    def test_sparse_matrix_formatter(self):
        """Test sparse matrix formatting."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        formatter = SparseMatrixFormatter()
        mat = sp.random(1000, 500, density=0.1, format="csr")

        assert formatter.can_format(mat)
        result = formatter.format(mat, FormatterContext())

        assert "csr" in result.type_name.lower()
        assert "sparse" in result.type_name.lower()
        assert "stored" in result.type_name.lower()

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

    def test_sparse_lil_formatter(self):
        """Test sparse formatter with LIL matrix."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        formatter = SparseMatrixFormatter()
        mat = sp.lil_matrix((10, 10))
        mat[0, 0] = 1
        mat[5, 5] = 2

        assert formatter.can_format(mat)
        result = formatter.format(mat, FormatterContext())
        assert "lil" in result.type_name.lower()

    def test_sparse_dok_formatter(self):
        """Test sparse formatter with DOK matrix."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        formatter = SparseMatrixFormatter()
        mat = sp.dok_matrix((10, 10))
        mat[0, 0] = 1
        mat[5, 5] = 2

        assert formatter.can_format(mat)
        result = formatter.format(mat, FormatterContext())
        assert "dok" in result.type_name.lower()

    def test_sparse_dia_formatter(self):
        """Test sparse formatter with DIA matrix."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        formatter = SparseMatrixFormatter()
        data = np.array([[1, 2, 3, 4], [5, 6, 7, 8]])
        offsets = np.array([0, 1])
        mat = sp.dia_matrix((data, offsets), shape=(4, 4))

        assert formatter.can_format(mat)
        result = formatter.format(mat, FormatterContext())
        assert "dia" in result.type_name.lower()

    def test_sparse_bsr_formatter(self):
        """Test sparse formatter with BSR matrix."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        formatter = SparseMatrixFormatter()
        mat = sp.bsr_matrix(
            np.array([[1, 0, 0, 0], [0, 0, 2, 0], [0, 0, 0, 3], [4, 0, 0, 0]])
        )

        assert formatter.can_format(mat)
        result = formatter.format(mat, FormatterContext())
        assert "bsr" in result.type_name.lower()

    def test_sparse_zero_elements(self):
        """Test sparse formatter with zero-element matrix."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        formatter = SparseMatrixFormatter()
        mat = sp.csr_matrix((0, 0))

        result = formatter.format(mat, FormatterContext())
        assert (
            "sparse" not in result.type_name.lower() or "stored" not in result.type_name
        )

    def test_sparse_formatter_duck_typing_fallback(self):
        """Test sparse formatter uses duck typing when scipy checks fail."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        class MockSparseArray:
            def __init__(self):
                self.nnz = 10
                self.shape = (5, 5)
                self.dtype = np.float64

            def tocsr(self):
                pass

        MockSparseArray.__module__ = "scipy.sparse._csr"

        formatter = SparseMatrixFormatter()
        mock_sparse = MockSparseArray()

        assert formatter.can_format(mock_sparse)
        result = formatter.format(mock_sparse, FormatterContext())
        assert "MockSparseArray" in result.type_name


class TestPandasFormatters:
    """Tests for pandas formatters."""

    def test_categorical_formatter(self):
        """Test categorical formatting."""
        from anndata._repr.formatters import CategoricalFormatter
        from anndata._repr.registry import FormatterContext

        formatter = CategoricalFormatter()
        cat_series = pd.Series(pd.Categorical(["A", "B", "C"] * 10))

        assert formatter.can_format(cat_series)
        result = formatter.format(cat_series, FormatterContext())

        assert "category" in result.type_name.lower()
        assert "(3)" in result.type_name

    def test_categorical_direct_object(self):
        """Test CategoricalFormatter with direct pd.Categorical object."""
        from anndata._repr.formatters import CategoricalFormatter
        from anndata._repr.registry import FormatterContext

        formatter = CategoricalFormatter()
        cat = pd.Categorical(["A", "B", "A", "C"])

        assert formatter.can_format(cat)
        result = formatter.format(cat, FormatterContext())
        assert "category" in result.type_name.lower()
        assert "(3)" in result.type_name

    def test_series_formatter_simple(self):
        """Test SeriesFormatter with simple numeric series."""
        from anndata._repr.formatters import SeriesFormatter
        from anndata._repr.registry import FormatterContext

        formatter = SeriesFormatter()
        series = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0])

        assert formatter.can_format(series)
        result = formatter.format(series, FormatterContext())
        assert "float64" in result.type_name

    def test_dataframe_formatter(self):
        """Test DataFrameFormatter."""
        from anndata._repr.formatters import DataFrameFormatter
        from anndata._repr.registry import FormatterContext

        formatter = DataFrameFormatter()
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})

        assert formatter.can_format(df)
        result = formatter.format(df, FormatterContext())
        assert "3 × 2" in result.type_name

        result_obsm = formatter.format(df, FormatterContext(section="obsm"))
        assert result_obsm.preview_html is not None
        assert "a" in result_obsm.preview_html
        assert "b" in result_obsm.preview_html

    def test_dataframe_formatter_expandable(self):
        """Test DataFrameFormatter with expandable to_html enabled."""
        import anndata
        from anndata._repr.formatters import DataFrameFormatter
        from anndata._repr.registry import FormatterContext

        formatter = DataFrameFormatter()
        ctx = FormatterContext()
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})

        result = formatter.format(df, ctx)
        assert result.expanded_html is None

        original = anndata.settings.repr_html_dataframe_expand
        try:
            anndata.settings.repr_html_dataframe_expand = True
            result_expanded = formatter.format(df, ctx)
            assert result_expanded.expanded_html is not None
            assert "<table" in result_expanded.expanded_html
        finally:
            anndata.settings.repr_html_dataframe_expand = original

    def test_dataframe_long_column_names_truncation(self):
        """Test DataFrame with long column names in preview."""
        from anndata._repr.formatters import DataFrameFormatter
        from anndata._repr.registry import FormatterContext

        formatter = DataFrameFormatter()
        long_names = {
            f"very_long_column_name_{i}_with_extra_text": [1] for i in range(3)
        }
        df = pd.DataFrame(long_names)

        result = formatter.format(df, FormatterContext(section="obsm"))
        assert result.preview_html is not None
        assert "anndata-columns" in result.preview_html
        assert "very_long_column_name_0" in result.preview_html


class TestBuiltinFormatters:
    """Tests for built-in type formatters."""

    def test_list_formatter(self):
        """Test list formatter."""
        from anndata._repr.formatters import ListFormatter
        from anndata._repr.registry import FormatterContext

        formatter = ListFormatter()
        obj = [1, 2, 3, 4, 5]

        assert formatter.can_format(obj)
        result = formatter.format(obj, FormatterContext())
        assert "list" in result.type_name.lower()
        assert result.preview is not None

    def test_dict_formatter(self):
        """Test dict formatter."""
        from anndata._repr.formatters import DictFormatter
        from anndata._repr.registry import FormatterContext

        formatter = DictFormatter()
        obj = {"a": 1, "b": 2, "c": 3}

        assert formatter.can_format(obj)
        result = formatter.format(obj, FormatterContext())
        assert "dict" in result.type_name.lower()
        assert result.preview is not None
        assert "a" in result.preview

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


class TestSpecialArrayFormatters:
    """Tests for special array type formatters (Dask, CuPy, Awkward, ArrayAPI)."""

    @pytest.mark.skipif(not HAS_DASK, reason="dask not installed")
    def test_dask_array_display(self):
        """Test Dask array shows chunk info."""
        import dask.array as da

        X = da.zeros((1000, 500), chunks=(100, 500))
        adata = AnnData(X)
        html = adata._repr_html_()
        assert "dask" in html.lower()

    @pytest.mark.skipif(not HAS_AWKWARD, reason="awkward not installed")
    def test_awkward_array_display(self):
        """Test Awkward array displays correctly."""
        import awkward as ak

        adata = AnnData(np.zeros((10, 5)))
        adata.obsm["ragged"] = ak.Array([[1, 2], [3, 4, 5], [6]] * 3 + [[7]])
        html = adata._repr_html_()
        assert "awkward" in html.lower()

    def test_cupy_array_formatter_with_mock(self):
        """Test CuPy array formatter with a mock CuPy array."""
        from anndata._repr.formatters import CuPyArrayFormatter
        from anndata._repr.registry import FormatterContext

        class MockDevice:
            id = 0

        class MockCuPyArray:
            def __init__(self):
                self.shape = (100, 50)
                self.dtype = np.float32
                self.device = MockDevice()

        MockCuPyArray.__module__ = "cupy._core.core"

        formatter = CuPyArrayFormatter()
        mock_arr = MockCuPyArray()

        assert formatter.can_format(mock_arr)
        result = formatter.format(mock_arr, FormatterContext())

        assert "cupy" in result.type_name.lower()
        assert "100" in result.type_name
        assert "50" in result.type_name
        assert result.css_class == "anndata-dtype--gpu"
        assert "GPU:0" in result.tooltip

    def test_awkward_array_formatter_with_mock(self):
        """Test Awkward array formatter with a mock object."""
        from anndata._repr.formatters import AwkwardArrayFormatter
        from anndata._repr.registry import FormatterContext

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
        assert result.css_class == "anndata-dtype--awkward"
        assert "var * int64" in result.tooltip

    def test_awkward_array_formatter_exception_handling(self):
        """Test Awkward array formatter handles exceptions."""
        from anndata._repr.formatters import AwkwardArrayFormatter
        from anndata._repr.registry import FormatterContext

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

        assert "awkward" in result.type_name.lower()
        assert "?" in result.type_name
        assert "unknown" in result.tooltip.lower()


class TestArrayAPIFormatter:
    """Tests for Array-API compatible array formatter."""

    def test_array_api_formatter_jax_like(self):
        """Test Array-API formatter with JAX-like array."""
        from anndata._repr.formatters import ArrayAPIFormatter
        from anndata._repr.registry import FormatterContext

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
        assert result.css_class == "anndata-dtype--array-api"
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
        """Test Array-API formatter with device_buffer attribute."""
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
        assert not formatter.can_format(arr)

    def test_array_api_formatter_excludes_pandas(self):
        """Test Array-API formatter excludes pandas objects."""
        from anndata._repr.formatters import ArrayAPIFormatter

        formatter = ArrayAPIFormatter()
        assert not formatter.can_format(pd.DataFrame({"a": [1, 2, 3]}))
        assert not formatter.can_format(pd.Series([1, 2, 3]))


class TestDtypeCSSClassHelpers:
    """Tests for dtype CSS class helper functions."""

    def test_get_dtype_css_class_complex(self):
        """Test CSS class for complex dtype."""
        from anndata._repr.formatters import _get_dtype_css_class

        complex_dtype = np.dtype(np.complex128)
        css_class = _get_dtype_css_class(complex_dtype)
        assert css_class == "anndata-dtype--float"

    def test_get_dtype_css_class_unknown(self):
        """Test CSS class for unknown dtype."""
        from anndata._repr.formatters import _get_dtype_css_class

        datetime_dtype = np.dtype("datetime64[ns]")
        css_class = _get_dtype_css_class(datetime_dtype)
        assert css_class == "anndata-dtype--object"

    def test_get_dtype_css_class_category(self):
        """Test CSS class for pandas category dtype."""
        from anndata._repr.formatters import _get_dtype_css_class

        cat_dtype = pd.CategoricalDtype(categories=["a", "b", "c"])
        css_class = _get_dtype_css_class(cat_dtype)
        assert css_class == "anndata-dtype--category"

    def test_get_dtype_css_class_timedelta(self):
        """Test CSS class for timedelta dtype."""
        from anndata._repr.formatters import _get_dtype_css_class

        timedelta_dtype = pd.to_timedelta([1, 2, 3], unit="s").dtype
        css_class = _get_dtype_css_class(timedelta_dtype)
        assert css_class == "anndata-dtype--object"


class TestAnnDataFormatter:
    """Tests for AnnData formatter."""

    def test_anndata_formatter_nested_depth_limit(self):
        """Test AnnData formatter respects max_depth for nested objects."""
        from anndata._repr.formatters import AnnDataFormatter
        from anndata._repr.registry import FormatterContext

        formatter = AnnDataFormatter()
        inner = AnnData(np.zeros((5, 3)))

        result_shallow = formatter.format(inner, FormatterContext(depth=0, max_depth=3))
        assert result_shallow.expanded_html is not None

        result_deep = formatter.format(inner, FormatterContext(depth=2, max_depth=3))
        assert result_deep.expanded_html is None


class TestFutureCompatibility:
    """Tests for future compatibility with upcoming array types."""

    def test_sparse_duck_typing_detection(self):
        """Test that sparse arrays can be detected via duck typing."""
        from anndata._repr.formatters import SparseMatrixFormatter
        from anndata._repr.registry import FormatterContext

        sparse_matrix = sp.csr_matrix([[1, 0, 2], [0, 0, 3], [4, 5, 6]])

        formatter = SparseMatrixFormatter()
        assert formatter.can_format(sparse_matrix)

        context = FormatterContext(depth=0)
        result = formatter.format(sparse_matrix, context)
        assert "csr" in result.type_name.lower()
        assert result.css_class == "anndata-dtype--sparse"
        assert "6 stored" in result.type_name

    def test_array_api_formatter_with_mock_jax_array(self):
        """Test ArrayAPIFormatter with a mock JAX-like array."""
        from anndata._repr.formatters import ArrayAPIFormatter
        from anndata._repr.registry import FormatterContext

        class MockJAXArray:
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

        formatter = ArrayAPIFormatter()
        can_format = formatter.can_format(mock_array)
        assert can_format

        context = FormatterContext(depth=0)
        result = formatter.format(mock_array, context)
        assert "100 × 50" in result.type_name
        assert "float32" in result.type_name
        assert result.css_class == "anndata-dtype--array-api"
        assert "jax" in result.tooltip.lower()  # Backend info in tooltip

    def test_unknown_array_type_graceful_fallback(self):
        """Test that completely unknown array types don't break HTML rendering."""
        FutureArray = type(
            "FutureArray",
            (),
            {
                "__module__": "future_lib.arrays",
                "__init__": lambda self: setattr(self, "_initialized", True),
                "shape": property(lambda self: (10, 5)),
                "dtype": property(lambda self: np.dtype("int64")),
                "ndim": property(lambda self: 2),
            },
        )

        adata = AnnData(np.zeros((10, 5)))
        future_array = FutureArray()
        adata.uns["future_data"] = future_array

        html = adata._repr_html_()
        assert html is not None
        assert len(html) > 0
        assert (
            "FutureArray" in html
            or "future" in html.lower()
            or "object" in html.lower()
        )

    def test_css_array_api_styling_exists(self):
        """Test that CSS styling for Array-API arrays is present."""
        from anndata._repr.css import get_css

        css = get_css()

        assert "anndata-dtype--array-api" in css
        pattern = r"\.anndata-dtype--array-api\s*\{[^}]*color:"
        assert re.search(pattern, css)


class TestCustomHtmlContent:
    """Tests for custom HTML content in Type Formatters."""

    def test_inline_html_content(self):
        """Test inline (non-expandable) custom HTML content in preview column."""
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
                    css_class="anndata-dtype--custom",
                    preview_html='<span class="test-inline">Inline Preview</span>',
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
                    css_class="anndata-dtype--tree",
                    expanded_html=tree_html,
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
            assert "expand" in html.lower() or "Expand" in html
            # Should have the tree content somewhere
            assert "test-tree" in html
        finally:
            formatter_registry.unregister_type_formatter(formatter)

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


class TestReadmeIcon:
    """Tests for README icon functionality."""

    def test_readme_icon_shown_with_content(self):
        """Test README icon appears when README is set."""
        adata = AnnData(np.zeros((5, 3)))
        adata.uns["README"] = "This is a test dataset with important documentation."
        html = adata._repr_html_()
        # Should show the README content or icon
        assert "README" in html or "readme" in html.lower()

    def test_readme_content_accessible(self):
        """Test README content is accessible."""
        adata = AnnData(np.zeros((5, 3)))
        readme_text = "Dataset contains single-cell RNA-seq data."
        adata.uns["README"] = readme_text
        html = adata._repr_html_()
        # The readme text should be somewhere in the HTML
        assert readme_text in html or "README" in html
