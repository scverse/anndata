"""
Built-in formatters for common types.

This module registers formatters for:
- NumPy arrays (dense and masked)
- SciPy sparse matrices
- Pandas DataFrames, Series, Categorical
- Dask arrays
- CuPy arrays (GPU)
- Awkward arrays
- AnnData objects (for recursive display in .uns)
- Python built-in types
- Color lists

The formatters are registered automatically when this module is imported.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from anndata._repr.registry import (
    FormattedOutput,
    TypeFormatter,
    formatter_registry,
)
from anndata._repr.utils import (
    format_number,
    is_serializable,
)

if TYPE_CHECKING:
    from typing import Any

    from anndata._repr.registry import FormatterContext


# =============================================================================
# NumPy Formatters
# =============================================================================


class NumpyArrayFormatter(TypeFormatter):
    """Formatter for numpy.ndarray."""

    priority = 100

    def can_format(self, obj: Any) -> bool:
        return isinstance(obj, np.ndarray)

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        arr: np.ndarray = obj
        shape_str = " × ".join(format_number(s) for s in arr.shape)
        dtype_str = str(arr.dtype)

        # Determine CSS class based on dtype
        css_class = _get_dtype_css_class(arr.dtype)

        details = {
            "shape": arr.shape,
            "dtype": dtype_str,
            "ndim": arr.ndim,
        }

        if arr.ndim == 2:
            type_name = f"ndarray ({shape_str}) {dtype_str}"
        elif arr.ndim == 1:
            type_name = f"ndarray ({shape_str},) {dtype_str}"
        else:
            type_name = f"ndarray {arr.shape} {dtype_str}"

        return FormattedOutput(
            type_name=type_name,
            css_class=css_class,
            details=details,
            is_serializable=True,
        )


class NumpyMaskedArrayFormatter(TypeFormatter):
    """Formatter for numpy.ma.MaskedArray."""

    priority = 110

    def can_format(self, obj: Any) -> bool:
        return isinstance(obj, np.ma.MaskedArray)

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        arr: np.ma.MaskedArray = obj
        shape_str = " × ".join(format_number(s) for s in arr.shape)
        dtype_str = str(arr.dtype)
        n_masked = int(np.sum(arr.mask)) if arr.mask is not np.ma.nomask else 0

        return FormattedOutput(
            type_name=f"MaskedArray ({shape_str}) {dtype_str}",
            css_class=_get_dtype_css_class(arr.dtype),
            tooltip=f"{n_masked} masked values" if n_masked > 0 else "",
            details={
                "shape": arr.shape,
                "dtype": dtype_str,
                "n_masked": n_masked,
            },
            is_serializable=True,
        )


# =============================================================================
# SciPy Sparse Formatters
# =============================================================================


class SparseMatrixFormatter(TypeFormatter):
    """
    Formatter for scipy.sparse matrices and arrays.

    Future-proofing notes:
    - PR #1927 (https://github.com/scverse/anndata/pull/1927) removes scipy sparse inheritance
    - Uses duck typing as fallback to detect sparse-like objects without relying on isinstance()
    - Handles both scipy.sparse (CPU) and cupyx.scipy.sparse (GPU) sparse arrays
    """

    priority = 100

    def can_format(self, obj: Any) -> bool:
        # First try scipy.sparse.issparse() if available (backward compatibility)
        try:
            import scipy.sparse as sp

            if sp.issparse(obj):
                return True
        except ImportError:
            pass

        # Fallback: Duck typing for sparse-like objects
        # Future-proof against PR #1927 removing scipy sparse inheritance
        # A sparse object should have: nnz (non-zero count), shape, dtype, and sparse conversion methods
        module = type(obj).__module__
        is_sparse_module = module.startswith(("scipy.sparse", "cupyx.scipy.sparse"))
        has_sparse_attrs = (
            hasattr(obj, "nnz")
            and hasattr(obj, "shape")
            and hasattr(obj, "dtype")
            and (hasattr(obj, "tocsr") or hasattr(obj, "tocsc"))
        )

        return is_sparse_module and has_sparse_attrs

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        shape_str = " × ".join(format_number(s) for s in obj.shape)
        dtype_str = str(obj.dtype)

        # Calculate sparsity
        n_elements = obj.shape[0] * obj.shape[1] if len(obj.shape) == 2 else 1
        if n_elements > 0:
            sparsity = 1 - (obj.nnz / n_elements)
            sparsity_str = f"{sparsity:.1%} sparse"
        else:
            sparsity_str = ""

        # Determine format name
        # Try scipy-specific checks first (backward compatibility)
        format_name = None
        try:
            import scipy.sparse as sp

            if sp.isspmatrix_csr(obj):
                format_name = "csr_matrix"
            elif sp.isspmatrix_csc(obj):
                format_name = "csc_matrix"
            elif sp.isspmatrix_coo(obj):
                format_name = "coo_matrix"
            elif sp.isspmatrix_lil(obj):
                format_name = "lil_matrix"
            elif sp.isspmatrix_dok(obj):
                format_name = "dok_matrix"
            elif sp.isspmatrix_dia(obj):
                format_name = "dia_matrix"
            elif sp.isspmatrix_bsr(obj):
                format_name = "bsr_matrix"
        except (ImportError, TypeError):
            # ImportError: scipy not available
            # TypeError: isspmatrix_* functions may fail on new sparse array types (PR #1927)
            pass

        # Fallback: Use type name (works for new sparse array classes like csr_array, csc_array)
        if format_name is None:
            format_name = type(obj).__name__

        return FormattedOutput(
            type_name=f"{format_name} ({shape_str}) {dtype_str}",
            css_class="dtype-sparse",
            tooltip=f"{format_number(obj.nnz)} stored elements, {sparsity_str}",
            details={
                "shape": obj.shape,
                "dtype": dtype_str,
                "nnz": obj.nnz,
                "format": format_name,
                "sparsity": sparsity if n_elements > 0 else None,
            },
            is_serializable=True,
        )


# =============================================================================
# Pandas Formatters
# =============================================================================


class DataFrameFormatter(TypeFormatter):
    """Formatter for pandas.DataFrame.

    Shows column names in the meta column. Can optionally show full DataFrame
    as expandable content via pandas ``_repr_html_()`` - controlled by setting
    ``anndata.settings.repr_html_dataframe_expand`` (default: False).

    When expanded, uses the rich Jupyter-style output from pandas. Configure
    the display with pandas options::

        pd.set_option("display.max_rows", 10)
        pd.set_option("display.max_columns", 5)
    """

    priority = 100

    def can_format(self, obj: Any) -> bool:
        return isinstance(obj, pd.DataFrame)

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        df: pd.DataFrame = obj
        n_rows, n_cols = len(df), len(df.columns)
        cols = list(df.columns)

        # Build compact column preview for meta column
        meta_preview = ""
        meta_preview_full = ""
        if n_cols > 0:
            max_cols_preview = 5
            max_total_len = 40
            meta_preview_full = f"columns: [{', '.join(str(c) for c in cols)}]"
            if n_cols <= max_cols_preview:
                col_str = ", ".join(str(c) for c in cols)
            else:
                col_str = ", ".join(str(c) for c in cols[:max_cols_preview]) + ", …"
            # Truncate if still too long
            if len(col_str) > max_total_len:
                col_str = col_str[: max_total_len - 1] + "…"
            meta_preview = f"[{col_str}]"

        # Check if expandable _repr_html_ is enabled
        try:
            from anndata import settings

            expand_dataframes = getattr(settings, "repr_html_dataframe_expand", False)
        except (ImportError, AttributeError):
            expand_dataframes = False

        html_content = None
        is_expandable = False
        if expand_dataframes and n_rows > 0 and n_cols > 0:
            # Use pandas _repr_html_() for native Jupyter-style output
            # Respects pd.options.display settings (max_rows, max_columns, etc.)
            try:
                html_content = df._repr_html_()
                is_expandable = True
            except Exception:  # noqa: BLE001
                # Intentional broad catch: _repr_html_() can fail in many ways
                # (memory, recursion, custom dtypes, etc.) - gracefully degrade
                pass

        return FormattedOutput(
            type_name=f"DataFrame ({format_number(n_rows)} × {format_number(n_cols)})",
            css_class="dtype-dataframe",
            html_content=html_content,
            is_expandable=is_expandable,
            details={
                "n_rows": n_rows,
                "n_cols": n_cols,
                "columns": cols,
                "meta_preview": meta_preview,
                "meta_preview_full": meta_preview_full,
                "has_columns_list": n_cols > 0,  # Flag for wrap button
            },
            is_serializable=True,
        )


class SeriesFormatter(TypeFormatter):
    """Formatter for pandas.Series."""

    priority = 100

    def can_format(self, obj: Any) -> bool:
        return isinstance(obj, pd.Series) and not isinstance(
            obj.dtype, pd.CategoricalDtype
        )

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        series: pd.Series = obj
        dtype_str = str(series.dtype)
        css_class = _get_pandas_dtype_css_class(series.dtype)

        return FormattedOutput(
            type_name=f"{dtype_str}",
            css_class=css_class,
            details={
                "length": len(series),
                "dtype": dtype_str,
            },
            is_serializable=True,
        )


class CategoricalFormatter(TypeFormatter):
    """Formatter for pandas.Categorical and categorical Series."""

    priority = 110

    def can_format(self, obj: Any) -> bool:
        return isinstance(obj, pd.Categorical) or (
            isinstance(obj, pd.Series) and hasattr(obj, "cat")
        )

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        if isinstance(obj, pd.Series):
            cat = obj.cat
            n_categories = len(cat.categories)
        else:
            cat = obj
            n_categories = len(obj.categories)

        return FormattedOutput(
            type_name=f"category ({n_categories})",
            css_class="dtype-category",
            details={
                "n_categories": n_categories,
                "ordered": getattr(cat, "ordered", False),
            },
            is_serializable=True,
        )


# =============================================================================
# Dask Array Formatter
# =============================================================================


class DaskArrayFormatter(TypeFormatter):
    """Formatter for dask.array.Array."""

    priority = 120

    def can_format(self, obj: Any) -> bool:
        try:
            import dask.array as da

            return isinstance(obj, da.Array)
        except ImportError:
            return False

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        shape_str = " × ".join(format_number(s) for s in obj.shape)
        dtype_str = str(obj.dtype)

        # Get chunk info
        chunks_str = str(obj.chunksize) if hasattr(obj, "chunksize") else "unknown"

        return FormattedOutput(
            type_name=f"dask.array ({shape_str}) {dtype_str}",
            css_class="dtype-dask",
            tooltip=f"Chunks: {chunks_str}, {obj.npartitions} partitions",
            details={
                "shape": obj.shape,
                "dtype": dtype_str,
                "chunks": obj.chunks,
                "npartitions": obj.npartitions,
            },
            is_serializable=True,
        )


# =============================================================================
# CuPy Array Formatter
# =============================================================================


class CuPyArrayFormatter(TypeFormatter):
    """Formatter for cupy.ndarray (GPU arrays)."""

    priority = 120

    def can_format(self, obj: Any) -> bool:
        return type(obj).__module__.startswith("cupy")

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        shape_str = " × ".join(format_number(s) for s in obj.shape)
        dtype_str = str(obj.dtype)

        device_info = ""
        if hasattr(obj, "device"):
            device_info = f"GPU:{obj.device.id}"

        return FormattedOutput(
            type_name=f"cupy.ndarray ({shape_str}) {dtype_str}",
            css_class="dtype-gpu",
            tooltip=device_info,
            details={
                "shape": obj.shape,
                "dtype": dtype_str,
                "device": getattr(obj, "device", None),
            },
            is_serializable=True,
        )


# =============================================================================
# Awkward Array Formatter
# =============================================================================


class AwkwardArrayFormatter(TypeFormatter):
    """Formatter for awkward.Array (ragged/jagged arrays)."""

    priority = 120

    def can_format(self, obj: Any) -> bool:
        return type(obj).__module__.startswith("awkward")

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        try:
            length = len(obj)
            type_str = str(obj.type) if hasattr(obj, "type") else "unknown"
        except Exception:  # noqa: BLE001
            # Intentional broad catch: awkward arrays can fail on len/type access
            # in edge cases (lazy evaluation, corrupt data) - show placeholder
            length = "?"
            type_str = "unknown"

        return FormattedOutput(
            type_name=f"awkward.Array ({length} records)",
            css_class="dtype-awkward",
            tooltip=f"Type: {type_str}",
            details={
                "length": length,
                "type": type_str,
            },
            is_serializable=True,
        )


# =============================================================================
# Array-API Compatible Array Formatter
# =============================================================================


class ArrayAPIFormatter(TypeFormatter):
    """
    Formatter for Array-API compatible arrays (JAX, PyTorch, TensorFlow, etc.).

    Future-proofing notes:
    - PR #2063 (https://github.com/scverse/anndata/pull/2063) adds Array-API compatibility
    - Handles JAX arrays, PyTorch tensors, TensorFlow tensors, and other array-like objects
    - Uses duck typing to detect array-like objects without isinstance checks
    - Low priority (50) ensures specific formatters (numpy, cupy, etc.) are tried first

    References:
    - Array API Standard: https://data-apis.org/array-api/latest/
    """

    priority = 50  # Lower than specific formatters (numpy=110, cupy=120) but higher than builtins

    def can_format(self, obj: Any) -> bool:
        # Duck typing: Check for array-like attributes
        # Must have shape, dtype, and ndim (array-api standard)
        has_array_attrs = (
            hasattr(obj, "shape") and hasattr(obj, "dtype") and hasattr(obj, "ndim")
        )

        if not has_array_attrs:
            return False

        # Exclude types that already have specific formatters
        # This prevents conflicts and ensures more specific formatters take precedence
        module = type(obj).__module__
        already_handled = isinstance(
            obj, (np.ndarray, pd.DataFrame, pd.Series)
        ) or module.startswith((
            "numpy",
            "pandas",
            "scipy.sparse",
            "cupy",  # Has CuPyArrayFormatter
            "cupyx",
            "awkward",  # Has AwkwardArrayFormatter
            "dask",  # Has DaskArrayFormatter
        ))

        return not already_handled

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        # Extract module and type information
        module_name = type(obj).__module__.split(".")[0]  # e.g., "jax" from "jax.numpy"
        type_name = type(obj).__name__

        shape_str = " × ".join(format_number(s) for s in obj.shape)
        dtype_str = str(obj.dtype)

        # Detect common array-api backends
        # Known backends: JAX, PyTorch, TensorFlow, MXNet, etc.
        known_backends = {
            "jax": "JAX",
            "jaxlib": "JAX",
            "torch": "PyTorch",
            "tensorflow": "TensorFlow",
            "tf": "TensorFlow",
            "mxnet": "MXNet",
        }
        backend_label = known_backends.get(module_name, module_name)

        # Try to get device information (GPU vs CPU)
        device_info = ""
        try:
            if hasattr(obj, "device"):
                device_info = f" on {obj.device}"
            elif hasattr(obj, "device_buffer"):  # JAX
                device_info = f" on {obj.device_buffer.device()}"
        except Exception:  # noqa: BLE001
            # Intentional broad catch: device access varies by backend and can fail
            pass

        return FormattedOutput(
            type_name=f"{type_name} ({shape_str}) {dtype_str}",
            css_class="dtype-array-api",
            tooltip=f"{backend_label} array{device_info}",
            details={
                "shape": obj.shape,
                "dtype": dtype_str,
                "backend": module_name,
                "type": type_name,
                "ndim": obj.ndim,
            },
            is_serializable=True,
        )


# =============================================================================
# AnnData Formatter (for nested AnnData in .uns)
# =============================================================================


class AnnDataFormatter(TypeFormatter):
    """Formatter for nested AnnData objects."""

    priority = 150

    def can_format(self, obj: Any) -> bool:
        # Check by class name to avoid circular imports
        return type(obj).__name__ == "AnnData" and hasattr(obj, "n_obs")

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        shape_str = f"{format_number(obj.n_obs)} × {format_number(obj.n_vars)}"

        return FormattedOutput(
            type_name=f"AnnData ({shape_str})",
            css_class="dtype-anndata",
            tooltip="Nested AnnData object",
            details={
                "n_obs": obj.n_obs,
                "n_vars": obj.n_vars,
                "is_view": getattr(obj, "is_view", False),
            },
            is_expandable=context.depth < context.max_depth,
            is_serializable=True,
        )


# =============================================================================
# Python Built-in Type Formatters
# =============================================================================


class NoneFormatter(TypeFormatter):
    """Formatter for None."""

    priority = 50

    def can_format(self, obj: Any) -> bool:
        return obj is None

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        return FormattedOutput(
            type_name="NoneType",
            css_class="dtype-object",
            is_serializable=True,
        )


class BoolFormatter(TypeFormatter):
    """Formatter for bool."""

    priority = 50

    def can_format(self, obj: Any) -> bool:
        return isinstance(obj, bool)

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        return FormattedOutput(
            type_name="bool",
            css_class="dtype-bool",
            is_serializable=True,
        )


class IntFormatter(TypeFormatter):
    """Formatter for int."""

    priority = 50

    def can_format(self, obj: Any) -> bool:
        return isinstance(obj, (int, np.integer)) and not isinstance(obj, bool)

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        return FormattedOutput(
            type_name="int",
            css_class="dtype-int",
            is_serializable=True,
        )


class FloatFormatter(TypeFormatter):
    """Formatter for float."""

    priority = 50

    def can_format(self, obj: Any) -> bool:
        return isinstance(obj, (float, np.floating))

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        return FormattedOutput(
            type_name="float",
            css_class="dtype-float",
            is_serializable=True,
        )


class StringFormatter(TypeFormatter):
    """Formatter for str."""

    priority = 50

    def can_format(self, obj: Any) -> bool:
        return isinstance(obj, str)

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        return FormattedOutput(
            type_name="str",
            css_class="dtype-string",
            is_serializable=True,
        )


class DictFormatter(TypeFormatter):
    """Formatter for dict."""

    priority = 50

    def can_format(self, obj: Any) -> bool:
        return isinstance(obj, dict)

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        n_items = len(obj)

        # Check serializability of contents
        is_serial, reason = is_serializable(obj)
        warnings = [] if is_serial else [reason]

        return FormattedOutput(
            type_name="dict",
            css_class="dtype-object",
            details={"n_items": n_items, "keys": list(obj.keys())[:10]},
            is_expandable=n_items > 0 and context.depth < context.max_depth,
            is_serializable=is_serial,
            warnings=warnings,
        )


class ListFormatter(TypeFormatter):
    """Formatter for list and tuple."""

    priority = 50

    def can_format(self, obj: Any) -> bool:
        return isinstance(obj, (list, tuple))

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        type_name = "list" if isinstance(obj, list) else "tuple"
        n_items = len(obj)

        # Check serializability
        is_serial, reason = is_serializable(obj)
        warnings = [] if is_serial else [reason]

        return FormattedOutput(
            type_name=type_name,
            css_class="dtype-object",
            details={"n_items": n_items},
            is_serializable=is_serial,
            warnings=warnings,
        )


# =============================================================================
# Color List Formatter
# =============================================================================


class ColorListFormatter(TypeFormatter):
    """Special formatter for color lists (*_colors in .uns)."""

    priority = 200  # High priority to catch before generic list

    def can_format(self, obj: Any) -> bool:
        # This is context-dependent, handled in section formatter
        return False

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        # Not used directly - handled specially in UnsSection
        raise NotImplementedError


# =============================================================================
# Helper Functions
# =============================================================================


def _get_dtype_css_class(dtype: np.dtype) -> str:
    """Get CSS class for a numpy dtype."""
    kind = dtype.kind
    if kind in ("i", "u"):
        return "dtype-int"
    elif kind == "f":
        return "dtype-float"
    elif kind == "b":
        return "dtype-bool"
    elif kind in ("U", "S", "O"):
        return "dtype-string"
    elif kind == "c":
        return "dtype-float"  # complex
    else:
        return "dtype-object"


def _get_pandas_dtype_css_class(dtype) -> str:
    """Get CSS class for a pandas dtype."""
    dtype_name = str(dtype)
    if "int" in dtype_name:
        return "dtype-int"
    elif "float" in dtype_name:
        return "dtype-float"
    elif "bool" in dtype_name:
        return "dtype-bool"
    elif dtype_name == "category":
        return "dtype-category"
    elif "object" in dtype_name or "string" in dtype_name:
        return "dtype-string"
    else:
        return "dtype-object"


# =============================================================================
# Register all formatters
# =============================================================================


def _register_builtin_formatters() -> None:
    """Register all built-in formatters with the global registry."""
    formatters = [
        # High priority (specific types)
        AnnDataFormatter(),
        DaskArrayFormatter(),
        CuPyArrayFormatter(),
        AwkwardArrayFormatter(),
        NumpyMaskedArrayFormatter(),
        CategoricalFormatter(),
        # Medium priority
        NumpyArrayFormatter(),
        SparseMatrixFormatter(),
        DataFrameFormatter(),
        SeriesFormatter(),
        # Low-medium priority (Array-API compatible arrays)
        # Must come after specific array formatters (numpy, cupy, etc.) but before builtins
        # Handles JAX, PyTorch, TensorFlow arrays added in PR #2063
        ArrayAPIFormatter(),
        # Low priority (builtins)
        NoneFormatter(),
        BoolFormatter(),
        IntFormatter(),
        FloatFormatter(),
        StringFormatter(),
        DictFormatter(),
        ListFormatter(),
    ]

    for formatter in formatters:
        formatter_registry.register_type_formatter(formatter)


# Auto-register on import
_register_builtin_formatters()
