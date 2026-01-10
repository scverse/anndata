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

import contextlib
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from .._repr_constants import (
    COLOR_PREVIEW_LIMIT,
    CSS_DTYPE_ANNDATA,
    CSS_DTYPE_ARRAY_API,
    CSS_DTYPE_AWKWARD,
    CSS_DTYPE_BOOL,
    CSS_DTYPE_CATEGORY,
    CSS_DTYPE_DASK,
    CSS_DTYPE_DATAFRAME,
    CSS_DTYPE_FLOAT,
    CSS_DTYPE_GPU,
    CSS_DTYPE_INT,
    CSS_DTYPE_OBJECT,
    CSS_DTYPE_SPARSE,
    CSS_DTYPE_STRING,
    CSS_DTYPE_UNKNOWN,
    CSS_NESTED_ANNDATA,
    CSS_TEXT_MUTED,
)
from .components import render_category_list
from .lazy import get_lazy_categorical_info, is_lazy_column
from .registry import (
    FormattedOutput,
    TypeFormatter,
    formatter_registry,
)
from .utils import (
    check_color_category_mismatch,
    escape_html,
    format_number,
    get_categories_for_display,
    get_matching_column_colors,
    get_setting,
    is_color_list,
    is_serializable,
    preview_dict,
    preview_number,
    preview_sequence,
    preview_string,
    sanitize_css_color,
    should_warn_string_column,
)

if TYPE_CHECKING:
    from typing import Any

    from .registry import FormatterContext


def _check_array_has_writer(array: Any) -> bool:
    """Check if an array type has a registered IO writer.

    This uses the actual IO registry, making it future-proof: if a writer
    is registered for a new type (e.g., datetime64), this will detect it.
    """
    try:
        from .._io.specs.registry import _REGISTRY

        _REGISTRY.get_spec(array)
        return True
    except (KeyError, TypeError):
        return False


def _check_series_backing_array(series: pd.Series) -> tuple[bool, str]:
    """Check if a Series' backing array type can be serialized.

    Uses the IO registry to check the underlying array. This is future-proof:
    if anndata adds support for datetime64/timedelta64/etc, this will detect it.

    Returns (is_serializable, reason_if_not).
    """
    # Standard numpy dtypes are always serializable (no registry check needed)
    # This covers: float16/32/64, int8/16/32/64, uint*, bool, complex*, bytes, str
    if series.dtype.kind in ("f", "i", "u", "b", "c", "S", "U"):
        return True, ""

    # Get the backing array for extension dtypes
    backing_array = series.array

    # NumpyExtensionArray wraps numpy arrays - check the underlying numpy array
    if type(backing_array).__name__ == "NumpyExtensionArray":
        # The underlying numpy array is serializable
        return True, ""

    # For other extension arrays (DatetimeArray, ArrowStringArray, etc.),
    # check the IO registry. This is future-proof: if anndata adds support
    # for datetime64, the registry will have a writer and this returns True.
    if _check_array_has_writer(backing_array):
        return True, ""

    # No writer registered - provide a helpful message
    dtype_name = str(series.dtype)
    return False, f"{dtype_name} not serializable"


def _check_series_serializability(series: pd.Series) -> tuple[bool, str]:
    """
    Check if an object-dtype Series contains serializable values.

    For object dtype columns, checks the first non-null value to determine
    if the column can be written to H5AD/Zarr. Uses anndata's actual IO
    mechanism to test serializability.

    Parameters
    ----------
    series
        Pandas Series with object dtype

    Returns
    -------
    tuple of (is_serializable, reason_if_not)
    """
    if len(series) == 0:
        return True, ""

    # Get first non-null value
    first_valid_idx = series.first_valid_index()
    if first_valid_idx is None:
        return True, ""  # All null

    value = series.loc[first_valid_idx]

    # Object dtype columns with non-string/numeric values are problematic
    # Check if value is a type that anndata can serialize in a DataFrame column
    if isinstance(value, (list, tuple)):
        # Lists/tuples in DataFrame columns are not directly serializable
        # (they work in uns as arrays, but not as DataFrame cell values)
        # NOTE: If https://github.com/scverse/anndata/issues/1923 is resolved,
        # lists of strings may become serializable - update this check accordingly
        return False, f"Contains {type(value).__name__}"
    elif isinstance(value, dict):
        return False, "Contains dict"
    elif not isinstance(value, str | bytes | np.generic | int | float | bool):
        # Custom objects are not serializable
        return False, f"Contains {type(value).__name__}"

    return True, ""


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

        if arr.ndim == 2:
            type_name = f"ndarray ({shape_str}) {dtype_str}"
        elif arr.ndim == 1:
            type_name = f"ndarray ({shape_str},) {dtype_str}"
        else:
            type_name = f"ndarray {arr.shape} {dtype_str}"

        # For obsm/varm sections, show number of columns in preview
        preview = None
        if context.section in ("obsm", "varm") and arr.ndim == 2:
            n_cols = arr.shape[1]
            preview = f"({format_number(n_cols)} columns)"

        return FormattedOutput(
            type_name=type_name,
            css_class=css_class,
            preview=preview,
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
            is_serializable=True,
        )


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

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:  # noqa: PLR0912
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

        # Build type_name with sparsity info inline
        nnz_formatted = format_number(obj.nnz)
        if sparsity_str:
            type_name = (
                f"{format_name} ({shape_str}) {dtype_str} · "
                f"{sparsity_str} ({nnz_formatted} stored)"
            )
        else:
            type_name = f"{format_name} ({shape_str}) {dtype_str}"

        return FormattedOutput(
            type_name=type_name,
            css_class=CSS_DTYPE_SPARSE,
            tooltip=f"{nnz_formatted} stored elements",
            is_serializable=True,
        )


class BackedSparseDatasetFormatter(TypeFormatter):
    """Formatter for anndata's backed sparse datasets (_CSRDataset, _CSCDataset).

    These are HDF5/Zarr-backed sparse matrices that stay on disk.
    Only metadata (shape, dtype, format) is read — no data is loaded.
    """

    priority = 110  # Higher than SparseMatrixFormatter to check first

    def can_format(self, obj: Any) -> bool:
        # Check for anndata's backed sparse dataset classes
        module = type(obj).__module__
        return module.startswith("anndata._core.sparse_dataset") and hasattr(
            obj, "format"
        )

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        shape_str = " × ".join(format_number(s) for s in obj.shape)
        dtype_str = str(obj.dtype)
        format_name = getattr(obj, "format", "sparse")

        return FormattedOutput(
            type_name=f"{format_name}_matrix ({shape_str}) {dtype_str} · on disk",
            css_class=CSS_DTYPE_SPARSE,
            tooltip="Backed sparse matrix (data stays on disk)",
            is_serializable=True,
        )


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

        # Build preview_html with column list for obsm/varm sections
        # Uses anndata-columns class for CSS truncation and JS wrap button
        preview_html = None
        if n_cols > 0 and context.section in ("obsm", "varm"):
            col_str = ", ".join(escape_html(str(c)) for c in cols)
            preview_html = f'<span class="anndata-columns">[{col_str}]</span>'

        # Check if expandable _repr_html_ is enabled
        expand_dataframes = get_setting("repr_html_dataframe_expand", default=False)

        expanded_html = None
        if expand_dataframes and n_rows > 0 and n_cols > 0:
            # Use pandas _repr_html_() for native Jupyter-style output
            # Respects pd.options.display settings (max_rows, max_columns, etc.)
            # Intentional broad catch: _repr_html_() can fail in many ways
            # (memory, recursion, custom dtypes, etc.) - gracefully degrade
            with contextlib.suppress(Exception):
                expanded_html = df._repr_html_()

        return FormattedOutput(
            type_name=f"DataFrame ({format_number(n_rows)} × {format_number(n_cols)})",
            css_class=CSS_DTYPE_DATAFRAME,
            expanded_html=expanded_html,
            preview_html=preview_html,
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
        css_class = _get_dtype_css_class(series.dtype)

        # Check serializability using the IO registry (future-proof)
        is_serial = True
        warnings = []

        # For non-object dtypes, check if the backing array has a registered writer
        # This is future-proof: if anndata adds datetime64 support, this will detect it
        if series.dtype != np.dtype("object"):
            is_serial, reason = _check_series_backing_array(series)
            if not is_serial:
                warnings.append(reason)

        # Object dtype columns need value-level checking
        elif len(series) > 0:
            is_serial, reason = _check_series_serializability(series)
            if not is_serial:
                warnings.append(reason)

        # Compute unique count for preview column (only for obs/var sections)
        preview = None
        n_unique = None
        if context.section in ("obs", "var"):
            if context.unique_limit > 0 and len(series) <= context.unique_limit:
                # nunique() fails on unhashable types (e.g., lists/dicts)
                with contextlib.suppress(TypeError):
                    n_unique = series.nunique()
            if n_unique is not None:
                preview = f"({n_unique} unique)"

            # Check for string->category conversion warning
            should_warn, warn_msg = should_warn_string_column(series, n_unique)
            if should_warn:
                warnings.append(warn_msg)

        return FormattedOutput(
            type_name=f"{dtype_str}",
            css_class=css_class,
            preview=preview,
            is_serializable=is_serial,
            warnings=warnings,
        )


class CategoricalFormatter(TypeFormatter):
    """Formatter for pandas.Categorical, categorical Series, and xarray DataArrays."""

    priority = 110

    def can_format(self, obj: Any) -> bool:
        # pandas Categorical
        if isinstance(obj, pd.Categorical):
            return True
        # pandas Series with categorical dtype
        if isinstance(obj, pd.Series) and hasattr(obj, "cat"):
            return True
        # Check for lazy categorical (CategoricalArray) without accessing dtype
        # which would trigger loading
        try:
            from anndata.experimental.backed._lazy_arrays import CategoricalArray

            if hasattr(obj, "variable") and hasattr(obj.variable, "_data"):
                lazy_indexed = obj.variable._data
                if hasattr(lazy_indexed, "array") and isinstance(
                    lazy_indexed.array, CategoricalArray
                ):
                    return True
        except ImportError:
            pass
        # Fallback: xarray DataArray with categorical dtype (will load data)
        return (
            hasattr(obj, "dtype")
            and isinstance(obj.dtype, pd.CategoricalDtype)
            and not isinstance(obj, pd.Series | pd.Categorical)
        )

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:  # noqa: PLR0912
        # Determine if this is a lazy (xarray DataArray) categorical
        is_lazy = is_lazy_column(obj)
        n_categories = 0

        # Get number of categories based on object type
        if isinstance(obj, pd.Series):
            n_categories = len(obj.cat.categories)
        elif isinstance(obj, pd.Categorical):
            n_categories = len(obj.categories)
        else:
            # Try to get info from lazy categorical without loading
            lazy_count, _lazy_ordered = get_lazy_categorical_info(obj)
            if lazy_count is not None:
                n_categories = lazy_count
            elif hasattr(obj, "dtype") and hasattr(obj.dtype, "categories"):
                # Fallback: access dtype.categories (will load data)
                n_categories = len(obj.dtype.categories)

        # Format type name - indicate lazy for xarray DataArrays
        type_name = (
            f"category ({n_categories}, lazy)"
            if is_lazy
            else f"category ({n_categories})"
        )

        # Build preview_html with category list and colors
        preview_html = None
        if context.section in ("obs", "var") and context.column_name is not None:
            try:
                # Get categories (respecting lazy loading limits)
                categories, was_truncated, n_total = get_categories_for_display(
                    obj, context, is_lazy=is_lazy
                )

                if len(categories) == 0:
                    # Metadata-only mode or no categories: show just count
                    if n_total is not None:
                        preview_html = f'<span class="{CSS_TEXT_MUTED}">({n_total} categories)</span>'
                    else:
                        preview_html = (
                            f'<span class="{CSS_TEXT_MUTED}">(categories)</span>'
                        )
                else:
                    # Get colors for categories
                    colors = None
                    if context.adata_ref is not None:
                        max_cats = context.max_categories
                        n_cats_to_show = min(len(categories), max_cats)
                        # For lazy with truncation, only load colors we need
                        color_limit = (
                            n_cats_to_show if (is_lazy and was_truncated) else None
                        )
                        colors = get_matching_column_colors(
                            context.adata_ref, context.column_name, limit=color_limit
                        )

                    # Render category list with colors
                    n_hidden = (
                        (n_total - len(categories))
                        if (n_total and was_truncated)
                        else 0
                    )
                    preview_html = render_category_list(
                        categories, colors, context.max_categories, n_hidden=n_hidden
                    )
            except Exception as e:  # noqa: BLE001
                # Never let preview generation crash the repr
                preview_html = (
                    f'<span class="{CSS_TEXT_MUTED}">(error: {type(e).__name__})</span>'
                )

        # Check for color mismatch warning
        warnings = []
        if (
            n_categories > 0
            and context.adata_ref is not None
            and context.column_name is not None
        ):
            color_warning = check_color_category_mismatch(
                context.adata_ref, context.column_name, n_categories
            )
            if color_warning:
                warnings.append(color_warning)

        return FormattedOutput(
            type_name=type_name,
            css_class=CSS_DTYPE_CATEGORY,
            preview_html=preview_html,
            is_serializable=True,
            warnings=warnings,
        )


class LazyColumnFormatter(TypeFormatter):
    """
    Formatter for lazy obs/var columns (xarray DataArray) from read_lazy().

    For lazy AnnData, obs/var columns are xarray DataArrays instead of pandas Series.
    This formatter shows the dtype with "(lazy)" indicator, without the shape since
    all columns in obs/var have the same length (n_obs or n_var).

    Note: Categorical columns are handled by CategoricalFormatter (higher priority).
    """

    priority = (
        60  # Higher than ArrayAPIFormatter (50), lower than CategoricalFormatter (110)
    )
    sections = ("obs", "var")  # Only apply to obs/var sections

    def can_format(self, obj: Any) -> bool:
        # xarray DataArray (categoricals already handled by higher-priority CategoricalFormatter)
        if not (
            hasattr(obj, "dtype") and hasattr(obj, "shape") and hasattr(obj, "ndim")
        ):
            return False

        # Exclude already-handled types
        if isinstance(obj, np.ndarray | pd.DataFrame | pd.Series | pd.Categorical):
            return False

        # Check if it looks like an xarray DataArray (has .data attribute)
        # This is a good heuristic for lazy obs/var columns
        return hasattr(obj, "data")

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        dtype_str = str(obj.dtype)

        # Map common dtypes to CSS classes
        if "int" in dtype_str:
            css_class = CSS_DTYPE_INT
        elif "float" in dtype_str:
            css_class = CSS_DTYPE_FLOAT
        elif "bool" in dtype_str:
            css_class = CSS_DTYPE_BOOL
        elif "str" in dtype_str or dtype_str == "object":
            css_class = CSS_DTYPE_OBJECT
        else:
            css_class = CSS_DTYPE_UNKNOWN

        # For lazy non-categorical columns, we can't compute unique count
        # without loading data, so we just indicate it's lazy
        preview = None
        if context.section in ("obs", "var"):
            preview = "(lazy)"

        return FormattedOutput(
            type_name=f"{dtype_str} (lazy)",
            css_class=css_class,
            preview=preview,
            is_serializable=True,
        )


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
        dtype_str = str(obj.dtype)

        # Get chunk info
        chunks_str = str(obj.chunksize) if hasattr(obj, "chunksize") else "unknown"

        # In obsm/varm/obsp/varp sections, don't show shape in type (redundant)
        # - obsp/varp: always n_obs × n_obs or n_var × n_var
        # - obsm/varm: preview column shows number of columns
        if context.section in ("obsm", "varm", "obsp", "varp"):
            type_name = f"dask.array {dtype_str} · chunks={chunks_str}"
        else:
            shape_str = " × ".join(format_number(s) for s in obj.shape)
            type_name = f"dask.array ({shape_str}) {dtype_str} · chunks={chunks_str}"

        # For obsm/varm sections, show number of columns in preview
        preview = None
        if context.section in ("obsm", "varm") and len(obj.shape) == 2:
            n_cols = obj.shape[1]
            preview = f"({format_number(n_cols)} columns)"

        return FormattedOutput(
            type_name=type_name,
            css_class=CSS_DTYPE_DASK,
            tooltip=f"{obj.npartitions} partitions",
            preview=preview,
            is_serializable=True,
        )


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

        # For obsm/varm sections, show number of columns in preview
        preview = None
        if context.section in ("obsm", "varm") and len(obj.shape) == 2:
            n_cols = obj.shape[1]
            preview = f"({format_number(n_cols)} columns)"

        return FormattedOutput(
            type_name=f"cupy.ndarray ({shape_str}) {dtype_str}",
            css_class=CSS_DTYPE_GPU,
            tooltip=device_info,
            preview=preview,
            is_serializable=True,
        )


class AwkwardArrayFormatter(TypeFormatter):
    """Formatter for awkward.Array (ragged/jagged arrays)."""

    priority = 120

    def can_format(self, obj: Any) -> bool:
        return type(obj).__module__.startswith("awkward")

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        length: int | None = None
        try:
            length = len(obj)
            type_str = str(obj.type) if hasattr(obj, "type") else "unknown"
        except Exception:  # noqa: BLE001
            # Intentional broad catch: awkward arrays can fail on len/type access
            # in edge cases (lazy evaluation, corrupt data) - show placeholder
            type_str = "unknown"

        length_str = str(length) if length is not None else "?"
        return FormattedOutput(
            type_name=f"awkward.Array ({length_str} records)",
            css_class=CSS_DTYPE_AWKWARD,
            tooltip=f"Type: {type_str}",
            is_serializable=True,
        )


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

        # For obsm/varm sections, show number of columns in preview
        preview = None
        if context.section in ("obsm", "varm") and obj.ndim == 2:
            n_cols = obj.shape[1]
            preview = f"({format_number(n_cols)} columns)"

        return FormattedOutput(
            type_name=f"{type_name} ({shape_str}) {dtype_str}",
            css_class=CSS_DTYPE_ARRAY_API,
            tooltip=f"{backend_label} array{device_info}",
            preview=preview,
            is_serializable=True,
        )


class AnnDataFormatter(TypeFormatter):
    """Formatter for nested AnnData objects."""

    priority = 150

    def can_format(self, obj: Any) -> bool:
        # Check by class name to avoid circular imports
        return type(obj).__name__ == "AnnData" and hasattr(obj, "n_obs")

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        shape_str = f"{format_number(obj.n_obs)} × {format_number(obj.n_vars)}"

        # Generate expanded HTML if within depth limit
        expanded_html = None
        if context.depth < context.max_depth - 1:
            # Lazy import to avoid circular dependency
            from .html import generate_repr_html

            nested_html = generate_repr_html(
                obj,
                depth=context.depth + 1,
                max_depth=context.max_depth,
                show_header=True,
                show_search=False,
            )
            expanded_html = f'<div class="{CSS_NESTED_ANNDATA}">{nested_html}</div>'

        return FormattedOutput(
            type_name=f"AnnData ({shape_str})",
            css_class=CSS_DTYPE_ANNDATA,
            tooltip="Nested AnnData object",
            expanded_html=expanded_html,
            is_serializable=True,
        )


class NoneFormatter(TypeFormatter):
    """Formatter for None."""

    priority = 50

    def can_format(self, obj: Any) -> bool:
        return obj is None

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        return FormattedOutput(
            type_name="NoneType",
            css_class=CSS_DTYPE_OBJECT,
            preview="None",
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
            css_class=CSS_DTYPE_BOOL,
            preview=preview_number(obj),
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
            css_class=CSS_DTYPE_INT,
            preview=preview_number(obj),
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
            css_class=CSS_DTYPE_FLOAT,
            preview=preview_number(obj),
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
            css_class=CSS_DTYPE_STRING,
            preview=preview_string(obj, context.max_string_length),
            is_serializable=True,
        )


class DictFormatter(TypeFormatter):
    """Formatter for dict."""

    priority = 50

    def can_format(self, obj: Any) -> bool:
        return isinstance(obj, dict)

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        # Check serializability of contents
        is_serial, reason = is_serializable(obj)
        warnings = [] if is_serial else [reason]

        return FormattedOutput(
            type_name="dict",
            css_class=CSS_DTYPE_OBJECT,
            preview=preview_dict(obj),
            is_serializable=is_serial,
            warnings=warnings,
        )


class ColorListFormatter(TypeFormatter):
    """Formatter for color lists (uns entries ending in _colors)."""

    priority = 60  # Higher than ListFormatter to check first

    def can_format(self, obj: Any) -> bool:
        # Requires context.column_name to be set (from uns key)
        return False  # Will be checked with context in format()

    def can_format_with_context(self, obj: Any, context: FormatterContext) -> bool:
        """Check if this is a color list based on key name and value."""
        key = context.column_name
        return key is not None and is_color_list(key, obj)

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        colors = list(obj) if hasattr(obj, "__iter__") else []
        n_colors = len(colors)

        # Build color swatch HTML with sanitized colors
        swatches = []
        for color in colors[:COLOR_PREVIEW_LIMIT]:
            # Sanitize color to prevent CSS injection
            safe_color = sanitize_css_color(str(color))
            if safe_color:
                swatches.append(
                    f'<span class="anndata-colors__swatch" '
                    f'style="background:{safe_color}" title="{escape_html(str(color))}"></span>'
                )
            else:
                # Invalid/unsafe color - show as text only, no style
                swatches.append(
                    f'<span class="anndata-colors__swatch anndata-colors__swatch--invalid" '
                    f'title="{escape_html(str(color))}">?</span>'
                )
        if n_colors > COLOR_PREVIEW_LIMIT:
            swatches.append(
                f'<span class="{CSS_TEXT_MUTED}">+{n_colors - COLOR_PREVIEW_LIMIT}</span>'
            )

        preview_html = f'<span class="anndata-colors">{"".join(swatches)}</span>'

        return FormattedOutput(
            type_name=f"colors ({n_colors})",
            css_class=CSS_DTYPE_OBJECT,
            preview_html=preview_html,
            is_serializable=True,
        )


class ListFormatter(TypeFormatter):
    """Formatter for list and tuple."""

    priority = 50

    def can_format(self, obj: Any) -> bool:
        return isinstance(obj, (list, tuple))

    def format(self, obj: Any, context: FormatterContext) -> FormattedOutput:
        type_name = "list" if isinstance(obj, list) else "tuple"

        # Check serializability
        is_serial, reason = is_serializable(obj)
        warnings = [] if is_serial else [reason]

        return FormattedOutput(
            type_name=type_name,
            css_class=CSS_DTYPE_OBJECT,
            preview=preview_sequence(obj),
            is_serializable=is_serial,
            warnings=warnings,
        )


def _get_dtype_css_class(dtype: np.dtype | pd.api.types.CategoricalDtype) -> str:  # noqa: PLR0911
    """Get CSS class for a numpy or pandas dtype."""
    # Check for pandas CategoricalDtype first (has kind="O" but is special)
    dtype_name = str(dtype)
    if dtype_name == "category":
        return CSS_DTYPE_CATEGORY

    # Try numpy dtype.kind (most reliable for standard dtypes)
    kind = getattr(dtype, "kind", None)
    if kind is not None:
        if kind in ("i", "u"):
            return CSS_DTYPE_INT
        if kind == "f":
            return CSS_DTYPE_FLOAT
        if kind == "b":
            return CSS_DTYPE_BOOL
        if kind in ("U", "S", "O"):
            return CSS_DTYPE_STRING
        if kind == "c":
            return CSS_DTYPE_FLOAT  # complex

    # Fallback to string matching for pandas extension dtypes
    if "int" in dtype_name:
        return CSS_DTYPE_INT
    if "float" in dtype_name:
        return CSS_DTYPE_FLOAT
    if "bool" in dtype_name:
        return CSS_DTYPE_BOOL
    if "object" in dtype_name or "string" in dtype_name:
        return CSS_DTYPE_STRING
    return CSS_DTYPE_OBJECT


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
        BackedSparseDatasetFormatter(),  # Before SparseMatrixFormatter (backed sparse)
        # Medium priority
        NumpyArrayFormatter(),
        SparseMatrixFormatter(),
        DataFrameFormatter(),
        SeriesFormatter(),
        # Lazy obs/var columns (xarray DataArray) - must come before ArrayAPIFormatter
        LazyColumnFormatter(),
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
        ColorListFormatter(),  # Before ListFormatter (higher priority for *_colors keys)
        ListFormatter(),
    ]

    for formatter in formatters:
        formatter_registry.register_type_formatter(formatter)


# Auto-register on import
_register_builtin_formatters()
