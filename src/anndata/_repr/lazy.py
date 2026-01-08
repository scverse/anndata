"""
Lazy loading utilities for AnnData HTML representation.

This module consolidates all logic related to detecting and handling lazy AnnData
objects (from read_lazy()). Lazy AnnData uses xarray-backed storage and requires
special handling to avoid triggering data loading during repr generation.

Key concepts:
- Lazy AnnData: Created by read_lazy(), obs/var are Dataset2D (xarray-backed)
- Lazy series: Individual columns from Dataset2D, implemented as xarray DataArrays
- CategoricalArray: anndata's lazy categorical implementation for zarr/h5 storage

Usage:
    from .lazy import (
        is_lazy_adata,
        is_lazy_column,
        get_lazy_category_count,
        get_lazy_categories,
        get_lazy_categorical_info,
    )
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any

    from .registry import FormatterContext


def _get_categorical_array(col: Any):
    """
    Get the underlying CategoricalArray from a lazy xarray DataArray.

    Navigates through the xarray structure:
    DataArray -> Variable -> LazilyIndexedArray -> CategoricalArray

    Parameters
    ----------
    col
        The column (potentially an xarray DataArray) to extract from

    Returns
    -------
    CategoricalArray if found, None otherwise
    """
    try:
        from anndata.experimental.backed._lazy_arrays import CategoricalArray

        # Navigate through xarray structure to find CategoricalArray
        # DataArray -> Variable -> LazilyIndexedArray -> CategoricalArray
        if hasattr(col, "variable") and hasattr(col.variable, "_data"):
            lazy_indexed = col.variable._data
            if hasattr(lazy_indexed, "array"):
                arr = lazy_indexed.array
                if isinstance(arr, CategoricalArray):
                    return arr
    except ImportError:
        pass
    return None


def is_lazy_adata(obj: Any) -> bool:
    """Check if an AnnData uses lazy loading (experimental read_lazy).

    Lazy AnnData has Dataset2D (xarray-backed) obs/var instead of regular DataFrames.

    Parameters
    ----------
    obj
        Object to check (typically an AnnData)

    Returns
    -------
    True if obj is a lazy AnnData
    """
    obs = getattr(obj, "obs", None)
    if obs is None:
        return False
    # Dataset2D has a different class name than DataFrame
    return obs.__class__.__name__ == "Dataset2D"


def _extract_path_from_lazy_array(arr: Any) -> dict[str, Any] | None:
    """Extract file path and format from a lazy array (CategoricalArray/MaskedArray)."""
    from pathlib import Path

    base_path = arr.base_path_or_zarr_group
    file_format = getattr(arr, "file_format", "")

    # H5AD files have a Path as base_path
    if isinstance(base_path, Path):
        fmt = "H5AD" if file_format == "h5" else "Zarr"
        return {"filename": str(base_path), "format": fmt}

    # For zarr groups, extract the store path (v2 uses .path, v3 uses .root)
    if hasattr(base_path, "store"):
        store = base_path.store
        store_path = getattr(store, "path", None) or getattr(store, "root", None)
        if store_path is not None:
            return {"filename": str(store_path), "format": "Zarr"}
        return {"filename": "", "format": "Zarr"}

    return None


def get_lazy_backing_info(obj: Any) -> dict[str, Any]:
    """Get backing file information from a lazy AnnData.

    Extracts the file path and format from the underlying lazy arrays
    (CategoricalArray or MaskedArray) in obs/var columns.

    Parameters
    ----------
    obj
        A lazy AnnData object (from read_lazy())

    Returns
    -------
    Dictionary with:
        - 'filename': str - path to the backing file (empty if not found)
        - 'format': str - 'H5AD' or 'Zarr' (empty if not found)
    """
    empty_result: dict[str, Any] = {"filename": "", "format": ""}

    if not is_lazy_adata(obj):
        return empty_result

    # Try to get path from adata.file (set for H5AD files opened via path)
    file_obj = getattr(obj, "file", None)
    if file_obj is not None:
        filename = getattr(file_obj, "filename", None)
        if filename is not None:
            filename_str = str(filename)
            fmt = "H5AD" if filename_str.endswith(".h5ad") else "Zarr"
            return {"filename": filename_str, "format": fmt}

    # Try to extract from underlying lazy arrays in obs/var
    obs = getattr(obj, "obs", None)
    ds = getattr(obs, "ds", None) if obs is not None and hasattr(obs, "ds") else None
    if ds is None:
        return empty_result

    # Search through columns for a backing array with path info
    try:
        from anndata.experimental.backed._lazy_arrays import (
            CategoricalArray,
            MaskedArray,
        )

        for col_name in ds.data_vars:
            col = ds[col_name]
            # Navigate: DataArray -> Variable -> LazilyIndexedArray -> BackingArray
            if not (hasattr(col, "variable") and hasattr(col.variable, "_data")):
                continue
            lazy_indexed = col.variable._data
            if not hasattr(lazy_indexed, "array"):
                continue
            arr = lazy_indexed.array
            if isinstance(arr, (CategoricalArray, MaskedArray)):
                result = _extract_path_from_lazy_array(arr)
                if result is not None:
                    return result
    except ImportError:
        pass

    return empty_result


def is_lazy_column(series: Any) -> bool:
    """
    Check if a Series-like object is lazy (backed by remote/lazy storage).

    This detects Series from Dataset2D (xarray-backed DataFrames used in
    lazy AnnData) to prevent operations that would trigger data loading.

    Note: We avoid accessing .data as that triggers loading for lazy
    CategoricalArrays. Instead we check for xarray-specific attributes.

    Parameters
    ----------
    series
        The column/series to check

    Returns
    -------
    True if series is a lazy column (xarray DataArray)
    """
    # Check for xarray DataArray structure without triggering data loading
    # xarray DataArrays have 'variable' and 'dims' attributes
    if hasattr(series, "variable") and hasattr(series, "dims"):
        return True
    # Check for xarray Variable backing
    return hasattr(series, "_variable")


def get_lazy_category_count(col: Any) -> int | None:
    """
    Get the number of categories for a lazy categorical without loading them.

    For lazy categoricals, we access the underlying CategoricalArray directly
    and read the category count from the zarr/h5 storage metadata, avoiding
    any data loading.

    Parameters
    ----------
    col
        The lazy categorical column (xarray DataArray)

    Returns
    -------
    Number of categories, or None if cannot be determined
    """
    # Try to get category count from CategoricalArray without loading
    cat_arr = _get_categorical_array(col)
    if cat_arr is not None:
        try:
            # Access the raw _categories group/array shape
            # For zarr: _categories is a Group with 'values' array
            # For h5: similar structure
            cats = cat_arr._categories
            if hasattr(cats, "keys"):  # It's a group
                values = cats["values"]
                return values.shape[0]
            elif hasattr(cats, "shape"):  # It's an array directly
                return cats.shape[0]
        except Exception:  # noqa: BLE001
            pass
    return None


def get_lazy_categorical_info(obj: Any) -> tuple[int | None, bool]:
    """
    Get category count and ordered flag from a lazy categorical without loading data.

    For lazy categoricals (xarray DataArray backed by CategoricalArray),
    this accesses the underlying storage metadata directly to get the count
    without loading the actual category values.

    Parameters
    ----------
    obj
        The object to check (xarray DataArray backed by CategoricalArray)

    Returns
    -------
    tuple of (n_categories, ordered)
        n_categories: Number of categories, or None if cannot be determined
        ordered: Whether the categorical is ordered
    """
    try:
        from anndata.experimental.backed._lazy_arrays import CategoricalArray

        # Navigate through xarray structure to find CategoricalArray
        if hasattr(obj, "variable") and hasattr(obj.variable, "_data"):
            lazy_indexed = obj.variable._data
            if hasattr(lazy_indexed, "array"):
                arr = lazy_indexed.array
                if isinstance(arr, CategoricalArray):
                    # Get count from storage metadata without loading
                    cats = arr._categories
                    if hasattr(cats, "keys"):  # It's a group (zarr)
                        values = cats["values"]  # type: ignore[index]
                        return values.shape[0], arr._ordered  # type: ignore[union-attr]
                    elif hasattr(cats, "shape"):  # It's an array directly
                        return cats.shape[0], arr._ordered  # type: ignore[union-attr]
    except (ImportError, Exception):  # noqa: BLE001
        pass
    return None, False


def get_lazy_categories(
    col: Any, context: FormatterContext
) -> tuple[list, bool, int | None]:
    """
    Get categories for a lazy categorical column, respecting limits.

    For lazy AnnData (from read_lazy()), this accesses the underlying
    CategoricalArray directly and reads only the needed categories from
    storage, avoiding loading the full categorical data.

    Parameters
    ----------
    col
        Column (lazy xarray DataArray) to get categories from
    context
        FormatterContext with max_lazy_categories limit

    Returns
    -------
    tuple of (categories_list, was_truncated, n_categories)
        categories_list: List of category values (empty if skipped)
        was_truncated: True if categories were truncated due to limit
        n_categories: Total number of categories (if known)
    """
    # Import here to avoid circular imports
    from .utils import _get_categories_from_column

    # Try to get category count without loading
    n_cats = get_lazy_category_count(col)

    # If max_lazy_categories is 0, skip loading entirely (metadata-only mode)
    if context.max_lazy_categories == 0:
        return [], True, n_cats

    # Determine if we need to truncate
    should_truncate = n_cats is not None and n_cats > context.max_lazy_categories
    n_to_read = context.max_lazy_categories if should_truncate else n_cats

    # Try to read categories directly from CategoricalArray storage.
    # We access _categories (private) to bypass the @cached_property which loads
    # ALL categories. Instead, we use read_elem_partial (official API) to read
    # only the first N categories. This is intentional - for large categoricals,
    # loading everything defeats the purpose of lazy loading.
    cat_arr = _get_categorical_array(col)
    if cat_arr is not None:
        try:
            from anndata._io.specs.registry import read_elem, read_elem_partial

            cats = cat_arr._categories
            # Get values array: zarr uses group with "values" key, h5 uses array directly
            values = cats["values"] if hasattr(cats, "keys") else cats
            if n_to_read is not None and n_to_read < (n_cats or float("inf")):
                categories = list(
                    read_elem_partial(values, indices=slice(0, n_to_read))
                )
            else:
                categories = list(read_elem(values))  # type: ignore[arg-type]
            return categories, should_truncate, n_cats
        except Exception:  # noqa: BLE001
            pass

    # Fallback to unified accessor (will trigger loading)
    try:
        return _get_categories_from_column(col), False, n_cats
    except Exception:  # noqa: BLE001
        return [], True, n_cats
