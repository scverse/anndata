"""
Lazy loading utilities for AnnData HTML representation.

This module consolidates all logic related to detecting and handling lazy AnnData
objects (from read_lazy()). Lazy AnnData uses xarray-backed storage and requires
special handling to avoid triggering data loading during repr generation.

Key concepts:
- Lazy AnnData: Created by read_lazy(), obs/var are Dataset2D (xarray-backed)
- Lazy series: Individual columns from Dataset2D, implemented as xarray DataArrays
- CategoricalArray: anndata's lazy categorical implementation for zarr/h5 storage

``CategoricalArray.categories``/``.dtype`` are cached properties that load *all*
categories, so this module reads the category count and the first few labels
from storage directly instead.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from .._core.xarray import Dataset2D
from ..experimental.backed._lazy_arrays import CategoricalArray, MaskedArray

if TYPE_CHECKING:
    from .registry import FormatterContext


def _get_backing_array(col: object) -> object | None:
    """Get the storage-backed array behind a lazy column without loading it.

    Navigates DataArray -> Variable -> LazilyIndexedArray -> backing array.
    """
    try:
        return col.variable._data.array  # type: ignore[attr-defined]
    except Exception:  # noqa: BLE001
        # Not a lazy column (AttributeError), or a broken object
        return None


def _get_categorical_array(col: object) -> CategoricalArray | None:
    """Get the underlying CategoricalArray of a lazy column, if it is one."""
    arr = _get_backing_array(col)
    return arr if isinstance(arr, CategoricalArray) else None


def _category_values(cat_arr: CategoricalArray) -> object:
    """Get the on-disk array of category labels (zarr groups nest it under "values")."""
    cats = cat_arr._categories
    return cats["values"] if hasattr(cats, "keys") else cats


def is_lazy_adata(obj: object) -> bool:
    """Check if an AnnData uses lazy loading (experimental read_lazy).

    Lazy AnnData has Dataset2D (xarray-backed) obs/var instead of regular DataFrames.

    Parameters
    ----------
    obj
        Object to check (typically an AnnData)

    Returns
    -------
    True if obj is a lazy AnnData

    Notes
    -----
    This function accesses the .obs attribute which may trigger I/O for some
    objects. If .obs raises an exception, returns False.
    """
    try:
        return isinstance(getattr(obj, "obs", None), Dataset2D)
    except Exception:  # noqa: BLE001
        # Intentional broad catch: .obs access may raise anything
        return False


def _extract_path_from_lazy_array(
    arr: CategoricalArray | MaskedArray,
) -> dict[str, str] | None:
    """Extract file path and format from a lazy array (CategoricalArray/MaskedArray)."""
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


def get_lazy_backing_info(obj: object) -> dict[str, str]:
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
    empty_result: dict[str, str] = {"filename": "", "format": ""}

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
    for col_name in ds.data_vars:
        arr = _get_backing_array(ds[col_name])
        if isinstance(arr, CategoricalArray | MaskedArray):
            result = _extract_path_from_lazy_array(arr)
            if result is not None:
                return result

    return empty_result


def is_lazy_column(series: object) -> bool:
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


def get_lazy_categorical_info(obj: object) -> tuple[int | None, bool]:
    """
    Get category count and ordered flag from a lazy categorical without loading data.

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
    cat_arr = _get_categorical_array(obj)
    if cat_arr is None:
        return None, False
    try:
        return _category_values(cat_arr).shape[0], cat_arr._ordered  # type: ignore[attr-defined]
    except Exception:  # noqa: BLE001
        return None, False


def get_lazy_category_count(col: object) -> int | None:
    """Get the number of categories of a lazy categorical without loading them."""
    return get_lazy_categorical_info(col)[0]


def get_lazy_categories(
    col: object, context: FormatterContext
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

            values = _category_values(cat_arr)
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
