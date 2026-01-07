from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from anndata._core.index import _subset
from anndata._core.views import as_view
from anndata._io.specs.lazy_methods import get_chunksize

from ..._settings import settings
from ...compat import (
    H5Array,
    XBackendArray,
    XDataArray,
    XZarrArrayWrapper,
    ZarrArray,
)

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal

    from pandas._libs.missing import NAType
    from pandas.core.dtypes.base import ExtensionDtype

    from anndata.compat import ZarrGroup

    from ...compat import Index1DNorm

    if TYPE_CHECKING:  # Double nesting so Sphinx can import the parent block
        from xarray.core.extension_array import PandasExtensionArray
        from xarray.core.indexing import ExplicitIndexer


class ZarrOrHDF5Wrapper[K: (H5Array, ZarrArray)](XZarrArrayWrapper):
    def __init__(self, array: K) -> None:
        self.chunks = array.chunks
        if isinstance(array, ZarrArray):
            super().__init__(array)
            return
        self._array = array
        self.shape = self._array.shape
        self.dtype = self._array.dtype

    def __getitem__(self, key: ExplicitIndexer):
        from xarray.core.indexing import IndexingSupport, explicit_indexing_adapter

        if isinstance(self._array, ZarrArray):
            return super().__getitem__(key)
        res = explicit_indexing_adapter(
            key, self.shape, IndexingSupport.OUTER_1VECTOR, self._getitem
        )
        return res

    def _getitem(self, key: tuple[int | np.integer | slice | np.ndarray]):
        if not isinstance(key, tuple):
            msg = f"`xr.core.indexing.explicit_indexing_adapter` should have produced a tuple, got {type(key)} instead"
            raise ValueError(msg)
        if (n_key_dims := len(key)) != 1:
            msg = f"Backed arrays currently only supported in 1d, got {n_key_dims} dims"
            raise ValueError(msg)
        key = key[0]
        # See https://github.com/h5py/h5py/issues/293 for why we need to convert.
        # See https://github.com/pydata/xarray/blob/fa03b5b4ae95a366f6de5b60f5cc4eb801cd51ec/xarray/core/indexing.py#L1259-L1263
        # for why we can expect sorted/deduped indexers (which are needed for hdf5).
        if (
            isinstance(key, np.ndarray)
            and np.issubdtype(key.dtype, np.integer)
            and isinstance(self._array, H5Array)
        ):
            key_mask = np.zeros(self._array.shape).astype("bool")
            key_mask[key] = True
            return self._array[key_mask]
        return self._array[key]


class LazyCategories:
    """Lazy accessor for category values supporting efficient slicing.

    Supports partial reads without loading all categories into memory.

    Examples
    --------
    .. code-block:: python

        cats = lazy_adata.obs["cell_type"].cat.categories
        len(cats)  # cheap (metadata only)
        cats[:10]  # partial read (head)
        cats[-5:]  # partial read (tail)
        np.array(cats)  # full load when needed
    """

    def __init__(self, cat_array: CategoricalArray):
        self._cat_array = cat_array

    def __len__(self) -> int:
        """Number of categories (cheap, metadata only)."""
        return self._cat_array.n_categories

    @property
    def shape(self) -> tuple[int]:
        """Shape of categories array (cheap, metadata only)."""
        return (len(self),)

    def __getitem__(self, key: int | slice) -> np.ndarray | str:
        """Get categories by index with efficient partial reads."""
        # If already cached, slice from cache
        if "categories" in self._cat_array.__dict__:
            return self._cat_array.categories[key]

        from anndata._io.specs.registry import read_elem_partial

        n_cats = len(self)
        single = isinstance(key, int)

        # Normalize to slice
        if single:
            if key < -n_cats or key >= n_cats:
                msg = f"index {key} is out of bounds for categories with size {n_cats}"
                raise IndexError(msg)
            idx = key if key >= 0 else n_cats + key
            key = slice(idx, idx + 1)

        start, stop, step = key.indices(n_cats)
        arr = read_elem_partial(
            self._cat_array._categories["values"], indices=slice(start, stop)
        )
        if step != 1:
            arr = arr[::step]
        return arr[0] if single else arr  # scalar for int key, array for slice

    def __iter__(self):
        """Iterate over all categories (loads all data)."""
        return iter(self._get_all())

    def __array__(self, dtype=None, copy=None):
        """Convert to numpy array (loads all data)."""
        arr = np.asarray(self._get_all())
        if dtype is not None:
            return arr.astype(dtype)
        return arr

    def __repr__(self) -> str:
        return f"LazyCategories(n={len(self)})"

    def _get_all(self) -> np.ndarray | pd.api.extensions.ExtensionArray:
        """Get all categories (uses CategoricalArray.categories cached_property)."""
        return self._cat_array.categories


class CategoricalArray[K: (H5Array, ZarrArray)](XBackendArray):
    """
    A wrapper class meant to enable working with lazy categorical data.
    We do not guarantee the stability of this API beyond that guaranteed
    by :class:`xarray.backends.BackendArray`.
    """

    _codes: ZarrOrHDF5Wrapper[K]
    _categories: ZarrArray | H5Array
    shape: tuple[int, ...]
    base_path_or_zarr_group: Path | ZarrGroup
    elem_name: str

    def __init__(
        self,
        codes: K,
        categories: ZarrArray | H5Array,
        base_path_or_zarr_group: Path | ZarrGroup,
        elem_name: str,
        *args,
        ordered: bool,
        **kwargs,
    ):
        self._categories = categories
        self._ordered = ordered
        self._codes = ZarrOrHDF5Wrapper(codes)
        self.shape = self._codes.shape
        self.base_path_or_zarr_group = base_path_or_zarr_group
        self.file_format = "zarr" if isinstance(codes, ZarrArray) else "h5"
        self.elem_name = elem_name

    @property
    def n_categories(self) -> int:
        """Number of categories (cheap, from metadata only)."""
        return self._categories["values"].shape[0]

    @cached_property
    def categories(self) -> np.ndarray | pd.api.extensions.ExtensionArray:
        """All categories (cached, loads all on first access)."""
        if isinstance(self._categories, ZarrArray):
            return self._categories[...]
        from anndata.io import read_elem

        return read_elem(self._categories)

    def __getitem__(self, key: ExplicitIndexer) -> PandasExtensionArray:
        from xarray.core.extension_array import PandasExtensionArray

        codes = self._codes[key]
        categorical_array = pd.Categorical.from_codes(
            codes=codes,
            categories=self.categories,
            ordered=self._ordered,
        )
        if settings.remove_unused_categories:
            categorical_array = categorical_array.remove_unused_categories()
        return PandasExtensionArray(categorical_array)

    @cached_property
    def dtype(self):
        return pd.CategoricalDtype(categories=self.categories, ordered=self._ordered)


# circumvent https://github.com/tox-dev/sphinx-autodoc-typehints/issues/580
type K = H5Array | ZarrArray


class MaskedArray[K: (H5Array, ZarrArray)](XBackendArray):
    """
    A wrapper class meant to enable working with lazy masked data.
    We do not guarantee the stability of this API beyond that guaranteed
    by :class:`xarray.backends.BackendArray`.
    """

    _mask: ZarrOrHDF5Wrapper[K]
    _values: ZarrOrHDF5Wrapper[K]
    _dtype_str: Literal["nullable-integer", "nullable-boolean", "nullable-string-array"]
    shape: tuple[int, ...]
    base_path_or_zarr_group: Path | ZarrGroup
    elem_name: str

    def __init__(
        self,
        values: ZarrArray | H5Array,
        dtype_str: Literal[
            "nullable-integer", "nullable-boolean", "nullable-string-array"
        ],
        mask: ZarrArray | H5Array,
        base_path_or_zarr_group: Path | ZarrGroup,
        elem_name: str,
    ):
        self._mask = ZarrOrHDF5Wrapper(mask)
        self._values = ZarrOrHDF5Wrapper(values)
        self._dtype_str = dtype_str
        self.shape = self._values.shape
        self.base_path_or_zarr_group = base_path_or_zarr_group
        self.file_format = "zarr" if isinstance(mask, ZarrArray) else "h5"
        self.elem_name = elem_name

    def __getitem__(self, key: ExplicitIndexer) -> PandasExtensionArray | np.ndarray:
        from xarray.core.extension_array import PandasExtensionArray

        values = self._values[key]
        mask = self._mask[key]
        if self._dtype_str == "nullable-integer":
            # numpy does not support nan ints
            extension_array = pd.arrays.IntegerArray(values, mask=mask)
        elif self._dtype_str == "nullable-boolean":
            extension_array = pd.arrays.BooleanArray(values, mask=mask)
        elif self._dtype_str == "nullable-string-array":
            # https://github.com/pydata/xarray/issues/10419
            values = values.astype(self.dtype)
            values[mask] = pd.NA
            return values
        else:
            msg = f"Invalid dtype_str {self._dtype_str}"
            raise RuntimeError(msg)
        return PandasExtensionArray(extension_array)

    @cached_property
    def dtype(self) -> np.dtypes.StringDType[NAType] | ExtensionDtype:
        if self._dtype_str == "nullable-integer":
            return pd.array(
                [],
                dtype=str(pd.api.types.pandas_dtype(self._values.dtype)).capitalize(),
            ).dtype
        elif self._dtype_str == "nullable-boolean":
            return pd.BooleanDtype()
        elif self._dtype_str == "nullable-string-array":
            # https://github.com/pydata/xarray/issues/10419
            return np.dtypes.StringDType(na_object=pd.NA)
        msg = f"Invalid dtype_str {self._dtype_str}"
        raise RuntimeError(msg)


@_subset.register(XDataArray)
def _subset_masked(
    a: XDataArray, subset_idx: tuple[Index1DNorm] | tuple[Index1DNorm, Index1DNorm]
):
    return a[subset_idx]


@as_view.register(XDataArray)
def _view_pd_boolean_array(a: XDataArray, view_args):
    return a


@get_chunksize.register(MaskedArray)
def _(a: MaskedArray):
    return get_chunksize(a._values)


@get_chunksize.register(CategoricalArray)
def _(a: CategoricalArray):
    return get_chunksize(a._codes)


def _get_categorical_array(xarray_obj: XDataArray) -> CategoricalArray | None:
    """Extract CategoricalArray from an xarray DataArray if present."""
    try:
        # Navigate: DataArray -> Variable -> LazilyIndexedArray -> CategoricalArray
        lazy_indexed = xarray_obj.variable._data
        if hasattr(lazy_indexed, "array") and isinstance(
            lazy_indexed.array, CategoricalArray
        ):
            return lazy_indexed.array
    except (AttributeError, TypeError):
        pass
    return None


def _register_cat_accessor():
    """Register the cat accessor on xarray DataArray."""
    try:
        import xarray as xr

        @xr.register_dataarray_accessor("cat")
        class CatAccessor:
            """Accessor for categorical operations on xarray DataArrays.

            Provides efficient access to category information. For lazy categorical
            columns backed by CategoricalArray, uses partial reads to avoid loading
            all data. For other categorical columns, falls back to dtype access.
            Returns None for non-categorical columns.

            Examples
            --------
            .. code-block:: python

                lazy_adata = ad.experimental.read_lazy("dataset.zarr")
                cats = lazy_adata.obs["cell_type"].cat.categories
                len(cats)  # cheap, metadata only
                cats[:5]  # efficient partial read
                cats[-3:]  # efficient partial read
                lazy_adata.obs["numeric_col"].cat.categories  # returns None
            """

            def __init__(self, xarray_obj: XDataArray):
                self._obj = xarray_obj
                self._cat_array = _get_categorical_array(xarray_obj)

            def _is_categorical(self) -> bool:
                """Check if the underlying data is categorical."""
                return isinstance(self._obj.dtype, pd.CategoricalDtype)

            @property
            def categories(self) -> LazyCategories | pd.Index | None:
                """Categories with support for efficient slicing.

                For lazy CategoricalArray, returns LazyCategories supporting
                efficient partial reads via slicing (e.g., categories[:10]).
                For other categoricals, returns the pandas Index from dtype.
                Returns None for non-categorical columns.

                Examples
                --------
                .. code-block:: python

                    cats = col.cat.categories
                    len(cats)  # cheap for lazy data
                    cats[:10]  # efficient partial read
                    cats[-5:]  # efficient partial read
                    np.array(cats)  # full load
                """
                if self._cat_array is not None:
                    return LazyCategories(self._cat_array)
                if self._is_categorical():
                    return self._obj.dtype.categories
                return None

            @property
            def codes(self) -> ZarrArray | H5Array | None:
                """Integer codes for the categorical values.

                For lazy CategoricalArray, returns the underlying zarr/h5 array
                containing the integer codes. Supports lazy slicing.
                For non-categorical columns, returns None.
                """
                if self._cat_array is not None:
                    return self._cat_array._codes._array
                return None

            @property
            def ordered(self) -> bool | None:
                """Whether the categorical is ordered.

                Returns True/False for categorical columns, None for non-categorical.
                """
                if self._cat_array is not None:
                    return bool(self._cat_array._ordered)
                if self._is_categorical():
                    return bool(self._obj.dtype.ordered)
                return None

        return CatAccessor
    except ImportError:
        return None


# Register the accessor when xarray is available
_CatAccessor = _register_cat_accessor()
