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


class LazyCategoricalDtype(pd.CategoricalDtype):
    """A CategoricalDtype that lazily loads categories from zarr/h5 storage.

    This dtype provides efficient access to categorical metadata without loading
    all categories into memory. Categories are loaded lazily on first full access
    and cached thereafter.

    Parameters
    ----------
    categories_array
        The underlying zarr or h5 array containing category values.
    ordered
        Whether the categorical is ordered.

    Examples
    --------
    >>> dtype = col.dtype  # LazyCategoricalDtype
    >>> dtype.n_categories  # cheap, metadata only
    100000
    >>> dtype.head_categories(5)  # partial read, first 5
    array(['cat_0', 'cat_1', 'cat_2', 'cat_3', 'cat_4'], dtype='<U6')
    >>> dtype.tail_categories(3)  # partial read, last 3
    array(['cat_99997', 'cat_99998', 'cat_99999'], dtype='<U6')
    >>> dtype.categories  # full load, cached
    Index(['cat_0', 'cat_1', ...], dtype='str')
    """

    # Attributes that should be preserved during copying/pickling
    _metadata = ("_categories_array", "_ordered_flag", "_cached_n_categories")

    def __new__(
        cls,
        categories_array: ZarrArray | H5Array | None = None,
        *,
        ordered: bool = False,
    ):
        # Create instance without calling parent __init__ with categories
        instance = object.__new__(cls)
        return instance

    def __init__(
        self,
        categories_array: ZarrArray | H5Array | None = None,
        *,
        ordered: bool = False,
    ):
        self._categories_array = categories_array
        self._ordered_flag = bool(ordered)
        self._cached_n_categories: int | None = None
        self.__categories: pd.Index | None = (
            None  # Double underscore to avoid conflicts
        )

    def _get_categories_array(self) -> ZarrArray | H5Array:
        """Get the underlying categories array (handles both encodings).

        For string-array encoding: _categories_array is directly the array.
        For nullable-string-array encoding: _categories_array is a Group with "values" key.
        """
        if isinstance(self._categories_array, (ZarrArray, H5Array)):
            return self._categories_array
        # nullable-string-array encoding: Group with "values" and "mask"
        return self._categories_array["values"]

    @property
    def categories(self) -> pd.Index | None:
        """Categories index. Loads all categories on first access and caches."""
        if self.__categories is None and self._categories_array is not None:
            arr = self._get_categories_array()
            if isinstance(arr, ZarrArray):
                values = arr[...]
            else:
                from anndata.io import read_elem

                values = read_elem(self._categories_array)
            self.__categories = pd.Index(values)
        return self.__categories

    @property
    def ordered(self) -> bool:
        """Whether the categorical is ordered."""
        return self._ordered_flag

    @property
    def n_categories(self) -> int:
        """Number of categories (cheap, metadata only)."""
        if self._cached_n_categories is not None:
            return self._cached_n_categories
        if self.__categories is not None:
            return len(self.__categories)
        if self._categories_array is not None:
            n = self._get_categories_array().shape[0]
            self._cached_n_categories = n
            return n
        return 0

    def head_categories(
        self, n: int = 5
    ) -> np.ndarray | pd.api.extensions.ExtensionArray:
        """Return first n categories without loading all into memory.

        Parameters
        ----------
        n
            Number of categories to return. Default 5.

        Returns
        -------
        np.ndarray or ExtensionArray
            The first n categories.

        Examples
        --------
        >>> dtype.head_categories()  # first 5
        >>> dtype.head_categories(10)  # first 10
        """
        # If already fully loaded, slice from cache
        if self.__categories is not None:
            return np.asarray(self.__categories[:n])

        if self._categories_array is None:
            return np.array([])

        from anndata._io.specs.registry import read_elem_partial

        arr = self._get_categories_array()
        total = self.n_categories
        return read_elem_partial(arr, indices=slice(0, min(n, total)))

    def tail_categories(
        self, n: int = 5
    ) -> np.ndarray | pd.api.extensions.ExtensionArray:
        """Return last n categories without loading all into memory.

        Parameters
        ----------
        n
            Number of categories to return. Default 5.

        Returns
        -------
        np.ndarray or ExtensionArray
            The last n categories.

        Examples
        --------
        >>> dtype.tail_categories()  # last 5
        >>> dtype.tail_categories(10)  # last 10
        """
        # If already fully loaded, slice from cache
        if self.__categories is not None:
            return np.asarray(self.__categories[-n:])

        if self._categories_array is None:
            return np.array([])

        from anndata._io.specs.registry import read_elem_partial

        arr = self._get_categories_array()
        total = self.n_categories
        start = max(total - n, 0)
        return read_elem_partial(arr, indices=slice(start, total))

    def __repr__(self) -> str:
        if self.__categories is not None:
            # Fully loaded - use standard repr
            return f"CategoricalDtype(categories={self.__categories!r}, ordered={self.ordered})"
        return f"LazyCategoricalDtype(n_categories={self.n_categories}, ordered={self.ordered})"

    @property
    def name(self) -> str:
        """String identifier for this dtype."""
        return "category"

    def __hash__(self) -> int:
        # Need to be hashable for pandas internals
        return hash((id(self._categories_array), self._ordered_flag))

    def __eq__(self, other) -> bool:
        if isinstance(other, LazyCategoricalDtype):
            return (
                self._categories_array is other._categories_array
                and self._ordered_flag == other._ordered_flag
            )
        if isinstance(other, pd.CategoricalDtype):
            # Compare with regular CategoricalDtype - need to load categories
            if self.ordered != other.ordered:
                return False
            if other.categories is None:
                return self.categories is None
            if self.categories is None:
                return False
            return self.categories.equals(other.categories)
        return False


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


class CategoricalArray[K: (H5Array, ZarrArray)](XBackendArray):
    """
    A wrapper class meant to enable working with lazy categorical data.
    We do not guarantee the stability of this API beyond that guaranteed
    by :class:`xarray.backends.BackendArray`.
    """

    _codes: ZarrOrHDF5Wrapper[K]
    _categories_array: ZarrArray | H5Array
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
        self._categories_array = categories
        self._ordered = ordered
        self._codes = ZarrOrHDF5Wrapper(codes)
        self.shape = self._codes.shape
        self.base_path_or_zarr_group = base_path_or_zarr_group
        self.file_format = "zarr" if isinstance(codes, ZarrArray) else "h5"
        self.elem_name = elem_name
        # Create the lazy dtype - this is where categories are cached
        self._lazy_dtype = LazyCategoricalDtype(
            categories_array=categories, ordered=ordered
        )

    @property
    def categories(self) -> pd.Index | None:
        """All categories. Loads and caches on first access."""
        return self._lazy_dtype.categories

    def __getitem__(self, key: ExplicitIndexer) -> PandasExtensionArray:
        from xarray.core.extension_array import PandasExtensionArray

        codes = self._codes[key]
        categorical_array = pd.Categorical.from_codes(
            codes=codes, categories=self.categories, ordered=self._ordered
        )
        if settings.remove_unused_categories:
            categorical_array = categorical_array.remove_unused_categories()
        return PandasExtensionArray(categorical_array)

    @property
    def dtype(self) -> LazyCategoricalDtype:
        """The dtype with lazy category loading support."""
        return self._lazy_dtype


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
