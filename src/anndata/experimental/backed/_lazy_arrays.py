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
    H5AsTypeView,
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

    from anndata.compat import H5Group, ZarrGroup

    from ...compat import Index1DNorm

    if TYPE_CHECKING:  # Double nesting so Sphinx can import the parent block
        from xarray.core.extension_array import PandasExtensionArray
        from xarray.core.indexing import ExplicitIndexer


class LazyCategoricalDtype(pd.CategoricalDtype):
    """A CategoricalDtype that lazily loads categories from zarr/h5 storage.

    This dtype provides efficient access to categorical metadata without loading
    all categories into memory via :meth:`head_categories`, :meth:`tail_categories`,
    and :attr:`n_categories`. Accessing :attr:`categories` will load all categories
    into memory.

    Parameters
    ----------
    categories_elem
        The underlying zarr or h5 array (or group for nullable-string-array
        encoding) containing category values. Can be None for empty dtype.
    ordered
        Whether the categorical is ordered.
    """

    # Attributes that should be preserved during copying/pickling
    _metadata = ("_categories_elem", "_ordered_flag")

    def __new__(
        cls,
        categories_elem: ZarrArray | H5Array | ZarrGroup | H5Group | None = None,
        *,
        ordered: bool = False,
    ):
        # Create instance without calling parent __init__ with categories
        instance = object.__new__(cls)
        return instance

    def __init__(
        self,
        categories_elem: ZarrArray | H5Array | ZarrGroup | H5Group | None = None,
        *,
        ordered: bool = False,
    ):
        # Can be None for edge cases (empty dtype). See test_lazy_categorical_dtype_empty_array.
        self._categories_elem = categories_elem
        self._ordered_flag = bool(ordered)

    def _get_categories_array(self) -> ZarrArray | H5Array:
        """Get the underlying categories array.

        For string-array encoding: _categories_elem is directly the array.
        For nullable-string-array encoding: _categories_elem would be a Group
        with "values" key (not currently used for categories in anndata, but
        handled defensively).
        """
        if isinstance(self._categories_elem, (ZarrArray, H5Array)):
            return self._categories_elem
        # nullable-string-array encoding: Group with "values" and "mask"
        return self._categories_elem["values"]

    @cached_property
    def categories(self) -> pd.Index | None:
        """Categories index. Loads all categories on first access and caches."""
        if self._categories_elem is None:
            return None
        arr = self._get_categories_array()
        if isinstance(arr, ZarrArray):
            values = arr[...]
        else:
            from anndata.io import read_elem

            values = read_elem(self._categories_elem)
        return pd.Index(values)

    @property
    def ordered(self) -> bool:
        """Whether the categorical is ordered."""
        return self._ordered_flag

    @property
    def n_categories(self) -> int:
        """Number of categories (cheap, metadata only)."""
        if self._categories_elem is None:
            return 0
        if "categories" in self.__dict__:
            return len(self.categories)
        return self._get_categories_array().shape[0]

    def _read_partial_categories(
        self, start: int, stop: int
    ) -> np.ndarray | pd.api.extensions.ExtensionArray:
        """Read a slice of categories from disk.

        Uses read_elem_partial for proper HDF5 string decoding.
        """
        from anndata._io.specs.registry import read_elem_partial

        arr = self._get_categories_array()
        return read_elem_partial(arr, indices=slice(start, stop))

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
        """
        if self._categories_elem is None:
            return np.array([])

        # If already fully loaded, slice from cache
        if "categories" in self.__dict__:
            return np.asarray(self.categories[:n])

        total = self.n_categories
        return self._read_partial_categories(0, min(n, total))

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
        """
        if self._categories_elem is None:
            return np.array([])

        # If already fully loaded, slice from cache
        if "categories" in self.__dict__:
            return np.asarray(self.categories[-n:])

        total = self.n_categories
        start = max(total - n, 0)
        return self._read_partial_categories(start, total)

    def __repr__(self) -> str:
        if "categories" in self.__dict__ and self.categories is not None:
            # Fully loaded - use standard repr
            return f"CategoricalDtype(categories={self.categories!r}, ordered={self.ordered})"
        return f"LazyCategoricalDtype(n_categories={self.n_categories}, ordered={self.ordered})"

    @property
    def name(self) -> str:
        """String identifier for this dtype.

        Required for string comparison (e.g., dtype == "category") used in
        anndata merge operations.
        """
        return "category"

    def __hash__(self) -> int:
        """Hash based on identity of underlying array and ordered flag.

        Required for use in sets and as dictionary keys (e.g., collecting
        unique dtypes across AnnData objects).
        """
        return hash((id(self._categories_elem), self._ordered_flag))

    def __eq__(self, other) -> bool:
        # Handle string comparison (e.g., dtype == "category")
        if isinstance(other, str):
            return other == self.name
        if isinstance(other, LazyCategoricalDtype):
            return (
                self._categories_elem is other._categories_elem
                and self._ordered_flag == other._ordered_flag
            )
        if not isinstance(other, pd.CategoricalDtype):
            return False
        # Compare with regular CategoricalDtype - need to load categories
        if self.ordered != other.ordered:
            return False
        if other.categories is None or self.categories is None:
            return other.categories is None and self.categories is None
        return self.categories.equals(other.categories)


class ZarrOrHDF5Wrapper[K: (H5Array | H5AsTypeView, ZarrArray)](XZarrArrayWrapper):
    def __init__(self, array: K) -> None:
        # AstypeView from h5py .astype() lacks chunks attribute
        self.chunks = getattr(array, "chunks", None)
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
            and isinstance(self._array, H5Array | H5AsTypeView)
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
    _categories_elem: K
    shape: tuple[int, ...]
    base_path_or_zarr_group: Path | ZarrGroup
    elem_name: str

    def __init__(
        self,
        codes: K,
        categories: K,
        base_path_or_zarr_group: Path | ZarrGroup,
        elem_name: str,
        *args,
        ordered: bool,
        **kwargs,
    ):
        self._categories_elem = categories
        self._ordered = ordered
        self._codes = ZarrOrHDF5Wrapper(codes)
        self.shape = self._codes.shape
        self.base_path_or_zarr_group = base_path_or_zarr_group
        self.file_format = "zarr" if isinstance(codes, ZarrArray) else "h5"
        self.elem_name = elem_name
        # Create the lazy dtype - this is where categories are cached
        self._lazy_dtype = LazyCategoricalDtype(
            categories_elem=categories, ordered=ordered
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
type K = H5Array | H5AsTypeView | ZarrArray


class MaskedArray[K: (H5Array | H5AsTypeView, ZarrArray)](XBackendArray):
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
        values: K,
        dtype_str: Literal[
            "nullable-integer", "nullable-boolean", "nullable-string-array"
        ],
        mask: K,
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
