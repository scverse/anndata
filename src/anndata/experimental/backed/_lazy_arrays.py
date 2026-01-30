from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from anndata._core.index import _subset
from anndata._core.views import as_view
from anndata._io.specs.lazy_methods import get_chunksize

from ..._io.utils import pandas_nullable_dtype
from ..._settings import settings
from ...compat import (
    H5Array,
    H5AsTypeView,
    XBackendArray,
    XDataArray,
    XZarrArrayWrapper,
    ZarrArray,
)

# Number of categories to show at head/tail in LazyCategoricalDtype repr
_N_CATEGORIES_REPR_SHOW = 10

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal

    from numpy.typing import NDArray
    from pandas._libs.missing import NAType
    from pandas.core.dtypes.dtypes import BaseMaskedDtype

    from anndata.compat import H5Group, ZarrGroup

    from ...compat import Index1DNorm

    if TYPE_CHECKING:  # Double nesting so Sphinx can import the parent block
        from xarray.core.extension_array import PandasExtensionArray
        from xarray.core.indexing import ExplicitIndexer
else:  # https://github.com/tox-dev/sphinx-autodoc-typehints/issues/580
    type K = H5Array | ZarrArray


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
        encoding) containing category values.
    ordered
        Whether the categorical is ordered.
    """

    # Attributes that should be preserved during copying/pickling
    _metadata = ("_categories_elem", "_ordered_flag")

    def __new__(
        cls,
        categories_elem: ZarrArray | H5Array | ZarrGroup | H5Group,
        *,
        ordered: bool = False,
    ):
        # Create instance without calling parent __init__ with categories
        instance = object.__new__(cls)
        return instance

    def __init__(
        self,
        categories_elem: ZarrArray | H5Array | ZarrGroup | H5Group,
        *,
        ordered: bool = False,
    ):
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
    def categories(self) -> pd.Index:
        """Categories index. Loads all categories on first access and caches."""
        from anndata.io import read_elem

        return pd.Index(read_elem(self._categories_elem))

    @property
    def ordered(self) -> bool:
        """Whether the categorical is ordered."""
        return self._ordered_flag

    @property
    def n_categories(self) -> int:
        """Number of categories (cheap, metadata only)."""
        if "categories" in self.__dict__:
            return len(self.categories)
        return self._get_categories_array().shape[0]

    def _get_categories_slice(
        self, n: int, *, from_end: bool = False
    ) -> np.ndarray | pd.api.extensions.ExtensionArray:
        """Get n categories from start or end.

        Parameters
        ----------
        n
            Number of categories to return.
        from_end
            If True, return last n categories. If False, return first n.

        Returns
        -------
        np.ndarray or ExtensionArray
            The requested categories.
        """
        # If already fully loaded, slice from cache
        if "categories" in self.__dict__:
            sliced = self.categories[-n:] if from_end else self.categories[:n]
            return np.asarray(sliced)

        # Read partial from disk
        from anndata._io.specs.registry import read_elem_partial

        arr = self._get_categories_array()
        total = arr.shape[0]
        if from_end:
            start, stop = max(total - n, 0), total
        else:
            start, stop = 0, min(n, total)
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
        return self._get_categories_slice(n, from_end=False)

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
        return self._get_categories_slice(n, from_end=True)

    def __repr__(self) -> str:
        n_total = self.n_categories
        ordered_str = ", ordered=True" if self.ordered else ""

        if n_total <= _N_CATEGORIES_REPR_SHOW * 2:
            # Small enough to show all categories
            if "categories" in self.__dict__:
                cats = list(self.categories)
            else:
                cats = list(self.head_categories(n_total))
            return f"LazyCategoricalDtype(categories={cats!r}{ordered_str})"

        # Show truncated: first n ... last n
        head = list(self.head_categories(_N_CATEGORIES_REPR_SHOW))
        tail = list(self.tail_categories(_N_CATEGORIES_REPR_SHOW))
        cats_display = [*head, "...", *tail]
        return f"LazyCategoricalDtype(categories={cats_display!r}, n={n_total}{ordered_str})"

    def __hash__(self) -> int:
        """Hash based on identity of underlying array and ordered flag.

        Required for use in sets and as dictionary keys (e.g., collecting
        unique dtypes across AnnData objects).
        """
        return hash((id(self._categories_elem), self._ordered_flag))

    def __eq__(self, other) -> bool:
        if isinstance(other, LazyCategoricalDtype):
            has_same_ordering = self._ordered_flag == other._ordered_flag
            are_arrays_equal = (self._categories_elem is other._categories_elem) or (
                self._get_categories_array() == other._get_categories_array()
            )
            return has_same_ordering and are_arrays_equal
        # Defer to pandas base implementation for all other comparisons
        # This handles string comparison ("category"), CategoricalDtype comparisons,
        # and all edge cases (None categories, ordered vs unordered, etc.)
        return super().__eq__(other)


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

    def __getitem__(
        self, key: ExplicitIndexer
    ) -> PandasExtensionArray | NDArray[np.str_]:
        values = self._values[key]
        mask = self._mask[key]

        if isinstance(self.dtype, np.dtypes.StringDType):
            # https://github.com/pydata/xarray/issues/10419
            values = values.astype(self.dtype)
            values[mask] = pd.NA
            return values

        from xarray.core.extension_array import PandasExtensionArray

        cls = self.dtype.construct_array_type()
        return PandasExtensionArray(cls(values, mask))

    @cached_property
    def dtype(self) -> BaseMaskedDtype | np.dtypes.StringDType[NAType]:
        if self._dtype_str == "nullable-string-array":
            # https://github.com/pydata/xarray/issues/10419
            return np.dtypes.StringDType(na_object=pd.NA)
        try:
            return pandas_nullable_dtype(self._values.dtype)
        except NotImplementedError:
            msg = f"Invalid dtype_str {self._dtype_str}"
            raise RuntimeError(msg) from None


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
