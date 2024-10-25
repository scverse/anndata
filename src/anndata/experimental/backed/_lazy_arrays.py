from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING, Generic, TypeVar

import pandas as pd

from anndata._core.index import _subset
from anndata._core.views import as_view
from anndata.compat import H5Array, ZarrArray

from ..._settings import settings
from ._compat import BackendArray, DataArray, ZarrArrayWrapper
from ._compat import xarray as xr

if TYPE_CHECKING:
    from typing import Literal

    from anndata._core.index import Index


K = TypeVar("K", H5Array, ZarrArray)


class ZarrOrHDF5Wrapper(ZarrArrayWrapper, Generic[K]):
    def __init__(self, array: K):
        if isinstance(array, ZarrArray):
            return super().__init__(array)
        if set(self.__slots__) != {"dtype", "shape", "_array"}:
            msg = "Expected attributes of xarray ZarrArrayWrapper have changed - "
            "please file an issue with anndata and consider downgrading xarray"
            raise AssertionError(msg)
        self._array = array
        self.shape = self._array.shape
        self.dtype = self._array.dtype

    def __getitem__(self, key):
        if isinstance(self._array, ZarrArray):
            return super().__getitem__(key)
        # adapted from https://github.com/pydata/xarray/blob/main/xarray/backends/h5netcdf_.py#L50-L58C13
        # TODO: locks?
        return xr.core.indexing.explicit_indexing_adapter(
            key,
            self.shape,
            xr.core.indexing.IndexingSupport.OUTER_1VECTOR,
            lambda key: self._array[key],
        )


class CategoricalArray(BackendArray, Generic[K]):
    _codes: ZarrOrHDF5Wrapper[K]
    _categories: ZarrArray | H5Array
    shape: tuple[int, ...]

    def __init__(
        self,
        codes: K,
        categories: ZarrArray | H5Array,
        ordered: bool,
        *args,
        **kwargs,
    ):
        self._categories = categories
        self._ordered = ordered
        self._codes = ZarrOrHDF5Wrapper(codes)
        self.shape = self._codes.shape

    @cached_property
    def categories(self):
        if isinstance(self._categories, ZarrArray):
            return self._categories[...]
        from ..._io.h5ad import read_dataset

        return read_dataset(self._categories)

    def __getitem__(
        self, key: xr.core.indexing.ExplicitIndexer
    ) -> xr.core.extension_array.PandasExtensionArray:
        codes = self._codes[key]
        categorical_array = pd.Categorical.from_codes(
            codes=codes, categories=self.categories, ordered=self._ordered
        )
        if settings.remove_unused_categories:
            categorical_array = categorical_array.remove_unused_categories()
        return xr.core.extension_array.PandasExtensionArray(categorical_array)

    @cached_property
    def dtype(self):
        return pd.CategoricalDtype(categories=self._categories, ordered=self._ordered)


class MaskedArray(BackendArray, Generic[K]):
    _mask: ZarrOrHDF5Wrapper[K]
    _values: ZarrOrHDF5Wrapper[K]
    _dtype_str: Literal["nullable-integer", "nullable-boolean"]
    shape: tuple[int, ...]

    def __init__(
        self,
        values: ZarrArray | H5Array,
        dtype_str: Literal["nullable-integer", "nullable-boolean"],
        mask: ZarrArray | H5Array | None = None,
    ):
        self._mask = ZarrOrHDF5Wrapper(mask)
        self._values = ZarrOrHDF5Wrapper(values)
        self._dtype_str = dtype_str
        self.shape = self._values.shape

    def __getitem__(self, key) -> xr.core.extension_array.PandasExtensionArray:
        values = self._values[key]
        if self._mask is not None:
            mask = self._mask[key]
            if self._dtype_str == "nullable-integer":
                # numpy does not support nan ints
                extension_array = pd.arrays.IntegerArray(values, mask=mask)
            elif self._dtype_str == "nullable-boolean":
                extension_array = pd.arrays.BooleanArray(values, mask=mask)
            else:
                raise ValueError(f"Invalid dtype_str {self._dtype_str}")
            return xr.core.extension_array.PandasExtensionArray(extension_array)
        return xr.core.extension_array.PandasExtensionArray(pd.array(values))

    @cached_property
    def dtype(self):
        if self._dtype_str == "nullable-integer":
            return pd.array(
                [],
                dtype=str(pd.api.types.pandas_dtype(self._values.dtype)).capitalize(),
            ).dtype
        elif self._dtype_str == "nullable-boolean":
            return pd.BooleanDtype()
        else:
            raise ValueError(f"Invalid dtype_str {self._dtype_str}")


@_subset.register(DataArray)
def _subset_masked(a: DataArray, subset_idx: Index):
    return a[subset_idx]


@as_view.register(DataArray)
def _view_pd_boolean_array(a: DataArray, view_args):
    return a
