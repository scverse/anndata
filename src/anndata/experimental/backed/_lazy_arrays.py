from __future__ import annotations

from functools import singledispatchmethod
from typing import Generic, TypeVar, Union

import pandas as pd
import xarray as xr

from anndata._core.index import Index, _subset
from anndata._core.views import as_view
from anndata.compat import H5Array, ZarrArray

K = TypeVar("K", bound=Union[H5Array, ZarrArray])


class ZarrOrHDF5Wrapper(xr.backends.zarr.ZarrArrayWrapper, Generic[K]):
    @singledispatchmethod  # type: ignore
    def __init__(self, array: ZarrArray):
        return super().__init__(array)

    @__init__.register
    def _(self, array: H5Array):
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


class CategoricalArray(xr.backends.BackendArray):
    def __init__(
        self,
        codes: ZarrArray | H5Array,
        categories: ZarrArray | H5Array,
        ordered: bool,
        drop_unused_cats: bool = False,
        *args,
        **kwargs,
    ):
        self._categories = categories
        self._ordered = ordered
        self._drop_unused_cats = drop_unused_cats
        self._categories_cache = None
        self._codes = ZarrOrHDF5Wrapper[type(codes)](codes)
        self.shape = self._codes.shape
        self.dtype = pd.CategoricalDtype(
            categories=self._categories, ordered=self._ordered
        )

    @property
    def categories(self):  # __slots__ and cached_property are incompatible
        if self._categories_cache is None:
            if isinstance(self._categories, ZarrArray):
                self._categories_cache = self._categories[...]
            else:
                if (
                    "read_dataset" not in dir()
                ):  # avoid circular dependency, not sure what caused this all of a sudden after merging https://github.com/scverse/anndata/pull/949/commits/dc9f12fcbca977841e967c8414b9f1032e069250
                    from ..._io.h5ad import read_dataset
                self._categories_cache = read_dataset(self._categories)
        return self._categories_cache

    def __getitem__(
        self, key: xr.core.indexing.ExplicitIndexer
    ) -> xr.core.extension_array.PandasExtensionArray:
        codes = self._codes[key]
        categorical_array = pd.Categorical.from_codes(
            codes=codes, categories=self.categories, ordered=self._ordered
        )
        if self._drop_unused_cats:
            return xr.core.extension_array.PandasExtensionArray(
                categorical_array.remove_unused_categories()
            )
        return xr.core.extension_array.PandasExtensionArray(categorical_array)


class MaskedArray(xr.backends.BackendArray):
    def __init__(
        self,
        values: ZarrArray | H5Array,
        dtype_str: str,
        mask: ZarrArray | H5Array | None = None,
    ):
        self._mask = ZarrOrHDF5Wrapper[type(mask)](mask)
        self._values = ZarrOrHDF5Wrapper[type(values)](values)
        self._dtype_str = dtype_str
        self.shape = self._values.shape
        self.dtype = pd.api.types.pandas_dtype(self._values.dtype)

    def __getitem__(self, key) -> xr.core.extension_array.PandasExtensionArray:
        # HACK! TODO(ilan-gold): open issue about hdf5 compat that doesn't allow initialization!
        values = self._values[key]
        if self._mask is not None:
            mask = self._mask[key]
            if self._dtype_str == "nullable-integer":
                # numpy does not support nan ints
                return xr.core.extension_array.PandasExtensionArray(
                    pd.arrays.IntegerArray(values, mask=mask)
                )
            elif self._dtype_str == "nullable-boolean":
                return xr.core.extension_array.PandasExtensionArray(
                    pd.arrays.BooleanArray(values, mask=mask)
                )
        return xr.core.extension_array.PandasExtensionArray(pd.array(values))


@_subset.register(xr.DataArray)
def _subset_masked(a: xr.DataArray, subset_idx: Index):
    return a[subset_idx]


@as_view.register(xr.DataArray)
def _view_pd_boolean_array(a: xr.DataArray, view_args):
    return a
