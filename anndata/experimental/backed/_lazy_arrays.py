from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd
import xarray as xr

from anndata._core.index import Index, _subset
from anndata._core.views import as_view
from anndata.compat import H5Array, ZarrArray

if TYPE_CHECKING:
    import numpy as np


class CategoricalArray(
    xr.backends.zarr.ZarrArrayWrapper
):  # Zarr works for hdf5, xarray only supports integration hdf5 in the netcdf context
    def __init__(
        self,
        codes: ZarrArray | H5Array,
        categories: ZarrArray | H5Array,
        ordered: bool,
        drop_unused_cats: bool,
        *args,
        **kwargs,
    ):
        self._categories = categories
        self._ordered = ordered
        self._drop_unused_cats = drop_unused_cats
        self._categories_cache = None
        self._array = codes
        self.shape = self._array.shape
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

    def __getitem__(self, key: xr.core.indexing.ExplicitIndexer) -> np.typing.ArrayLike:
        codes = super().__getitem__(key)
        categorical_array = pd.Categorical.from_codes(
            codes=codes, categories=self.categories, ordered=self._ordered
        )
        if self._drop_unused_cats:
            return xr.core.indexing.ExtensionDuckArray(
                categorical_array.remove_unused_categories()
            )
        return xr.core.indexing.ExtensionDuckArray(categorical_array)


class MaskedArray(xr.backends.zarr.ZarrArrayWrapper):
    def __init__(
        self,
        values: ZarrArray | H5Array,
        dtype_str: str,
        *args,
        mask: ZarrArray | H5Array | None = None,
        **kwargs,
    ):
        self._mask = mask
        self._values = values
        self._dtype_str = dtype_str
        self._array = values
        self.shape = self._array.shape
        self.dtype = pd.api.types.pandas_dtype(self._array.dtype)

    def __getitem__(self, key):
        # HACK! TODO(ilan-gold): open issue about hdf5 compat that doesn't allow initialization!
        self._array = self._values
        values = super().__getitem__(key)
        if self._mask is not None:
            self._array = self._mask
            mask = super().__getitem__(key)
            if self._dtype_str == "nullable-integer":
                # numpy does not support nan ints
                return xr.core.indexing.ExtensionDuckArray(
                    pd.arrays.IntegerArray(values, mask=mask)
                )
            elif self._dtype_str == "nullable-boolean":
                return xr.core.indexing.ExtensionDuckArray(
                    pd.arrays.BooleanArray(values, mask=mask)
                )
        return xr.core.indexing.ExtensionDuckArray(pd.array(values))


@_subset.register(xr.DataArray)
def _subset_masked(a: xr.DataArray, subset_idx: Index):
    return a[subset_idx]


@as_view.register(xr.DataArray)
def _view_pd_boolean_array(a: xr.DataArray, view_args):
    return a
