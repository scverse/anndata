from typing import Tuple
from anndata._core.index import Index, _subset
from anndata._core.views import _resolve_idx, as_view
from anndata._io.h5ad import read_dataset
from anndata.compat import ZarrArray

import pandas as pd
import numpy as np
from xarray.core.indexing import ExplicitlyIndexedNDArrayMixin, BasicIndexer, OuterIndexer
import xarray as xr


class MaskedArrayMixIn(ExplicitlyIndexedNDArrayMixin):

    def __eq__(self, __o) -> np.ndarray:
        return self[...] == __o

    def __ne__(self, __o) -> np.ndarray:
        return ~(self == __o)
    
    @property
    def shape(self) -> Tuple[int, ...]:
        """Shape of this array

        Returns:
            Tuple[int, ...]: A shape that looks like a 1-d shape i.e., (#, )
        """
        return self.values.shape


class LazyCategoricalArray(MaskedArrayMixIn):
    __slots__ = (
        "values",
        "attrs",
        "_categories",
        "_categories_cache",
        "group",
    )

    def __init__(self, codes, categories, attrs, *args, **kwargs):
        """Class for lazily reading categorical data from formatted zarr group.   Used as base for `LazilyIndexedArray`.

        Args:
            codes (Union[zarr.Array, h5py.Dataset]): values (integers) of the array, one for each element
            categories (Union[zarr.Array, h5py.Dataset]): mappings from values to strings
            attrs (Union[zarr.Array, h5py.Dataset]): attrs containing boolean "ordered"
        """
        self.values = codes
        self._categories = categories
        self._categories_cache = None
        self.attrs = dict(attrs)

    @property
    def categories(self):  # __slots__ and cached_property are incompatible
        if self._categories_cache is None:
            if isinstance(self._categories, ZarrArray):
                self._categories_cache = self._categories[...]
            else:
                self._categories_cache = read_dataset(self._categories)
        return self._categories_cache

    @property
    def dtype(self) -> pd.CategoricalDtype:
        return pd.CategoricalDtype(self.categories, self.ordered)

    @property
    def ordered(self):
        return bool(self.attrs["ordered"])

    def __getitem__(self, selection) -> pd.Categorical:
        idx = selection
        if isinstance(selection, BasicIndexer) or isinstance(selection, OuterIndexer):
            idx = selection.tuple[0] # need to better understand this
        if isinstance(self.values, ZarrArray):
            codes = self.values.oindex[idx]
        else:
            codes = self.values[idx]
        if codes.shape == ():  # handle 0d case
            codes = np.array([codes])
        return pd.Categorical.from_codes(
            codes=codes,
            categories=self.categories,
            ordered=self.ordered,
        ).remove_unused_categories()

    def __repr__(self) -> str:
        return f"LazyCategoricalArray(codes=..., categories={self.categories}, ordered={self.ordered})"

    def copy(self) -> "LazyCategoricalArray":
        """Returns a copy of this array which can then be safely edited

        Returns:
            LazyCategoricalArray: copied LazyCategoricalArray
        """
        arr = LazyCategoricalArray(
            self.values, self._categories, self.attrs
        )  # self.categories reads in data
        return arr


class LazyMaskedArray(MaskedArrayMixIn):
    __slots__ = ("mask", "values", "_dtype_str")

    def __init__(self, values, mask, dtype_str, *args, **kwargs):
        """Class for lazily reading categorical data from formatted zarr group.  Used as base for `LazilyIndexedArray`.

        Args:
            values (Union[zarr.Array, h5py.Dataset]): Integer/Boolean array of values
            mask (Union[zarr.Array, h5py.Dataset]): mask indicating which values are non-null
            dtype_str (Nullable): one of `nullable-integer` or `nullable-boolean`
        """
        self.values = values
        self.mask = mask
        self._dtype_str = dtype_str

    @property
    def dtype(self) -> pd.CategoricalDtype:
        if self.mask is not None:
            if self._dtype_str == "nullable-integer":
                return pd.arrays.IntegerArray
            elif self._dtype_str == "nullable-boolean":
                return pd.arrays.BooleanArray
        return pd.array

    def __getitem__(self, selection) -> pd.Categorical:
        idx = selection
        if isinstance(selection, BasicIndexer) or isinstance(selection, OuterIndexer):
            idx = selection.tuple[0] # need to understand this better
        if type(idx) == int:
            idx = slice(idx, idx + 1)
        values = np.array(self.values[idx])
        if self.mask is not None:
            mask = np.array(self.mask[idx])
            if self._dtype_str == "nullable-integer":
                return pd.arrays.IntegerArray(values, mask=mask)
            elif self._dtype_str == "nullable-boolean":
                return pd.arrays.BooleanArray(values, mask=mask)
        return pd.array(values)

    def __repr__(self) -> str:
        if self._dtype_str == "nullable-integer":
            return "LazyNullableIntegerArray"
        elif self._dtype_str == "nullable-boolean":
            return "LazyNullableBooleanArray"

    def copy(self) -> "LazyMaskedArray":
        """Returns a copy of this array which can then be safely edited

        Returns:
            LazyMaskedArray: copied LazyMaskedArray
        """
        arr = LazyMaskedArray(self.values, self.mask, self._dtype_str)
        return arr


@_subset.register(xr.DataArray)
def _subset_masked(a: xr.DataArray, subset_idx: Index):
    return a[subset_idx]

@as_view.register(pd.Categorical)
def _view_pd_categorical(a: pd.Categorical, view_args):
    return a


@as_view.register(pd.api.extensions.ExtensionArray)
def _view_pd_array(a: pd.api.extensions.ExtensionArray, view_args):
    return a


@as_view.register(pd.arrays.IntegerArray)
def _view_pd_integer_array(a: pd.arrays.IntegerArray, view_args):
    return a


@as_view.register(pd.arrays.BooleanArray)
def _view_pd_boolean_array(a: pd.arrays.BooleanArray, view_args):
    return a

@as_view.register(xr.DataArray)
def _view_pd_boolean_array(a: xr.DataArray, view_args):
    return a
