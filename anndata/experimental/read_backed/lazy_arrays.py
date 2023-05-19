from typing import Tuple
from anndata._core.index import Index, _subset
from anndata._core.views import _resolve_idx, as_view

import pandas as pd
import numpy as np
from xarray.core.indexing import ExplicitlyIndexedNDArrayMixin


class MaskedArrayMixIn(ExplicitlyIndexedNDArrayMixin):
    def _resolve_idx(self, new_idx):
        return (
            new_idx
            if self.subset_idx is None
            else _resolve_idx(self.subset_idx, new_idx, self.shape[0])
        )

    @property
    def subset_idx(self):
        return self._subset_idx

    @subset_idx.setter
    def subset_idx(self, new_idx):
        self._subset_idx = self._resolve_idx(new_idx)

    @property
    def shape(self) -> Tuple[int, ...]:
        if self.subset_idx is None:
            return self.values.shape
        if isinstance(self.subset_idx, slice):
            if self.subset_idx == slice(None, None, None):
                return self.values.shape
            return (self.subset_idx.stop - self.subset_idx.start,)
        else:
            return (len(self.subset_idx),)

    def __eq__(self, __o) -> np.ndarray:
        return self[...] == __o

    def __ne__(self, __o) -> np.ndarray:
        return ~(self == __o)


class LazyCategoricalArray(MaskedArrayMixIn):
    __slots__ = (
        "values",
        "attrs",
        "_categories",
        "_categories_cache",
        "_subset_idx",
        "group",
    )

    def __init__(self, codes, categories, attrs, *args, **kwargs):
        """Class for lazily reading categorical data from formatted zarr group

        Args:
            group (zarr.Group): group containing "codes" and "categories" key as well as "ordered" attr
        """
        self.values = codes
        self._categories = categories
        self._categories_cache = None
        self._subset_idx = None
        self.attrs = dict(attrs)

    @property
    def categories(self):  # __slots__ and cached_property are incompatible
        if self._categories_cache is None:
            self._categories_cache = self._categories[...]
        return self._categories_cache

    @property
    def dtype(self) -> pd.CategoricalDtype:
        return pd.CategoricalDtype(self.categories, self.ordered)

    @property
    def ordered(self):
        return bool(self.attrs["ordered"])

    def __getitem__(self, selection) -> pd.Categorical:
        idx = self._resolve_idx(selection)
        codes = self.values.oindex[idx]
        if codes.shape == ():  # handle 0d case
            codes = np.array([codes])
        return pd.Categorical.from_codes(
            codes=codes,
            categories=self.categories,
            ordered=self.ordered,
        ).remove_unused_categories()

    def __repr__(self) -> str:
        return f"LazyCategoricalArray(codes=..., categories={self.categories}, ordered={self.ordered})"

    def copy(self):
        arr = LazyCategoricalArray(self.values, self.categories, self.attrs)
        arr.subset_idx = self.subset_idx
        return arr


class LazyMaskedArray(MaskedArrayMixIn):
    __slots__ = ("mask", "values", "_subset_idx", "_dtype_str")

    def __init__(self, values, mask, dtype_str, *args, **kwargs):
        """Class for lazily reading categorical data from formatted zarr group

        Args:
            group (zarr.Group): group containing "codes" and "categories" key as well as "ordered" attr
            dtype_str (Nullable): one of `nullable-integer` or `nullable-boolean`
        """
        self.values = values
        self.mask = mask
        self._subset_idx = None
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
        idx = self._resolve_idx(selection)
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

    def copy(self):
        arr = LazyMaskedArray(self.values, self.mask, self._dtype_str)
        arr.subset_idx = self.subset_idx
        return arr


@_subset.register(MaskedArrayMixIn)
def _subset_masked(a: MaskedArrayMixIn, subset_idx: Index):
    a_copy = a.copy()
    a_copy.subset_idx = subset_idx
    return a_copy


@as_view.register(MaskedArrayMixIn)
def _view_masked(a: MaskedArrayMixIn, view_args):
    return a


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
