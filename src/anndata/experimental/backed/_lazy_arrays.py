from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING, Generic, TypeVar

import pandas as pd

from anndata._core.index import _subset
from anndata._core.views import as_view
from anndata._io.specs.lazy_methods import get_chunksize
from anndata.compat import H5Array, ZarrArray

from ..._settings import settings
from ._compat import BackendArray, DataArray, ZarrArrayWrapper
from ._compat import xarray as xr

if TYPE_CHECKING:
    from pathlib import Path
    from typing import Literal

    import numpy as np

    from anndata._core.index import Index
    from anndata.compat import ZarrGroup


K = TypeVar("K", H5Array, ZarrArray)


class ZarrOrHDF5Wrapper(ZarrArrayWrapper, Generic[K]):
    def __init__(self, array: K):
        self.chunks = array.chunks
        if isinstance(array, ZarrArray):
            return super().__init__(array)
        self._array = array
        self.shape = self._array.shape
        self.dtype = self._array.dtype

    def __getitem__(self, key: xr.core.indexing.ExplicitIndexer):
        if isinstance(self._array, ZarrArray):
            return super().__getitem__(key)
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
    base_path_or_zarr_group: Path | ZarrGroup
    elem_name: str

    def __init__(
        self,
        codes: K,
        categories: ZarrArray | H5Array,
        ordered: bool,
        base_path_or_zarr_group: Path | ZarrGroup,
        elem_name: str,
        *args,
        **kwargs,
    ):
        self._categories = categories
        self._ordered = ordered
        self._codes = ZarrOrHDF5Wrapper(codes)
        self.shape = self._codes.shape
        self.base_path_or_zarr_group = base_path_or_zarr_group
        self.file_format = "zarr" if isinstance(codes, ZarrArray) else "h5"
        self.elem_name = elem_name

    @cached_property
    def categories(self) -> np.ndarray:
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
    base_path_or_zarr_group: Path | ZarrGroup
    elem_name: str

    def __init__(
        self,
        values: ZarrArray | H5Array,
        dtype_str: Literal["nullable-integer", "nullable-boolean"],
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

    def __getitem__(
        self, key: xr.core.indexing.ExplicitIndexer
    ) -> xr.core.extension_array.PandasExtensionArray:
        values = self._values[key]
        mask = self._mask[key]
        if self._dtype_str == "nullable-integer":
            # numpy does not support nan ints
            extension_array = pd.arrays.IntegerArray(values, mask=mask)
        elif self._dtype_str == "nullable-boolean":
            extension_array = pd.arrays.BooleanArray(values, mask=mask)
        else:
            raise RuntimeError(f"Invalid dtype_str {self._dtype_str}")
        return xr.core.extension_array.PandasExtensionArray(extension_array)

    @cached_property
    def dtype(self):
        if self._dtype_str == "nullable-integer":
            return pd.array(
                [],
                dtype=str(pd.api.types.pandas_dtype(self._values.dtype)).capitalize(),
            ).dtype
        elif self._dtype_str == "nullable-boolean":
            return pd.BooleanDtype()
        raise RuntimeError(f"Invalid dtype_str {self._dtype_str}")


@_subset.register(DataArray)
def _subset_masked(a: DataArray, subset_idx: Index):
    return a[subset_idx]


@as_view.register(DataArray)
def _view_pd_boolean_array(a: DataArray, view_args):
    return a


@get_chunksize.register(MaskedArray)
def _(a: MaskedArray):
    return get_chunksize(a._values)


@get_chunksize.register(CategoricalArray)
def _(a: CategoricalArray):
    return get_chunksize(a._codes)
