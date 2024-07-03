from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal

if TYPE_CHECKING:
    from collections.abc import Iterable


import pandas as pd
import xarray as xr

from ..._core.anndata import _gen_dataframe, _remove_unused_categories
from ..._core.file_backing import to_memory
from ..._core.index import Index, _subset
from ..._core.views import as_view


def get_index_dim(ds):
    assert (
        len(ds.sizes) == 1
    ), f"xarray Dataset should not have more than 1 dims, found {len(ds)}"
    return list(ds.sizes.keys())[0]


class Dataset2D(xr.Dataset):
    @property
    def index(self) -> pd.Index:
        coord = list(self.coords.keys())[0]
        return pd.Index(self.coords[coord].data)

    @property
    def shape(
        self,
    ):  # aligned mapping classes look for this for DataFrames so this ensures usability with e.g., obsm
        return [self.dims[get_index_dim(self)], len(self)]

    @property
    def iloc(self):
        class IlocGetter:
            def __init__(self, ds):
                self._ds = ds

            def __getitem__(self, idx):
                coord = list(self._ds.coords.keys())[0]
                return self._ds.isel(**{coord: idx})

        return IlocGetter(self)


@_subset.register(Dataset2D)
def _(a: xr.DataArray, subset_idx: Index):
    key = get_index_dim(a)
    if (
        isinstance(subset_idx, tuple) and len(subset_idx) == 1
    ):  # xarray seems to have some code looking for a second entry in tuples
        return a.isel(**{key: subset_idx[0]})
    return a.isel(**{key: subset_idx})


@as_view.register(Dataset2D)
def _(a: Dataset2D, view_args):
    return a


@_gen_dataframe.register(Dataset2D)
def _gen_dataframe_xr(
    anno: Dataset2D,
    index_names: Iterable[str],
    *,
    source: Literal["X", "shape"],
    attr: Literal["obs", "var"],
    length: int | None = None,
):
    return anno


@_remove_unused_categories.register(Dataset2D)
def _remove_unused_categories_xr(
    df_full: Dataset2D, df_sub: Dataset2D, uns: dict[str, Any]
):
    pass  # this is handled automatically by the categorical arrays themselves i.e., they dedup upon access.


@to_memory.register(Dataset2D)
def to_memory(ds: Dataset2D, copy=False):
    df = ds.to_dataframe()
    df.index.name = None  # matches old AnnData object
    return df
