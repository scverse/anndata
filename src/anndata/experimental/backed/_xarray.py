from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

from ..._core.anndata import AnnData, _gen_dataframe
from ..._core.file_backing import to_memory
from ..._core.index import _subset
from ..._core.views import as_view

try:
    from xarray import Dataset
except ImportError:

    class Dataset:
        def __repr__(self) -> str:
            return "mock Dataset"


if TYPE_CHECKING:
    from collections.abc import Hashable, Iterable
    from typing import Any, Literal

    from ..._core.index import Index
    from ._compat import xarray as xr


def get_index_dim(ds: xr.DataArray) -> Hashable:
    if len(ds.sizes) != 1:
        msg = f"xarray Dataset should not have more than 1 dims, found {len(ds)}"
        raise ValueError(msg)
    return list(ds.indexes.keys())[0]


class Dataset2D(Dataset):
    __slots__ = ()

    @property
    def index(self) -> pd.Index:
        """:attr:`~anndata.AnnData` internally looks for :attr:`~pandas.DataFrame.index` so this ensures usability

        Returns
        -------
        The index of the of the dataframe as resolved from :attr:`~xarray.Dataset.coords`.
        """
        coord = get_index_dim(self)
        return self.indexes[coord]

    @index.setter
    def index(self, val) -> None:
        coord = get_index_dim(self)
        self.coords[coord] = val

    @property
    def shape(self) -> tuple[int, int]:
        """:attr:`~anndata.AnnData` internally looks for :attr:`~pandas.DataFrame.shape` so this ensures usability

        Returns
        -------
        The (2D) shape of the dataframe resolved from :attr:`~xarray.Dataset.sizes`.
        """
        return (self.sizes[get_index_dim(self)], len(self))

    @property
    def iloc(self):
        """:attr:`~anndata.AnnData` internally looks for :attr:`~pandas.DataFrame.iloc` so this ensures usability

        Returns
        -------
        Handler class for doing the iloc-style indexing using :meth:`~xarray.Dataset.isel`.
        """

        class IlocGetter:
            def __init__(self, ds):
                self._ds = ds

            def __getitem__(self, idx):
                coord = get_index_dim(self._ds)
                return self._ds.isel(**{coord: idx})

        return IlocGetter(self)

    @property
    def columns(self) -> pd.Index:
        """
        :class:`~anndata.AnnData` internally looks for :attr:`~pandas.DataFrame.columns` so this ensures usability

        Returns
        -------
        :class:`pandas.Index` that represents the "columns."
        """
        columns_list = list(self.keys())
        return pd.Index(columns_list)


@_subset.register(Dataset2D)
def _(a: Dataset2D, subset_idx: Index):
    key = a.attrs["indexing_key"]
    # xarray seems to have some code looking for a second entry in tuples
    if isinstance(subset_idx, tuple) and len(subset_idx) == 1:
        subset_idx = subset_idx[0]
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


@AnnData._remove_unused_categories.register(Dataset2D)
def _remove_unused_categories_xr(
    df_full: Dataset2D, df_sub: Dataset2D, uns: dict[str, Any]
):
    pass  # this is handled automatically by the categorical arrays themselves i.e., they dedup upon access.


@to_memory.register(Dataset2D)
def to_memory(ds: Dataset2D, *, copy=False):
    df = ds.to_dataframe()
    index_key = ds.attrs.get("indexing_key", None)
    if df.index.name != index_key and index_key is not None:
        df = df.set_index(index_key)
    df.index.name = None  # matches old AnnData object
    return df
