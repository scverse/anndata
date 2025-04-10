from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

from ..compat import XDataset

if TYPE_CHECKING:
    from collections.abc import Hashable, Iterable
    from typing import Any, Literal

    from .index import Index
    from ..compat import XArray


class Dataset2D(XDataset):
    """
    A wrapper class meant to enable working with lazy dataframe data.
    We do not guarantee the stability of this API beyond that guaranteed
    by :class:`xarray.Dataset` and the `to_memory` function, a thin wrapper
    around :meth:`xarray.Dataset.to_dataframe` to ensure roundtrip
    compatibility here.
    """

    __slots__ = ()

    @property
    def index_dim(self) -> str:
        if len(self.sizes) != 1:
            msg = f"xarray Dataset should not have more than 1 dims, found {len(self.sizes)} {self.sizes}, {self}"
            raise ValueError(msg)
        return next(iter(self.coords.keys()))

    @property
    def xr_index(self) -> XArray:
        return self[self.index_dim]

    @property
    def index(self) -> pd.Index:
        """:attr:`~anndata.AnnData` internally looks for :attr:`~pandas.DataFrame.index` so this ensures usability

        Returns
        -------
        The index of the of the dataframe as resolved from :attr:`~xarray.Dataset.coords`.
        """
        return self.indexes[self.index_dim]

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
        return (self.sizes[self.index_dim], len(self))

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
                coord = self._ds.index_dim
                return self._ds.isel(**{coord: idx})

        return IlocGetter(self)

    def to_memory(self, *, copy=False) -> pd.DataFrame:
        df = self.to_dataframe()
        index_key = self.attrs.get("indexing_key", None)
        if df.index.name != index_key and index_key is not None:
            df = df.set_index(index_key)
        df.index.name = None  # matches old AnnData object
        return df

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
