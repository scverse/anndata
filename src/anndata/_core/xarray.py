from __future__ import annotations

from functools import wraps
from typing import TYPE_CHECKING
import warnings

import pandas as pd

from ..compat import XDataArray, XDataset, XVariable


def requires_xarray(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            import xarray  # noqa: F401
        except ImportError as e:
            msg = "xarray is required to read dataframes lazily. Please install xarray."
            raise ImportError(msg) from e
        return func(*args, **kwargs)

    return wrapper


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
    def is_backed(self) -> bool:
        return self.attrs.get("is_backed", False)

    @is_backed.setter
    def is_backed(self, isbacked: bool):
        if not isbacked and "is_backed" in self.attrs:
            del self.attrs["is_backed"]
        else:
            self.attrs["is_backed"] = isbacked

    @property
    def index_dim(self) -> str:
        if len(self.sizes) != 1:
            msg = f"xarray Dataset should not have more than 1 dims, found {len(self.sizes)} {self.sizes}, {self}"
            raise ValueError(msg)
        return next(iter(self.coords.keys()))

    @property
    def true_index_dim(self) -> str:
        return self.attrs.get("indexing_key", self.index_dim)

    @true_index_dim.setter
    def true_index_dim(self, val: str):
        if val is None or (val == self.index_dim and "indexing_key" in self.attrs):
            del self.attrs["indexing_key"]
        elif val not in self.dims:
            if val not in self.data_vars:
                msg = f"Unknown variable `{val}`."
                raise ValueError(msg)
            self.attrs["indexing_key"] = val

    @property
    def xr_index(self) -> XDataArray:
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
        index_dim = self.index_dim
        self.coords[index_dim] = (index_dim, val)
        if isinstance(val, pd.Index) and val.name is not None and val.name != index_dim:
            self.update(self.rename({self.index_dim: val.name}))
            del self.coords[index_dim]
        # without `indexing_key` explicitly set on `self.attrs`, `self.true_index_dim` will use the `self.index_dim`
        if "indexing_key" in self.attrs:
            del self.attrs["indexing_key"]

    @property
    def true_xr_index(self) -> XDataArray:
        return self[self.true_index_dim]

    @property
    def true_index(self) -> pd.Index:
        return self.true_xr_index.to_index()

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

    def __getitem__(self, idx) -> Dataset2D:
        ret = super().__getitem__(idx)
        if len(idx) == 0 and not isinstance(idx, tuple):  # empty XDataset
            ret.coords[self.index_dim] = self.xr_index
        return ret

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
        columns = set(self.keys())
        index_key = self.attrs.get("indexing_key", None)
        if index_key is not None:
            columns.discard(index_key)
        return pd.Index(columns)

    def __setitem__(self, key, value):
        """
        Setting can only be performed when the incoming value is “standalone” like :class:`nump.ndarray` to mimic pandas.
        One can also use the tuple setting style like `ds["foo"] = (ds.index_dim, value)` to set the value, although the index name must match.
        Similarly, one can use the :class:`xarray.DataArray` but it must have the same (one and only one) dim name/coord name as `self.index_dim`.

        For supported setter values see :meth:`xarray.Dataset.__setitem__`.
        """
        if key == self.index_dim:
            msg = f"Cannot set {self.index_dim} as a variable. Use `index` instead."
            raise KeyError(msg)
        if isinstance(value, tuple):
            if isinstance(value[0], tuple):
                if value[0][0] != self.index_dim:
                    msg = f"Dimension tuple should have only {self.index_dim} as its dimension, found {value[0][0]}"
                    raise ValueError(msg)
                if len(value[0]) > 1:
                    msg = "Dimension tuple is too long."
                    raise ValueError(msg)
            elif value[0] != self.index_dim:
                msg = f"Setting value tuple should have first entry {self.index_dim}, found {value[0]}"
                raise ValueError(msg)
        elif isinstance(value, XDataArray | XDataset | XVariable):
            value_typ = type(value).__name__
            # https://docs.xarray.dev/en/stable/generated/xarray.Dataset.dims.html#xarray.Dataset.dims
            # Unfortunately `dims` not the same across data structures.
            with warnings.catch_warnings(action="ignore"):
                dims = (
                    list(value.dims.keys())
                    if isinstance(value, XDataset)
                    else value.dims
                )
            if (
                isinstance(value, XDataArray)
                and value.name is not None
                and value.name != key
            ):
                msg = f"{value_typ} should have name {key}, found {value.name}"
                raise ValueError(msg)
            if len(dims) != 1:
                msg = f"{value_typ} should have only one dimension, found {len(dims)}"
                raise ValueError(msg)
            if dims[0] != self.index_dim:
                msg = f"{value_typ} should have dimension {self.index_dim}, found {dims[0]}"
                raise ValueError(msg)
            if not isinstance(value, XVariable) and (
                self.index_dim not in value.coords
                or value.coords[self.index_dim].name != self.index_dim
            ):
                msg = f"{value_typ} should have coordinate {self.index_dim} with same name, found {value.coords} with name {value.coords[next(iter(value.coords.keys()))].name}"
                raise ValueError(msg)
        else:
            # maintain setting behavior of a 2D dataframe i.e., one dim
            value = (self.index_dim, value)
        super().__setitem__(key, value)
