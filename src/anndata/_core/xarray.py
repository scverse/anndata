from __future__ import annotations

import warnings
from dataclasses import dataclass
from functools import wraps
from typing import TYPE_CHECKING, overload

import numpy as np
import pandas as pd

from ..compat import XDataArray, XDataset, XVariable

if TYPE_CHECKING:
    from collections.abc import Hashable, Iterable, Iterator, Mapping
    from typing import Any, Literal

    from .._types import Dataset2DIlocIndexer


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


class Dataset2D:
    r"""
    Bases :class:`~collections.abc.Mapping`\ [:class:`~collections.abc.Hashable`, :class:`~xarray.DataArray` | :class:`~anndata.experimental.backed.Dataset2D`\ ]

    A wrapper class meant to enable working with lazy dataframe data according to
    :class:`~anndata.AnnData`'s internal API.  This class ensures that "dataframe-invariants"
    are respected, namely that there is only one 1d dim and coord with the same name i.e.,
    like a :class:`pandas.DataFrame`.

    You should not have to initiate this class yourself.  Setting an :class:`xarray.Dataset`
    into a relevant part of the :class:`~anndata.AnnData` object will attempt to wrap that
    object in this object, trying to enforce the "dataframe-invariants."

    Because xarray requires :attr:`xarray.Dataset.coords` to be in-memory, this class provides
    handling for an out-of-memory index via :attr:`~anndata.experimental.backed.Dataset2D.true_index`.
    This feature is helpful for loading remote data faster where the index itself may not be initially useful
    for constructing the object e.g., cell ids.
    """

    @staticmethod
    def _validate_shape_invariants(ds: XDataset):
        """
        Validate that the dataset has only one dimension, which is the index dimension.
        This is a requirement for 2D datasets.
        """
        if not isinstance(ds, XDataset):
            msg = f"Expected an xarray Dataset, found {type(ds)}"
            raise TypeError(msg)
        if (is_coords_too_long := (len(ds.coords) != 1)) or len(ds.dims) != 1:
            string, length, rep = (
                ("coordinate", len(ds.coords), ds.coords)
                if is_coords_too_long
                else ("dimension", len(ds.dims), ds.dims)
            )
            msg = f"Dataset should have exactly one {string}, found {length}: {rep}"
            raise ValueError(msg)
        if next(iter(ds.dims)) != next(iter(ds.coords)):
            msg = f"Dataset dimension {next(iter(ds.dims))} does not match coordinate {next(iter(ds.coords))}."
            raise ValueError(msg)

    def __init__(self, ds: XDataset):
        Dataset2D._validate_shape_invariants(ds)
        self._ds = ds

    @property
    def ds(self) -> XDataset:
        """The underlying :class:`xarray.Dataset`."""
        return self._ds

    def keys(self) -> list[Hashable]:
        return list(iter(self.ds))

    @property
    def is_backed(self) -> bool:
        """
        Check whether or not the object is backed, used to indicate if there are any in-memory objects.
        Must be externally set, defaults false.
        """
        return self.ds.attrs.get("is_backed", False)

    @is_backed.setter
    def is_backed(self, isbacked: bool) -> bool:
        if not isbacked and "is_backed" in self.ds.attrs:
            del self.ds.attrs["is_backed"]
        else:
            self.ds.attrs["is_backed"] = isbacked

    @property
    def index_dim(self) -> str:
        """The underlying computational index i.e., the lone coordinate dimension."""
        if len(self.ds.sizes) != 1:
            msg = f"xarray Dataset should not have more than 1 dims, found {len(self.ds.sizes)} {self.ds.sizes}, {self}"
            raise ValueError(msg)
        return next(iter(self.ds.coords.keys()))

    @property
    def true_index_dim(self) -> str:
        """
        Because xarray loads its coordinates/indexes in memory,
        we allow for signaling that a given variable, which is not a coordinate, is the "true" index.

        For example, the true index may be cell names but loading these over an internet connection may not be
        desirable or necessary for most use cases such as getting a quick preview of the columns or loading only
        one column that isn't the index.

        This property is the key of said variable. The default is `index_dim` if this variable has not been set.
        """
        return self.ds.attrs.get("indexing_key", self.index_dim)

    @true_index_dim.setter
    def true_index_dim(self, val: str):
        if val is None or (val == self.index_dim and "indexing_key" in self.ds.attrs):
            del self.ds.attrs["indexing_key"]
        elif val not in self.ds.dims:
            if val not in self.ds.data_vars:
                msg = f"Unknown variable `{val}`."
                raise ValueError(msg)
            self.ds.attrs["indexing_key"] = val

    @property
    def xr_index(self) -> XDataArray:
        """The coordinate of :attr:`anndata.experimental.backed.Dataset2D.index_dim`"""
        return self.ds[self.index_dim]

    @property
    def index(self) -> pd.Index:
        """:attr:`~anndata.AnnData` internally looks for :attr:`~pandas.DataFrame.index` so this ensures usability
        A :class:`pandas.Index` object corresponding to :attr:`anndata.experimental.backed.Dataset2D.index_dim`
        Returns
        -------
        The index of the of the dataframe as resolved from :attr:`~xarray.Dataset.coords`.
        """
        return self.ds.indexes[self.index_dim]

    @index.setter
    def index(self, val) -> None:
        index_dim = self.index_dim
        self.ds.coords[index_dim] = (index_dim, val)
        if isinstance(val, pd.Index) and val.name is not None and val.name != index_dim:
            self.ds.update(self.ds.rename({self.index_dim: val.name}))
            del self.ds.coords[index_dim]
        # without `indexing_key` explicitly set on `self.ds.attrs`, `self.true_index_dim` will use the `self.index_dim`
        if "indexing_key" in self.ds.attrs:
            del self.ds.attrs["indexing_key"]

    @property
    def true_xr_index(self) -> XDataArray:
        """The index :class:`~anndata.AnnData` is actually interested in e.g., cell names, for verification."""
        return self.ds[self.true_index_dim]

    @property
    def true_index(self) -> pd.Index:
        """:attr:`~anndata.experimental.backed.Dataset2D.true_xr_index` as a :class:`pandas.Index`"""
        return self.true_xr_index.to_index()

    @property
    def shape(self) -> tuple[int, int]:
        """:attr:`~anndata.AnnData` internally looks for :attr:`~pandas.DataFrame.shape` so this ensures usability

        Returns
        -------
        The (2D) shape of the dataframe resolved from :attr:`~xarray.Dataset.sizes`.
        """
        return (self.ds.sizes[self.index_dim], len(self.ds))

    @property
    def iloc(self) -> Dataset2DIlocIndexer:
        """:attr:`~anndata.AnnData` internally looks for :attr:`~pandas.DataFrame.iloc` so this ensures usability

        Returns
        -------
        Handler class for doing the iloc-style indexing using :meth:`~xarray.Dataset.isel`.
        """

        return IlocGetter(self.ds, self.index_dim)

    # See https://github.com/pydata/xarray/blob/568f3c1638d2d34373408ce2869028faa3949446/xarray/core/dataset.py#L1239-L1248
    # for typing
    @overload
    def __getitem__(self, key: Hashable) -> XDataArray: ...
    @overload
    def __getitem__(self, key: Iterable[Hashable]) -> Dataset2D: ...
    def __getitem__(
        self, key: Mapping[Any, Any] | Hashable | Iterable[Hashable]
    ) -> Dataset2D | XDataArray:
        ret = self.ds.__getitem__(key)
        if len(key) == 0 and not isinstance(key, tuple):  # empty XDataset
            ret.coords[self.index_dim] = self.xr_index
        if isinstance(ret, XDataset):
            # If we get an xarray Dataset, we return a Dataset2D
            as_2d = Dataset2D(ret)

            as_2d.true_index_dim = self.true_index_dim
            as_2d.is_backed = self.is_backed
            return as_2d
        return ret

    def to_memory(self, *, copy: bool = False) -> pd.DataFrame:
        """
        Converts to :class:`pandas.DataFrame`.
        The index of the dataframe comes from :attr:`~anndata.experimental.backed.Dataset2D.true_index_dim`
        if it differs from :attr:`~anndata.experimental.backed.Dataset2D.index_dim`.

        Parameters
        ----------
        copy
            Unused argument

        Returns
        -------
            :class:`pandas.DataFrame` with index set accordingly.
        """
        # https://github.com/pydata/xarray/issues/10419
        non_nullable_string_cols = {
            col
            for col in self.columns
            if not self[col].attrs.get("is_nullable_string", False)
        }
        df = self.ds.to_dataframe()
        index_key = self.ds.attrs.get("indexing_key", None)
        if df.index.name != index_key and index_key is not None:
            df = df.set_index(index_key)
        for col in set(self.columns) - non_nullable_string_cols:
            df[col] = df[col].astype(dtype="string")
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
        columns = set(self.ds.keys())
        index_key = self.ds.attrs.get("indexing_key", None)
        if index_key is not None:
            columns.discard(index_key)
        return pd.Index(columns)

    def __setitem__(
        self, key: Hashable | Iterable[Hashable] | Mapping, value: Any
    ) -> None:
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
        self.ds.__setitem__(key, value)

    def copy(
        self,
        data: Mapping | None = None,
        *,
        deep: bool = False,
    ) -> Dataset2D:
        """
        Return a copy of the Dataset2D object.
        See :meth:`xarray.Dataset.copy` for more information.
        """
        as_2d = Dataset2D(self.ds.copy(deep=deep, data=data))
        as_2d.true_index_dim = self.true_index_dim
        as_2d.is_backed = self.is_backed
        return as_2d

    def __iter__(self) -> Iterator[Hashable]:
        return iter(self.ds)

    def __len__(self) -> int:
        return len(self.ds)

    @property
    def dtypes(self) -> pd.Series:
        """
        Return a Series with the dtypes of the variables in the Dataset2D.
        """
        return self.ds.dtypes

    def equals(self, b: object) -> bool:
        """Thin wrapper around :meth:`xarray.Dataset.equals`"""
        if isinstance(b, Dataset2D):
            b = b.ds
        return self.ds.equals(b)

    def reindex(
        self,
        index: pd.Index | None = None,
        axis: Literal[0] = 0,
        fill_value: Any | None = np.nan,
    ) -> Dataset2D:
        """Reindex the current object against a new index.

        Parameters
        ----------
        index
            The new index for reindexing, by default None
        axis
            Provided for API consistency, should not be called over axis!=0, by default 0
        fill_value
            The value with which to fill in via :meth:`pandas.Series.reindex`, by default np.nan

        Returns
        -------
            Reindexed dataset.
        """
        index_dim = self.index_dim
        if axis != 0:  # pragma: no cover
            msg = f"Only axis 0 is supported, got axis: {axis}"
            raise ValueError(msg)
        # Dataset.reindex() can't handle ExtensionArrays
        extension_arrays = {
            col: data
            for col, data in self._items()
            if pd.api.types.is_extension_array_dtype(data.dtype)
        }
        el = self.ds.drop_vars(extension_arrays.keys())
        el = el.reindex({index_dim: index}, method=None, fill_value=fill_value)
        for col, data in extension_arrays.items():
            el[col] = XDataArray.from_series(
                pd.Series(data.data, index=self.index).reindex(
                    index.rename(self.index.name) if index is not None else index,
                    fill_value=fill_value,
                )
            )
        return Dataset2D(el)

    # Used "publicly" in src/anndata/_core/merge.py but not intended for public use.
    def _items(self):
        for col in self:
            yield col, self[col]


@dataclass(frozen=True)
class IlocGetter:
    _ds: XDataset
    _coord: str

    def __getitem__(self, idx) -> Dataset2D:
        # xarray seems to have some code looking for a second entry in tuples,
        # so we unpack the tuple
        if isinstance(idx, tuple) and len(idx) == 1:
            idx = idx[0]
        return Dataset2D(self._ds.isel(**{self._coord: idx}))
