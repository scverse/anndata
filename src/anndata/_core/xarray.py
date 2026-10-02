"""
:class:`Dataset2D`, the lazy dataframe AnnData stores as ``obs``/``var``.

Alongside it sit the conversions to pandas (:meth:`Dataset2D.to_memory`) and to dask
(:meth:`Dataset2D.to_dask_dataframe`), the inverse of the latter, and the narwhals plugin
handing narwhals whichever of the two applies.
"""

from __future__ import annotations

import warnings
from collections.abc import Hashable, Mapping, Sized
from dataclasses import dataclass
from functools import partial, wraps
from itertools import accumulate, pairwise
from typing import TYPE_CHECKING, overload

import narwhals as nw
import numpy as np
import pandas as pd

from anndata._warnings import warn

from ..compat import DaskArray, XDataArray, XDataset, XVariable, pandas_as_str

if TYPE_CHECKING:
    from collections.abc import Callable, Collection, Iterable, Iterator, KeysView
    from typing import Any, Literal, Self

    from dask.dataframe import DataFrame as DaskDataFrame
    from narwhals._dask.dataframe import DaskLazyFrame
    from narwhals._pandas_like.dataframe import PandasLikeDataFrame
    from narwhals.utils import Version
    from pandas.api.typing.aliases import Scalar

    from .._types import Dataset2DIlocIndexer


NULLABLE_STRING_ATTR = "is_nullable_string"
"""Variable attribute marking values standing in for a pandas nullable-string array.

xarray cannot hold such an array, so whoever builds the variable stores the values some
other way and sets this flag; :meth:`Dataset2D.to_memory` reads it to convert them back.
"""


def requires_xarray[R, **P](func: Callable[P, R]) -> Callable[P, R]:
    @wraps(func)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
        try:
            import xarray  # noqa: F401
        except ImportError as e:
            msg = "xarray is required to read dataframes lazily. Please install xarray."
            raise ImportError(msg) from e
        return func(*args, **kwargs)

    return wrapper


class Dataset2D(Mapping[Hashable, "XDataArray | Dataset2D"]):
    r"""
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
    def _validate_shape_invariants(ds: XDataset) -> None:
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

    def keys(self) -> KeysView[Hashable]:
        return self.ds.keys()

    @property
    def is_backed(self) -> bool:
        """
        Check whether or not the object is backed, used to indicate if there are any in-memory objects.
        Must be externally set, defaults false.
        """
        return self.ds.attrs.get("is_backed", False)

    @is_backed.setter
    def is_backed(self, isbacked: bool) -> None:
        if not isbacked and "is_backed" in self.ds.attrs:
            del self.ds.attrs["is_backed"]
        else:
            self.ds.attrs["is_backed"] = isbacked

    @property
    def is_lazy(self) -> bool:
        """Whether any column still defers its data to dask.

        Unlike :attr:`is_backed`, which records where the data came from, this reports
        whether anything is left to compute: a store-backed dataset whose columns were all
        read eagerly (as categoricals always are) is not lazy.
        """
        return any(isinstance(column.data, DaskArray) for _, column in self._items())

    @property
    def index_dim(self) -> str:
        """The underlying computational index i.e., the lone coordinate dimension."""
        if len(self.ds.sizes) != 1:
            msg = f"xarray Dataset should not have more than 1 dims, found {len(self.ds.sizes)} {self.ds.sizes}, {self}"
            raise ValueError(msg)
        if not isinstance(dim := next(iter(self.ds.coords.keys())), str):
            msg = f"Index dimension should be a string, found {dim!r}"
            raise ValueError(msg)
        return dim

    @property
    def true_index_dim(self) -> str:
        """Key of the “true” index.

        Because xarray loads its coordinates/indexes in memory,
        we allow for signaling that a given variable, which is not a coordinate, is the "true" index.

        For example, the true index may be cell names but loading these over an internet connection may not be
        desirable or necessary for most use cases such as getting a quick preview of the columns or loading only
        one column that isn't the index.

        This property is the key of said variable. The default is `index_dim` if this variable has not been set.
        """
        return self.ds.attrs.get("indexing_key", self.index_dim)

    @true_index_dim.setter
    def true_index_dim(self, val: str | None) -> None:
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
        """A :class:`pandas.Index` object corresponding to :attr:`anndata.experimental.backed.Dataset2D.index_dim`.

        :attr:`~anndata.AnnData` internally looks for :attr:`~pandas.DataFrame.index` so this ensures usability.

        Returns
        -------
        The index of the of the dataframe as resolved from :attr:`~xarray.Dataset.coords`.
        """
        return self.ds.indexes[self.index_dim]

    @index.setter
    def index(self, val: object | pd.Index | XDataArray) -> None:
        index_dim = self.index_dim
        if (
            isinstance(val, pd.Index | XDataArray)
            and val.name is not None
            and val.name != index_dim
        ):
            # swap the names of the dimensions out and drop the old index variable, setting `coords` in the process if `val` came from this dataset.
            self._ds = self.ds.swap_dims({index_dim: val.name}).drop_vars(index_dim)
            # swapping dims only changes the name, but not the underlying value i.e., the coordinate, if the underlying value was not present in the dataset.
            # If we were to `__setitem__` on `.coords` without checking, `val` could have the old `index_dim` as its `name` because it was present in the dataset.
            if val.name not in self.ds.coords:
                self.ds.coords[val.name] = val
            self._validate_shape_invariants(self._ds)
        else:
            self.ds.coords[index_dim] = (index_dim, val)
        # without `indexing_key` explicitly set on `self.ds.attrs`, `self.true_index_dim` will use the `self.index_dim`
        if "indexing_key" in self.ds.attrs and (
            hasattr(val, "name") and val.name == self.ds.attrs["indexing_key"]
        ):
            del self.ds.attrs["indexing_key"]

    @property
    def true_xr_index(self) -> XDataArray:
        """The index :class:`~anndata.AnnData` is actually interested in e.g., cell names, for verification."""
        return self.ds[self.true_index_dim]

    @property
    def true_index(self) -> pd.Index:
        """:attr:`~anndata.experimental.backed.Dataset2D.true_xr_index` as a :class:`pandas.Index`."""
        idx = self.true_xr_index.to_index()
        idx.name = self.true_xr_index.name
        return idx

    @property
    def shape(self) -> tuple[int, int]:
        """:attr:`~anndata.AnnData` internally looks for :attr:`~pandas.DataFrame.shape` so this ensures usability.

        Returns
        -------
        The (2D) shape of the dataframe resolved from :attr:`~xarray.Dataset.sizes`.
        """
        return (self.ds.sizes[self.index_dim], len(self.ds))

    @property
    def iloc(self) -> Dataset2DIlocIndexer:
        """:attr:`~anndata.AnnData` internally looks for :attr:`~pandas.DataFrame.iloc` so this ensures usability.

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
    def __getitem__(self, key: Collection[Hashable]) -> Dataset2D: ...
    def __getitem__(
        self, key: Mapping[Any, Any] | Hashable | Iterable[Hashable]
    ) -> Dataset2D | XDataArray:
        ret = self.ds.__getitem__(key)
        if is_empty := (
            isinstance(key, Sized) and len(key) == 0 and not isinstance(key, tuple)
        ):
            ret.coords[self.index_dim] = self.xr_index
        if isinstance(ret, XDataset):
            # If we get an xarray Dataset, we return a Dataset2D
            as_2d = Dataset2D(ret)
            if not is_empty and self.true_index_dim not in [
                *as_2d.columns,
                as_2d.index_dim,
            ]:
                as_2d[self.true_index_dim] = self.true_index
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
        from anndata._core.merge import (
            DS_CONCAT_DUMMY_INDEX_NAME,
            DS_MERGE_DUMMY_INDEX_NAME,
        )

        index_key = self.ds.attrs.get("indexing_key", None)
        all_columns = {*self.columns, *([] if index_key is None else [index_key])}
        # https://github.com/pydata/xarray/issues/10419
        non_nullable_string_cols = {
            col
            for col in all_columns
            if not self[col].attrs.get(NULLABLE_STRING_ATTR, False)
        }
        df = self.ds.to_dataframe()
        for col in all_columns - non_nullable_string_cols:
            df[col] = (
                pandas_as_str(df[col]) if col == index_key else df[col].astype("string")
            )
        if df.index.name != index_key and index_key is not None:
            df = df.set_index(index_key)
        if df.index.name in {
            "_index",
            DS_CONCAT_DUMMY_INDEX_NAME,
            DS_MERGE_DUMMY_INDEX_NAME,
        }:
            df.index.name = None  # matches old AnnData object
        return df

    def to_dask_dataframe(self) -> DaskDataFrame:
        """Convert to a :class:`dask.dataframe.DataFrame` without reading anything.

        Partitions follow the chunking of the lazy columns and are read through
        :meth:`to_memory`, so computing the result gives back what reading the whole thing
        at once would, index and dtypes included.

        Returns
        -------
            :class:`dask.dataframe.DataFrame` deferring to the same storage this does.
        """
        import dask.dataframe as dd

        # An empty slice reports every categorical column as having no categories, which
        # would not match the partitions dask validates against this metadata.
        meta = _partition_to_memory(self, (0, 0)).astype({
            name: dtype
            for name in self.columns
            if isinstance(dtype := self[name].dtype, pd.CategoricalDtype)
        })
        return dd.from_map(
            partial(_partition_to_memory, self), _partition_bounds(self), meta=meta
        )

    @classmethod
    def from_dask_dataframe(cls, frame: DaskDataFrame, *, dim_name: str) -> Self:
        """Wrap a :class:`dask.dataframe.DataFrame`, reading as little as it takes.

        The inverse of :meth:`to_dask_dataframe`. The axis names have to be in memory,
        since xarray keeps coordinates there: they come from the ``dim_name`` column when
        ``frame`` has one and from its index otherwise. Extension-dtype columns are read as
        well, as xarray can only hold those eagerly and rendering them as dask arrays would
        flatten categories and nulls into ``object``. Everything else stays deferred.

        Parameters
        ----------
        frame
            The dataframe to wrap.
        dim_name
            Name of the index dimension, and the column to take the names from if present.

        Returns
        -------
            :class:`Dataset2D` deferring to the same storage ``frame`` does.
        """
        if dim_name in frame.columns:
            index = pd.Index(frame[dim_name].compute().to_numpy(), name=dim_name)
            frame = frame.drop(columns=[dim_name])
        else:
            index = frame.index.compute().rename(dim_name)
        lengths = tuple(frame.map_partitions(len).compute())
        variables = {
            str(column): _dask_variable(
                frame[column], dim_name=dim_name, lengths=lengths
            )
            for column in frame.columns
        }
        return cls(XDataset(variables, coords={dim_name: index}))

    @property
    def columns(self) -> pd.Index:
        """
        :class:`~anndata.AnnData` internally looks for :attr:`~pandas.DataFrame.columns` so this ensures usability

        Returns
        -------
        :class:`pandas.Index` that represents the "columns."
        """
        columns = list(self.ds.keys())
        index_key = self.ds.attrs.get("indexing_key", None)
        if index_key in columns:
            columns.remove(index_key)
        return pd.Index(columns)

    @columns.setter
    def columns(self, val) -> None:
        if len(self.columns.symmetric_difference(val)) > 0:
            msg = "Trying to rename the keys of the mapping with new names - please use a different API to rename the keys of the underlying dataset mapping."
            raise ValueError(msg)
        warn(
            "Renaming or reordering columns on `Dataset2D` has no effect because the underlying data structure has no apparent ordering on its keys",
            UserWarning,
        )

    def __setitem__(
        self, key: Hashable | Iterable[Hashable] | Mapping, value: object
    ) -> None:
        """
        Setting can only be performed when the incoming value is “standalone” like :class:`nump.ndarray` to mimic pandas.
        One can also use the tuple setting style like `ds["foo"] = (ds.index_dim, value)` to set the value, although the index name must match.
        Similarly, one can use the :class:`xarray.DataArray` but it must have the same (one and only one) dim name/coord name as `self.index_dim`.

        For supported setter values see :meth:`xarray.Dataset.__setitem__`.
        """
        if key == self.index_dim:
            msg = f"Cannot set the index dimension {self.index_dim} as if it were a variable. Use `ds.index = ...` instead."
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
    def dtypes(self) -> Mapping[Hashable, np.dtype]:
        """
        Return a Mapping with the dtypes of the variables in the Dataset2D.
        """
        return self.ds.dtypes

    def equals(self, b: object) -> bool:
        """Thin wrapper around :meth:`xarray.Dataset.equals`"""
        if isinstance(b, Dataset2D):
            b = b.ds
        if not isinstance(b, XDataset):
            msg = f"Cannot compare a Dataset2D to {type(b).__name__}"
            raise TypeError(msg)
        return self.ds.equals(b)

    def reindex(
        self,
        index: pd.Index | None = None,
        axis: Literal[0] = 0,
        fill_value: Scalar | None = np.nan,
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


def _dask_variable(
    column: Any,  # noqa: ANN401
    *,
    dim_name: str,
    lengths: tuple[int, ...],
) -> XVariable:
    """One column of a dask dataframe, read eagerly only where xarray cannot defer it."""
    dtype = column.dtype
    if isinstance(dtype, pd.StringDtype):
        # xarray rejects the Arrow-backed array pandas gives a string column, so hold it
        # as NumPy strings and mark it for `to_memory` to convert back.
        return XVariable(
            [dim_name],
            column.compute().to_numpy(dtype=np.dtypes.StringDType(na_object=pd.NA)),
            attrs={NULLABLE_STRING_ATTR: True},
        )
    if pd.api.types.is_extension_array_dtype(dtype):
        return XVariable([dim_name], column.compute().array)
    return XVariable([dim_name], column.to_dask_array(lengths=lengths))


def _partition_bounds(ds: Dataset2D) -> list[tuple[int, int]]:
    """Row ranges to read one at a time, following the chunking of the lazy columns."""
    chunks = next(
        (
            column.data.chunks[0]
            for _, column in ds._items()
            if isinstance(column.data, DaskArray)
            and column.data.ndim == 1
            and not any(np.isnan(chunk) for chunk in column.data.chunks[0])
        ),
        (ds.shape[0],),
    )
    return list(pairwise(accumulate(chunks, initial=0)))


def _partition_to_memory(ds: Dataset2D, bounds: tuple[int, int]) -> pd.DataFrame:
    return ds.iloc[slice(*bounds)].to_memory()


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


###################
# Narwhals plugin
###################

# Narwhals finds this module through the `narwhals.plugins` entry point in `pyproject.toml`
# and looks up the three names below on it. Nothing in anndata calls them.

NATIVE_PACKAGE = "anndata"


def is_native(native_object: object, /) -> bool:
    """Return whether ``native_object`` is a :class:`Dataset2D`."""
    return isinstance(native_object, Dataset2D)


class Dataset2DNamespace:
    """Routes a :class:`Dataset2D` to a compliant frame over one of its conversions.

    Narwhals already supports both dataframes a ``Dataset2D`` converts to, so we hand it a
    compliant frame over one of those rather than reimplementing the protocol.
    :attr:`~Dataset2D.is_lazy` picks which: a lazy dataset goes through
    :meth:`~Dataset2D.to_dask_dataframe` and becomes a :class:`narwhals.LazyFrame`, so
    wrapping it reads nothing, while an in-memory one goes through
    :meth:`~Dataset2D.to_memory` and becomes a :class:`narwhals.DataFrame` keeping its row
    labels as the index, recoverable via :func:`narwhals.maybe_get_index`.

    Both ``from_native`` constructors read ``_implementation`` and ``_version`` off their
    ``context``, so this namespace doubles as that context.
    """

    _implementation: nw.Implementation = nw.Implementation.PANDAS

    def __init__(self, *, version: Version) -> None:
        self._version = version

    def from_native(
        self, native_object: Dataset2D, /
    ) -> PandasLikeDataFrame | DaskLazyFrame:
        if native_object.is_lazy:
            from narwhals._dask.dataframe import DaskLazyFrame

            return DaskLazyFrame.from_native(
                native_object.to_dask_dataframe(), context=self
            )
        from narwhals._pandas_like.dataframe import PandasLikeDataFrame

        return PandasLikeDataFrame.from_native(native_object.to_memory(), context=self)


def __narwhals_namespace__(version: Version) -> Dataset2DNamespace:
    """Return the compliant namespace narwhals uses to wrap a :class:`Dataset2D`."""
    return Dataset2DNamespace(version=version)
