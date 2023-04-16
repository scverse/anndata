from functools import singledispatch
from os import PathLike
from pathlib import Path
from typing import Optional, Union, Iterator, Literal
from collections.abc import Mapping

import h5py
import pandas as pd

from . import anndata
from .sparse_dataset import SparseDataset
from ..compat import ZarrArray, DaskArray, AwkArray


class AnnDataFileManager:
    """Backing file manager for AnnData."""

    def __init__(
        self,
        adata: "anndata.AnnData",
        filename: Optional[PathLike] = None,
        filemode: Optional[Literal["r", "r+"]] = None,
    ):
        self._adata = adata
        self.filename = filename
        self._filemode = filemode
        self._file = None
        if filename:
            self.open()

    def __repr__(self) -> str:
        if self.filename is None:
            return "Backing file manager: no file is set."
        else:
            return f"Backing file manager of file {self.filename}."

    def __contains__(self, x) -> bool:
        return x in self._file

    def __iter__(self) -> Iterator[str]:
        return iter(self._file)

    def __getitem__(self, key: str) -> Union[h5py.Group, h5py.Dataset, SparseDataset]:
        return self._file[key]

    def __setitem__(
        self, key: str, value: Union[h5py.Group, h5py.Dataset, SparseDataset]
    ):
        self._file[key] = value

    def __delitem__(self, key: str):
        del self._file[key]

    @property
    def filename(self) -> Path:
        return self._filename

    @filename.setter
    def filename(self, filename: Optional[PathLike]):
        self._filename = None if filename is None else Path(filename)

    def open(
        self,
        filename: Optional[PathLike] = None,
        filemode: Optional[Literal["r", "r+"]] = None,
    ):
        if filename is not None:
            self.filename = filename
        if filemode is not None:
            self._filemode = filemode
        if self.filename is None:
            raise ValueError("Cannot open backing file if backing not initialized.")
        self._file = h5py.File(self.filename, self._filemode)

    def close(self):
        """Close the backing file, remember filename, do *not* change to memory mode."""
        if self._file is not None:
            self._file.close()

    def _to_memory_mode(self):
        """Close the backing file, forget filename, *do* change to memory mode."""
        self._adata.__X = self._adata.X[()]
        self._file.close()
        self._file = None
        self._filename = None

    @property
    def is_open(self) -> bool:
        """State of backing file."""
        if self._file is None:
            return False
        # try accessing the id attribute to see if the file is open
        return bool(self._file.id)


@singledispatch
def to_memory(x, copy=False):
    """Permissivley convert objects to in-memory representation.

    If they already are in-memory, (or are just unrecognized) pass a copy through.
    """
    if copy and hasattr(x, "copy"):
        return x.copy()
    else:
        return x


@to_memory.register(ZarrArray)
@to_memory.register(h5py.Dataset)
def _(x, copy=False):
    return x[...]


@to_memory.register(SparseDataset)
def _(x: SparseDataset, copy=False):
    return x.to_memory()


@to_memory.register(DaskArray)
def _(x, copy=False):
    return x.compute()


@to_memory.register(Mapping)
def _(x: Mapping, copy=False):
    return {k: to_memory(v, copy=copy) for k, v in x.items()}


@to_memory.register(AwkArray)
def _(x, copy=False):
    from copy import copy as _copy

    if copy:
        return _copy(x)
    else:
        return x


class HDF5DataFrame:
    def __init__(self, group: h5py.Group):
        assert (
            group.attrs["encoding-type"] == "dataframe"
        ), "HDF5 group at path '{group.name}' is not encoded as a dataframe"

        self._group = group
        self._attrs = self._group.attrs

        self._index = self._group[self._attrs["_index"]].asstr()
        self.columns = pd.Index(self._attrs["column-order"])

        for column in self.columns:
            setattr(self, column, self._group[column])

    @property
    def index(self):
        return pd.Index(self._index[:])

    # def __getitem__(self, key):
    #     if isinstance(key, str) and key in self.columns:
    #         return self._group[key]

    #     elif isinstance(key, slice):
    #         return self._group[key]

    #     if isinstance(index, tuple) and self.attr in ("obs", "obsm"):
    #         oidx = index[0]
    #         if len(index) > 1:
    #             vidx = index[1]

    #     if oidx is None:
    #         view = self.adset[index]
    #     else:
    #         view = self.adset[oidx]
    #     attr_arr = getattr(view, self.attr)
    #     if self.key is not None:
    #         attr_arr = attr_arr[self.key]
    #     return attr_arr if vidx is None else attr_arr[:, vidx]

    # @property
    # def shape(self):
    #     shape = self.adset.shape
    #     if self.attr in ["X", "layers"]:
    #         return shape
    #     elif self.attr == "obs":
    #         return (shape[0],)
    #     elif self.attr == "obsm" and self.key is not None:
    #         return shape[0], self[:1].shape[1]
    #     else:
    #         return None

    # @property
    # def ndim(self):
    #     return len(self.shape) if self.shape is not None else 0

    # @property
    # def dtype(self):
    #     _dtypes = self.adset._dtypes
    #     if _dtypes is not None and self.attr in _dtypes:
    #         return _dtypes[self.attr][self.key]

    #     attr = self[:1]
    #     if hasattr(attr, "dtype"):
    #         return attr.dtype
    #     else:
    #         return None
