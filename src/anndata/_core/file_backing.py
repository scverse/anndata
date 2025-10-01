from __future__ import annotations

import weakref
from collections.abc import Mapping
from functools import singledispatch
from pathlib import Path, PurePosixPath
from typing import TYPE_CHECKING

import h5py

from ..compat import AwkArray, DaskArray, ZarrArray, ZarrGroup
from .sparse_dataset import BaseCompressedSparseDataset
from .xarray import Dataset2D

if TYPE_CHECKING:
    from collections.abc import Iterator
    from os import PathLike
    from typing import Literal

    from .._types import ArrayStorageType
    from . import anndata


class AnnDataFileManager:
    """Backing file manager for AnnData."""

    def __init__(
        self,
        adata: anndata.AnnData,
        filename: PathLike[str] | str | None = None,
        filemode: Literal["r", "r+"] | None = None,
    ):
        self._adata_ref = weakref.ref(adata)
        self.filename = filename
        self._filemode = filemode
        self._file = None
        if filename:
            self.open()

    def __getstate__(self):
        state = self.__dict__.copy()
        state["_adata_ref"] = state["_adata_ref"]()
        return state

    def __setstate__(self, state):
        self.__dict__ = state.copy()
        self.__dict__["_adata_ref"] = weakref.ref(state["_adata_ref"])

    @property
    def _adata(self):
        return self._adata_ref()

    def __repr__(self) -> str:
        if self.filename is None:
            return "Backing file manager: no file is set."
        else:
            return f"Backing file manager of file {self.filename}."

    def __contains__(self, x) -> bool:
        return x in self._file

    def __iter__(self) -> Iterator[str]:
        return iter(self._file)

    def __getitem__(
        self, key: str
    ) -> h5py.Group | h5py.Dataset | BaseCompressedSparseDataset:
        return self._file[key]

    def __setitem__(
        self,
        key: str,
        value: h5py.Group | h5py.Dataset | BaseCompressedSparseDataset,
    ):
        self._file[key] = value

    def __delitem__(self, key: str):
        del self._file[key]

    @property
    def filename(self) -> Path:
        return self._filename

    @filename.setter
    def filename(self, filename: PathLike[str] | str | None):
        self._filename = None if filename is None else Path(filename)

    def open(
        self,
        filename: PathLike[str] | str | None = None,
        filemode: Literal["r", "r+"] | None = None,
    ):
        if filename is not None:
            self.filename = filename
        if filemode is not None:
            self._filemode = filemode
        if self.filename is None:
            msg = "Cannot open backing file if backing not initialized."
            raise ValueError(msg)
        self._file = h5py.File(self.filename, self._filemode)

    def close(self):
        """Close the backing file, remember filename, do *not* change to memory mode."""
        if self._file is not None:
            self._file.close()

    def _to_memory_mode(self):
        """Close the backing file, forget filename, *do* change to memory mode."""
        self._adata._X = self._adata.X[()]
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
def to_memory(x, *, copy: bool = False):
    """Permissivley convert objects to in-memory representation.

    If they already are in-memory, (or are just unrecognized) pass a copy through.
    """
    if copy and hasattr(x, "copy"):
        return x.copy()
    else:
        return x


@to_memory.register(ZarrArray)
@to_memory.register(h5py.Dataset)
def _(x: ArrayStorageType, *, copy: bool = False):
    return x[...]


@to_memory.register(BaseCompressedSparseDataset)
def _(x: BaseCompressedSparseDataset, *, copy: bool = False):
    return x.to_memory()


@to_memory.register(DaskArray)
def _(x: DaskArray, *, copy: bool = False):
    return x.compute()


@to_memory.register(Mapping)
def _(x: Mapping, *, copy: bool = False):
    return {k: to_memory(v, copy=copy) for k, v in x.items()}


@to_memory.register(AwkArray)
def _(x: AwkArray, *, copy: bool = False):
    from copy import copy as _copy

    if copy:
        return _copy(x)
    else:
        return x


@to_memory.register(Dataset2D)
def _(x: Dataset2D, *, copy: bool = False):
    return x.to_memory(copy=copy)


@singledispatch
def filename(x):
    msg = f"Not implemented for {type(x)}"
    raise NotImplementedError(msg)


@filename.register(h5py.Group)
@filename.register(h5py.Dataset)
def _(x):
    return x.file.filename


@filename.register(ZarrArray)
@filename.register(ZarrGroup)
def _(x):
    return x.store.path


@singledispatch
def get_elem_name(x):
    msg = f"Not implemented for {type(x)}"
    raise NotImplementedError(msg)


@get_elem_name.register(h5py.Group)
def _(x):
    return x.name


@get_elem_name.register(ZarrGroup)
def _(x):
    return PurePosixPath(x.path).name
