from __future__ import annotations

from collections.abc import Iterator, Mapping
from functools import singledispatch
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import h5py

from ..compat import AwkArray, DaskArray, ZarrArray, ZarrGroup
from .sparse_dataset import BaseCompressedSparseDataset

if TYPE_CHECKING:
    from os import PathLike

    from . import anndata


class AnnDataFileManager:
    """Backing file manager for AnnData."""

    def __init__(
        self,
        adata: anndata.AnnData,
        filename: PathLike | None = None,
        filemode: Literal["r", "r+"] | None = None,
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
    def filename(self, filename: PathLike | None):
        self._filename = None if filename is None else Path(filename)

    def open(
        self,
        filename: PathLike | None = None,
        filemode: Literal["r", "r+"] | None = None,
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


@to_memory.register(BaseCompressedSparseDataset)
def _(x: BaseCompressedSparseDataset, copy=True):
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


@singledispatch
def filename(x):
    raise NotImplementedError(f"Not implemented for {type(x)}")


@filename.register(h5py.Group)
@filename.register(h5py.Dataset)
def _(x):
    return x.file.filename


@filename.register(ZarrArray)
@filename.register(ZarrGroup)
def _(x):
    return x.store.path
