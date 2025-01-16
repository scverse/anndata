from __future__ import annotations

from importlib.util import find_spec
from typing import TYPE_CHECKING

from ._core.sparse_dataset import sparse_dataset
from ._io.h5ad import read_h5ad, write_h5ad
from ._io.read import (
    read_csv,
    read_excel,
    read_hdf,
    read_loom,
    read_mtx,
    read_text,
    read_umi_tools,
)
from ._io.specs import read_elem, write_elem
from ._io.write import write_csvs, write_loom

if find_spec("zarr") or TYPE_CHECKING:
    from ._io.zarr import read_zarr, write_zarr
else:  # pragma: no cover

    def read_zarr(*args, **kw):
        msg = "zarr is not installed"
        raise ImportError(msg)

    def write_zarr(*args, **kw):
        msg = "zarr is not installed"
        raise ImportError(msg)


__all__ = [
    "read_csv",
    "read_excel",
    "read_h5ad",
    "read_hdf",
    "read_loom",
    "read_mtx",
    "read_text",
    "read_umi_tools",
    "read_zarr",
    "write_csvs",
    "write_h5ad",
    "write_loom",
    "write_zarr",
    "write_elem",
    "read_elem",
    "sparse_dataset",
]
