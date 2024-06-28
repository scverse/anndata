from __future__ import annotations

from .h5ad import read_h5ad, write_h5ad
from .read import (
    read_csv,
    read_excel,
    read_hdf,
    read_loom,
    read_mtx,
    read_text,
    read_umi_tools,
    read_zarr,
)
from .write import write_csvs, write_loom


def write_zarr(*args, **kw):
    from .zarr import write_zarr

    return write_zarr(*args, **kw)


# We use this in test by attribute access
from . import specs  # noqa: F401, E402

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
]
