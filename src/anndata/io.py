from __future__ import annotations

import os

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


def read_zarr(*args, **kw):
    from ._io.zarr import read_zarr

    return read_zarr(*args, **kw)


def write_zarr(*args, **kw):
    from ._io.zarr import write_zarr

    return write_zarr(*args, **kw)


if os.environ.get("READTHEDOCS", "") != "True":
    # This makes the lazy import above into a direct one. Importing zarr is slow.
    from ._io.zarr import read_zarr, write_zarr  # noqa: TCH004


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
