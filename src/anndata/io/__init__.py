from __future__ import annotations

import warnings

with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore", message=r"Importing read_.* from `anndata` is deprecated"
    )
    from .._io import (
        read_csv,
        read_elem,
        read_excel,
        read_h5ad,
        read_hdf,
        read_loom,
        read_mtx,
        read_text,
        read_umi_tools,
        read_zarr,
        sparse_dataset,
        write_csvs,
        write_elem,
        write_h5ad,
        write_loom,
        write_zarr,
    )

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
