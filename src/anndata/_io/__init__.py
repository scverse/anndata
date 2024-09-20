from __future__ import annotations

import warnings

__all__: list[str] = []


def __getattr__(key: str):
    from .._core.sparse_dataset import sparse_dataset
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
    from .specs import read_elem, write_elem
    from .write import write_csvs, write_loom

    def write_zarr(*args, **kw):
        from .zarr import write_zarr

        return write_zarr(*args, **kw)

    # We use this in test by attribute access
    from . import specs  # noqa: F401, E402

    name_to_implementation = {
        "sparse_dataset": sparse_dataset,
        "read_h5ad": read_h5ad,
        "write_h5ad": write_h5ad,
        "read_csv": read_csv,
        "read_excel": read_excel,
        "read_hdf": read_hdf,
        "read_loom": read_loom,
        "read_mtx": read_mtx,
        "read_text": read_text,
        "read_umi_tools": read_umi_tools,
        "read_zarr": read_zarr,
        "read_elem": read_elem,
        "write_elem": write_elem,
        "write_csvs": write_csvs,
        "write_loom": write_loom,
        "write_zarr": write_zarr,
        "specs": specs,
    }
    if key not in name_to_implementation:
        raise AttributeError(f"module {__name__!r} has no attribute {key!r}")
    warnings.warn(
        f"Importing {key} from `anndata._io` is deprecated. "
        "Please use anndata.io instead.",
        FutureWarning,
    )
    return name_to_implementation[key]
