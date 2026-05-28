from __future__ import annotations

from typing import Annotated, Literal

import scverse_misc
from pydantic import Field


class Settings(
    scverse_misc.Settings, exported_object_name="settings", docstring_style="scverse"
):
    remove_unused_categories: bool = True
    """Whether or not to remove unused categories with :class:`~pandas.Categorical`."""

    check_uniqueness: bool = True
    """Whether or not to check uniqueness of the `obs` indices on `__init__` of :class:`~anndata.AnnData`."""

    copy_on_write_X: Annotated[
        bool,
        Field(
            deprecated="This will be removed in 0.14 (deprecated in 0.13) and copy-on-write will be default",
        ),
    ] = True
    """Whether to copy-on-write X. Currently `my_adata_view[subset].X = value` will write back to the original AnnData object at the `subset` location. `X` is the only element where this behavior is implemented though."""

    allow_write_nullable_strings: bool | None = None
    """Whether or not to allow writing of `pd.arrays.[Arrow]StringArray`. When set to `None`, it will be inferred from `pd.options.future.infer_string`. When set to `False` explicitly, we will try writing `string` arrays in the old, non-nullable format."""

    zarr_write_format: Literal[2, 3] = 3
    """Which version of zarr to write to when anndata must internally open a write-able zarr group."""

    use_sparse_array_on_read: bool = False
    """Whether or not to use :class:`scipy.sparse.sparray` as the default class when reading in data"""

    min_rows_for_chunked_h5_copy: Annotated[int, Field(gt=0)] = 1000
    """Minimum number of rows at a time to copy when writing out an H5 Dataset to a new location"""

    disallow_forward_slash_in_h5ad: bool = True
    """Whether or not to disallow the `/` character in keys for h5ad files"""

    write_csr_csc_indices_with_min_possible_dtype: bool = False
    """Write a csr or csc matrix with the minimum possible data type for `indices`, always unsigned integer."""

    auto_shard_zarr_v3: bool | None = True
    """Whether or not to use zarr's auto computation of sharding for v3.  For v2 this setting will be ignored. The setting will apply to all calls to anndata's writing mechanism (write_zarr / write_elem) and will **not** override any user-defined kwargs for shards."""

    restrict_index_types: bool = True
    """
    Whether to force coercion to a string index upon declaration of the `AnnData` object or setting `obs`/`var`.
    "Setting this to `False` will e.g. also allow `MultiIndex` indexes upon declaration/setting.
    "Only integer indices i.e., those caught by :func:`pandas.api.types.is_integer_dtype` will always be converted to strings.
    """

    parallel_h5_read: bool = True
    """
    Enables substantially faster reads of h5ad files via multi-threading. Currently supports
    these compression codecs: HDF5's gzip (zlib), zstd, and blosc with any inner codec. Will
    fall back to standard single-threading with a warning if h5ad file compression method
    is not supported or an error occurs. Set this parameter to False to restrict h5ad file reads
    to a single thread.
    """

    parallel_h5_read_min_mb: Annotated[int, Field(ge=0)] = 64
    """
    Minimum uncompressed array size (MiB) below which the serial HDF5 read
    path is used. Default 64 MiB: avoids thread-pool overhead on the many
    small ``obs``/``var`` arrays a normal ``read_h5ad`` traverses; gains from
    multi-threaded decompression come from ``X``/``layers`` arrays, which are
    typically far above this threshold.
    """

    parallel_h5_read_workers: Annotated[int | None, Field(ge=1)] = None
    """
    Thread-pool size for parallel HDF5 reads. ``None`` (default) chooses
    ``min(os.cpu_count(), 16)``. For long-running jobs on machines with
    more than 16 logical cores and high disk bandwidth, consider raising this
    to the logical core count and reducing until performance peaks.
    """


settings = Settings()
