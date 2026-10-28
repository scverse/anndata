from __future__ import annotations

from typing import Annotated

import scverse_misc
from pydantic import Field


class Settings(scverse_misc.Settings):
    remove_unused_categories: bool = True
    """Whether or not to remove unused categories with :class:`~pandas.Categorical`."""

    check_uniqueness: bool = True
    """Whether or not to check uniqueness of the `obs` indices on `__init__` of :class:`~anndata.AnnData`."""

    allow_write_nullable_strings: bool | None = None
    """Whether or not to allow writing of `pd.arrays.[Arrow]StringArray`. When set to `None`, it will be inferred from `pd.options.future.infer_string`. When set to `False` explicitly, we will try writing `string` arrays in the old, non-nullable format."""

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

    validate_before_write: bool = False
    """\
    Whether :func:`~anndata.io.write_h5ad` and :func:`~anndata.io.write_zarr` check every element
    for a writer (see :meth:`~anndata.AnnData.unwriteable`) before creating the file or store.

    When `True`, nothing is written if a problem is found: a single problem is raised as the
    exception writing it would have raised, several are raised together as an
    :exc:`ExceptionGroup` whose message lists their paths (catch with ``except*``).
    The check is best effort: it catches elements and columns without an on-disk encoding,
    unconvertible object columns, non-2D ``X``/layers and dataframe structure problems,
    but backend, codec and IO failures and writers that reject particular values can still fail
    mid-write. It reads no data from backed files and computes nothing, but object-dtype
    dataframe columns are scanned and converted once more (allocating) to find out whether
    they are writable.
    When `False`, writing stops at the first unwritable element, reports only that one,
    and may leave a partially written file.
    """

    restrict_index_types: bool = True
    """
    Whether to force coercion to a string index upon declaration of the `AnnData` object or setting `obs`/`var`.
    "Setting this to `False` will e.g. also allow `MultiIndex` indexes upon declaration/setting.
    "Only integer indices i.e., those caught by :func:`pandas.api.types.is_integer_dtype` will always be converted to strings.
    """


settings = Settings()
