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
    """
    Whether :func:`~anndata.io.write_h5ad` and :func:`~anndata.io.write_zarr` check every element like :meth:`~anndata.AnnData.unwriteable` before writing anything (with `use=`, only the elements written), and raise all problems at once.
    One problem raises its own exception, of the type the writer would raise: :exc:`TypeError` for an unconvertible column or an invalid key (the writer raises :exc:`AttributeError` for a non-`str` key), :exc:`ValueError` for a structural problem or a reference cycle (the writer recurses until :exc:`RecursionError`), :exc:`RuntimeError` for nullable strings that need :attr:`anndata.settings.allow_write_nullable_strings`, and anndata's registry error for a missing writer.
    Several raise an :exc:`ExceptionGroup` whose members keep their types, so `except* TypeError` and `except* ValueError` can handle some of them; the rest, such as missing writers, propagate.
    """

    restrict_index_types: bool = True
    """
    Whether to force coercion to a string index upon declaration of the `AnnData` object or setting `obs`/`var`.
    "Setting this to `False` will e.g. also allow `MultiIndex` indexes upon declaration/setting.
    "Only integer indices i.e., those caught by :func:`pandas.api.types.is_integer_dtype` will always be converted to strings.
    """


settings = Settings()
