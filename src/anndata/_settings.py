from __future__ import annotations

from enum import Enum
from typing import Annotated

import scverse_misc
from packaging.version import Version
from pydantic import Field
from scverse_misc import Deprecation, deprecated


class WriteCompat(Version, Enum):
    """Which anndata version’s write behavior to target.

    Each member is the oldest anndata version whose I/O behavior we promise to be
    compatible with: encodings it cannot read are disallowed,
    and restrictions added after it are relaxed.
    See :doc:`/write-compat` for the whole table.
    """

    V0_12 = "0.12"
    """anndata 0.12 (2025-07-16).

    Relaxes the restrictions added in 0.13, i.e. allows writing
    3-dimensional :attr:`~anndata.AnnData.X`/:attr:`~anndata.AnnData.layers` and
    `/` in `h5ad` keys of `obs`, `var` and `uns`.
    """

    V0_13 = "0.13"
    """anndata 0.13 (2026-07-07).

    No relaxations. This is the default.
    """


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

    write_compat: WriteCompat = WriteCompat.V0_13
    """Which anndata version’s write behavior to target, see :class:`~anndata.WriteCompat`.

    Write functions accept a `compat` argument to override this for a single call."""

    disallow_forward_slash_in_h5ad: Annotated[
        bool | None,
        Field(
            deprecated=deprecated(
                Deprecation(
                    "0.14", "This will be removed in 0.15, use `write_compat` instead."
                )
            ),
        ),
    ] = None
    """Whether or not to disallow the `/` character in keys for h5ad files.

    `None` derives it from :attr:`write_compat`, i.e. disallows it for `"0.13"` and newer."""

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


settings = Settings()


def forward_slash_disallowed() -> bool:
    """Whether `/` is disallowed in `h5ad` keys, honoring the deprecated setting."""
    # bypass the descriptor, as our own reads shouldn’t emit its deprecation warning
    if (explicit := settings.__dict__["disallow_forward_slash_in_h5ad"]) is not None:
        return explicit
    return settings.write_compat >= WriteCompat.V0_13
