from __future__ import annotations

from typing import Annotated, Literal

import scverse_misc
from pydantic import Field

from ._repr_constants import (
    DEFAULT_FOLD_THRESHOLD,
    DEFAULT_MAX_CATEGORIES,
    DEFAULT_MAX_DEPTH,
    DEFAULT_MAX_FIELD_WIDTH,
    DEFAULT_MAX_ITEMS,
    DEFAULT_MAX_LAZY_CATEGORIES,
    DEFAULT_MAX_README_SIZE,
    DEFAULT_TYPE_WIDTH,
    DEFAULT_UNIQUE_LIMIT,
)


class Settings(scverse_misc.Settings):
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

    # HTML representation settings

    repr_html_enabled: bool = True
    """Whether to use rich HTML representation in Jupyter notebooks. Set to False to use plain text repr."""

    repr_html_fold_threshold: int = DEFAULT_FOLD_THRESHOLD
    """Auto-fold sections in HTML repr when they have more than this many entries."""

    repr_html_max_depth: int = DEFAULT_MAX_DEPTH
    """Maximum recursion depth for nested AnnData objects in HTML repr."""

    repr_html_max_items: int = DEFAULT_MAX_ITEMS
    """Maximum number of items to show per section in HTML repr."""

    repr_html_max_categories: int = DEFAULT_MAX_CATEGORIES
    """Maximum number of category values to display inline in HTML repr."""

    repr_html_max_lazy_categories: int = DEFAULT_MAX_LAZY_CATEGORIES
    """Maximum categories to load for lazy categoricals in HTML repr. For lazy AnnData (from :func:`~anndata.experimental.read_lazy`), loading categories requires reading from disk. This limit prevents loading too many categories. Set to 0 to disable loading categories entirely (metadata-only mode)."""

    repr_html_unique_limit: int = DEFAULT_UNIQUE_LIMIT
    """Maximum number of rows to compute unique counts for in HTML repr. Set to 0 to disable."""

    repr_html_dataframe_expand: bool = False
    """Whether to show expandable pandas DataFrame previews in HTML repr. When enabled, DataFrames in `obsm`/`varm` can be expanded to show their content using pandas `_repr_html_()` (rich Jupyter-style output). Configure pandas display options to control output: `pd.set_option('display.max_rows', 10)`"""

    repr_html_max_field_width: int = DEFAULT_MAX_FIELD_WIDTH
    """Maximum width in pixels for the field name column in HTML repr."""

    repr_html_type_width: int = DEFAULT_TYPE_WIDTH
    """Width in pixels for the type column in HTML repr."""

    repr_html_max_readme_size: int = DEFAULT_MAX_README_SIZE
    """Maximum size in characters for README content in HTML repr. READMEs larger than this will be truncated with a note. Set to 0 to disable truncation (not recommended for very large READMEs)."""


settings = Settings()
