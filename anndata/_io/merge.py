from .specs import read_elem, write_elem
import zarr
from .._core.merge import (
    _resolve_dim,
    StrategiesLiteral,
    intersect_keys,
    unify_categorical_dtypes,
    merge_indices,
    merge_dataframes,
    resolve_merge_strategy,
)
from .._core.sparse_dataset import SparseDataset
from ..experimental import write_dispatched, read_dispatched
from scipy import sparse
from scipy.sparse import spmatrix
from collections import OrderedDict
from collections.abc import Mapping, MutableSet
from functools import reduce, singledispatch
import pandas as pd
from pathlib import Path
from typing import (
    Any,
    Callable,
    Collection,
    Iterable,
    Optional,
    Tuple,
    Set,
    Sequence,
    TypeVar,
    Union,
    Literal,
    MutableMapping,
)
import typing

import numpy as np
from scipy.sparse import spmatrix

from ..compat import AwkArray, DaskArray, ZarrGroup, ZarrArray


def _df_index(df: ZarrGroup) -> pd.Index:
    index_key = df.attrs["_index"]
    return pd.Index(read_elem(df[index_key]))


def _has_same_attrs(groups: list[ZarrGroup], path: str, attrs: Set[str]) -> bool:
    if len(groups) == 0:
        return True

    init_g = groups[0][path]

    return all(
        all(g[path].attrs[attr] == init_g.attrs[attr] for attr in attrs) for g in groups
    )


def _attrs_equal(
    groups: list[ZarrGroup], path: str, attrs_map: Mapping[str, str]
) -> bool:
    if len(groups) == 0:
        raise ValueError("List should not be empty")

    init_g = groups[0][path]

    return all(
        init_g.attrs[attr] == val for attr, val in attrs_map.items()
    ) and _has_same_attrs(groups, path, attrs_map.keys())


def _index_equal(groups: list[ZarrGroup], path) -> bool:
    names = _df_index(groups[0][path])
    for g in groups[1:]:
        curr_names = _df_index(g[path])
        if not np.array_equal(names, curr_names):
            return False
    return True


def write_concat_mappings_aligned(mappings, output_group: ZarrGroup, keys, path, axis=0, index=None):
    mapping_group = output_group.create_group(path)
    mapping_group.attrs.update(
        {
            "encoding-type": "dict",
            "encoding-version": "0.1.0",
        }
    )
    for k in keys:
        elems = [m[k] for m in mappings]
        write_concat_sequence_aligned(
            elems, output_group=mapping_group, out_path=k, axis=axis, index=index)


ARRAY_LIKE = {
    "array",
    "dataframe",
    "csc_matrix",
    "csr_matrix",
    "awkward-array",
}

SPARSE_MATRIX = {"csc_matrix", "csr_matrix"}

EAGER_TYPES = {"dataframe", "awkward-array"}


def read_group(group: ZarrGroup):

    def callback(func, elem_name: str, elem, iospec):
        if iospec.encoding_type in SPARSE_MATRIX:
            return SparseDataset(elem)
        elif iospec.encoding_type in EAGER_TYPES:
            return read_elem(elem)
        elif iospec.encoding_type == "array":
            return elem
        else:
            return func(elem)

    return read_dispatched(group, callback=callback)


def read_only_dict(group: ZarrGroup):
    def callback(func, elem_name: str, elem, iospec):
        if iospec.encoding_type == "dict":
            return {
                k: v for k, v in elem.items()
            }
        else:
            return elem
    return read_dispatched(group, callback=callback)


def dim_indices(group: ZarrGroup, *, axis=None, dim=None) -> pd.Index:
    """Helper function to get adata.{dim}_names."""
    _, dim = _resolve_dim(axis=axis, dim=dim)
    return group[dim].attrs["_index"]


def get_encoding_type(group: ZarrGroup) -> str:
    return group.attrs.get("encoding-type")


def get_shape(group: ZarrGroup) -> str:
    return group.attrs.get("shape")


def default_fill_value_group(groups: Sequence[ZarrGroup]):
    # TODO refer to original one
    """Given some arrays, returns what the default fill value should be.

    This is largely due to backwards compat, and might not be the ideal solution.
    """
    if any(get_encoding_type(g) in SPARSE_MATRIX for g in groups):
        return 0
    else:
        return np.nan


def write_concat_sequence_aligned(groups: Sequence[ZarrGroup], output_group,
                                  out_path, axis=0, index=None):
    """
    array, dataframe, csc_matrix, csc_matrix, awkward-array
    """

    if any(get_encoding_type(g) == "dataframe" for g in groups):
        # TODO: This is hacky, 0 is a sentinel for outer_concat_aligned_mapping
        # TODO: what is missing val thing?
        if not all(
            get_encoding_type(g) == "dataframe" or 0 in get_shape(g)
            for g in groups
        ):
            raise NotImplementedError(
                "Cannot concatenate a dataframe with other array types."
            )
        # TODO: behaviour here should be chosen through a merge strategy
        df = pd.concat(
            unify_categorical_dtypes([read_group(g) for g in groups]),
            ignore_index=True,
            axis=axis,
        )
        df.index = index
        write_elem(output_group, out_path, df)

    elif any(get_encoding_type(g) == "awkard-array" for g in groups):
        from ..compat import awkward as ak

        if not all(
            get_encoding_type(g) == "awkard-array" or 0 in get_shape(g) for g in groups
        ):
            raise NotImplementedError(
                "Cannot concatenate an AwkwardArray with other array types."
            )

        res = ak.concatenate([read_group(g) for g in groups], axis=axis)
        write_elem(output_group, out_path, res)
    # If all are compatible sparse matrices
    elif (all(get_encoding_type(g) == "csc_matrix" for g in groups) and axis == 1) or \
            (all(get_encoding_type(g) == "csr_matrix" for g in groups) and axis == 0):

        datasets: Sequence[SparseDataset] = [read_group(g) for g in groups]
        write_elem(output_group, out_path, datasets[0])
        out_dataset: SparseDataset = read_group(output_group[out_path])
        for ds in datasets[1:]:
            out_dataset.append(ds)
    # If all are arrays
    elif all(get_encoding_type(g) == "array" for g in groups):
        arrays: Sequence[ZarrArray] = [g for g in groups]
        output_group.create_dataset(out_path, data=arrays[0])
        out_array: ZarrArray = output_group[out_path]
        for arr in arrays[1:]:
            out_array.append(arr)
    else:
        raise NotImplementedError(
            f"Concatenation of these types is not yet implemented: {[get_encoding_type(g) for g in groups],axis}."
        )


def concat_on_disk(
    in_files: Union[Collection[str], typing.MutableMapping],
    out_file: Union[str, typing.MutableMapping],
    overwrite: bool = False,
    *,
    axis: Literal[0, 1] = 0,
    join: Literal["inner", "outer"] = "inner",
    merge: Union[StrategiesLiteral, Callable, None] = None,
    uns_merge: Union[StrategiesLiteral, Callable, None] = None,
    label: Optional[str] = None,
    keys: Optional[Collection] = None,
    index_unique: Optional[str] = None,
    fill_value: Optional[Any] = None,
    pairwise: bool = False,
):
    """Concatenates multiple AnnData objects along a specified axis using their
    corresponding stores or paths, and writes the resulting AnnData object
    to a target location on disk.

    Unlike the `concat` function, this method does not require
    loading the input AnnData objects into memory,
    making it a memory-efficient alternative for large datasets.
    The resulting object written to disk should be equivalent
    to the concatenation of the loaded AnnData objects using
    the `concat` function.



    Params
    ------
    in_files
        The corresponding stores or paths of AnnData objects to
        be concatenated. If a Mapping is passed, keys are used for the `keys`
        argument and values are concatenated.
    out_file
        The target path or store to write the result in.
    overwrite
        If False the and if a file already exists it will raise an error,
        otherwise it will overwrite.
    axis
        Which axis to concatenate along.
    join
        How to align values when concatenating. If "outer", the union of the other axis
        is taken. If "inner", the intersection. See :doc:`concatenation <../concatenation>`
        for more.
    merge
        How elements not aligned to the axis being concatenated along are selected.
        Currently implemented strategies include:

        * `None`: No elements are kept.
        * `"same"`: Elements that are the same in each of the objects.
        * `"unique"`: Elements for which there is only one possible value.
        * `"first"`: The first element seen at each from each position.
        * `"only"`: Elements that show up in only one of the objects.
    uns_merge
        How the elements of `.uns` are selected. Uses the same set of strategies as
        the `merge` argument, except applied recursively.
    label
        Column in axis annotation (i.e. `.obs` or `.var`) to place batch information in.
        If it's None, no column is added.
    keys
        Names for each object being added. These values are used for column values for
        `label` or appended to the index if `index_unique` is not `None`. Defaults to
        incrementing integer labels.
    index_unique
        Whether to make the index unique by using the keys. If provided, this
        is the delimeter between "{orig_idx}{index_unique}{key}". When `None`,
        the original indices are kept.
    fill_value
        When `join="outer"`, this is the value that will be used to fill the introduced
        indices. By default, sparse arrays are padded with zeros, while dense arrays and
        DataFrames are padded with missing values.
    pairwise
        Whether pairwise elements along the concatenated dimension should be included.
        This is False by default, since the resulting arrays are often not meaningful.

    Notes
    -----

    .. warning::

        If you use `join='outer'` this fills 0s for sparse data when
        variables are absent in a batch. Use this with care. Dense data is
        filled with `NaN`.
    """
    # Argument normalization
    merge = resolve_merge_strategy(merge)
    uns_merge = resolve_merge_strategy(uns_merge)

    if isinstance(in_files, Mapping):
        if keys is not None:
            raise TypeError(
                "Cannot specify categories in both mapping keys and using `keys`. "
                "Only specify this once."
            )
        keys, in_files = list(in_files.keys()), list(in_files.values())
    else:
        in_files = list(in_files)

    if keys is None:
        keys = np.arange(len(in_files)).astype(str)

    _, dim = _resolve_dim(axis=axis)
    alt_axis, alt_dim = _resolve_dim(axis=1 - axis)

    # TODO: Generalize this to more than zarr
    groups = [zarr.open(store=p) for p in in_files]
    output_group = zarr.open(store=out_file)
    assert len(groups) > 1

    # TODO: This is just a temporary assertion
    # All dim_names must be equal
    if not _index_equal(groups, alt_dim):
        raise ValueError(f"{alt_dim}_names must be equal")

    # All groups must be anndata
    if not _attrs_equal(groups, path="", attrs_map={"encoding-type": "anndata"}):
        raise ValueError("All groups must be anndata")

    # Write metadata
    output_group.attrs.update(
        {"encoding-type": "anndata", "encoding-version": "0.1.0"})

    # Label column
    label_col = pd.Categorical.from_codes(
        np.repeat(np.arange(len(groups)), [
                  g["X"].attrs["shape"][axis] for g in groups]),
        categories=keys,
    )

    # Combining indexes
    concat_indices = pd.concat(
        [pd.Series(_df_index(g[dim])) for g in groups],
        ignore_index=True
    )
    if index_unique is not None:
        concat_indices = concat_indices.str.cat(
            label_col.map(str), sep=index_unique)
    concat_indices = pd.Index(concat_indices)

    alt_indices = merge_indices([_df_index(g[alt_dim])
                                for g in groups], join=join)

    Xs = [g["X"] for g in groups]
    write_concat_sequence_aligned(
        groups=Xs, output_group=output_group, out_path="X", axis=axis)

    # Annotation for concatenation axis
    concat_annot = pd.concat(
        unify_categorical_dtypes([read_elem(g[dim]) for g in groups]),
        join=join,
        ignore_index=True,
    )
    concat_annot.index = concat_indices
    if label is not None:
        concat_annot[label] = label_col
    write_elem(output_group, dim, concat_annot)

    # Annotation for other axis
    alt_annot = merge_dataframes(
        [read_elem(g[alt_dim]) for g in groups], alt_indices, merge
    )
    write_elem(output_group, alt_dim, alt_annot)

    mapping_names = [
        ("layers", None, axis),
        (f"{dim}m", concat_indices, 0),
    ]
    for m, m_index, m_axis in mapping_names:
        maps = [read_only_dict(g[m]) for g in groups]
        write_concat_mappings_aligned(
            maps, output_group, intersect_keys(maps), m, axis=m_axis, index=m_index)

    alt_mapping = merge(
        [read_group(g[alt_dim]) for g in groups]
    )
    write_elem(output_group, alt_dim, alt_mapping)
    if not alt_mapping:
        alt_df = pd.DataFrame(index=alt_indices)
        write_elem(output_group, alt_dim, alt_df)