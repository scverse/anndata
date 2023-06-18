from functools import singledispatch
from .._io.specs import read_elem, write_elem
from .._core.merge import (
    _resolve_dim,
    StrategiesLiteral,
    intersect_keys,
    unify_dtypes,
    merge_indices,
    merge_dataframes,
    resolve_merge_strategy,
    gen_reindexer,
    gen_inner_reindexers,
    reduce,
    concat_arrays,
    MissingVal,
    Reindexer,
)
from .._core.sparse_dataset import SparseDataset
from .._core.file_backing import to_memory
from . import read_dispatched
import pandas as pd
from typing import (
    Any,
    Callable,
    Collection,
    Optional,
    Sequence,
    Union,
    Literal,
    MutableMapping,
    Mapping,
)
import typing
from pathlib import Path
import numpy as np

from ..compat import ZarrGroup, ZarrArray, H5Group, H5Array

from scipy.sparse import csr_matrix, csc_matrix

MAX_LOAD_SIZE = 60_000_000

SPARSE_MATRIX = {"csc_matrix", "csr_matrix"}

EAGER_TYPES = {"dataframe", "awkward-array"}

###################
# Utilities
###################

# Wrapper to reindexer that stores if there is a change
# and won't do anything if there is
class IdentityReindexer:
    def __init__(self):
        self.no_change = True

    def __call__(self, x, *args, **kwargs):
        return x


def _df_index(df: Union[ZarrGroup, H5Group]) -> pd.Index:
    index_key = df.attrs["_index"]
    return pd.Index(read_elem(df[index_key]))


def _requires_reindexing(indices) -> bool:
    init_elem = indices[0]
    return any(not np.array_equal(init_elem, elem) for elem in indices[1:])


def write_concat_mappings(
    mappings,
    output_group: Union[ZarrGroup, H5Group],
    keys,
    path,
    axis=0,
    index=None,
    reindexers=None,
    fill_value=None,
):
    """
    Write a list of mappings to a zarr/h5 group.
    """
    mapping_group = output_group.create_group(path)
    mapping_group.attrs.update(
        {
            "encoding-type": "dict",
            "encoding-version": "0.1.0",
        }
    )
    for k in keys:
        elems = [m[k] for m in mappings]
        write_concat_sequence(
            elems,
            output_group=mapping_group,
            output_path=k,
            axis=axis,
            index=index,
            reindexers=reindexers,
            fill_value=fill_value,
        )


def read_as_backed(group: Union[ZarrGroup, H5Group]):
    """Read the group until
    SparseDataset, Array or EAGER_TYPES are encountered.
    """

    def callback(func, elem_name: str, elem, iospec):
        if iospec.encoding_type in SPARSE_MATRIX:
            return SparseDataset(elem)
        elif iospec.encoding_type in EAGER_TYPES:
            return read_elem(elem)
        elif iospec.encoding_type == "array":
            return elem
        elif iospec.encoding_type == "dict":
            return dict(elem.items())
        else:
            return func(elem)

    return read_dispatched(group, callback=callback)


def get_encoding_type(group: Union[ZarrGroup, H5Group]) -> str:
    """Get the encoding type of a group."""
    return group.attrs.get("encoding-type")


def get_shape(group: Union[ZarrGroup, H5Group]) -> str:
    """Get the shape of a group (array, sparse matrix, etc.)"""
    if group.attrs["encoding-type"] == "array":
        return group.shape
    if group.attrs["encoding-type"] in SPARSE_MATRIX:
        return group.attrs.get("shape")
    raise NotImplementedError(f"Cannot get shape of {group.attrs['encoding-type']}")


def _gen_elem_to_append(elems, axis=0, reindexers=None, fill_value=None):
    for elem, ri in zip(elems, reindexers):
        if ri.no_change:
            yield elem
        else:
            yield ri(elem, axis=1 - axis, fill_value=fill_value)



def _gen_slice_to_append(
    datasets: Sequence[SparseDataset],
    reindexers=None,
    axis=0,
    fill_value=None,
):
    for ds, ri in zip(datasets, reindexers):
        n_slices = ds.shape[axis] * ds.shape[1 - axis] // MAX_LOAD_SIZE
        if n_slices < 2:
            yield ri(to_memory(ds), axis=1 - axis, fill_value=fill_value)
        else:
            slice_size = MAX_LOAD_SIZE // ds.shape[1 - axis]
            if slice_size == 0:
                slice_size = 1
            rem_slices = ds.shape[axis]
            idx = 0
            while rem_slices > 0:
                ds_part = None
                if axis == 0:
                    ds_part = ds[idx : idx + slice_size, :]
                elif axis == 1:
                    ds_part = ds[:, idx : idx + slice_size]

                yield ri(ds_part, axis=1 - axis, fill_value=fill_value)
                rem_slices -= slice_size
                idx += slice_size


def _write_concat_sparse(
    datasets: Sequence[SparseDataset],
    output_group,
    output_path,
    axis=0,
    reindexers=None,
    fill_value=None,
):
    elems = _gen_slice_to_append(datasets, reindexers, axis, fill_value)
    init_elem = (csr_matrix, csc_matrix)[axis](next(elems))
    write_elem(output_group, output_path, init_elem)
    del init_elem
    out_dataset: SparseDataset = read_as_backed(output_group[output_path])
    for temp_elem in elems:
        out_dataset.append((csr_matrix, csc_matrix)[axis](temp_elem))
        del temp_elem


def write_concat_arrays(
    arrays: Sequence[Union[ZarrArray, H5Array, SparseDataset]],
    output_group,
    output_path,
    axis=0,
    reindexers=None,
    fill_value=None,
    join="inner",
):
    init_elem = arrays[0]
    init_type = type(init_elem)
    if not all(isinstance(a, init_type) for a in arrays):
        raise NotImplementedError(
            f"All elements must be the same type instead got types: {[type(a) for a in arrays]}"
        )
    if reindexers is None:
        if join == "inner":
            reindexers = gen_inner_reindexers(arrays, None, axis=axis)
        else:
            raise NotImplementedError("Cannot reindex arrays with outer join.")
    if isinstance(init_elem, SparseDataset):
        expected_sparse_fmt = ["csr", "csc"][axis]
        if all(a.format_str == expected_sparse_fmt for a in arrays):
            _write_concat_sparse(
                arrays, output_group, output_path, axis, reindexers, fill_value
            )
        else:
            raise NotImplementedError(
                f"Concat of following not supported: {[a.format_str for a in arrays]}"
            )
    elif isinstance(init_elem, H5Array):
        dim_shapes = [a.shape[axis] for a in arrays]
        dim_res_len = sum(dim_shapes)

        if all(r.no_change for r in reindexers):
            alt_dim_res_len = arrays[0].shape[1 - axis]
        else:
            alt_dim_res_len = len(reindexers[0].new_idx)
            assert all(len(ri.new_idx) == alt_dim_res_len for ri in reindexers)

        new_shape = (dim_res_len, alt_dim_res_len)[:: 1 - 2 * axis]
        out_array: H5Array = output_group.create_dataset(output_path, shape=new_shape)
        idx = 0
        for arr in _gen_slice_to_append(
            arrays, axis=axis, reindexers=reindexers, fill_value=fill_value
        ):
            added_size = arr.shape[axis]
            if axis == 0:
                out_array[idx : idx + added_size, :] = arr
            else:
                out_array[:, idx : idx + added_size] = arr
            idx += added_size
        output_group[output_path].attrs.update(
            {"encoding-type": "array", "encoding-version": "0.2.0"}
        )

    elif isinstance(init_elem, ZarrArray):
        import dask.array as da

        darrays = (da.from_zarr(a) for a in arrays)
        res = da.concatenate(
            _gen_elem_to_append(
                darrays,
                axis=axis,
                reindexers=reindexers,
                fill_value=fill_value,
            ),
            axis=axis,
        )
        write_elem(output_group, output_path, res)
        output_group[output_path].attrs.update(
            {"encoding-type": "array", "encoding-version": "0.2.0"}
        )


def write_concat_sequence(
    groups: Sequence[Union[ZarrGroup, H5Group]],
    output_group,
    output_path,
    axis=0,
    index=None,
    reindexers=None,
    fill_value=None,
    join="inner",
):
    """
    array, dataframe, csc_matrix, csc_matrix
    """
    arrays: Union[pd.DataFrame, SparseDataset, H5Array, ZarrArray] = [
        read_as_backed(g) for g in groups
    ]
    if any(isinstance(a, pd.DataFrame) for a in arrays):
        if reindexers is None:
            if join == "inner":
                reindexers = gen_inner_reindexers(arrays, None, axis=axis)
            else:
                raise NotImplementedError("Cannot reindex dataframes with outer join.")
        if not all(
            isinstance(a, pd.DataFrame) or a is MissingVal or 0 in a.shape
            for a in arrays
        ):
            raise NotImplementedError(
                "Cannot concatenate a dataframe with other array types."
            )
        df = concat_arrays(
            arrays=arrays,
            reindexers=reindexers,
            axis=axis,
            index=index,
            fill_value=fill_value,
        )
        write_elem(output_group, output_path, df)
    # If all are compatible sparse matrices
    elif all(
        isinstance(a, (pd.DataFrame, SparseDataset, H5Array, ZarrArray)) for a in arrays
    ):
        write_concat_arrays(
            arrays, output_group, output_path, axis, reindexers, fill_value, join
        )
    else:
        raise NotImplementedError(
            f"Concatenation of these types is not yet implemented: {[get_encoding_type(g) for g in groups] } with axis={axis}."
        )


def _get_group_from_hdf5(store, *args, **kwargs) -> H5Group:
    import h5py

    return h5py.File(store, *args, **kwargs)


@singledispatch
def _get_group(store, *args, **kwargs) -> Union[ZarrGroup, H5Group]:
    raise NotImplementedError("This is not yet implemented.")


@_get_group.register
def _(store: Path, *args, **kwargs) -> Union[ZarrGroup, H5Group]:
    if store.suffix == ".h5ad":
        return _get_group_from_hdf5(store, *args, **kwargs)
    import zarr

    return zarr.open_group(store, *args, **kwargs)


@_get_group.register
def _(store: str, *args, **kwargs) -> Union[ZarrGroup, H5Group]:
    if store.endswith(".h5ad"):
        return _get_group_from_hdf5(store, *args, **kwargs)
    import zarr

    return zarr.open_group(store, *args, **kwargs)


@_get_group.register
def _(store: dict, *args, **kwargs) -> ZarrGroup:
    import zarr

    return zarr.open_group(store, *args, **kwargs)


@_get_group.register
def _(store: ZarrGroup, *args, **kwargs) -> ZarrGroup:
    return store


@_get_group.register
def _(store: H5Group, *args, **kwargs) -> H5Group:
    return store


def _get_groups_from_paths(
    in_files: Union[
        Collection[str],
        Collection[Path],
        MutableMapping,
        Collection[Union[ZarrGroup, H5Group]],
    ],
    out_file: Union[str, MutableMapping, ZarrGroup, Path],
    overwrite: bool,
) -> typing.Tuple[Sequence[ZarrGroup], ZarrGroup]:
    """Returns the groups to be concatenated and the output group."""
    mode = "w" if overwrite else "w-"

    res_out_file = _get_group(out_file, mode=mode)
    res_in_files = [_get_group(f) for f in in_files]
    return res_in_files, res_out_file


def _write_alt_mapping(groups, output_group, alt_dim, alt_indices, merge):
    alt_mapping = merge([read_as_backed(g[alt_dim]) for g in groups])
    # If its empty, we need to write an empty dataframe with the correct index
    if not alt_mapping:
        alt_df = pd.DataFrame(index=alt_indices)
        write_elem(output_group, alt_dim, alt_df)
    else:
        write_elem(output_group, alt_dim, alt_mapping)


def _write_alt_annot(groups, output_group, alt_dim, alt_indices, merge):
    # Annotation for other axis
    alt_annot = merge_dataframes(
        [read_elem(g[alt_dim]) for g in groups], alt_indices, merge
    )
    write_elem(output_group, alt_dim, alt_annot)


def _write_dim_annot(groups, output_group, dim, concat_indices, label, label_col, join):
    concat_annot = pd.concat(
        unify_dtypes([read_elem(g[dim]) for g in groups]),
        join=join,
        ignore_index=True,
    )
    concat_annot.index = concat_indices
    if label is not None:
        concat_annot[label] = label_col
    write_elem(output_group, dim, concat_annot)


def concat_on_disk(
    in_files: Union[
        Collection[str],
        typing.MutableMapping,
        Collection[ZarrGroup],
        Collection[H5Group],
    ],
    out_file: Union[str, typing.MutableMapping, ZarrGroup, H5Group],
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
    if pairwise:
        raise NotImplementedError("pairwise concatenation not yet implemented")
    if join != "inner":
        raise NotImplementedError("only inner join is currently supported")

    merge = resolve_merge_strategy(merge)
    uns_merge = resolve_merge_strategy(uns_merge)
    if len(in_files) <= 1:
        raise ValueError("Must pass at least two files to concatenate.")
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
    _, alt_dim = _resolve_dim(axis=1 - axis)

    groups, output_group = _get_groups_from_paths(
        in_files, out_file, overwrite=overwrite
    )

    use_reindexing = False

    alt_dims = [_df_index(g[alt_dim]) for g in groups]
    # All dim_names must be equal if reindexing not applied
    if _requires_reindexing(alt_dims):
        use_reindexing = True

    # All groups must be anndata
    if not all(get_encoding_type(g) == "anndata" for g in groups):
        raise ValueError("All groups must be anndata")

    # Write metadata
    output_group.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})

    # Label column
    label_col = pd.Categorical.from_codes(
        np.repeat(np.arange(len(groups)), [get_shape(g["X"])[axis] for g in groups]),
        categories=keys,
    )

    # Combining indexes
    concat_indices = pd.concat(
        [pd.Series(_df_index(g[dim])) for g in groups], ignore_index=True
    )
    if index_unique is not None:
        concat_indices = concat_indices.str.cat(label_col.map(str), sep=index_unique)

    # Resulting indices for {dim} and {alt_dim}
    concat_indices = pd.Index(concat_indices)

    alt_indices = merge_indices(alt_dims, join=join)

    reindexers = None
    if use_reindexing:
        reindexers = [
            gen_reindexer(alt_indices, alt_old_index) for alt_old_index in alt_dims
        ]
    else:
        reindexers = [IdentityReindexer()] * len(groups)

    # Write {dim}
    _write_dim_annot(groups, output_group, dim, concat_indices, label, label_col, join)

    # Write {alt_dim}
    _write_alt_annot(groups, output_group, alt_dim, alt_indices, merge)

    # Write {alt_dim}m
    _write_alt_mapping(groups, output_group, alt_dim, alt_indices, merge)

    # Write X
    Xs = [g["X"] for g in groups]
    write_concat_sequence(
        groups=Xs,
        output_group=output_group,
        output_path="X",
        axis=axis,
        reindexers=reindexers,
        fill_value=fill_value,
    )

    # Write Layers and {dim}m
    mapping_names = [
        (
            f"{dim}m",
            concat_indices,
            0,
            None if use_reindexing else [IdentityReindexer()] * len(groups),
        ),
        ("layers", None, axis, reindexers),
    ]
    for m, m_index, m_axis, m_reindexers in mapping_names:
        maps = [read_as_backed(g[m]) for g in groups]
        write_concat_mappings(
            maps,
            output_group,
            intersect_keys(maps),
            m,
            axis=m_axis,
            index=m_index,
            reindexers=m_reindexers,
            fill_value=fill_value,
        )
