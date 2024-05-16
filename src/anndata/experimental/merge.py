from __future__ import annotations

import os
import shutil
from collections.abc import Collection, Iterable, Mapping, Sequence
from functools import singledispatch
from pathlib import Path
from typing import Any, Callable, Literal

import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix, csr_matrix

from .._core.file_backing import to_memory
from .._core.merge import (
    MissingVal,
    Reindexer,
    StrategiesLiteral,
    _resolve_dim,
    concat_arrays,
    gen_inner_reindexers,
    gen_reindexer,
    intersect_keys,
    merge_dataframes,
    merge_indices,
    resolve_merge_strategy,
    unify_dtypes,
)
from .._core.sparse_dataset import BaseCompressedSparseDataset, sparse_dataset
from .._io.specs import read_elem, write_elem
from ..compat import H5Array, H5Group, ZarrArray, ZarrGroup, _map_cat_to_str
from . import read_dispatched

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


# Checks if given indices are equal to each other in the whole list.
def _indices_equal(indices: Iterable[pd.Index]) -> bool:
    init_elem = indices[0]
    return all(np.array_equal(init_elem, elem) for elem in indices[1:])


def _gen_slice_to_append(
    datasets: Sequence[BaseCompressedSparseDataset],
    reindexers,
    max_loaded_elems: int,
    axis=0,
    fill_value=None,
):
    for ds, ri in zip(datasets, reindexers):
        n_slices = ds.shape[axis] * ds.shape[1 - axis] // max_loaded_elems
        if n_slices < 2:
            yield (csr_matrix, csc_matrix)[axis](
                ri(to_memory(ds), axis=1 - axis, fill_value=fill_value)
            )
        else:
            slice_size = max_loaded_elems // ds.shape[1 - axis]
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

                yield (csr_matrix, csc_matrix)[axis](
                    ri(ds_part, axis=1 - axis, fill_value=fill_value)
                )
                rem_slices -= slice_size
                idx += slice_size


###################
# File Management
###################


@singledispatch
def as_group(store, *args, **kwargs) -> ZarrGroup | H5Group:
    raise NotImplementedError("This is not yet implemented.")


@as_group.register(os.PathLike)
@as_group.register(str)
def _(store: os.PathLike | str, *args, **kwargs) -> ZarrGroup | H5Group:
    store = Path(store)
    if store.suffix == ".h5ad":
        import h5py

        return h5py.File(store, *args, **kwargs)
    import zarr

    return zarr.open_group(store, *args, **kwargs)


@as_group.register(ZarrGroup)
@as_group.register(H5Group)
def _(store, *args, **kwargs):
    return store


###################
# Reading
###################


def read_as_backed(group: ZarrGroup | H5Group):
    """
    Read the group until
    BaseCompressedSparseDataset, Array or EAGER_TYPES are encountered.
    """

    def callback(func, elem_name: str, elem, iospec):
        if iospec.encoding_type in SPARSE_MATRIX:
            return sparse_dataset(elem)
        elif iospec.encoding_type in EAGER_TYPES:
            return read_elem(elem)
        elif iospec.encoding_type == "array":
            return elem
        elif iospec.encoding_type == "dict":
            return {k: read_as_backed(v) for k, v in elem.items()}
        else:
            return func(elem)

    return read_dispatched(group, callback=callback)


def _df_index(df: ZarrGroup | H5Group) -> pd.Index:
    index_key = df.attrs["_index"]
    return pd.Index(read_elem(df[index_key]))


###################
# Writing
###################


def write_concat_dense(
    arrays: Sequence[ZarrArray | H5Array],
    output_group: ZarrGroup | H5Group,
    output_path: ZarrGroup | H5Group,
    axis: Literal[0, 1] = 0,
    reindexers: Reindexer = None,
    fill_value=None,
):
    """
    Writes the concatenation of given dense arrays to disk using dask.
    """
    import dask.array as da

    darrays = (
        da.from_array(a, chunks="auto" if a.chunks is None else a.chunks)
        for a in arrays
    )

    res = da.concatenate(
        [
            ri(a, axis=1 - axis, fill_value=fill_value)
            for a, ri in zip(darrays, reindexers)
        ],
        axis=axis,
    )
    write_elem(output_group, output_path, res)
    output_group[output_path].attrs.update(
        {"encoding-type": "array", "encoding-version": "0.2.0"}
    )


def write_concat_sparse(
    datasets: Sequence[BaseCompressedSparseDataset],
    output_group: ZarrGroup | H5Group,
    output_path: ZarrGroup | H5Group,
    max_loaded_elems: int,
    axis: Literal[0, 1] = 0,
    reindexers: Reindexer = None,
    fill_value=None,
):
    """
    Writes and concatenates sparse datasets into a single output dataset.

    Args:
        datasets (Sequence[BaseCompressedSparseDataset]): A sequence of BaseCompressedSparseDataset objects to be concatenated.
        output_group (Union[ZarrGroup, H5Group]): The output group where the concatenated dataset will be written.
        output_path (Union[ZarrGroup, H5Group]): The output path where the concatenated dataset will be written.
        max_loaded_elems (int): The maximum number of sparse elements to load at once.
        axis (Literal[0, 1], optional): The axis along which the datasets should be concatenated.
            Defaults to 0.
        reindexers (Reindexer, optional): A reindexer object that defines the reindexing operation to be applied.
            Defaults to None.
        fill_value (Any, optional): The fill value to use for missing elements. Defaults to None.
    """
    elems = None
    if all(ri.no_change for ri in reindexers):
        elems = iter(datasets)
    else:
        elems = _gen_slice_to_append(
            datasets, reindexers, max_loaded_elems, axis, fill_value
        )
    number_non_zero = sum(len(d.group["indices"]) for d in datasets)
    init_elem = next(elems)
    indptr_dtype = "int64" if number_non_zero >= np.iinfo(np.int32).max else "int32"
    write_elem(
        output_group,
        output_path,
        init_elem,
        dataset_kwargs=dict(indptr_dtype=indptr_dtype),
    )
    del init_elem
    out_dataset: BaseCompressedSparseDataset = read_as_backed(output_group[output_path])
    for temp_elem in elems:
        out_dataset.append(temp_elem)
        del temp_elem


def _write_concat_mappings(
    mappings,
    output_group: ZarrGroup | H5Group,
    keys,
    path,
    max_loaded_elems,
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
        _write_concat_sequence(
            elems,
            output_group=mapping_group,
            output_path=k,
            axis=axis,
            index=index,
            reindexers=reindexers,
            fill_value=fill_value,
            max_loaded_elems=max_loaded_elems,
        )


def _write_concat_arrays(
    arrays: Sequence[ZarrArray | H5Array | BaseCompressedSparseDataset],
    output_group,
    output_path,
    max_loaded_elems,
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
            reindexers = gen_inner_reindexers(arrays, new_index=None, axis=axis)
        else:
            raise NotImplementedError("Cannot reindex arrays with outer join.")

    if isinstance(init_elem, BaseCompressedSparseDataset):
        expected_sparse_fmt = ["csr", "csc"][axis]
        if all(a.format == expected_sparse_fmt for a in arrays):
            write_concat_sparse(
                arrays,
                output_group,
                output_path,
                max_loaded_elems,
                axis,
                reindexers,
                fill_value,
            )
        else:
            raise NotImplementedError(
                f"Concat of following not supported: {[a.format for a in arrays]}"
            )
    else:
        write_concat_dense(
            arrays, output_group, output_path, axis, reindexers, fill_value
        )


def _write_concat_sequence(
    arrays: Sequence[pd.DataFrame | BaseCompressedSparseDataset | H5Array | ZarrArray],
    output_group,
    output_path,
    max_loaded_elems,
    axis=0,
    index=None,
    reindexers=None,
    fill_value=None,
    join="inner",
):
    """
    array, dataframe, csc_matrix, csc_matrix
    """
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
    elif all(
        isinstance(a, (pd.DataFrame, BaseCompressedSparseDataset, H5Array, ZarrArray))
        for a in arrays
    ):
        _write_concat_arrays(
            arrays,
            output_group,
            output_path,
            max_loaded_elems,
            axis,
            reindexers,
            fill_value,
            join,
        )
    else:
        raise NotImplementedError(
            f"Concatenation of these types is not yet implemented: {[type(a) for a in arrays] } with axis={axis}."
        )


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
        unify_dtypes(read_elem(g[dim]) for g in groups),
        join=join,
        ignore_index=True,
    )
    concat_annot.index = concat_indices
    if label is not None:
        concat_annot[label] = label_col
    write_elem(output_group, dim, concat_annot)


def concat_on_disk(
    in_files: Collection[str | os.PathLike] | Mapping[str, str | os.PathLike],
    out_file: str | os.PathLike,
    *,
    max_loaded_elems: int = 100_000_000,
    axis: Literal[0, 1] = 0,
    join: Literal["inner", "outer"] = "inner",
    merge: StrategiesLiteral | Callable[[Collection[Mapping]], Mapping] | None = None,
    uns_merge: (
        StrategiesLiteral | Callable[[Collection[Mapping]], Mapping] | None
    ) = None,
    label: str | None = None,
    keys: Collection[str] | None = None,
    index_unique: str | None = None,
    fill_value: Any | None = None,
    pairwise: bool = False,
) -> None:
    """\
    Concatenates multiple AnnData objects along a specified axis using their
    corresponding stores or paths, and writes the resulting AnnData object
    to a target location on disk.

    Unlike :func:`anndata.concat`, this method does not require
    loading the input AnnData objects into memory,
    making it a memory-efficient alternative for large datasets.
    The resulting object written to disk should be equivalent
    to the concatenation of the loaded AnnData objects using
    :func:`anndata.concat`.

    To adjust the maximum amount of data loaded in memory; for sparse
    arrays use the max_loaded_elems argument; for dense arrays
    see the Dask documentation, as the Dask concatenation function is used
    to concatenate dense arrays in this function

    Params
    ------
    in_files
        The corresponding stores or paths of AnnData objects to
        be concatenated. If a Mapping is passed, keys are used for the `keys`
        argument and values are concatenated.
    out_file
        The target path or store to write the result in.
    max_loaded_elems
        The maximum number of elements to load in memory when concatenating
        sparse arrays. Note that this number also includes the empty entries.
        Set to 100m by default meaning roughly 400mb will be loaded
        to memory simultaneously.
    axis
        Which axis to concatenate along.
    join
        How to align values when concatenating. If `"outer"`, the union of the other axis
        is taken. If `"inner"`, the intersection. See :doc:`concatenation <../concatenation>`
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
        is the delimiter between `"{orig_idx}{index_unique}{key}"`. When `None`,
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

    Examples
    --------

    See :func:`anndata.concat` for the semantics.
    The following examples highlight the differences this function has.

    First, let’s get some “big” datasets with a compatible ``var`` axis:

    >>> import httpx
    >>> import scanpy as sc
    >>> base_url = "https://datasets.cellxgene.cziscience.com"
    >>> def get_cellxgene_data(id_: str):
    ...     out_path = sc.settings.datasetdir / f'{id_}.h5ad'
    ...     if out_path.exists():
    ...         return out_path
    ...     file_url = f"{base_url}/{id_}.h5ad"
    ...     sc.settings.datasetdir.mkdir(parents=True, exist_ok=True)
    ...     with httpx.stream('GET', file_url) as r, out_path.open('wb') as f:
    ...         r.raise_for_status()
    ...         for data in r.iter_bytes():
    ...             f.write(data)
    ...     return out_path
    >>> path_b_cells = get_cellxgene_data('a93eab58-3d82-4b61-8a2f-d7666dcdb7c4')
    >>> path_fetal = get_cellxgene_data('d170ff04-6da0-4156-a719-f8e1bbefbf53')

    Now we can concatenate them on-disk:

    >>> import anndata as ad
    >>> ad.experimental.concat_on_disk(
    ...     dict(b_cells=path_b_cells, fetal=path_fetal),
    ...     'merged.h5ad',
    ...     label='dataset',
    ... )
    >>> adata = ad.read_h5ad('merged.h5ad', backed=True)
    >>> adata.X
    CSRDataset: backend hdf5, shape (490, 15585), data_dtype float32
    >>> adata.obs['dataset'].value_counts()  # doctest: +SKIP
    dataset
    fetal      344
    b_cells    146
    Name: count, dtype: int64
    """
    if len(in_files) == 0:
        raise ValueError("No objects to concatenate.")

    # Argument normalization
    if pairwise:
        raise NotImplementedError("pairwise concatenation not yet implemented")
    if join != "inner":
        raise NotImplementedError("only inner join is currently supported")

    merge = resolve_merge_strategy(merge)
    uns_merge = resolve_merge_strategy(uns_merge)

    out_file = Path(out_file)
    if not out_file.parent.exists():
        raise FileNotFoundError(f"Parent directory of {out_file} does not exist.")

    if isinstance(in_files, Mapping):
        if keys is not None:
            raise TypeError(
                "Cannot specify categories in both mapping keys and using `keys`. "
                "Only specify this once."
            )
        keys, in_files = list(in_files.keys()), list(in_files.values())
    else:
        in_files = list(in_files)

    if len(in_files) == 1:
        shutil.copy2(in_files[0], out_file)
        return

    if keys is None:
        keys = np.arange(len(in_files)).astype(str)

    _, dim = _resolve_dim(axis=axis)
    _, alt_dim = _resolve_dim(axis=1 - axis)

    output_group = as_group(out_file, mode="w")
    groups = [as_group(f) for f in in_files]

    use_reindexing = False

    alt_dims = [_df_index(g[alt_dim]) for g in groups]
    # All dim_names must be equal if reindexing not applied
    if not _indices_equal(alt_dims):
        use_reindexing = True

    # All groups must be anndata
    if not all(g.attrs.get("encoding-type") == "anndata" for g in groups):
        raise ValueError("All groups must be anndata")

    # Write metadata
    output_group.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})

    # Read the backed objects of Xs
    Xs = [read_as_backed(g["X"]) for g in groups]

    # Label column
    label_col = pd.Categorical.from_codes(
        np.repeat(np.arange(len(groups)), [x.shape[axis] for x in Xs]),
        categories=keys,
    )

    # Combining indexes
    concat_indices = pd.concat(
        [pd.Series(_df_index(g[dim])) for g in groups], ignore_index=True
    )
    if index_unique is not None:
        concat_indices = concat_indices.str.cat(
            _map_cat_to_str(label_col), sep=index_unique
        )

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

    _write_concat_arrays(
        arrays=Xs,
        output_group=output_group,
        output_path="X",
        axis=axis,
        reindexers=reindexers,
        fill_value=fill_value,
        max_loaded_elems=max_loaded_elems,
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
        _write_concat_mappings(
            maps,
            output_group,
            intersect_keys(maps),
            m,
            max_loaded_elems=max_loaded_elems,
            axis=m_axis,
            index=m_index,
            reindexers=m_reindexers,
            fill_value=fill_value,
        )
