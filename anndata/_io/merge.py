from .specs import read_elem, write_elem
import zarr
from .._core.merge import _resolve_dim
from .._core.sparse_dataset import SparseDataset
from .specs.registry import read_groups

from collections import OrderedDict
from collections.abc import Mapping, MutableSet
from functools import reduce, singledispatch

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
    MutableMapping
)
import typing

import numpy as np
from scipy.sparse import spmatrix

from ..compat import AwkArray, DaskArray


def _df_index(df: zarr.Group) -> np.ndarray:
    index_key = df.attrs["_index"]
    return read_elem(df[index_key])


def _has_same_attrs(groups: list[zarr.Group], path: str, attrs: Set[str]) -> bool:

    if len(groups) == 0:
        return True

    init_g = groups[0][path]

    return all(all(g[path].attrs[attr] == init_g.attrs[attr] for attr in attrs) for g in groups)


def _attrs_equal(groups: list[zarr.Group], path: str, attrs_map: Mapping[str, str]) -> bool:

    if len(groups) == 0:
        raise ValueError("List should not be empty")

    init_g = groups[0][path]

    return all(init_g.attrs[attr] == val for attr, val in attrs_map.items()) and \
        _has_same_attrs(groups, path, attrs_map.keys())


# def _attrs_is_in(groups: list[zarr.Group], path: str, attrs_map: Mapping[str, Set[str]]) -> bool:

#     if len(groups) == 0:
#         raise ValueError("List should not be empty")

#     return all(g[path].attrs[attr] in val for attr, val in attrs_map.items()
#                for g in groups)


# def _encoding_type_in(groups: list[zarr.Group], attr: str, ) -> bool:
#     init_encoding_type = groups[0][attr].attrs["encoding-type"]
#     init_encoding_version = groups[0][attr].attrs["encoding-version"]

#     return all(g[attr].attrs["encoding-type"] == init_encoding_type and
#                g[attr].attrs["encoding-version"] == init_encoding_version
#                for g in groups)


# def _assert_same_encoding_type_n_version(groups: list[zarr.Group], attr: str):
#     init_encoding_type = groups[0][attr].attrs["encoding-type"]
#     init_encoding_version = groups[0][attr].attrs["encoding-version"]
#     for g in groups:
#         assert g[attr].attrs["encoding-type"] == init_encoding_type
#         assert g[attr].attrs["encoding-version"] == init_encoding_version


def _index_equal(groups: list[zarr.Group], path) -> bool:
    names = _df_index(groups[0][path])
    for g in groups[1:]:
        curr_names = _df_index(g[path])
        if not np.array_equal(names, curr_names):
            return False
    return True


@singledispatch
def append_items(elem1, elem2, axis=0):
    raise NotImplementedError(f"Not implemented for type {type(elems),elems}")


@append_items.register
def _(elem1: SparseDataset, elem2: SparseDataset,  axis=0):
    supported_fmt = ["csr", "csc"][axis]

    if elem1.format_str != supported_fmt:
        raise ValueError(
            f"{elem1.format_str} not supported for axis={axis} concatenation.")
    # write_elem(output_group, path, elems[0])
    # written_elem = SparseDataset(output_group[path])
    # for e in elems[1:]:
        # written_elem.append(e)
    elem1.append(elem2)


@append_items.register
def _(elem1: zarr.Array, elem2: zarr.Array, axis=0):
    # write_elem(output_group, path, elems[0])
    # written_elem: zarr.Array = output_group[path]
    # for e in elems[1:]:
    elem1.append(elem2, axis=axis)


def _write_concat(elems, output_group: zarr.Group, path, axis=0):
    write_elem(output_group, path, elems[0])
    written_elem = SparseDataset(output_group[path])
    for g in elems[1:]:
        append_items(written_elem, g, axis=axis)


def _write_dim(groups, dim, output_group, same_names=False):

    dim_names_key = groups[0][dim].attrs["_index"]
    dim_group = output_group.create_group(dim)
    dim_group.attrs.update({
        "_index": dim_names_key,
        "column-order": [],
        "encoding-type": "dataframe",
        "encoding-version": "0.2.0",
    })
    dim_names = None
    if same_names:
        dim_names = _df_index(groups[0][dim])
    else:
        dim_names = np.concatenate([_df_index(g[dim]) for g in groups])

    write_elem(dim_group, dim_names_key, dim_names)



def concat_on_disk_zarr(
    groups: list[zarr.Group],
    output_group: zarr.Group,
    axis=0
):

    assert len(groups) > 1

    # var or obs
    _, dim = _resolve_dim(axis=axis)
    alt_axis, alt_dim = _resolve_dim(axis=1 - axis)

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

    # Write dim names
    _write_dim(groups, dim, output_group, same_names=False)
    _write_dim(groups, alt_dim, output_group, same_names=True)

    
    # dim_names_key = groups[0][dim].attrs["_index"]
    # dim_group = output_group.create_group(dim)
    # dim_names = _df_index(groups[0][dim])
    # dim_group.attrs.update({
    #     "_index": dim_names_key,
    #     "column-order": [],
    #     "encoding-type": "dataframe",
    #     "encoding-version": "0.2.0",
    # })
    # write_elem(dim_group, dim_names_key, dim_names)

    # # Write alt-dim names
    # alt_names_key = groups[0][alt_dim].attrs["_index"]
    # alt_dim_group = output_group.create_group(alt_dim)
    # alt_dim_group.attrs.update({
    #     "_index": alt_names_key,
    #     "column-order": [],
    #     "encoding-type": "dataframe",
    #     "encoding-version": "0.2.0",
    # })
    # write_elem(alt_dim_group, alt_names_key, np.concatenate(
    #     [_df_index(g[alt_dim]) for g in groups]))

    # Get on disk representation of X
    Xs = [read_groups(g["X"]) for g in groups]
    _write_concat(Xs, output_group, "X", axis=axis)



# def concat_on_disk(
#     in_files: Union[Collection[str], typing.MutableMapping],
#     out_file: Union[str, typing.MutableMapping],
#     overwrite: bool = False,
#     *,
#     axis: Literal[0, 1] = 0,
#     join: Literal["inner", "outer"] = "inner",
#     merge: Union[StrategiesLiteral, Callable, None] = None,
#     uns_merge: Union[StrategiesLiteral, Callable, None] = None,
#     label: Optional[str] = None,
#     keys: Optional[Collection] = None,
#     index_unique: Optional[str] = None,
#     fill_value: Optional[Any] = None,
#     pairwise: bool = False,
# ):
#     """Concatenates multiple AnnData objects along a specified axis using their
#     corresponding stores or paths, and writes the resulting AnnData object
#     to a target location on disk.

#     Unlike the `concat` function, this method does not require
#     loading the input AnnData objects into memory,
#     making it a memory-efficient alternative for large datasets.
#     The resulting object written to disk should be equivalent
#     to the concatenation of the loaded AnnData objects using
#     the `concat` function.



#     Params
#     ------
#     in_files
#         The corresponding stores or paths of AnnData objects to
#         be concatenated. If a Mapping is passed, keys are used for the `keys`
#         argument and values are concatenated.
#     out_file
#         The target path or store to write the result in.
#     overwrite
#         If False the and if a file already exists it will raise an error,
#         otherwise it will overwrite.
#     format
#         TODO
#     axis
#         Which axis to concatenate along.
#     join
#         How to align values when concatenating. If "outer", the union of the other axis
#         is taken. If "inner", the intersection. See :doc:`concatenation <../concatenation>`
#         for more.
#     merge
#         How elements not aligned to the axis being concatenated along are selected.
#         Currently implemented strategies include:

#         * `None`: No elements are kept.
#         * `"same"`: Elements that are the same in each of the objects.
#         * `"unique"`: Elements for which there is only one possible value.
#         * `"first"`: The first element seen at each from each position.
#         * `"only"`: Elements that show up in only one of the objects.
#     uns_merge
#         How the elements of `.uns` are selected. Uses the same set of strategies as
#         the `merge` argument, except applied recursively.
#     label
#         Column in axis annotation (i.e. `.obs` or `.var`) to place batch information in.
#         If it's None, no column is added.
#     keys
#         Names for each object being added. These values are used for column values for
#         `label` or appended to the index if `index_unique` is not `None`. Defaults to
#         incrementing integer labels.
#     index_unique
#         Whether to make the index unique by using the keys. If provided, this
#         is the delimeter between "{orig_idx}{index_unique}{key}". When `None`,
#         the original indices are kept.
#     fill_value
#         When `join="outer"`, this is the value that will be used to fill the introduced
#         indices. By default, sparse arrays are padded with zeros, while dense arrays and
#         DataFrames are padded with missing values.
#     pairwise
#         Whether pairwise elements along the concatenated dimension should be included.
#         This is False by default, since the resulting arrays are often not meaningful.

#     Notes
#     -----

#     .. warning::

#         If you use `join='outer'` this fills 0s for sparse data when
#         variables are absent in a batch. Use this with care. Dense data is
#         filled with `NaN`.
#     """
#     # Argument normalization
#     # merge = resolve_merge_strategy(merge)
#     # uns_merge = resolve_merge_strategy(uns_merge)

#     # if isinstance(in_files, Mapping):
#     #     if keys is not None:
#     #         raise TypeError(
#     #             "Cannot specify categories in both mapping keys and using `keys`. "
#     #             "Only specify this once."
#     #         )
#     #     keys, in_files = list(in_files.keys()), list(in_files.values())
#     # else:
#     #     in_files = list(in_files)

#     # if keys is None:
#     #     keys = np.arange(len(in_files)).astype(str)

#     # axis, dim = _resolve_dim(axis=axis)
#     # alt_axis, alt_dim = _resolve_dim(axis=1 - axis)

#     # # Label column
#     # label_col = pd.Categorical.from_codes(
#     #     np.repeat(np.arange(len(paths)), [a.shape[axis] for a in paths]),
#     #     categories=keys,
#     # )

#     # # Combining indexes
#     # concat_indices = pd.concat(
#     #     [pd.Series(dim_indices(a, axis=axis)) for a in paths], ignore_index=True
#     # )
#     # if index_unique is not None:
#     #     concat_indices = concat_indices.str.cat(label_col.map(str), sep=index_unique)
#     # concat_indices = pd.Index(concat_indices)

#     # alt_indices = merge_indices(
#     #     [dim_indices(a, axis=alt_axis) for a in paths], join=join
#     # )
#     # reindexers = [
#     #     gen_reindexer(alt_indices, dim_indices(a, axis=alt_axis)) for a in adatas
#     # ]

#     # # Annotation for concatenation axis
#     # concat_annot = pd.concat(
#     #     unify_categorical_dtypes([getattr(a, dim) for a in adatas]),
#     #     join=join,
#     #     ignore_index=True,
#     # )
#     # concat_annot.index = concat_indices
#     # if label is not None:
#     #     concat_annot[label] = label_col

#     # # Annotation for other axis
#     # alt_annot = merge_dataframes(
#     #     [getattr(a, alt_dim) for a in adatas], alt_indices, merge
#     # )

#     # X = concat_Xs(adatas, reindexers, axis=axis, fill_value=fill_value)

#     # if join == "inner":
#     #     layers = inner_concat_aligned_mapping(
#     #         [a.layers for a in adatas], axis=axis, reindexers=reindexers
#     #     )
#     #     concat_mapping = inner_concat_aligned_mapping(
#     #         [getattr(a, f"{dim}m") for a in adatas], index=concat_indices
#     #     )
#     #     if pairwise:
#     #         concat_pairwise = concat_pairwise_mapping(
#     #             mappings=[getattr(a, f"{dim}p") for a in adatas],
#     #             shapes=[a.shape[axis] for a in adatas],
#     #             join_keys=intersect_keys,
#     #         )
#     #     else:
#     #         concat_pairwise = {}
#     # elif join == "outer":
#     #     layers = outer_concat_aligned_mapping(
#     #         [a.layers for a in adatas], reindexers, axis=axis, fill_value=fill_value
#     #     )
#     #     concat_mapping = outer_concat_aligned_mapping(
#     #         [getattr(a, f"{dim}m") for a in adatas],
#     #         index=concat_indices,
#     #         fill_value=fill_value,
#     #     )
#     #     if pairwise:
#     #         concat_pairwise = concat_pairwise_mapping(
#     #             mappings=[getattr(a, f"{dim}p") for a in adatas],
#     #             shapes=[a.shape[axis] for a in adatas],
#     #             join_keys=union_keys,
#     #         )
#     #     else:
#     #         concat_pairwise = {}

#     # # TODO: Reindex lazily, so we don't have to make those copies until we're sure we need the element
#     # alt_mapping = merge(
#     #     [
#     #         {k: r(v, axis=0) for k, v in getattr(a, f"{alt_dim}m").items()}
#     #         for r, a in zip(reindexers, adatas)
#     #     ],
#     # )
#     # alt_pairwise = merge(
#     #     [
#     #         {k: r(r(v, axis=0), axis=1) for k, v in getattr(a, f"{alt_dim}p").items()}
#     #         for r, a in zip(reindexers, adatas)
#     #     ]
#     # )
#     # uns = uns_merge([a.uns for a in adatas])

#     # raw = None
#     # has_raw = [a.raw is not None for a in adatas]
#     # if all(has_raw):
#     #     raw = concat(
#     #         [
#     #             AnnData(
#     #                 X=a.raw.X,
#     #                 obs=pd.DataFrame(index=a.obs_names),
#     #                 var=a.raw.var,
#     #                 varm=a.raw.varm,
#     #             )
#     #             for a in adatas
#     #         ],
#     #         join=join,
#     #         label=label,
#     #         keys=keys,
#     #         index_unique=index_unique,
#     #         fill_value=fill_value,
#     #         axis=axis,
#     #     )
#     # elif any(has_raw):
#     #     warn(
#     #         "Only some AnnData objects have `.raw` attribute, "
#     #         "not concatenating `.raw` attributes.",
#     #         UserWarning,
#     #     )
#     # return AnnData(
#     #     **{
#     #         "X": X,
#     #         "layers": layers,
#     #         dim: concat_annot,
#     #         alt_dim: alt_annot,
#     #         f"{dim}m": concat_mapping,
#     #         f"{alt_dim}m": alt_mapping,
#     #         f"{dim}p": concat_pairwise,
#     #         f"{alt_dim}p": alt_pairwise,
#     #         "uns": uns,
#     #         "raw": raw,
#     #     }
#     # )
