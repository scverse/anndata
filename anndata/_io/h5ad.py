import re
import collections.abc as cabc
from functools import _find_impl, partial
from warnings import warn
from pathlib import Path
from types import MappingProxyType
from typing import Callable, Type, TypeVar, Union
from typing import Collection, Sequence, Mapping

import h5py
import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype
from scipy import sparse

from .._core.sparse_dataset import SparseDataset
from .._core.file_backing import AnnDataFileManager
from .._core.anndata import AnnData
from .._core.raw import Raw
from ..compat import _from_fixed_length_strings, _clean_uns, Literal
from .utils import (
    report_read_key_on_error,
    report_write_key_on_error,
    idx_chunks_along_axis,
    write_attribute,
    read_attribute,
    _read_legacy_raw,
    EncodingVersions,
)


H5Group = Union[h5py.Group, h5py.File]
H5Dataset = Union[h5py.Dataset]
T = TypeVar("T")


def _to_hdf5_vlen_strings(value: np.ndarray) -> np.ndarray:
    """This corrects compound dtypes to work with hdf5 files."""
    new_dtype = []
    for dt_name, (dt_type, _) in value.dtype.fields.items():
        if dt_type.kind in ("U", "O"):
            new_dtype.append((dt_name, h5py.special_dtype(vlen=str)))
        else:
            new_dtype.append((dt_name, dt_type))
    return value.astype(new_dtype)


def write_h5ad(
    filepath: Union[Path, str],
    adata: AnnData,
    *,
    force_dense: bool = None,
    as_dense: Sequence[str] = (),
    dataset_kwargs: Mapping = MappingProxyType({}),
    **kwargs,
) -> None:
    if force_dense is not None:
        warn(
            "The `force_dense` argument is deprecated. Use `as_dense` instead.",
            FutureWarning,
        )
    if force_dense is True:
        if adata.raw is not None:
            as_dense = ("X", "raw/X")
        else:
            as_dense = ("X",)
    if isinstance(as_dense, str):
        as_dense = [as_dense]
    if "raw.X" in as_dense:
        as_dense = list(as_dense)
        as_dense[as_dense.index("raw.X")] = "raw/X"
    if any(val not in {"X", "raw/X"} for val in as_dense):
        raise NotImplementedError(
            "Currently, only `X` and `raw/X` are supported values in `as_dense`"
        )
    if "raw/X" in as_dense and adata.raw is None:
        raise ValueError("Cannot specify writing `raw/X` to dense if it doesn’t exist.")

    adata.strings_to_categoricals()
    if adata.raw is not None:
        adata.strings_to_categoricals(adata.raw.var)
    dataset_kwargs = {**dataset_kwargs, **kwargs}
    filepath = Path(filepath)
    mode = "a" if adata.isbacked else "w"
    if adata.isbacked:  # close so that we can reopen below
        adata.file.close()
    with h5py.File(filepath, mode) as f:
        if "X" in as_dense and isinstance(adata.X, (sparse.spmatrix, SparseDataset)):
            write_sparse_as_dense(f, "X", adata.X, dataset_kwargs=dataset_kwargs)
        elif not (adata.isbacked and Path(adata.filename) == Path(filepath)):
            # If adata.isbacked, X should already be up to date
            write_attribute(f, "X", adata.X, dataset_kwargs=dataset_kwargs)
        if "raw/X" in as_dense and isinstance(
            adata.raw.X, (sparse.spmatrix, SparseDataset)
        ):
            write_sparse_as_dense(
                f, "raw/X", adata.raw.X, dataset_kwargs=dataset_kwargs
            )
            write_attribute(f, "raw/var", adata.raw.var, dataset_kwargs=dataset_kwargs)
            write_attribute(
                f, "raw/varm", adata.raw.varm, dataset_kwargs=dataset_kwargs
            )
        else:
            write_attribute(f, "raw", adata.raw, dataset_kwargs=dataset_kwargs)
        write_attribute(f, "obs", adata.obs, dataset_kwargs=dataset_kwargs)
        write_attribute(f, "var", adata.var, dataset_kwargs=dataset_kwargs)
        write_attribute(f, "obsm", adata.obsm, dataset_kwargs=dataset_kwargs)
        write_attribute(f, "varm", adata.varm, dataset_kwargs=dataset_kwargs)
        write_attribute(f, "obsp", adata.obsp, dataset_kwargs=dataset_kwargs)
        write_attribute(f, "varp", adata.varp, dataset_kwargs=dataset_kwargs)
        write_attribute(f, "layers", adata.layers, dataset_kwargs=dataset_kwargs)
        write_attribute(f, "uns", adata.uns, dataset_kwargs=dataset_kwargs)
    if adata.isbacked:
        adata.file.open(filepath, "r+")


def _write_method(cls: Type[T]) -> Callable[[H5Group, str, T], None]:
    return _find_impl(cls, H5AD_WRITE_REGISTRY)


@write_attribute.register(h5py.File)
@write_attribute.register(h5py.Group)
def write_attribute_h5ad(f: H5Group, key: str, value, *args, **kwargs):
    if key in f:
        del f[key]
    _write_method(type(value))(f, key, value, *args, **kwargs)


def write_raw(f, key, value, dataset_kwargs=MappingProxyType({})):
    group = f.create_group(key)
    group.attrs["encoding-type"] = "raw"
    group.attrs["encoding-version"] = EncodingVersions.raw.value
    group.attrs["shape"] = value.shape
    write_attribute(f, "raw/X", value.X, dataset_kwargs=dataset_kwargs)
    write_attribute(f, "raw/var", value.var, dataset_kwargs=dataset_kwargs)
    write_attribute(f, "raw/varm", value.varm, dataset_kwargs=dataset_kwargs)


@report_write_key_on_error
def write_not_implemented(f, key, value, dataset_kwargs=MappingProxyType({})):
    # If it’s not an array, try and make it an array. If that fails, pickle it.
    # Maybe rethink that, maybe this should just pickle,
    # and have explicit implementations for everything else
    raise NotImplementedError(
        f"Failed to write value for {key}, "
        f"since a writer for type {type(value)} has not been implemented yet."
    )


@report_write_key_on_error
def write_basic(f, key, value, dataset_kwargs=MappingProxyType({})):
    f.create_dataset(key, data=value, **dataset_kwargs)


@report_write_key_on_error
def write_list(f, key, value, dataset_kwargs=MappingProxyType({})):
    write_array(f, key, np.array(value), dataset_kwargs=dataset_kwargs)


@report_write_key_on_error
def write_none(f, key, value, dataset_kwargs=MappingProxyType({})):
    pass


@report_write_key_on_error
def write_scalar(f, key, value, dataset_kwargs=MappingProxyType({})):
    # Can’t compress scalars, error is thrown
    if "compression" in dataset_kwargs:
        dataset_kwargs = dict(dataset_kwargs)
        dataset_kwargs.pop("compression")
    write_array(f, key, np.array(value), dataset_kwargs=dataset_kwargs)


@report_write_key_on_error
def write_array(f, key, value, dataset_kwargs=MappingProxyType({})):
    # Convert unicode to fixed length strings
    if value.dtype.kind in {"U", "O"}:
        value = value.astype(h5py.special_dtype(vlen=str))
    elif value.dtype.names is not None:
        value = _to_hdf5_vlen_strings(value)
    f.create_dataset(key, data=value, **dataset_kwargs)


@report_write_key_on_error
def write_sparse_compressed(
    f, key, value, fmt: Literal["csr", "csc"], dataset_kwargs=MappingProxyType({})
):
    g = f.create_group(key)
    g.attrs["encoding-type"] = f"{fmt}_matrix"
    g.attrs["encoding-version"] = EncodingVersions[f"{fmt}_matrix"].value
    g.attrs["shape"] = value.shape

    # Allow resizing
    if "maxshape" not in dataset_kwargs:
        dataset_kwargs = dict(maxshape=(None,), **dataset_kwargs)

    g.create_dataset("data", data=value.data, **dataset_kwargs)
    g.create_dataset("indices", data=value.indices, **dataset_kwargs)
    g.create_dataset("indptr", data=value.indptr, **dataset_kwargs)


write_csr = partial(write_sparse_compressed, fmt="csr")
write_csc = partial(write_sparse_compressed, fmt="csc")


@report_write_key_on_error
def write_sparse_dataset(f, key, value, dataset_kwargs=MappingProxyType({})):
    write_sparse_compressed(
        f, key, value.to_backed(), fmt=value.format_str, dataset_kwargs=dataset_kwargs,
    )


@report_write_key_on_error
def write_sparse_as_dense(f, key, value, dataset_kwargs=MappingProxyType({})):
    real_key = None  # Flag for if temporary key was used
    if key in f:
        if (
            isinstance(value, (h5py.Group, h5py.Dataset, SparseDataset))
            and value.file.filename == f.filename
        ):  # Write to temporary key before overwriting
            real_key = key
            # Transform key to temporary, e.g. raw/X -> raw/_X, or X -> _X
            key = re.sub(r"(.*)(\w(?!.*/))", r"\1_\2", key.rstrip("/"))
        else:
            del f[key]  # Wipe before write
    dset = f.create_dataset(key, shape=value.shape, dtype=value.dtype, **dataset_kwargs)
    compressed_axis = int(isinstance(value, sparse.csc_matrix))
    for idx in idx_chunks_along_axis(value.shape, compressed_axis, 1000):
        dset[idx] = value[idx].toarray()
    if real_key is not None:
        del f[real_key]
        f[real_key] = f[key]
        del f[key]


@report_write_key_on_error
def write_dataframe(f, key, df, dataset_kwargs=MappingProxyType({})):
    # Check arguments
    for reserved in ("__categories", "_index"):
        if reserved in df.columns:
            raise ValueError(f"{reserved!r} is a reserved name for dataframe columns.")
    group = f.create_group(key)
    group.attrs["encoding-type"] = "dataframe"
    group.attrs["encoding-version"] = EncodingVersions.dataframe.value
    group.attrs["column-order"] = list(df.columns)

    if df.index.name is not None:
        index_name = df.index.name
    else:
        index_name = "_index"
    group.attrs["_index"] = index_name

    write_series(group, index_name, df.index, dataset_kwargs=dataset_kwargs)
    for colname, series in df.items():
        write_series(group, colname, series, dataset_kwargs=dataset_kwargs)


@report_write_key_on_error
def write_series(group, key, series, dataset_kwargs=MappingProxyType({})):
    # group here is an h5py type, otherwise categoricals won’t write
    if series.dtype == object:  # Assuming it’s string
        group.create_dataset(
            key,
            data=series.values,
            dtype=h5py.special_dtype(vlen=str),
            **dataset_kwargs,
        )
    elif is_categorical_dtype(series):
        # This should work for categorical Index and Series
        categorical: pd.Categorical = series.values
        categories: np.ndarray = categorical.categories.values
        codes: np.ndarray = categorical.codes
        category_key = f"__categories/{key}"

        write_array(group, category_key, categories, dataset_kwargs=dataset_kwargs)
        write_array(group, key, codes, dataset_kwargs=dataset_kwargs)

        group[key].attrs["categories"] = group[category_key].ref
        group[category_key].attrs["ordered"] = categorical.ordered
    else:
        group[key] = series.values


def write_mapping(f, key, value, dataset_kwargs=MappingProxyType({})):
    for sub_key, sub_value in value.items():
        write_attribute(f, f"{key}/{sub_key}", sub_value, dataset_kwargs=dataset_kwargs)


H5AD_WRITE_REGISTRY = {
    Raw: write_raw,
    object: write_not_implemented,
    h5py.Dataset: write_basic,
    list: write_list,
    type(None): write_none,
    str: write_scalar,
    float: write_scalar,
    np.floating: write_scalar,
    bool: write_scalar,
    np.bool_: write_scalar,
    int: write_scalar,
    np.integer: write_scalar,
    np.ndarray: write_array,
    sparse.csr_matrix: write_csr,
    sparse.csc_matrix: write_csc,
    SparseDataset: write_sparse_dataset,
    pd.DataFrame: write_dataframe,
    cabc.Mapping: write_mapping,
}


def read_h5ad_backed(filename: Union[str, Path], mode: Literal["r", "r+"]) -> AnnData:
    d = dict(filename=filename, filemode=mode)

    f = h5py.File(filename, mode)

    attributes = ["obsm", "varm", "obsp", "varp", "uns", "layers"]
    df_attributes = ["obs", "var"]

    d.update({k: read_attribute(f[k]) for k in attributes if k in f})
    for k in df_attributes:
        if k in f:  # Backwards compat
            d[k] = read_dataframe(f[k])

    d["raw"] = _read_raw(f, attrs={"var", "varm"})

    X_dset = f.get("X", None)
    if X_dset is None:
        pass
    elif isinstance(X_dset, h5py.Group):
        d["dtype"] = X_dset["data"].dtype
    elif hasattr(X_dset, "dtype"):
        d["dtype"] = f["X"].dtype
    else:
        raise ValueError()

    _clean_uns(d)

    return AnnData(**d)


def read_h5ad(
    filename: Union[str, Path],
    backed: Union[Literal["r", "r+"], bool, None] = None,
    *,
    as_sparse: Sequence[str] = (),
    as_sparse_fmt: Type[sparse.spmatrix] = sparse.csr_matrix,
    chunk_size: int = 6000,  # TODO, probably make this 2d chunks
) -> AnnData:
    """\
    Read `.h5ad`-formatted hdf5 file.

    Parameters
    ----------
    filename
        File name of data file.
    backed
        If `'r'`, load :class:`~anndata.AnnData` in `backed` mode
        instead of fully loading it into memory (`memory` mode).
        If you want to modify backed attributes of the AnnData object,
        you need to choose `'r+'`.
    as_sparse
        If an array was saved as dense, passing its name here will read it as
        a sparse_matrix, by chunk of size `chunk_size`.
    as_sparse_fmt
        Sparse format class to read elements from `as_sparse` in as.
    chunk_size
        Used only when loading sparse dataset that is stored as dense.
        Loading iterates through chunks of the dataset of this row size
        until it reads the whole dataset.
        Higher size means higher memory consumption and higher (to a point)
        loading speed.
    """
    if backed not in {None, False}:
        mode = backed
        if mode is True:
            mode = "r+"
        assert mode in {"r", "r+"}
        return read_h5ad_backed(filename, mode)

    if as_sparse_fmt not in (sparse.csr_matrix, sparse.csc_matrix):
        raise NotImplementedError(
            "Dense formats can only be read to CSR or CSC matrices at this time."
        )
    if isinstance(as_sparse, str):
        as_sparse = [as_sparse]
    else:
        as_sparse = list(as_sparse)
    for i in range(len(as_sparse)):
        if as_sparse[i] in {("raw", "X"), "raw.X"}:
            as_sparse[i] = "raw/X"
        elif as_sparse[i] not in {"raw/X", "X"}:
            raise NotImplementedError(
                "Currently only `X` and `raw/X` can be read as sparse."
            )

    rdasp = partial(
        read_dense_as_sparse, sparse_format=as_sparse_fmt, axis_chunk=chunk_size,
    )

    with h5py.File(filename, "r") as f:
        d = {}
        for k in f.keys():
            # Backwards compat for old raw
            if k == "raw" or k.startswith("raw."):
                continue
            if k == "X" and "X" in as_sparse:
                d[k] = rdasp(f[k])
            elif k == "raw":
                assert False, "unexpected raw format"
            elif k in {"obs", "var"}:
                d[k] = read_dataframe(f[k])
            else:  # Base case
                d[k] = read_attribute(f[k])

        d["raw"] = _read_raw(f, as_sparse, rdasp)

        X_dset = f.get("X", None)
        if X_dset is None:
            pass
        elif isinstance(X_dset, h5py.Group):
            d["dtype"] = X_dset["data"].dtype
        elif hasattr(X_dset, "dtype"):
            d["dtype"] = f["X"].dtype
        else:
            raise ValueError()

    _clean_uns(d)  # backwards compat

    return AnnData(**d)


def _read_raw(
    f: Union[h5py.File, AnnDataFileManager],
    as_sparse: Collection[str] = (),
    rdasp: Callable[[h5py.Dataset], sparse.spmatrix] = None,
    *,
    attrs: Collection[str] = ("X", "var", "varm"),
):
    if as_sparse:
        assert rdasp is not None, "must supply rdasp if as_sparse is supplied"
    raw = {}
    if "X" in attrs and "raw/X" in f:
        read_x = rdasp if "raw/X" in as_sparse else read_attribute
        raw["X"] = read_x(f["raw/X"])
    for v in ("var", "varm"):
        if v in attrs and f"raw/{v}" in f:
            raw[v] = read_attribute(f[f"raw/{v}"])
    return _read_legacy_raw(f, raw, read_dataframe, read_attribute, attrs=attrs)


@report_read_key_on_error
def read_dataframe_legacy(dataset) -> pd.DataFrame:
    """Read pre-anndata 0.7 dataframes."""
    df = pd.DataFrame(_from_fixed_length_strings(dataset[()]))
    df.set_index(df.columns[0], inplace=True)
    return df


@report_read_key_on_error
def read_dataframe(group) -> pd.DataFrame:
    if not isinstance(group, h5py.Group):
        return read_dataframe_legacy(group)
    columns = list(group.attrs["column-order"])
    idx_key = group.attrs["_index"]
    df = pd.DataFrame(
        {k: read_series(group[k]) for k in columns},
        index=read_series(group[idx_key]),
        columns=list(columns),
    )
    if idx_key != "_index":
        df.index.name = idx_key
    return df


@report_read_key_on_error
def read_series(dataset) -> Union[np.ndarray, pd.Categorical]:
    if "categories" in dataset.attrs:
        categories = dataset.attrs["categories"]
        if isinstance(categories, h5py.Reference):
            categories_dset = dataset.parent[dataset.attrs["categories"]]
            categories = categories_dset[...]
            ordered = bool(categories_dset.attrs.get("ordered", False))
        else:
            # TODO: remove this code at some point post 0.7
            # TODO: Add tests for this
            warn(
                f"Your file {str(dataset.file.name)!r} has invalid categorical "
                "encodings due to being written from a development version of "
                "AnnData. Rewrite the file ensure you can read it in the future.",
                FutureWarning,
            )
        return pd.Categorical.from_codes(dataset[...], categories, ordered=ordered)
    else:
        return dataset[...]


# @report_read_key_on_error
# def read_sparse_dataset_backed(group: h5py.Group) -> sparse.spmatrix:
#     return SparseDataset(group)


@read_attribute.register(h5py.Group)
@report_read_key_on_error
def read_group(group: h5py.Group) -> Union[dict, pd.DataFrame, sparse.spmatrix]:
    if "h5sparse_format" in group.attrs:  # Backwards compat
        return SparseDataset(group).to_memory()

    encoding_type = group.attrs.get("encoding-type")
    if encoding_type:
        EncodingVersions[encoding_type].check(
            group.name, group.attrs["encoding-version"]
        )
    if encoding_type in {None, "raw"}:
        pass
    elif encoding_type == "dataframe":
        return read_dataframe(group)
    elif encoding_type in {"csr_matrix", "csc_matrix"}:
        return SparseDataset(group).to_memory()
    else:
        raise ValueError(f"Unfamiliar `encoding-type`: {encoding_type}.")
    d = dict()
    for sub_key, sub_value in group.items():
        d[sub_key] = read_attribute(sub_value)
    return d


@read_attribute.register(h5py.Dataset)
@report_read_key_on_error
def read_dataset(dataset: h5py.Dataset):
    value = dataset[()]
    if not hasattr(value, "dtype"):
        return value
    elif isinstance(value.dtype, str):
        pass
    elif issubclass(value.dtype.type, np.string_):
        value = value.astype(str)
        # Backwards compat, old datasets have strings as one element 1d arrays
        if len(value) == 1:
            return value[0]
    elif len(value.dtype.descr) > 1:  # Compound dtype
        # For backwards compat, now strings are written as variable length
        value = _from_fixed_length_strings(value)
    if value.shape == ():
        value = value[()]
    return value


@report_read_key_on_error
def read_dense_as_sparse(
    dataset: h5py.Dataset, sparse_format: sparse.spmatrix, axis_chunk: int
):
    if sparse_format == sparse.csr_matrix:
        return read_dense_as_csr(dataset, axis_chunk)
    elif sparse_format == sparse.csc_matrix:
        return read_dense_as_csc(dataset, axis_chunk)
    else:
        raise ValueError(f"Cannot read dense array as type: {sparse_format}")


def read_dense_as_csr(dataset, axis_chunk=6000):
    sub_matrices = []
    for idx in idx_chunks_along_axis(dataset.shape, 0, axis_chunk):
        dense_chunk = dataset[idx]
        sub_matrix = sparse.csr_matrix(dense_chunk)
        sub_matrices.append(sub_matrix)
    return sparse.vstack(sub_matrices, format="csr")


def read_dense_as_csc(dataset, axis_chunk=6000):
    sub_matrices = []
    for idx in idx_chunks_along_axis(dataset.shape, 1, axis_chunk):
        sub_matrix = sparse.csc_matrix(dataset[idx])
        sub_matrices.append(sub_matrix)
    return sparse.hstack(sub_matrices, format="csc")
