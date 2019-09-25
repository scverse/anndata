from collections.abc import Mapping
from functools import _find_impl, singledispatch
from pathlib import Path
from types import MappingProxyType
from typing import Callable, Optional, Type, TypeVar, Union
from warnings import warn

import h5py as h5py
import numpy as np
import pandas as pd
from pandas.api.types import is_categorical_dtype
from scipy import sparse

from .. import h5py as adh5py
from ..core.anndata import AnnData, Raw
from ..compat import _from_fixed_length_strings, _clean_uns
from .utils import report_key_on_error


H5Group = Union[h5py.Group, h5py.File, adh5py.Group, adh5py.File]
H5Dataset = Union[h5py.Dataset, adh5py.Dataset]
T = TypeVar("T")


def _to_hdf5_vlen_strings(value: np.ndarray) -> np.ndarray:
    """
    This corrects compound dtypes to work with hdf5 files.
    """
    new_dtype = []
    for dt_name, dt_type in value.dtype.descr:
        if dt_type[1] == "U":
            new_dtype.append((dt_name, h5py.special_dtype(vlen=str)))
        else:
            new_dtype.append((dt_name, dt_type))
    return value.astype(new_dtype)


def write_h5ad(
    filepath: Union[Path, str],
    adata: AnnData,
    force_dense: bool = False,
    dataset_kwargs: Mapping = MappingProxyType({}),
    **kwargs,
) -> None:
    """Write ``.h5ad``-formatted hdf5 file.

    .. note::

        Setting compression to ``'gzip'`` can save disk space but
        will slow down writing and subsequent reading. Prior to
        v0.6.16, this was the default for parameter
        ``compression``.

    Generally, if you have sparse data that are stored as a dense
    matrix, you can dramatically improve performance and reduce
    disk space by converting to a :class:`~scipy.sparse.csr_matrix`::

        from scipy.sparse import csr_matrix
        adata.X = csr_matrix(adata.X)

    Parameters
    ----------
    adata
        AnnData object to write.
    filename
        Filename of data file. Defaults to backing file.
    compression : ``None``,  {``'gzip'``, ``'lzf'``} (default: ``None``)
        See the h5py :ref:`dataset_compression`.
    force_dense
        Write sparse data as a dense matrix.
    """
    adata.strings_to_categoricals()
    if adata.raw is not None:
        adata.strings_to_categoricals(adata.raw.var)
    dataset_kwargs = dataset_kwargs.copy()
    dataset_kwargs.update(kwargs)
    filepath = Path(filepath)
    mode = "a" if adata.isbacked else "w"
    if adata.isbacked:  # close so that we can reopen below
        adata.file.close()
    with adh5py.File(filepath, mode, force_dense=force_dense) as f:
        if not (adata.isbacked and Path(adata.filename) == Path(filepath)):
            # Otherwise, X should already be up to date
            write_attribute(f, "X", adata.X, dataset_kwargs)
        write_attribute(f, "obs", adata.obs, dataset_kwargs)
        write_attribute(f, "var", adata.var, dataset_kwargs)
        write_attribute(f, "obsm", adata.obsm, dataset_kwargs)
        write_attribute(f, "varm", adata.varm, dataset_kwargs)
        write_attribute(f, "obsp", adata.obsp, dataset_kwargs)
        write_attribute(f, "varp", adata.varp, dataset_kwargs)
        write_attribute(f, "layers", adata.layers, dataset_kwargs)
        write_attribute(f, "uns", adata.uns, dataset_kwargs)
        write_attribute(f, "raw", adata.raw, dataset_kwargs)
    if adata.isbacked:
        adata.file.open(filepath, 'r+')


def _write_method(cls: Type[T]) -> Callable[[H5Group, str, T], None]:
    return _find_impl(cls, H5AD_WRITE_REGISTRY)


def write_attribute(f: H5Group, key: str, value, dataset_kwargs: Mapping):
    if key in f:
        del f[key]
    _write_method(type(value))(f, key, value, dataset_kwargs)


def write_raw(f, key, value, dataset_kwargs=MappingProxyType({})):
    group = f.create_group(key)
    group.attrs["encoding-type"] = "raw"
    group.attrs["encoding-version"] = "0.1.0"
    group.attrs["shape"] = value.shape
    write_attribute(f, "raw/X", value.X, dataset_kwargs)
    write_attribute(f, "raw/var", value.var, dataset_kwargs)
    write_attribute(f, "raw/varm", value.varm, dataset_kwargs)


def write_not_implemented(f, key, value, dataset_kwargs=MappingProxyType({})):
    # If it's not an array, try and make it an array. If that fails, pickle it.
    # Maybe rethink that, maybe this should just pickle, and have explicit implementations for everything else
    raise NotImplementedError(
        f"Failed to write value for {key}, since a writer for type {type(value)}"
        f" has not been implemented yet."
    )


def write_basic(f, key, value, dataset_kwargs=MappingProxyType({})):
    f.create_dataset(key, value, **dataset_kwargs)


def write_list(f, key, value, dataset_kwargs=MappingProxyType({})):
    write_array(f, key, np.array(value), dataset_kwargs)


def write_none(f, key, value, dataset_kwargs=MappingProxyType({})):
    pass


def write_scalar(f, key, value, dataset_kwargs=MappingProxyType({})):
    if (
        "compression" in dataset_kwargs
    ):  # Can't compress scalars, error is thrown
        dataset_kwargs = dataset_kwargs.copy()
        dataset_kwargs.pop("compression")
    write_array(f, key, np.array(value), dataset_kwargs)


def write_array(f, key, value, dataset_kwargs=MappingProxyType({})):
    # Convert unicode to fixed length strings
    if value.dtype.kind in {'U', 'O'}:
        value = value.astype(adh5py.special_dtype(vlen=str))
    elif value.dtype.names is not None:
        value = _to_hdf5_vlen_strings(value)
    f.create_dataset(key, data=value, **dataset_kwargs)


def write_dataframe(f, key, df, dataset_kwargs=MappingProxyType({})):
    # Check arguments
    for reserved in ("__categories", "_index"):
        if reserved in df.columns:
            raise ValueError(
                f"'{reserved}' is a reserved name for dataframe columns."
            )
    group = f.h5f.create_group(key)
    group.attrs["encoding-type"] = "dataframe"
    group.attrs["encoding-version"] = "0.1.0"
    group.attrs["column-order"] = list(df.columns)

    if df.index.name is not None:
        index_name = df.index.name
    else:
        index_name = "_index"
    group.attrs["_index"] = index_name

    write_series(group, index_name, df.index, dataset_kwargs)
    for colname, series in df.items():
        write_series(group, colname, series, dataset_kwargs)


def write_series(group, key, series, dataset_kwargs=MappingProxyType({})):
    # group here is an h5py type, otherwise categoricals won't write
    if series.dtype == object:  # Assuming it's string
        group.create_dataset(
            key,
            data=series.values,
            dtype=h5py.special_dtype(vlen=str),
            **dataset_kwargs,
        )
    elif is_categorical_dtype(series):
        cats = series.cat.categories.values
        codes = series.cat.codes.values
        category_key = f"__categories/{key}"

        write_array(group, category_key, cats, dataset_kwargs)
        write_array(group, key, codes, dataset_kwargs)

        group[key].attrs["categories"] = group[category_key].ref
        group[category_key].attrs["ordered"] = series.cat.ordered
    else:
        group[key] = series.values


def write_mapping(f, key, value, dataset_kwargs=MappingProxyType({})):
    for sub_key, sub_value in value.items():
        write_attribute(f, f"{key}/{sub_key}", sub_value, dataset_kwargs)


H5AD_WRITE_REGISTRY = {
    Raw: write_raw,
    object: write_not_implemented,
    adh5py.Dataset: write_basic,
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
    sparse.spmatrix: write_basic,
    adh5py.SparseDataset: write_basic,
    pd.DataFrame: write_dataframe,
    Mapping: write_mapping,
}


def read_h5ad_backed(filename: Union[str, Path], mode: str) -> AnnData:
    d = {"filename": filename, "filemode": mode}
    raw = {}

    f = adh5py.File(filename, mode)

    attributes = ["obsm", "varm", "obsp", "varp", "uns", "layers"]
    df_attributes = ["obs", "var"]

    d.update({k: read_attribute(f[k]) for k in attributes if k in f})
    for k in df_attributes:
        if k in f:  # Backwards compat
            d[k] = read_dataframe(f[k])

    if "raw" in f:
        if "raw/var" in f:
            raw["var"] = read_attribute(f["raw/var"])
        if "raw/varm" in f:
            raw["varm"] = read_attribute(f["raw/varm"])
    else:  # Legacy case
        if "raw.var" in f:
            raw["var"] = read_dataframe(f["raw.var"])  # Backwards compat
        if "raw.varm" in f:
            raw["varm"] = read_attribute(f["raw.varm"])

    if len(raw) > 0:
        d["raw"] = raw
    if f.get("X", None) is not None:
        d["dtype"] = f["X"].dtype

    _clean_uns(d)

    return AnnData(**d)


# TODO: chunks, what are possible values for chunk_size?
def read_h5ad(
    filename: Union[str, Path],
    backed: Optional[Union[str, bool]] = None,
    chunk_size=None,
) -> AnnData:
    """Read ``.h5ad``-formatted hdf5 file.

    Parameters
    ----------
    filename
        File name of data file.
    backed : {``None``, ``'r'``, ``'r+'``}
        If ``'r'``, load :class:`~anndata.AnnData` in ``backed`` mode instead
        of fully loading it into memory (`memory` mode). If you want to modify
        backed attributes of the AnnData object, you need to choose ``'r+'``.
    chunk_size
        Used only when loading sparse dataset that is stored as dense.
        Loading iterates through chunks of the dataset of this row size
        until it reads the whole dataset.
        Higher size means higher memory consumption and higher loading speed.
    """
    if backed not in {None, False}:
        mode = backed
        if mode is True:
            mode = "r+"
        assert mode in {"r", "r+"}
        return read_h5ad_backed(filename, mode)

    with adh5py.File(filename, "r") as f:
        d = {}
        for k in f.keys():
            # Backwards compat for old raw
            if k.startswith("raw."):
                continue
            if k in {"obs", "var"}:
                d[k] = read_dataframe(f[k])
            else:  # Base case
                d[k] = read_attribute(f[k])

        # Backwards compat for reading legacy raw
        raw = {}
        if "raw.var" in f:
            raw["var"] = read_dataframe(f["raw.var"])  # Backwards compat
        if "raw.varm" in f:
            raw["varm"] = read_attribute(f["raw.varm"])
        if "raw.X" in f:
            raw["X"] = read_attribute(f["raw.X"])
        if len(raw) > 0:
            assert (
                "raw" not in d
            ), f"File {filename} has both legacy and current raw formats."
            d["raw"] = raw

    if d.get("X", None) is not None:
        d["dtype"] = d["X"].dtype

    _clean_uns(d)  # backwards compat

    return AnnData(**d)


@singledispatch
def read_attribute(value):
    raise NotImplementedError()


@report_key_on_error
def read_dataframe_legacy(dataset) -> pd.DataFrame:
    """
    Read pre-anndata 0.7 dataframes.
    """
    df = pd.DataFrame(_from_fixed_length_strings(dataset[()]))
    df.set_index(df.columns[0], inplace=True)
    return df


@report_key_on_error
def read_dataframe(group) -> pd.DataFrame:
    if not isinstance(group, adh5py.Group):
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


@report_key_on_error
def read_series(dataset) -> Union[np.ndarray, pd.Categorical]:
    if "categories" in dataset.attrs:
        categories = dataset.attrs["categories"]
        if isinstance(categories, h5py.Reference):
            categories_dset = dataset.parent[dataset.attrs["categories"]]
            categories = categories_dset[...]
            ordered = categories_dset.attrs.get("ordered", False)
        else:
            # TODO: remove this code at some point post 0.7
            # TODO: Add tests for this
            warn(
                f"Your file '{dataset.file.name}' has invalid categorical "
                "encodings due to being written from a development version of "
                "AnnData. Rewrite the file ensure you can read it in the future.",
                FutureWarning,
            )
        return pd.Categorical.from_codes(
            dataset[...], categories, ordered=ordered
        )
    else:
        return dataset[...]


@read_attribute.register(adh5py.Group)
@report_key_on_error
def read_group(group: adh5py.Group) -> Union[dict, pd.DataFrame]:
    if group.attrs.get("encoding-type", "") == "dataframe":
        return read_dataframe(group)
    d = dict()
    for sub_key, sub_value in group.items():
        d[sub_key] = read_attribute(sub_value)
    return d


@read_attribute.register(adh5py.Dataset)
@report_key_on_error
def read_dataset(dataset: adh5py.Dataset):
    value = dataset[()]
    if not hasattr(value, "dtype"):
        return value
    elif isinstance(value.dtype, str):
        pass
    elif issubclass(value.dtype.type, np.string_):
        value = value.astype(str)
        if (
            len(value) == 1
        ):  # Backwards compat, old datasets have strings written as one element 1d arrays
            return value[0]
    elif len(value.dtype.descr) > 1:  # Compound dtype
        # For backwards compat, now strings are written as variable length
        value = _from_fixed_length_strings(value)
    if value.shape == ():
        value = value[()]
    return value


@read_attribute.register(type(None))
def read_attribute_none(value) -> None:
    return None


@read_attribute.register(adh5py.SparseDataset)
@report_key_on_error
def read_sparse_dataset(value) -> sparse.spmatrix:
    return value.value
