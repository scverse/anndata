from __future__ import annotations

import re
from functools import partial
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, TypeVar
from warnings import warn

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from anndata._warnings import OldFormatWarning

from .._core.anndata import AnnData
from .._core.file_backing import filename
from .._core.sparse_dataset import BaseCompressedSparseDataset
from ..compat import (
    CSMatrix,
    _clean_uns,
    _decode_structured_array,
    _from_fixed_length_strings,
)
from ..experimental import read_dispatched
from .specs import read_elem, write_elem
from .specs.registry import IOSpec, write_spec
from .utils import (
    H5PY_V3,
    _read_legacy_raw,
    idx_chunks_along_axis,
    no_write_dataset_2d,
    report_read_key_on_error,
    report_write_key_on_error,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Collection, Mapping, Sequence
    from os import PathLike
    from typing import Any, Literal

    from .._core.file_backing import AnnDataFileManager

T = TypeVar("T")


@no_write_dataset_2d
def write_h5ad(
    filepath: PathLike[str] | str,
    adata: AnnData,
    *,
    as_dense: Sequence[str] = (),
    convert_strings_to_categoricals: bool = True,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
    **kwargs,
) -> None:
    """See :meth:`~anndata.AnnData.write_h5ad`."""
    if isinstance(as_dense, str):
        as_dense = [as_dense]
    if "raw.X" in as_dense:
        as_dense = list(as_dense)
        as_dense[as_dense.index("raw.X")] = "raw/X"
    if any(val not in {"X", "raw/X"} for val in as_dense):
        msg = "Currently, only `X` and `raw/X` are supported values in `as_dense`"
        raise NotImplementedError(msg)
    if "raw/X" in as_dense and adata.raw is None:
        msg = "Cannot specify writing `raw/X` to dense if it doesn’t exist."
        raise ValueError(msg)

    if convert_strings_to_categoricals:
        adata.strings_to_categoricals()
        if adata.raw is not None:
            adata.strings_to_categoricals(adata.raw.var)
    dataset_kwargs = {**dataset_kwargs, **kwargs}
    filepath = Path(filepath)
    mode = "a" if adata.isbacked else "w"
    if adata.isbacked:  # close so that we can reopen below
        adata.file.close()

    with h5py.File(filepath, mode) as f:
        # TODO: Use spec writing system for this
        # Currently can't use write_dispatched here because this function is also called to do an
        # inplace update of a backed object, which would delete "/"
        f = f["/"]
        f.attrs.setdefault("encoding-type", "anndata")
        f.attrs.setdefault("encoding-version", "0.1.0")

        if "X" in as_dense and isinstance(
            adata.X, CSMatrix | BaseCompressedSparseDataset
        ):
            write_sparse_as_dense(f, "X", adata.X, dataset_kwargs=dataset_kwargs)
        elif not (adata.isbacked and Path(adata.filename) == Path(filepath)):
            # If adata.isbacked, X should already be up to date
            write_elem(f, "X", adata.X, dataset_kwargs=dataset_kwargs)
        if "raw/X" in as_dense and isinstance(
            adata.raw.X, CSMatrix | BaseCompressedSparseDataset
        ):
            write_sparse_as_dense(
                f, "raw/X", adata.raw.X, dataset_kwargs=dataset_kwargs
            )
            write_elem(f, "raw/var", adata.raw.var, dataset_kwargs=dataset_kwargs)
            write_elem(
                f, "raw/varm", dict(adata.raw.varm), dataset_kwargs=dataset_kwargs
            )
        elif adata.raw is not None:
            write_elem(f, "raw", adata.raw, dataset_kwargs=dataset_kwargs)
        write_elem(f, "obs", adata.obs, dataset_kwargs=dataset_kwargs)
        write_elem(f, "var", adata.var, dataset_kwargs=dataset_kwargs)
        write_elem(f, "obsm", dict(adata.obsm), dataset_kwargs=dataset_kwargs)
        write_elem(f, "varm", dict(adata.varm), dataset_kwargs=dataset_kwargs)
        write_elem(f, "obsp", dict(adata.obsp), dataset_kwargs=dataset_kwargs)
        write_elem(f, "varp", dict(adata.varp), dataset_kwargs=dataset_kwargs)
        write_elem(f, "layers", dict(adata.layers), dataset_kwargs=dataset_kwargs)
        write_elem(f, "uns", dict(adata.uns), dataset_kwargs=dataset_kwargs)


@report_write_key_on_error
@write_spec(IOSpec("array", "0.2.0"))
def write_sparse_as_dense(
    f: h5py.Group,
    key: str,
    value: CSMatrix | BaseCompressedSparseDataset,
    *,
    dataset_kwargs: Mapping[str, Any] = MappingProxyType({}),
):
    real_key = None  # Flag for if temporary key was used
    if key in f:
        if isinstance(value, BaseCompressedSparseDataset) and (
            filename(value.group) == filename(f)
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


def read_h5ad_backed(
    filename: str | PathLike[str], mode: Literal["r", "r+"]
) -> AnnData:
    d = dict(filename=filename, filemode=mode)

    f = h5py.File(filename, mode)

    attributes = ["obsm", "varm", "obsp", "varp", "uns", "layers"]
    df_attributes = ["obs", "var"]

    if "encoding-type" in f.attrs:
        attributes.extend(df_attributes)
    else:
        for k in df_attributes:
            if k in f:  # Backwards compat
                d[k] = read_dataframe(f[k])

    d.update({k: read_elem(f[k]) for k in attributes if k in f})

    d["raw"] = _read_raw(f, attrs={"var", "varm"})

    adata = AnnData(**d)

    # Backwards compat to <0.7
    if isinstance(f["obs"], h5py.Dataset):
        _clean_uns(adata)

    return adata


def read_h5ad(
    filename: PathLike[str] | str,
    backed: Literal["r", "r+"] | bool | None = None,
    *,
    as_sparse: Sequence[str] = (),
    as_sparse_fmt: type[CSMatrix] = sparse.csr_matrix,
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

        Currently, `backed` only support updates to `X`. That means any
        changes to other slots like `obs` will not be written to disk in
        `backed` mode. If you would like save changes made to these slots
        of a `backed` :class:`~anndata.AnnData`, write them to a new file
        (see :meth:`~anndata.AnnData.write`). For an example, see
        :ref:`read-partial`.
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
        msg = "Dense formats can only be read to CSR or CSC matrices at this time."
        raise NotImplementedError(msg)
    as_sparse = [as_sparse] if isinstance(as_sparse, str) else list(as_sparse)
    for i in range(len(as_sparse)):
        if as_sparse[i] in {("raw", "X"), "raw.X"}:
            as_sparse[i] = "raw/X"
        elif as_sparse[i] not in {"raw/X", "X"}:
            msg = "Currently only `X` and `raw/X` can be read as sparse."
            raise NotImplementedError(msg)

    rdasp = partial(
        read_dense_as_sparse, sparse_format=as_sparse_fmt, axis_chunk=chunk_size
    )

    with h5py.File(filename, "r") as f:

        def callback(func, elem_name: str, elem, iospec):
            if iospec.encoding_type == "anndata" or elem_name.endswith("/"):
                return AnnData(
                    **{
                        # This is covering up backwards compat in the anndata initializer
                        # In most cases we should be able to call `func(elen[k])` instead
                        k: read_dispatched(elem[k], callback)
                        for k in elem
                        if not k.startswith("raw.")
                    }
                )
            elif elem_name.startswith("/raw."):
                return None
            elif elem_name == "/X" and "X" in as_sparse:
                return rdasp(elem)
            elif elem_name == "/raw":
                return _read_raw(f, as_sparse, rdasp)
            elif elem_name in {"/obs", "/var"}:
                # Backwards compat
                return read_dataframe(elem)
            return func(elem)

        adata = read_dispatched(f, callback=callback)

        # Backwards compat (should figure out which version)
        if "raw.X" in f:
            raw = AnnData(**_read_raw(f, as_sparse, rdasp))
            raw.obs_names = adata.obs_names
            adata.raw = raw

        # Backwards compat to <0.7
        if isinstance(f["obs"], h5py.Dataset):
            _clean_uns(adata)

    return adata


def _read_raw(
    f: h5py.File | AnnDataFileManager,
    as_sparse: Collection[str] = (),
    rdasp: Callable[[h5py.Dataset], CSMatrix] | None = None,
    *,
    attrs: Collection[str] = ("X", "var", "varm"),
) -> dict:
    if as_sparse:
        assert rdasp is not None, "must supply rdasp if as_sparse is supplied"
    raw = {}
    if "X" in attrs and "raw/X" in f:
        read_x = rdasp if "raw/X" in as_sparse else read_elem
        raw["X"] = read_x(f["raw/X"])
    for v in ("var", "varm"):
        if v in attrs and f"raw/{v}" in f:
            raw[v] = read_elem(f[f"raw/{v}"])
    return _read_legacy_raw(f, raw, read_dataframe, read_elem, attrs=attrs)


@report_read_key_on_error
def read_dataframe_legacy(dataset: h5py.Dataset) -> pd.DataFrame:
    """Read pre-anndata 0.7 dataframes."""
    msg = (
        f"{dataset.name!r} was written with a very old version of AnnData. "
        "Consider rewriting it."
    )
    warn(msg, OldFormatWarning, stacklevel=2)
    if H5PY_V3:
        df = pd.DataFrame(
            _decode_structured_array(
                _from_fixed_length_strings(dataset[()]), dtype=dataset.dtype
            )
        )
    else:
        df = pd.DataFrame(_from_fixed_length_strings(dataset[()]))
    df.set_index(df.columns[0], inplace=True)
    return df


def read_dataframe(group: h5py.Group | h5py.Dataset) -> pd.DataFrame:
    """Backwards compat function"""
    if not isinstance(group, h5py.Group):
        return read_dataframe_legacy(group)
    else:
        return read_elem(group)


@report_read_key_on_error
def read_dataset(dataset: h5py.Dataset):
    if H5PY_V3:
        string_dtype = h5py.check_string_dtype(dataset.dtype)
        if (string_dtype is not None) and (string_dtype.encoding == "utf-8"):
            dataset = dataset.asstr()
    value = dataset[()]
    if not hasattr(value, "dtype"):
        return value
    elif isinstance(value.dtype, str):
        pass
    elif issubclass(value.dtype.type, np.bytes_):
        value = value.astype(str)
        # Backwards compat, old datasets have strings as one element 1d arrays
        if len(value) == 1:
            return value[0]
    elif len(value.dtype.descr) > 1:  # Compound dtype
        # For backwards compat, now strings are written as variable length
        dtype = value.dtype
        value = _from_fixed_length_strings(value)
        if H5PY_V3:
            value = _decode_structured_array(value, dtype=dtype)
    if value.shape == ():
        value = value[()]
    return value


@report_read_key_on_error
def read_dense_as_sparse(
    dataset: h5py.Dataset, sparse_format: CSMatrix, axis_chunk: int
):
    if sparse_format == sparse.csr_matrix:
        return read_dense_as_csr(dataset, axis_chunk)
    elif sparse_format == sparse.csc_matrix:
        return read_dense_as_csc(dataset, axis_chunk)
    else:
        msg = f"Cannot read dense array as type: {sparse_format}"
        raise ValueError(msg)


def read_dense_as_csr(dataset: h5py.Dataset, axis_chunk: int = 6000):
    sub_matrices = []
    for idx in idx_chunks_along_axis(dataset.shape, 0, axis_chunk):
        dense_chunk = dataset[idx]
        sub_matrix = sparse.csr_matrix(dense_chunk)
        sub_matrices.append(sub_matrix)
    return sparse.vstack(sub_matrices, format="csr")


def read_dense_as_csc(dataset: h5py.Dataset, axis_chunk: int = 6000):
    sub_matrices = []
    for idx in idx_chunks_along_axis(dataset.shape, 1, axis_chunk):
        sub_matrix = sparse.csc_matrix(dataset[idx])
        sub_matrices.append(sub_matrix)
    return sparse.hstack(sub_matrices, format="csc")
