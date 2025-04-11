from __future__ import annotations

import bz2
import gzip
from collections import OrderedDict
from os import PathLike, fspath
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING
from warnings import warn

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from .. import AnnData
from ..compat import old_positionals
from .utils import is_float

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable, Iterator, Mapping


def read_csv(
    filename: PathLike[str] | str | Iterator[str],
    delimiter: str | None = ",",
    first_column_names: bool | None = None,
    dtype: str = "float32",
) -> AnnData:
    """\
    Read `.csv` file.

    Same as :func:`~anndata.io.read_text` but with default delimiter `','`.

    Parameters
    ----------
    filename
        Data file.
    delimiter
        Delimiter that separates data within text file.
        If `None`, will split at arbitrary number of white spaces,
        which is different from enforcing splitting at single white space `' '`.
    first_column_names
        Assume the first column stores row names.
    dtype
        Numpy data type.
    """
    return read_text(filename, delimiter, first_column_names, dtype)


def read_excel(
    filename: PathLike[str] | str, sheet: str | int, dtype: str = "float32"
) -> AnnData:
    """\
    Read `.xlsx` (Excel) file.

    Assumes that the first columns stores the row names and the first row the
    column names.

    Parameters
    ----------
    filename
        File name to read from.
    sheet
        Name of sheet in Excel file.
    """
    # rely on pandas for reading an excel file
    from pandas import read_excel

    df = read_excel(fspath(filename), sheet)
    X = df.values[:, 1:]
    row = dict(row_names=df.iloc[:, 0].values.astype(str))
    col = dict(col_names=np.array(df.columns[1:], dtype=str))
    return AnnData(X, row, col)


def read_umi_tools(filename: PathLike[str] | str, dtype=None) -> AnnData:
    """\
    Read a gzipped condensed count matrix from umi_tools.

    Parameters
    ----------
    filename
        File name to read from.
    """
    # import pandas for conversion of a dict of dicts into a matrix
    # import gzip to read a gzipped file :-)
    table = pd.read_table(filename, dtype={"gene": "category", "cell": "category"})

    X = sparse.csr_matrix(
        (table["count"], (table["cell"].cat.codes, table["gene"].cat.codes)),
        dtype=dtype,
    )
    obs = pd.DataFrame(index=pd.Index(table["cell"].cat.categories, name="cell"))
    var = pd.DataFrame(index=pd.Index(table["gene"].cat.categories, name="gene"))

    return AnnData(X=X, obs=obs, var=var)


def read_hdf(filename: PathLike[str] | str, key: str) -> AnnData:
    """\
    Read `.h5` (hdf5) file.

    Note: Also looks for fields `row_names` and `col_names`.

    Parameters
    ----------
    filename
        Filename of data file.
    key
        Name of dataset in the file.
    """
    with h5py.File(filename, "r") as f:
        # the following is necessary in Python 3, because only
        # a view and not a list is returned
        keys = list(f)
        if key == "":
            msg = (
                f"The file {filename} stores the following sheets:\n{keys}\n"
                f"Call read/read_hdf5 with one of them."
            )
            raise ValueError(msg)
        # read array
        X = f[key][()]
        # try to find row and column names
        rows_cols = [{}, {}]
        for iname, name in enumerate(["row_names", "col_names"]):
            if name in keys:
                rows_cols[iname][name] = f[name][()]
    adata = AnnData(X, rows_cols[0], rows_cols[1])
    return adata


def _fmt_loom_axis_attrs(
    input: Mapping, idx_name: str, dimm_mapping: Mapping[str, Iterable[str]]
) -> tuple[pd.DataFrame, Mapping[str, np.ndarray]]:
    axis_df = pd.DataFrame()
    axis_mapping = {}
    for key, names in dimm_mapping.items():
        axis_mapping[key] = np.array([input.pop(name) for name in names]).T

    for k, v in input.items():
        if v.ndim > 1 and v.shape[1] > 1:
            axis_mapping[k] = v
        else:
            axis_df[k] = v

    if idx_name in axis_df:
        axis_df.set_index(idx_name, drop=True, inplace=True)

    return axis_df, axis_mapping


@old_positionals(
    "sparse",
    "cleanup",
    "X_name",
    "obs_names",
    "obsm_names",
    "var_names",
    "varm_names",
    "dtype",
    "obsm_mapping",
    "varm_mapping",
)
def read_loom(  # noqa: PLR0912, PLR0913
    filename: PathLike[str] | str,
    *,
    sparse: bool = True,
    cleanup: bool = False,
    X_name: str = "spliced",
    obs_names: str = "CellID",
    obsm_names: Mapping[str, Iterable[str]] | None = None,
    var_names: str = "Gene",
    varm_names: Mapping[str, Iterable[str]] | None = None,
    dtype: str = "float32",
    obsm_mapping: Mapping[str, Iterable[str]] = MappingProxyType({}),
    varm_mapping: Mapping[str, Iterable[str]] = MappingProxyType({}),
    **kwargs,
) -> AnnData:
    """\
    Read `.loom`-formatted hdf5 file.

    This reads the whole file into memory.

    Beware that you have to explicitly state when you want to read the file as
    sparse data.

    Parameters
    ----------
    filename
        The filename.
    sparse
        Whether to read the data matrix as sparse.
    cleanup
        Whether to collapse all obs/var fields that only store
        one unique value into `.uns['loom-.']`.
    X_name
        Loompy key with which the data matrix :attr:`~anndata.AnnData.X` is initialized.
    obs_names
        Loompy key where the observation/cell names are stored.
    obsm_mapping
        Loompy keys which will be constructed into observation matrices
    var_names
        Loompy key where the variable/gene names are stored.
    varm_mapping
        Loompy keys which will be constructed into variable matrices
    **kwargs:
        Arguments to loompy.connect

    Example
    -------

    .. code:: python

        pbmc = anndata.io.read_loom(
            "pbmc.loom",
            sparse=True,
            X_name="lognorm",
            obs_names="cell_names",
            var_names="gene_names",
            obsm_mapping={
                "X_umap": ["umap_1", "umap_2"]
            }
        )
    """
    # Deprecations
    if obsm_names is not None:
        msg = (
            "Argument obsm_names has been deprecated in favour of `obsm_mapping`. "
            "In 0.9 this will be an error."
        )
        warn(msg, FutureWarning, stacklevel=2)
        if obsm_mapping != {}:
            msg = (
                "Received values for both `obsm_names` and `obsm_mapping`. This is "
                "ambiguous, only pass `obsm_mapping`."
            )
            raise ValueError(msg)
        obsm_mapping = obsm_names
    if varm_names is not None:
        msg = (
            "Argument varm_names has been deprecated in favour of `varm_mapping`. "
            "In 0.9 this will be an error."
        )
        warn(msg, FutureWarning, stacklevel=2)
        if varm_mapping != {}:
            msg = (
                "Received values for both `varm_names` and `varm_mapping`. This is "
                "ambiguous, only pass `varm_mapping`."
            )
            raise ValueError(msg)
        varm_mapping = varm_names

    filename = fspath(filename)  # allow passing pathlib.Path objects
    from loompy import connect

    if TYPE_CHECKING:
        from loompy import LoomConnection

    lc: LoomConnection
    with connect(filename, "r", **kwargs) as lc:
        assert lc.layers is not None

        if X_name not in lc.layers:
            X_name = ""
        X = lc.layers[X_name].sparse().T.tocsr() if sparse else lc.layers[X_name][()].T
        X = X.astype(dtype, copy=False)

        layers = OrderedDict()
        if X_name != "":
            layers["matrix"] = (
                lc.layers[""].sparse().T.tocsr() if sparse else lc.layers[""][()].T
            )
        for key, layer in lc.layers.items():
            if key != "":
                layers[key] = layer.sparse().T.tocsr() if sparse else layer[()].T

        # TODO: Figure out the singleton obs elements
        obs, obsm = _fmt_loom_axis_attrs(dict(lc.col_attrs), obs_names, obsm_mapping)
        var, varm = _fmt_loom_axis_attrs(dict(lc.row_attrs), var_names, varm_mapping)

        uns = {}
        if cleanup:
            uns_obs = {}
            for key in obs.columns:
                if len(obs[key].unique()) == 1:
                    uns_obs[key] = obs[key].iloc[0]
                    del obs[key]
            if uns_obs:
                uns["loom-obs"] = uns_obs
            uns_var = {}
            for key in var.columns:
                if len(var[key].unique()) == 1:
                    uns_var[key] = var[key].iloc[0]
                    del var[key]
            if uns_var:
                uns["loom-var"] = uns_var

        adata = AnnData(
            X,
            obs=obs,
            var=var,
            layers=layers,
            obsm=obsm if obsm else None,
            varm=varm if varm else None,
            uns=uns,
        )
    return adata


def read_mtx(filename: PathLike[str] | str, dtype: str = "float32") -> AnnData:
    """\
    Read `.mtx` file.

    Parameters
    ----------
    filename
        The filename.
    dtype
        Numpy data type.
    """
    from scipy.io import mmread

    # could be rewritten accounting for dtype to be more performant
    X = mmread(fspath(filename)).astype(dtype)
    from scipy.sparse import csr_matrix

    X = csr_matrix(X)
    return AnnData(X)


def read_text(
    filename: PathLike[str] | str | Iterator[str],
    delimiter: str | None = None,
    first_column_names: bool | None = None,
    dtype: str = "float32",
) -> AnnData:
    """\
    Read `.txt`, `.tab`, `.data` (text) file.

    Same as :func:`~anndata.io.read_csv` but with default delimiter `None`.

    Parameters
    ----------
    filename
        Data file, filename or stream.
    delimiter
        Delimiter that separates data within text file. If `None`, will split at
        arbitrary number of white spaces, which is different from enforcing
        splitting at single white space `' '`.
    first_column_names
        Assume the first column stores row names.
    dtype
        Numpy data type.
    """
    if not isinstance(filename, PathLike | str | bytes):
        return _read_text(filename, delimiter, first_column_names, dtype)

    filename = Path(filename)
    if filename.suffix == ".gz":
        with gzip.open(str(filename), mode="rt") as f:
            return _read_text(f, delimiter, first_column_names, dtype)
    elif filename.suffix == ".bz2":
        with bz2.open(str(filename), mode="rt") as f:
            return _read_text(f, delimiter, first_column_names, dtype)
    else:
        with filename.open() as f:
            return _read_text(f, delimiter, first_column_names, dtype)


def _iter_lines(file_like: Iterable[str]) -> Generator[str, None, None]:
    """Helper for iterating only nonempty lines without line breaks"""
    for line in file_like:
        line = line.rstrip("\r\n")
        if line:
            yield line


def _read_text(  # noqa: PLR0912, PLR0915
    f: Iterator[str],
    delimiter: str | None,
    first_column_names: bool | None,
    dtype: str,
) -> AnnData:
    comments = []
    data = []
    lines = _iter_lines(f)
    col_names = []
    row_names = []
    # read header and column names
    for line in lines:
        if line.startswith("#"):
            comment = line.lstrip("# ")
            if comment:
                comments.append(comment)
        else:
            if delimiter is not None and delimiter not in line:
                msg = f"Did not find delimiter {delimiter!r} in first line."
                raise ValueError(msg)
            line_list = line.split(delimiter)
            # the first column might be row names, so check the last
            if not is_float(line_list[-1]):
                col_names = line_list
                # logg.msg("    assuming first line in file stores column names", v=4)
            elif not is_float(line_list[0]) or first_column_names:
                first_column_names = True
                row_names.append(line_list[0])
                data.append(np.array(line_list[1:], dtype=dtype))
            else:
                data.append(np.array(line_list, dtype=dtype))
            break
    if not col_names:
        # try reading col_names from the last comment line
        if len(comments) > 0:
            # logg.msg("    assuming last comment line stores variable names", v=4)
            col_names = np.array(comments[-1].split())
        # just numbers as col_names
        else:
            # logg.msg("    did not find column names in file", v=4)
            col_names = np.arange(len(data[0])).astype(str)
    col_names = np.array(col_names, dtype=str)
    # read another line to check if first column contains row names or not
    if first_column_names is None:
        first_column_names = False
    for line in lines:
        line_list = line.split(delimiter)
        if first_column_names or not is_float(line_list[0]):
            # logg.msg("    assuming first column in file stores row names", v=4)
            first_column_names = True
            row_names.append(line_list[0])
            data.append(np.array(line_list[1:], dtype=dtype))
        else:
            data.append(np.array(line_list, dtype=dtype))
        break
    # if row names are just integers
    if len(data) > 1 and data[0].size != data[1].size:
        # logg.msg(
        #     "    assuming first row stores column names and first column row names",
        #     v=4,
        # )
        first_column_names = True
        col_names = np.array(data[0]).astype(int).astype(str)
        row_names.append(data[1][0].astype(int).astype(str))
        data = [data[1][1:]]
    # parse the file
    for line in lines:
        line_list = line.split(delimiter)
        if first_column_names:
            row_names.append(line_list[0])
            data.append(np.array(line_list[1:], dtype=dtype))
        else:
            data.append(np.array(line_list, dtype=dtype))
    # logg.msg("    read data into list of lists", t=True, v=4)
    # transform to array, this takes a long time and a lot of memory
    # but it’s actually the same thing as np.genfromtxt does
    # - we don’t use the latter as it would involve another slicing step
    #   in the end, to separate row_names from float data, slicing takes
    #   a lot of memory and CPU time
    if data[0].size != data[-1].size:
        msg = (
            f"Length of first line ({data[0].size}) is different "
            f"from length of last line ({data[-1].size})."
        )
        raise ValueError(msg)
    data = np.array(data, dtype=dtype)
    # logg.msg("    constructed array from list of list", t=True, v=4)
    # transform row_names
    if not row_names:
        row_names = np.arange(len(data)).astype(str)
        # logg.msg("    did not find row names in file", v=4)
    else:
        row_names = np.array(row_names)
        for iname, name in enumerate(row_names):
            row_names[iname] = name.strip('"')
    # adapt col_names if necessary
    if col_names.size > data.shape[1]:
        col_names = col_names[1:]
    for iname, name in enumerate(col_names):
        col_names[iname] = name.strip('"')
    return AnnData(
        data,
        obs=dict(obs_names=row_names),
        var=dict(var_names=col_names),
    )
