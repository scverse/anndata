from __future__ import annotations

import bz2
import gzip
from importlib.metadata import version
from os import PathLike, fspath
from pathlib import Path
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
from packaging.version import Version
from scipy import sparse

from .. import AnnData
from .._settings import settings
from ..compat import old_positionals, pandas_as_str
from .utils import is_float

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable, Iterator, Mapping, MutableMapping


@old_positionals("first_column_names", "dtype")
def read_csv(
    filename: PathLike[str] | str | Iterator[str],
    delimiter: str | None = ",",
    *,
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
    return read_text(
        filename, delimiter, first_column_names=first_column_names, dtype=dtype
    )


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
    row = dict(
        row_names=pandas_as_str(df.iloc[:, 0])
        if settings.restrict_index_types
        else df.iloc[:, 0]
    )
    col = dict(
        col_names=pandas_as_str(df.columns[1:])
        if settings.restrict_index_types
        else df.columns[1:]
    )
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
        rows_cols: list[dict[str, np.ndarray]] = [{}, {}]
        for iname, name in enumerate(["row_names", "col_names"]):
            if name in keys:
                rows_cols[iname][name] = f[name][()]
    adata = AnnData(X, rows_cols[0], rows_cols[1])
    return adata


def _fmt_loom_axis_attrs(
    input: MutableMapping, idx_name: str, dimm_mapping: Mapping[str, Iterable[str]]
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
        axis_df = axis_df.set_index(idx_name, drop=True)

    return axis_df, axis_mapping


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
    # https://github.com/scverse/anndata/issues/2477 for spmatrix
    path = fspath(filename)
    if Version(version("scipy")) >= Version("1.18.0rc1"):
        X = mmread(path, spmatrix=True)
    else:
        X = mmread(path)
    from scipy.sparse import csr_matrix

    return AnnData(csr_matrix(X.astype(dtype)))


@old_positionals("first_column_names", "dtype")
def read_text(
    filename: PathLike[str] | str | Iterator[str],
    delimiter: str | None = None,
    *,
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
        return _read_text(
            filename, delimiter, first_column_names=first_column_names, dtype=dtype
        )

    filename = Path(filename)
    if filename.suffix == ".gz":
        with gzip.open(str(filename), mode="rt") as f:
            return _read_text(
                f, delimiter, first_column_names=first_column_names, dtype=dtype
            )
    elif filename.suffix == ".bz2":
        with bz2.open(str(filename), mode="rt") as f:
            return _read_text(
                f, delimiter, first_column_names=first_column_names, dtype=dtype
            )
    else:
        with filename.open() as f:
            return _read_text(
                f, delimiter, first_column_names=first_column_names, dtype=dtype
            )


def _iter_lines(file_like: Iterable[str]) -> Generator[str, None, None]:
    """Helper for iterating only nonempty lines without line breaks"""
    for line in file_like:
        line = line.rstrip("\r\n")
        if line:
            yield line


def _read_text(  # noqa: PLR0912, PLR0915
    f: Iterator[str],
    delimiter: str | None,
    *,
    first_column_names: bool | None,
    dtype: str,
) -> AnnData:
    comments: list[str] = []
    rows: list[np.ndarray] = []
    lines = _iter_lines(f)
    col_names: list[str] = []
    row_names: list[str] = []
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
                rows.append(np.array(line_list[1:], dtype=dtype))
            else:
                rows.append(np.array(line_list, dtype=dtype))
            break
    if col_names:
        cols = np.array(col_names, dtype=str)
    # try reading col_names from the last comment line
    elif len(comments) > 0:
        # logg.msg("    assuming last comment line stores variable names", v=4)
        cols = np.array(comments[-1].split(), dtype=str)
    # just numbers as col_names
    else:
        # logg.msg("    did not find column names in file", v=4)
        cols = np.arange(len(rows[0])).astype(str)
    # read another line to check if first column contains row names or not
    if first_column_names is None:
        first_column_names = False
    for line in lines:
        line_list = line.split(delimiter)
        if first_column_names or not is_float(line_list[0]):
            # logg.msg("    assuming first column in file stores row names", v=4)
            first_column_names = True
            row_names.append(line_list[0])
            rows.append(np.array(line_list[1:], dtype=dtype))
        else:
            rows.append(np.array(line_list, dtype=dtype))
        break
    # if row names are just integers
    if len(rows) > 1 and rows[0].size != rows[1].size:
        # logg.msg(
        #     "    assuming first row stores column names and first column row names",
        #     v=4,
        # )
        first_column_names = True
        cols = np.array(rows[0]).astype(int).astype(str)
        row_names.append(rows[1][0].astype(int).astype(str))
        rows = [rows[1][1:]]
    # parse the file
    for line in lines:
        line_list = line.split(delimiter)
        if first_column_names:
            row_names.append(line_list[0])
            rows.append(np.array(line_list[1:], dtype=dtype))
        else:
            rows.append(np.array(line_list, dtype=dtype))
    # logg.msg("    read data into list of lists", t=True, v=4)
    # transform to array, this takes a long time and a lot of memory
    # but it’s actually the same thing as np.genfromtxt does
    # - we don’t use the latter as it would involve another slicing step
    #   in the end, to separate row_names from float data, slicing takes
    #   a lot of memory and CPU time
    if rows[0].size != rows[-1].size:
        msg = (
            f"Length of first line ({rows[0].size}) is different "
            f"from length of last line ({rows[-1].size})."
        )
        raise ValueError(msg)
    data = np.array(rows, dtype=dtype)
    # logg.msg("    constructed array from list of list", t=True, v=4)
    # transform row_names
    if not row_names:
        obs_names = np.arange(len(data)).astype(str)
        # logg.msg("    did not find row names in file", v=4)
    else:
        obs_names = np.array(row_names)
        for iname, name in enumerate(obs_names):
            obs_names[iname] = name.strip('"')
    # adapt col_names if necessary
    if cols.size > data.shape[1]:
        cols = cols[1:]
    for iname, name in enumerate(cols):
        cols[iname] = name.strip('"')
    return AnnData(
        data,
        obs=dict(obs_names=obs_names),
        var=dict(var_names=cols),
    )
