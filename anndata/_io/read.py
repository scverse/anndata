from pathlib import Path
from os import PathLike, fspath
from types import MappingProxyType
from typing import Union, Optional, Mapping, Tuple
from typing import Iterable, Iterator, Generator
from collections import OrderedDict
import gzip
import bz2
from warnings import warn

import h5py
import numpy as np
import pandas as pd

from .. import AnnData
from ..compat import _deprecate_positional_args
from .utils import is_float
from .h5ad import read_h5ad

try:
    from .zarr import read_zarr
except ImportError as e:  # noqa: F841

    def read_zarr(*_, **__):
        raise e


def read_csv(
    filename: Union[PathLike, Iterator[str]],
    delimiter: Optional[str] = ",",
    first_column_names: Optional[bool] = None,
    dtype: str = "float32",
) -> AnnData:
    """\
    Read `.csv` file.

    Same as :func:`~anndata.read_text` but with default delimiter `','`.

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
    filename: PathLike, sheet: Union[str, int], dtype: str = "float32"
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
    return AnnData(X, row, col, dtype=dtype)


def read_umi_tools(filename: PathLike, dtype: str = "float32") -> AnnData:
    """\
    Read a gzipped condensed count matrix from umi_tools.

    Parameters
    ----------
    filename
        File name to read from.
    """
    # import pandas for conversion of a dict of dicts into a matrix
    # import gzip to read a gzipped file :-)
    import gzip
    from pandas import DataFrame

    dod = {}  # this will contain basically everything
    fh = gzip.open(fspath(filename))
    _ = fh.readline()  # read the first line

    for line in fh:
        # gzip read bytes, hence the decoding
        t = line.decode("ascii").split("\t")
        try:
            dod[t[1]].update({t[0]: int(t[2])})
        except KeyError:
            dod[t[1]] = {t[0]: int(t[2])}

    df = DataFrame.from_dict(dod, orient="index")  # build the matrix
    df.fillna(value=0.0, inplace=True)  # many NaN, replace with zeros
    return AnnData(
        np.array(df),
        dict(obs_names=df.index),
        dict(var_names=df.columns),
        dtype=dtype,
    )


def read_hdf(filename: PathLike, key: str) -> AnnData:
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
        keys = [k for k in f.keys()]
        if key == "":
            raise ValueError(
                f"The file {filename} stores the following sheets:\n{keys}\n"
                f"Call read/read_hdf5 with one of them."
            )
        # read array
        X = f[key][()]
        # try to find row and column names
        rows_cols = [{}, {}]
        for iname, name in enumerate(["row_names", "col_names"]):
            if name in keys:
                rows_cols[iname][name] = f[name][()]
    adata = AnnData(X, rows_cols[0], rows_cols[1], dtype=X.dtype.name)
    return adata


def _fmt_loom_axis_attrs(
    input: Mapping, idx_name: str, dimm_mapping: Mapping[str, Iterable[str]]
) -> Tuple[pd.DataFrame, Mapping[str, np.ndarray]]:
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


@_deprecate_positional_args(version="0.9")
def read_loom(
    filename: PathLike,
    *,
    sparse: bool = True,
    cleanup: bool = False,
    X_name: str = "spliced",
    obs_names: str = "CellID",
    obsm_names: Optional[Mapping[str, Iterable[str]]] = None,
    var_names: str = "Gene",
    varm_names: Optional[Mapping[str, Iterable[str]]] = None,
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

        pbmc = anndata.read_loom(
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
        warn(
            "Argument obsm_names has been deprecated in favour of `obsm_mapping`. "
            "In 0.9 this will be an error.",
            FutureWarning,
        )
        if obsm_mapping != {}:
            raise ValueError(
                "Recieved values for both `obsm_names` and `obsm_mapping`. This is "
                "ambiguous, only pass `obsm_mapping`."
            )
        obsm_mapping = obsm_names
    if varm_names is not None:
        warn(
            "Argument varm_names has been deprecated in favour of `varm_mapping`. "
            "In 0.9 this will be an error.",
            FutureWarning,
        )
        if varm_mapping != {}:
            raise ValueError(
                "Recieved values for both `varm_names` and `varm_mapping`. This is "
                "ambiguous, only pass `varm_mapping`."
            )
        varm_mapping = varm_names

    filename = fspath(filename)  # allow passing pathlib.Path objects
    from loompy import connect

    with connect(filename, "r", **kwargs) as lc:
        if X_name not in lc.layers.keys():
            X_name = ""
        X = lc.layers[X_name].sparse().T.tocsr() if sparse else lc.layers[X_name][()].T

        layers = OrderedDict()
        if X_name != "":
            layers["matrix"] = (
                lc.layers[""].sparse().T.tocsr() if sparse else lc.layers[""][()].T
            )
        for key in lc.layers.keys():
            if key != "":
                layers[key] = (
                    lc.layers[key].sparse().T.tocsr()
                    if sparse
                    else lc.layers[key][()].T
                )

        # TODO: Figure out the singleton obs elements
        obs, obsm = _fmt_loom_axis_attrs(dict(lc.col_attrs), obs_names, obsm_mapping)
        var, varm = _fmt_loom_axis_attrs(dict(lc.row_attrs), var_names, varm_mapping)

        uns = {}
        if cleanup:
            uns_obs = {}
            for key in list(obs.keys()):
                if len(set(obs[key])) == 1:
                    uns_obs[f"{key}"] = obs[key][0]
                    del obs[key]
            if uns_obs:
                uns["loom-obs"] = uns_obs
            uns_var = {}
            for key in list(var.keys()):
                if len(set(var[key])) == 1:
                    uns_var[f"{key}"] = var[key][0]
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
            dtype=dtype,
        )
    return adata


def read_mtx(filename: PathLike, dtype: str = "float32") -> AnnData:
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
    return AnnData(X, dtype=dtype)


def read_text(
    filename: Union[PathLike, Iterator[str]],
    delimiter: Optional[str] = None,
    first_column_names: Optional[bool] = None,
    dtype: str = "float32",
) -> AnnData:
    """\
    Read `.txt`, `.tab`, `.data` (text) file.

    Same as :func:`~anndata.read_csv` but with default delimiter `None`.

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
    if not isinstance(filename, (PathLike, str, bytes)):
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


def iter_lines(file_like: Iterable[str]) -> Generator[str, None, None]:
    """Helper for iterating only nonempty lines without line breaks"""
    for line in file_like:
        line = line.rstrip("\r\n")
        if line:
            yield line


def _read_text(
    f: Iterator[str],
    delimiter: Optional[str],
    first_column_names: Optional[bool],
    dtype: str,
) -> AnnData:
    comments = []
    data = []
    lines = iter_lines(f)
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
                raise ValueError(f"Did not find delimiter {delimiter!r} in first line.")
            line_list = line.split(delimiter)
            # the first column might be row names, so check the last
            if not is_float(line_list[-1]):
                col_names = line_list
                # logg.msg("    assuming first line in file stores column names", v=4)
            else:
                if not is_float(line_list[0]) or first_column_names:
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
    # transfrom to array, this takes a long time and a lot of memory
    # but it’s actually the same thing as np.genfromtxt does
    # - we don’t use the latter as it would involve another slicing step
    #   in the end, to separate row_names from float data, slicing takes
    #   a lot of memory and CPU time
    if data[0].size != data[-1].size:
        raise ValueError(
            f"Length of first line ({data[0].size}) is different "
            f"from length of last line ({data[-1].size})."
        )
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
        dtype=dtype,
    )


def load_sparse_csr(d, key="X"):
    from scipy.sparse.csr import csr_matrix

    key_csr = f"{key}_csr"
    d[key] = csr_matrix(
        (d[f"{key_csr}_data"], d[f"{key_csr}_indices"], d[f"{key_csr}_indptr"]),
        shape=d[f"{key_csr}_shape"],
    )
    del_sparse_matrix_keys(d, key_csr)
    return d


def del_sparse_matrix_keys(mapping, key_csr):
    del mapping[f"{key_csr}_data"]
    del mapping[f"{key_csr}_indices"]
    del mapping[f"{key_csr}_indptr"]
    del mapping[f"{key_csr}_shape"]
