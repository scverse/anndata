from __future__ import annotations

from pathlib import Path
from os import PathLike, fspath
from types import MappingProxyType
from typing import Union, Optional, Literal, Mapping, Tuple
from typing import Iterable, Iterator, Generator
from collections import OrderedDict
import gzip
import bz2
from warnings import warn

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from .. import AnnData
from ..compat import _deprecate_positional_args
from .utils import is_float
from .h5ad import read_h5ad

try:
    from .zarr import read_zarr
except ImportError as _e:
    e = _e

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
    return AnnData(X, row, col)


def read_umi_tools(filename: PathLike, dtype=None) -> AnnData:
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
    adata = AnnData(X, rows_cols[0], rows_cols[1])
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
                "Received values for both `obsm_names` and `obsm_mapping`. This is "
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
                "Received values for both `varm_names` and `varm_mapping`. This is "
                "ambiguous, only pass `varm_mapping`."
            )
        varm_mapping = varm_names

    filename = fspath(filename)  # allow passing pathlib.Path objects
    from loompy import connect

    with connect(filename, "r", **kwargs) as lc:
        if X_name not in lc.layers.keys():
            X_name = ""
        X = lc.layers[X_name].sparse().T.tocsr() if sparse else lc.layers[X_name][()].T
        X = X.astype(dtype, copy=False)

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
    return AnnData(X)


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


# Reading 10x formats


def read_10x_h5(
    filename: Union[str, Path],
    *,
    genome: Optional[str] = None,
    feature_types: str | list[str] | None = None,
) -> AnnData:
    """\
    Read 10x-Genomics-formatted hdf5 file.

    Parameters
    ----------
    filename
        Path to a 10x hdf5 file.
    genome
        Filter expression to genes within this genome. For legacy 10x h5
        files, this must be provided if the data contains more than one genome.
    Feature types
        Which feature types to read in. Reads all in by default.

    Returns
    -------
    Annotated data matrix, where observations/cells are named by their
    barcode and variables/genes by gene name. Stores the following information:

    :attr:`~anndata.AnnData.X`
        The data matrix is stored
    :attr:`~anndata.AnnData.obs_names`
        Cell names
    :attr:`~anndata.AnnData.var_names`
        Gene names
    :attr:`~anndata.AnnData.var`\\ `['gene_ids']`
        Gene IDs
    :attr:`~anndata.AnnData.var`\\ `['feature_types']`
        Feature types
    """
    # Coerce feature_type to be either list of string or None
    if feature_types is not None:
        if isinstance(feature_types, str):
            feature_types = [feature_types]

    with h5py.File(str(filename), "r") as f:
        v3 = "/matrix" in f
    if v3:
        adata = _read_v3_10x_h5(filename)
        if genome:
            if genome not in adata.var["genome"].values:
                raise ValueError(
                    f"Could not find data corresponding to genome '{genome}' in '{filename}'. "
                    f'Available genomes are: {list(adata.var["genome"].unique())}.'
                )
            adata = adata[:, adata.var["genome"] == genome]
        if feature_types is not None:
            adata = adata[:, adata.var["feature_types"].isin(feature_types)]
        if adata.is_view:
            adata = adata.copy()
    else:
        adata = _read_legacy_10x_h5(filename, genome=genome)
    return adata


def _read_v3_10x_h5(filename):
    """
    Read hdf5 file from Cell Ranger v3 or later versions.
    """
    with h5py.File(str(filename), "r") as f:
        dsets = {}
        _collect_datasets(dsets, f["matrix"])

        from scipy.sparse import csr_matrix

        M, N = dsets["shape"]
        data = dsets["data"]
        if dsets["data"].dtype == np.dtype("int32"):
            data = dsets["data"].view("float32")
            data[:] = dsets["data"]
        matrix = csr_matrix(
            (data, dsets["indices"], dsets["indptr"]),
            shape=(N, M),
        )
        adata = AnnData(
            matrix,
            obs=dict(obs_names=dsets["barcodes"].astype(str)),
            var=dict(
                var_names=dsets["name"].astype(str),
                gene_ids=dsets["id"].astype(str),
                feature_types=dsets["feature_type"].astype(str),
                genome=dsets["genome"].astype(str),
            ),
        )
    return adata


def _read_legacy_10x_h5(filename, *, genome=None):
    """
    Read hdf5 file from Cell Ranger v2 or earlier versions.
    """
    with h5py.File(str(filename), "r") as f:
        children = list(f.keys())
        if not genome:
            if len(children) > 1:
                raise ValueError(
                    f"'{filename}' contains more than one genome. For legacy 10x h5 "
                    "files you must specify the genome if more than one is present. "
                    f"Available genomes are: {children}"
                )
            genome = children[0]
        elif genome not in children:
            raise ValueError(
                f"Could not find genome '{genome}' in '{filename}'. "
                f"Available genomes are: {children}"
            )

        dsets = {}
        _collect_datasets(dsets, f[genome])

        # AnnData works with csr matrices
        # 10x stores the transposed data, so we do the transposition right away
        from scipy.sparse import csr_matrix

        M, N = dsets["shape"]
        data = dsets["data"]
        if dsets["data"].dtype == np.dtype("int32"):
            data = dsets["data"].view("float32")
            data[:] = dsets["data"]
        matrix = csr_matrix(
            (data, dsets["indices"], dsets["indptr"]),
            shape=(N, M),
        )
        # the csc matrix is automatically the transposed csr matrix
        # as scanpy expects it, so, no need for a further transpostion
        adata = AnnData(
            matrix,
            obs=dict(obs_names=dsets["barcodes"].astype(str)),
            var=dict(
                var_names=dsets["gene_names"].astype(str),
                gene_ids=dsets["genes"].astype(str),
            ),
        )
        return adata


def _collect_datasets(dsets: dict, group: h5py.Group):
    for k, v in group.items():
        if isinstance(v, h5py.Dataset):
            dsets[k] = v[()]
        else:
            _collect_datasets(dsets, v)


def read_10x_mtx(
    path: Union[Path, str],
    *,
    var_names: Literal["gene_symbols", "gene_ids"] = "gene_symbols",
    make_unique: bool = True,
    gex_only: bool = True,
    prefix: str = None,
) -> AnnData:
    """\
    Read 10x-Genomics-formatted mtx directory.

    Parameters
    ----------
    path
        Path to directory for `.mtx` and `.tsv` files,
        e.g. './filtered_gene_bc_matrices/hg19/'.
    var_names
        The variables index.
    make_unique
        Whether to make the variables index unique by appending '-1',
        '-2' etc. or not.
    cache
        If `False`, read from source, if `True`, read from fast 'h5ad' cache.
    cache_compression
        See the h5py :ref:`dataset_compression`.
        (Default: `settings.cache_compression`)
    gex_only
        Only keep 'Gene Expression' data and ignore other feature types,
        e.g. 'Antibody Capture', 'CRISPR Guide Capture', or 'Custom'
    prefix
        Any prefix before `matrix.mtx`, `genes.tsv` and `barcodes.tsv`. For instance,
        if the files are named `patientA_matrix.mtx`, `patientA_genes.tsv` and
        `patientA_barcodes.tsv` the prefix is `patientA_`.
        (Default: no prefix)

    Returns
    -------
    An :class:`~anndata.AnnData` object
    """
    path = Path(path)
    prefix = "" if prefix is None else prefix
    genefile_exists = (path / f"{prefix}genes.tsv").is_file()
    read = _read_legacy_10x_mtx if genefile_exists else _read_v3_10x_mtx
    adata = read(
        str(path),
        var_names=var_names,
        make_unique=make_unique,
        prefix=prefix,
    )
    if genefile_exists or not gex_only:
        return adata
    else:
        gex_rows = list(
            map(lambda x: x == "Gene Expression", adata.var["feature_types"])
        )
        return adata[:, gex_rows].copy()


def _read_legacy_10x_mtx(
    path,
    var_names="gene_symbols",
    make_unique=True,
    *,
    prefix="",
):
    """
    Read mex from output from Cell Ranger v2 or earlier versions
    """
    from anndata.utils import make_index_unique

    path = Path(path)
    adata = read_mtx(
        path / f"{prefix}matrix.mtx",
    ).T  # transpose the data
    genes = pd.read_csv(path / f"{prefix}genes.tsv", header=None, sep="\t")
    if var_names == "gene_symbols":
        var_names = genes[1].values
        if make_unique:
            var_names = make_index_unique(pd.Index(var_names))
        adata.var_names = var_names
        adata.var["gene_ids"] = genes[0].values
    elif var_names == "gene_ids":
        adata.var_names = genes[0].values
        adata.var["gene_symbols"] = genes[1].values
    else:
        raise ValueError("`var_names` needs to be 'gene_symbols' or 'gene_ids'")
    adata.obs_names = pd.read_csv(path / f"{prefix}barcodes.tsv", header=None)[0].values
    return adata


def _read_v3_10x_mtx(
    path,
    *,
    var_names="gene_symbols",
    make_unique=True,
    prefix="",
):
    """
    Read mtx from output from Cell Ranger v3 or later versions
    """
    from anndata.utils import make_index_unique

    path = Path(path)
    adata = read_mtx(
        path / f"{prefix}matrix.mtx.gz",
    ).T  # transpose the data
    genes = pd.read_csv(path / f"{prefix}features.tsv.gz", header=None, sep="\t")
    if var_names == "gene_symbols":
        var_names = genes[1].values
        if make_unique:
            var_names = make_index_unique(pd.Index(var_names))
        adata.var_names = var_names
        adata.var["gene_ids"] = genes[0].values
    elif var_names == "gene_ids":
        adata.var_names = genes[0].values
        adata.var["gene_symbols"] = genes[1].values
    else:
        raise ValueError("`var_names` needs to be 'gene_symbols' or 'gene_ids'")
    adata.var["feature_types"] = genes[2].values
    adata.obs_names = pd.read_csv(path / f"{prefix}barcodes.tsv.gz", header=None)[
        0
    ].values
    return adata
