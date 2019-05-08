from pathlib import Path
from typing import Union, Optional
from typing import Iterable, Iterator, Generator
from collections import OrderedDict
import gzip
import bz2
import numpy as np
import warnings

from .. import AnnData
from .. import h5py
from ..compat import PathLike, fspath
from .utils import *


def read_csv(
    filename: Union[PathLike, Iterator[str]],
    delimiter: Optional[str]=',',
    first_column_names: Optional[bool]=None,
    dtype: str='float32',
) -> AnnData:
    """Read ``.csv`` file.

    Same as :func:`~anndata.read_text` but with default delimiter ``','``.

    Parameters
    ----------
    filename
        Data file.
    delimiter
        Delimiter that separates data within text file. If ``None``, will split at
        arbitrary number of white spaces, which is different from enforcing
        splitting at single white space ``' '``.
    first_column_names
        Assume the first column stores row names.
    dtype
        Numpy data type.
    """
    return read_text(filename, delimiter, first_column_names, dtype)


def read_excel(
    filename: PathLike,
    sheet: Union[str, int],
    dtype: str='float32',
) -> AnnData:
    """Read ``.xlsx`` (Excel) file.

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
    row = {'row_names': df.iloc[:, 0].values.astype(str)}
    col = {'col_names': np.array(df.columns[1:], dtype=str)}
    return AnnData(X, row, col, dtype=dtype)


def read_umi_tools(filename: PathLike, dtype: str='float32') -> AnnData:
    """Read a gzipped condensed count matrix from umi_tools.

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
    header = fh.readline()  # read the first line

    for line in fh:
        t = line.decode('ascii').split('\t')  # gzip read bytes, hence the decoding
        try:
            dod[t[1]].update({t[0]:int(t[2])})
        except KeyError:
            dod[t[1]] = {t[0]:int(t[2])}

    df = DataFrame.from_dict(dod, orient='index')  # build the matrix
    df.fillna(value=0., inplace=True)  # many NaN, replace with zeros
    return AnnData(np.array(df), {'obs_names': df.index}, {'var_names': df.columns}, dtype=dtype)


def read_hdf(filename: PathLike, key: str) -> AnnData:
    """Read ``.h5`` (hdf5) file.

    Note: Also looks for fields ``row_names`` and ``col_names``.

    Parameters
    ----------
    filename
        Filename of data file.
    key
        Name of dataset in the file.
    """
    with h5py.File(filename, 'r') as f:
        # the following is necessary in Python 3, because only
        # a view and not a list is returned
        keys = [k for k in f.keys()]
        if key == '':
            raise ValueError((
                'The file {} stores the following sheets:\n{}\n'
                'Call read/read_hdf5 with one of them.'
            ).format(filename, keys))
        # read array
        X = f[key][()]
        # try to find row and column names
        rows_cols = [{}, {}]
        for iname, name in enumerate(['row_names', 'col_names']):
            if name in keys:
                rows_cols[iname][name] = f[name][()]
    adata = AnnData(X, rows_cols[0], rows_cols[1], dtype=X.dtype.name)
    return adata


def read_loom(filename: PathLike, sparse: bool = True, cleanup: bool = False, X_name: str = 'spliced',
              obs_names: str = 'CellID', var_names: str = 'Gene', dtype: str='float32', **kwargs) -> AnnData:
    """Read `.loom`-formatted hdf5 file.

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
        Whether to collapse all obs/var fields that only store one unique value into `.uns['loom-.']`.
    X_name
        Loompy key with which the data matrix `.X` is initialized.
    obs_names
        Loompy key where the observation/cell names are stored.
    var_names
        Loompy key where the variable/gene names are stored.
    **kwargs
        Arguments to loompy.connect
    """
    filename = fspath(filename)  # allow passing pathlib.Path objects
    from loompy import connect
    with connect(filename, 'r', **kwargs) as lc:

        if X_name not in lc.layers.keys(): X_name = ''
        X = lc.layers[X_name].sparse().T.tocsr() if sparse else lc.layers[X_name][()].T

        layers = OrderedDict()
        if X_name != '': layers['matrix'] = lc.layers[''].sparse().T.tocsr() if sparse else lc.layers[''][()].T
        for key in lc.layers.keys():
            if key != '': layers[key] = lc.layers[key].sparse().T.tocsr() if sparse else lc.layers[key][()].T

        obs = dict(lc.col_attrs)
        if obs_names in obs.keys(): obs['obs_names'] = obs.pop(obs_names)
        obsm_attrs = [k for k, v in obs.items() if v.ndim > 1 and v.shape[1] > 1]

        obsm = {}
        for key in obsm_attrs:
            obsm[key] = obs.pop(key)

        var = dict(lc.row_attrs)
        if var_names in var.keys(): var['var_names'] = var.pop(var_names)
        varm_attrs = [k for k, v in var.items() if v.ndim > 1 and v.shape[1] > 1]

        varm = {}
        for key in varm_attrs:
            varm[key] = var.pop(key)

        uns = {}
        if cleanup:
            uns_obs = {}
            for key in list(obs.keys()):
                if len(set(obs[key])) == 1:
                    uns_obs[f'{key}'] = obs[key][0]
                    del obs[key]
            if uns_obs:
                uns['loom-obs'] = uns_obs
            uns_var = {}
            for key in list(var.keys()):
                if len(set(var[key])) == 1:
                    uns_var[f'{key}'] = var[key][0]
                    del var[key]
            if uns_var:
                uns['loom-var'] = uns_var

        adata = AnnData(
            X,
            obs=obs,
            var=var,
            layers=layers,
            obsm=obsm if obsm else None,
            varm=varm if varm else None,
            uns=uns,
            dtype=dtype)
    return adata


def read_mtx(filename: PathLike, dtype: str='float32') -> AnnData:
    """Read ``.mtx`` file.

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
    delimiter: Optional[str]=None,
    first_column_names: Optional[bool]=None,
    dtype: str='float32',
) -> AnnData:
    """Read ``.txt``, ``.tab``, ``.data`` (text) file.

    Same as :func:`~anndata.read_csv` but with default delimiter ``None``.

    Parameters
    ----------
    filename
        Data file, filename or stream.
    delimiter
        Delimiter that separates data within text file. If ``None``, will split at
        arbitrary number of white spaces, which is different from enforcing
        splitting at single white space ``' '``.
    first_column_names
        Assume the first column stores row names.
    dtype
        Numpy data type.
    """
    if not isinstance(filename, (PathLike, str, bytes)):
        return _read_text(filename, delimiter, first_column_names, dtype)

    filename = Path(filename)
    if filename.suffix == '.gz':
        with gzip.open(str(filename), mode='rt') as f:
            return _read_text(f, delimiter, first_column_names, dtype)
    elif filename.suffix == '.bz2':
        with bz2.open(str(filename), mode='rt') as f:
            return _read_text(f, delimiter, first_column_names, dtype)
    else:
        with filename.open() as f:
            return _read_text(f, delimiter, first_column_names, dtype)


def iter_lines(file_like: Iterable[str]) -> Generator[str, None, None]:
    """ Helper for iterating only nonempty lines without line breaks"""
    for line in file_like:
        line = line.rstrip('\r\n')
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
        if line.startswith('#'):
            comment = line.lstrip('# ')
            if comment:
                comments.append(comment)
        else:
            if delimiter is not None and delimiter not in line:
                raise ValueError('Did not find delimiter "{}" in first line.'
                                 .format(delimiter))
            line_list = line.split(delimiter)
            # the first column might be row names, so check the last
            if not is_float(line_list[-1]):
                col_names = line_list
                # logg.msg('    assuming first line in file stores column names', v=4)
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
            # logg.msg('    assuming last comment line stores variable names', v=4)
            col_names = np.array(comments[-1].split())
        # just numbers as col_names
        else:
            # logg.msg('    did not find column names in file', v=4)
            col_names = np.arange(len(data[0])).astype(str)
    col_names = np.array(col_names, dtype=str)
    # read another line to check if first column contains row names or not
    if first_column_names is None: first_column_names = False
    for line in lines:
        line_list = line.split(delimiter)
        if first_column_names or not is_float(line_list[0]):
            # logg.msg('    assuming first column in file stores row names', v=4)
            first_column_names = True
            row_names.append(line_list[0])
            data.append(np.array(line_list[1:], dtype=dtype))
        else:
            data.append(np.array(line_list, dtype=dtype))
        break
    # if row names are just integers
    if len(data) > 1 and data[0].size != data[1].size:
        # logg.msg('    assuming first row stores column names and '
        #          'first column row names', v=4)
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
    # logg.msg('    read data into list of lists', t=True, v=4)
    # transfrom to array, this takes a long time and a lot of memory
    # but it's actually the same thing as np.genfromtxt does
    # - we don't use the latter as it would involve another slicing step
    #   in the end, to separate row_names from float data, slicing takes
    #   a lot of memory and CPU time
    if data[0].size != data[-1].size:
        raise ValueError(
            'length of first line ({}) is different from length of last line ({})'
            .format(data[0].size, data[-1].size))
    data = np.array(data, dtype=dtype)
    # logg.msg('    constructed array from list of list', t=True, v=4)
    # transform row_names
    if not row_names:
        row_names = np.arange(len(data)).astype(str)
        # logg.msg('    did not find row names in file', v=4)
    else:
        row_names = np.array(row_names)
        for iname, name in enumerate(row_names):
            row_names[iname] = name.strip('"')
    # adapt col_names if necessary
    if col_names.size > data.shape[1]:
        col_names = col_names[1:]
    for iname, name in enumerate(col_names):
        col_names[iname] = name.strip('"')
    return AnnData(data,
                   obs={'obs_names': row_names},
                   var={'var_names': col_names},
                   dtype=dtype)


def read_zarr(store):
    """Read from a hierarchical Zarr array store.

    Parameters
    ----------
    store
        The filename, a :class:`~typing.MutableMapping`, or a Zarr storage class.
    """
    if isinstance(store, Path):
        store = str(store)
    import zarr
    f = zarr.open(store, mode='r')
    d = {}
    for key in f.keys():
        _read_key_value_from_zarr(f, d, key)
    return AnnData(*AnnData._args_from_dict(d))


def _read_key_value_from_zarr(f, d, key, key_write=None):
    import zarr
    if key_write is None: key_write = key
    if isinstance(f[key], zarr.hierarchy.Group):
        d[key_write] = OrderedDict() if key == 'uns' else {}
        for k in f[key].keys():
            _read_key_value_from_zarr(f, d[key_write], key + '/' + k, k)
        return
    # the '()' means 'load everything into memory' (by contrast, ':'
    # only works if not reading a scalar type)
    if key != 'X':
        value = f[key][()]
    else:
        value = f[key] # don't load X into memory

    d[key_write] = value
    return


def read_h5ad(filename, backed: Optional[str] = None, chunk_size: int = 6000):
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
    if isinstance(backed, bool):
        # We pass `None`s through to h5py.File, and its default is “a”
        # (=“r+”, but create the file if it doesn’t exist)
        backed = 'r+' if backed else None
        warnings.warn(
            "In a future version, read_h5ad will no longer explicitly support "
            "boolean arguments. Specify the read mode, or leave `backed=None`.",
            DeprecationWarning,
        )
    if backed:
        # open in backed-mode
        return AnnData(filename=filename, filemode=backed)
    else:
        # load everything into memory
        constructor_args = _read_args_from_h5ad(filename=filename, chunk_size=chunk_size)
        X = constructor_args[0]
        dtype = None
        if X is not None:
            dtype = X.dtype.name  # maintain dtype, since 0.7
        return AnnData(*_read_args_from_h5ad(filename=filename, chunk_size=chunk_size), dtype=dtype)


def _read_args_from_h5ad(
    adata: AnnData = None,
    filename: Optional[PathLike] = None,
    mode: Optional[str] = None,
    chunk_size: int = 6000
):
    """Return a tuple with the parameters for initializing AnnData.

    Parameters
    ----------
    filename
        Defaults to the objects filename if ``None``.
    """
    if filename is None and (adata is None or adata.filename is None):
        raise ValueError('Need either a filename or an AnnData object with file backing')

    # we need to be able to call the function without reference to self when
    # not reading in backed mode
    backed = mode is not None
    if filename is None and not backed:
        filename = adata.filename

    d = {}
    if backed:
        f = adata.file._file
    else:
        f = h5py.File(filename, 'r')
    for key in f.keys():
        if backed and key in AnnData._BACKED_ATTRS:
            d[key] = None
        else:
            _read_key_value_from_h5(f, d, key, chunk_size=chunk_size)
    # backwards compat: save X with the correct name
    if 'X' not in d:
        if backed == 'r+':
            for key in AnnData._H5_ALIASES['X']:
                if key in d:
                    del f[key]
                    f.create_dataset('X', data=d[key])
                    break
    # backwards compat: store sparse matrices properly
    csr_keys = [key.replace('_csr_data', '')
                for key in d if '_csr_data' in key]
    for key in csr_keys:
        d = load_sparse_csr(d, key=key)
    if not backed:
        f.close()
    return AnnData._args_from_dict(d)


def _read_key_value_from_h5(f, d, key, key_write=None, chunk_size=6000):
    if key_write is None: key_write = key
    if isinstance(f[key], h5py.Group):
        d[key_write] = OrderedDict() if key == 'uns' else {}
        for k in f[key].keys():
            _read_key_value_from_h5(f, d[key_write], key + '/' + k, k, chunk_size)
        return

    ds = f[key]

    if isinstance(ds, h5py.Dataset) and 'sparse_format' in ds.attrs:
        value = h5py._load_h5_dataset_as_sparse(ds, chunk_size)
    elif isinstance(ds, h5py.Dataset):
        value = np.empty(ds.shape, ds.dtype)
        if 0 not in ds.shape:
            ds.read_direct(value)
    elif isinstance(ds, h5py.SparseDataset):
        value = ds.value
    else:
        value = ds[()]

    def postprocess_reading(key, value):
        # record arrays should stay record arrays and not become scalars
        if value.ndim == 1 and len(value) == 1 and value.dtype.names is None:
            value = value[0]
        if hasattr(value, 'dtype') and value.dtype.kind == 'S':
            value = value.astype(str)
        # transform byte strings in recarrays to unicode strings
        # TODO: come up with a better way of solving this, see also below
        if (key not in AnnData._H5_ALIASES['obs']
            and key not in AnnData._H5_ALIASES['var']
            and key != 'raw.var'
            and not isinstance(value, dict)
            and hasattr(value, 'dtype') and value.dtype.names is not None):
            new_dtype = [((dt[0], 'U{}'.format(int(int(dt[1][2:])/4)))
                          if dt[1][1] == 'S' else dt) for dt in value.dtype.descr]
            value = value.astype(new_dtype)
        return key, value

    key, value = postprocess_reading(key, value)
    d[key_write] = value
    return


def load_sparse_csr(d, key='X'):
    from scipy.sparse.csr import csr_matrix
    key_csr = key + '_csr'
    d[key] = csr_matrix((d[key_csr + '_data'],
                         d[key_csr + '_indices'],
                         d[key_csr + '_indptr']),
                        shape=d[key_csr + '_shape'])
    del_sparse_matrix_keys(d, key_csr)
    return d


def del_sparse_matrix_keys(mapping, key_csr):
    del mapping[key_csr + '_data']
    del mapping[key_csr + '_indices']
    del mapping[key_csr + '_indptr']
    del mapping[key_csr + '_shape']
