from pathlib import Path
from typing import Union, Optional, Iterable, Generator, Iterator
from scipy.sparse import issparse
import numpy as np

from ..base import AnnData
from .. import h5py
from .utils import *


def read_csv(
    filename: Union[Path, str, Iterator[str]],
    delimiter: Optional[str]=',',
    first_column_names: Optional[bool]=None,
    dtype: str='float32',
) -> AnnData:
    """Read `.csv` file.

    Same as :func:`~anndata.read_text` but with default delimiter ','.

    Parameters
    ----------
    filename : `str` or `pathlib.Path` or `file`-like object
        Data file.
    delimiter : `str`, optional (default: `','`)
        Delimiter that separates data within text file. If `None`, will split at
        arbitrary number of white spaces, which is different from enforcing
        splitting at single white space ' '.
    first_column_names : `bool` or `None`, optional (default: `None`)
        Assume the first column stores row names.
    dtype : `str`, optional (default: 'float32')

    Returns
    -------
    An :class:`~scanpy.api.AnnData` object.
    """
    return read_text(filename, delimiter, first_column_names, dtype)


def read_excel(
    filename: Union[Path, str],
    sheet: Union[str, int],
) -> AnnData:
    """Read `.xlsx` (Excel) file.

    Assumes that the first columns stores the row names and the first row the
    column names.

    Parameters
    ----------
    filename : `str`
        File name to read from.
    sheet : `str`
        Name of sheet in Excel file.

    Returns
    -------
    An :class:`~anndata.AnnData` object.
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    # rely on pandas for reading an excel file
    from pandas import read_excel
    df = read_excel(filename, sheet)
    X = df.values[:, 1:].astype(float)
    row = {'row_names': df.iloc[:, 0].values.astype(str)}
    col = {'col_names': np.array(df.columns[1:], dtype=str)}
    return AnnData(X, row, col)


def read_umi_tools(filename: Union[Path, str]) -> AnnData:
    """Read a gzipped condensed count matrix from umi_tools.

    Parameters
    ----------
    filename : `str`
        File name to read from.

    Returns
    -------
    An :class:`~anndata.AnnData` object.
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    # import pandas for conversion of a dict of dicts into a matrix
    # import gzip to read a gzipped file :-)
    import gzip
    from pandas import DataFrame
    
    dod = {}  # this will contain basically everything
    fh = gzip.open(filename)
    header = fh.readline()  # read the first line
    
    for line in fh:
        t = line.decode('ascii').split('\t')  # gzip read bytes, hence the decoding
        try:
            dod[t[1]].update({t[0]:int(t[2])})
        except KeyError:
            dod[t[1]] = {t[0]:int(t[2])}
    
    df = DataFrame.from_dict(dod, orient='index')  # build the matrix
    df.fillna(value = 0., inplace=True)  # many NaN, replace with zeros
    return AnnData(np.array(df), {'obs_names': df.index}, {'var_names': df.columns})


def read_hdf(filename: Union[Path, str], key: str) -> AnnData:
    """Read `.h5` (hdf5) file.

    Note: Also looks for fields 'row_names' and 'col_names'.

    Parameters
    ----------
    filename : `str`
        Filename of data file.
    key : `str`
        Name of dataset in the file.

    Returns
    -------
    adata : :class:`~anndata.AnnData`
        An annotated data matrix.
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    with h5py.File(filename, 'r') as f:
        # the following is necessary in Python 3, because only
        # a view and not a list is returned
        keys = [k for k in f.keys()]
        if key == '':
            raise ValueError('The file ' + filename +
                             ' stores the following sheets:\n' + str(keys) +
                             '\n Call read/read_hdf5 with one of them.')
        # read array
        X = f[key][()]
        # try to find row and column names
        rows_cols = [{}, {}]
        for iname, name in enumerate(['row_names', 'col_names']):
            if name in keys:
                rows_cols[iname][name] = f[name][()]
    adata = AnnData(X, rows_cols[0], rows_cols[1])
    return adata


def read_loom(filename: Union[Path, str]) -> AnnData:
    """Read `.loom`-formatted hdf5 file.

    Parameters
    ----------
    filename : `str`
        The filename.

    Returns
    -------
    An :class:`~anndata.AnnData` object.
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    from loompy import connect
    lc = connect(filename, 'r')
    with h5py.File(filename, 'r') as f:
        X = f['matrix'][()]
    adata = AnnData(
        X.T,
        obs=dict(lc.col_attrs),  # not ideal: make the generator a dict...
        var=dict(lc.row_attrs))
    lc.close()
    return adata


def read_mtx(filename: Union[Path, str], dtype: str='float32') -> AnnData:
    """Read `.mtx` file.

    Returns
    -------
    An :class:`~anndata.AnnData` object.
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    from scipy.io import mmread
    # could be rewritten accounting for dtype to be more performant
    X = mmread(filename).astype(dtype)
    from scipy.sparse import csr_matrix
    X = csr_matrix(X)
    return AnnData(X)


def read_text(
    filename: Union[Path, str, Iterator[str]],
    delimiter: Optional[str]=None,
    first_column_names: Optional[bool]=None,
    dtype: str='float32',
) -> AnnData:
    """Read `.txt`, `.tab`, `.data` (text) file.

    Same as :func:`~anndata.read_csv` but with default delimiter `None`.

    Parameters
    ----------
    filename : `str`, `pathlib.Path`, or `file`-like object
        Data file.
    delimiter : `str`, optional (default: `None`)
        Delimiter that separates data within text file. If `None`, will split at
        arbitrary number of white spaces, which is different from enforcing
        splitting at single white space ' '.
    first_column_names : `bool` or `None`, optional (default: `None`)
        Assume the first column stores row names.
    dtype : `str`, optional (default: 'float32')

    Returns
    -------
    An :class:`~scanpy.api.AnnData` object.
    """
    if isinstance(filename, (str, Path)):
        with Path(filename).open() as f:
            return _read_text(f, delimiter, first_column_names, dtype)
    else:
        return _read_text(filename, delimiter, first_column_names, dtype)


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
                   var={'var_names': col_names})


def read_h5ad(filename, backed=False):
    """Read `.h5ad`-formatted hdf5 file.

    Parameters
    ----------
    filename : `str`
        File name of data file.
    backed : {`False`, `True`, 'r', 'r+'}, optional (default: `False`)
        Load :class:`~scanpy.api.AnnData` in `backed` mode instead of fully
        loading it into memory (`memory` mode). `True` and 'r' are
        equivalent. If you want to modify backed attributes of the AnnData
        object, you need to choose 'r+'.

    Returns
    -------
    An :class:`~anndata.AnnData` object.
    """

    if backed:
        # open in backed-mode
        return AnnData(filename=filename, filemode=backed)
    else:
        # load everything into memory
        d = _read_h5ad(filename=filename)
        return AnnData(d)


def _read_h5ad(adata=None, filename=None, mode=None):
    """Return a dict with arrays for initializing AnnData.

    Parameters
    ----------
    filename : `str` or `None`, optional (default: `None`)
        Defaults to the objects filename if `None`.
    """
    if filename is None and (adata is None or adata.filename is None):
        raise ValueError('Need either a filename or an AnnData object with file backing')

    # we need to be able to call the function without reference to self when
    # not reading in backed mode
    backed = False
    if filename is None:
        backed = True if mode is None else mode
        filename = adata.filename

    filename = str(filename)  # allow passing pathlib.Path objects
    d = {}
    if backed:
        f = adata.file._file
    else:
        f = h5py.File(filename, 'r')
    for key in f.keys():
        if backed and key in AnnData._BACKED_ATTRS:
            d[key] = None
        else:
            _read_key_value_from_h5(f, d, key)
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
    return d


def _read_key_value_from_h5(f, d, key, key_write=None):
    if key_write is None: key_write = key
    if isinstance(f[key], h5py.Group):
        d[key_write] = {}
        for k in f[key].keys():
            _read_key_value_from_h5(f, d[key_write], key + '/' + k, k)
        return
    # the '()' means 'load everything into memory' (by contrast, ':'
    # only works if not reading a scalar type)
    value = f[key][()]

    def postprocess_reading(key, value):
        if value.ndim == 1 and len(value) == 1:
            value = value[0]
        if value.dtype.kind == 'S':
            value = value.astype(str)
            # backwards compat:
            # recover a dictionary that has been stored as a string
            if len(value) > 0:
                if value[0] == '{' and value[-1] == '}': value = eval(value)
        # transform byte strings in recarrays to unicode strings
        # TODO: come up with a better way of solving this, see also below
        if (key not in AnnData._H5_ALIASES['obs']
            and key not in AnnData._H5_ALIASES['var']
            and key != 'raw.var'
            and not isinstance(value, dict) and value.dtype.names is not None):
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
