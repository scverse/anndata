import h5py
import numpy as np
from ..base import AnnData
from .utils import *


def read_csv(filename, delimiter=',', first_column_names=None, dtype='float32'):
    """Read `.csv` file.

    Same as :func:`~anndata.read_txt` but with default delimiter ','.

    Parameters
    ----------
    filename : `str`
        Filename of data file.
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
    return read_txt(filename, delimiter, first_column_names, dtype)


def read_excel(filename, sheet):
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


def read_hdf(filename, key):
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


def read_loom(filename):
    """Read `.loom`-formatted hdf5 file.

    Parameters
    ----------
    filename : `str`
        The filename.

    Returns
    -------
    An :class:`~anndata.AnnData` object.
    """
    from loompy import connect
    lc = connect(filename, 'r')
    with h5py.File(filename, 'r') as f:
        X = f['matrix'][()]
    adata = AnnData(
        X.T,
        obs=lc.col_attrs,
        var=lc.row_attrs)
    lc.close()
    return adata


def read_mtx(filename, dtype='float32'):
    """Read `.mtx` file.

    Returns
    -------
    An :class:`~anndata.AnnData` object.
    """
    from scipy.io import mmread
    # could be rewritten accounting for dtype to be more performant
    X = mmread(filename).astype(dtype)
    from scipy.sparse import csr_matrix
    X = csr_matrix(X)
    return AnnData(X)


def read_text(filename, delimiter=None, first_column_names=None, dtype='float32'):
    """Read `.txt`, `.tab`, `.data` (text) file.

    Same as :func:`~anndata.read_csv` but with default delimiter `None`.

    Parameters
    ----------
    filename : `str`
        Filename of data file.
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
    filename = str(filename)  # allow passing pathlib.Path objects
    header = ''
    data = []
    length = -1
    f = open(filename)
    col_names = []
    row_names = []
    # read header and column names
    for line in f:
        if line.startswith('#'):
            header += line
        else:
            if delimiter is not None and delimiter not in line:
                raise ValueError('Did not find delimiter "{}" in first line.'
                                 .format(delimiter))
            line_list = line.split(delimiter)
            if not is_float(line_list[0]):
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
        if len(header) > 0:
            # logg.msg('    assuming last comment line stores variable names', v=4)
            col_names = np.array(header.split('\n')[-2].strip('#').split())
        # just numbers as col_names
        else:
            # logg.msg('    did not find column names in file', v=4)
            col_names = np.arange(len(data[0])).astype(str)
    col_names = np.array(col_names, dtype=str)
    # read another line to check if first column contains row names or not
    if first_column_names is None: first_column_names = False
    for line in f:
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
    for line in f:
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
            'length of first line {} is different from length of last line {}'
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


