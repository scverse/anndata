import os
import warnings
from collections import Mapping
from pathlib import Path
from typing import Union
import pandas as pd
import numpy as np
from scipy.sparse import issparse
import logging as logg

from ..base import AnnData
from .. import h5py


def write_csvs(dirname: Union[Path, str], adata: AnnData, skip_data: bool = True, sep: str = ','):
    """See :meth:`~anndata.AnnData.write_csvs`.
    """
    dirname = Path(dirname)
    if dirname.suffix == '.csv':
        dirname = dirname.with_suffix('')
    # write the following at warning level, it's very important for the users
    # logg.info('writing \'.csv\' files to', dirname)
    if not dirname.is_dir():
        dirname.mkdir(parents=True, exist_ok=True)
    dir_uns = dirname / 'uns'
    if not dir_uns.is_dir():
        dir_uns.mkdir(parents=True, exist_ok=True)
    d = dict(
        obs=adata._obs,
        var=adata._var,
        obsm=adata._obsm.to_df(),
        varm=adata._varm.to_df(),
    )
    if not skip_data:
        d['X'] = pd.DataFrame(
            adata._X.toarray() if issparse(adata._X) else adata._X)
    d_write = {**d, **adata._uns}
    not_yet_raised_sparse_warning = True
    for key, value in d_write.items():
        if issparse(value):
            if not_yet_raised_sparse_warning:
                warnings.warn('Omitting to write sparse annotation.')
                not_yet_raised_sparse_warning = False
            continue
        filename = dirname
        if key not in {'X', 'var', 'obs', 'obsm', 'varm'}:
            filename = dir_uns
        filename /= '{}.csv'.format(key)
        df = value
        if not isinstance(value, pd.DataFrame):
            value = np.array(value)
            if np.ndim(value) == 0:
                value = value[None]
            try:
                df = pd.DataFrame(value)
            except Exception as e:
                warnings.warn('Omitting to write {!r}.'.format(key), type(e))
                continue
        df.to_csv(
            filename, sep=sep,
            header=True if key in {'obs', 'var', 'obsm', 'varm'} else False,
            index=True if key in {'obs', 'var'} else False,
        )


def write_loom(filename: Union[Path, str], adata: AnnData):
    filename = str(filename)  # allow passing Path object
    row_attrs = adata.var.to_dict('list')
    row_attrs['var_names'] = adata.var_names.values
    col_attrs = adata.obs.to_dict('list')
    col_attrs['obs_names'] = adata.obs_names.values
    X = adata.X.T
    if issparse(X):
        logg.info(
            '... writing to \'.loom\' file densifies sparse matrix')
        X = X.tocoo()
    from loompy import create
    if os.path.exists(filename):
        os.remove(filename)
    create(filename, X, row_attrs=row_attrs, col_attrs=col_attrs)


def _write_h5ad(filename: Union[Path, str], adata: AnnData, **kwargs):
    filename = str(filename)  # allow passing pathlib.Path objects
    if not filename.endswith(('.h5', '.h5ad')):
        raise ValueError('Filename needs to end with \'.h5ad\'.')
    if adata.isbacked:
        # close so that we can reopen below
        adata.file.close()
    # create directory if it doesn't exist
    dirname = os.path.dirname(filename)
    if dirname != '' and not os.path.exists(dirname):
        os.makedirs(dirname)
    d = adata._to_dict_fixed_width_arrays()
    # we're writing to a different location than the backing file
    # - load the matrix into the memory...
    if adata.isbacked and filename != adata.filename:
        d['X'] = adata.X[:]
    # need to use 'a' if backed, otherwise we loose the backed objects
    with h5py.File(filename, 'a' if adata.isbacked else 'w') as f:
        for key, value in d.items():
            _write_key_value_to_h5(f, key, value, **kwargs)
    if adata.isbacked:
        adata.file.open(filename, 'r+')


def _write_key_value_to_h5(f, key, value, **kwargs):
    if isinstance(value, Mapping):
        for k, v in value.items():
            if not isinstance(k, str):
                warnings.warn('dict key {} transformed to str upon writing to h5,'
                              'using string keys is recommended'
                              .format(k))
            _write_key_value_to_h5(f, key + '/' + str(k), v, **kwargs)
        return

    def preprocess_writing(value):
        if value is None:
            return value
        elif issparse(value):
            return value
        elif isinstance(value, dict):
            # old hack for storing dicts, is never reached
            # in the current implementation, can be removed in the future
            value = np.array([str(value)])
        else:
            # make sure value is an array
            value = np.array(value)
            # hm, why that?
            if value.ndim == 0: value = np.array([value])
        # make sure string format is chosen correctly
        if value.dtype.kind == 'U': value = value.astype(np.string_)
        return value

    value = preprocess_writing(value)

    # for some reason, we need the following for writing string arrays
    if key in f.keys() and value is not None: del f[key]

    # ignore arrays with empty dtypes
    if value is None or not value.dtype.descr:
        return
    try:
        if key in set(f.keys()):
            is_valid_group = isinstance(f[key], h5py.Group) \
                and f[key].shape == value.shape \
                and f[key].dtype == value.dtype \
                and not isinstance(f[key], h5py.SparseDataset)
            if not is_valid_group and not isinstance(value, sparse.spmatrix):
                f[key][()] = value
                return
            else:
                del f[key]
        f.create_dataset(key, data=value, **kwargs)
    except TypeError:
        # try writing as byte strings
        try:
            if value.dtype.names is None:
                if key in set(f.keys()):
                    if (f[key].shape == value.shape
                            and f[key].dtype == value.dtype):
                        f[key][()] = value.astype('S')
                        return
                    else:
                        del f[key]
                f.create_dataset(key, data=value.astype('S'), **kwargs)
            else:
                new_dtype = [(dt[0], 'S{}'.format(int(dt[1][2:])*4))
                             for dt in value.dtype.descr]
                if key in set(f.keys()):
                    if (f[key].shape == value.shape
                            and f[key].dtype == value.dtype):
                        f[key][()] = value.astype(new_dtype)
                        return
                    else:
                        del f[key]
                f.create_dataset(
                    key, data=value.astype(new_dtype), **kwargs)
        except Exception as e:
            warnings.warn('Could not save field with key = "{}" '
                          'to hdf5 file.'.format(key))
