import os
import h5py
import warnings
import pandas as pd
import numpy as np
from scipy.sparse import issparse
from . import utils


def write_anndata(filename, adata,
                  compression='gzip', compression_opts=None):
    def preprocess_writing(value):
        if isinstance(value, dict):
            # hack for storing dicts
            value = np.array([str(value)])
        else:
            value = np.array(value)
            if value.ndim == 0: value = np.array([value])
        # make sure string format is chosen correctly
        if value.dtype.kind == 'U': value = value.astype(np.string_)
        return value
    filename = str(filename)  # allow passing pathlib.Path objects
    if not filename.endswith(('.h5', '.anndata')):
        raise ValueError('Filename needs to end with \'.anndata\'.')
    if not os.path.exists(os.path.dirname(filename)):
        # logg.info('creating directory', directory + '/', 'for saving output files')
        os.makedirs(directory)
    # output the following at warning level, it's very important for the users
    d_write = {}
    for key, value in adata._to_dict_fixed_width_arrays().items():
        if issparse(value):
            for k, v in utils.save_sparse_csr(value, key=key).items():
                d_write[k] = v
        else:
            d_write[key] = preprocess_writing(value)
            # some output about the data to write
            # print(type(value), value.dtype, value.dtype.kind, value.shape)
    with h5py.File(filename, 'w') as f:
        for key, value in d_write.items():
            try:
                # ignore arrays with empty dtypes
                if value.dtype.descr:
                    f.create_dataset(key, data=value,
                                     compression=compression,
                                     compression_opts=compression_opts)
            except TypeError:
                # try writing it as byte strings
                try:
                    if value.dtype.names is None:
                        f.create_dataset(key, data=value.astype('S'),
                                         compression=compression,
                                         compression_opts=compression_opts)
                    else:
                        new_dtype = [(dt[0], 'S{}'.format(int(dt[1][2:])*4))
                                     for dt in value.dtype.descr]
                        f.create_dataset(key, data=value.astype(new_dtype),
                                         compression=compression,
                                         compression_opts=compression_opts)
                except Exception as e:
                    # logg.info(str(e))
                    warnings.warn('Could not save field with key = "{}" to hdf5 file.'
                                  .format(key))
                    

def write_csvs(dirname, adata, skip_data=True):
    if dirname.endswith('.csv'):
        dirname = dirname.replace('.csv', '/')
    if not dirname.endswith('/'): dirname += '/'
    # write the following at warning level, it's very important for the users
    # logg.info('writing \'.csv\' files to', dirname)
    if not os.path.exists(dirname): os.makedirs(dirname)
    if not os.path.exists(dirname + 'uns'): os.makedirs(dirname + 'uns')
    d = {'obs': adata._obs,
         'var': adata._var,
         'obsm': adata._obsm.to_df(),
         'varm': adata._varm.to_df()}
    if not skip_data:
        d['X'] = pd.DataFrame(
            adata._X.toarray() if sp.issparse(adata._X) else adata._X)
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
            filename += 'uns/'
        filename += key + '.csv'
        df = value
        if not isinstance(value, pd.DataFrame):
            if np.ndim(value) == 0: value = value[None]
            try:
                df = pd.DataFrame(value)
            except:
                warnings.warn('Omitting to write \'{}\'.'.format(key))
                continue
        df.to_csv(filename,
                  header=True if key in {'obs', 'var', 'obsm', 'varm'} else False,
                  index=True if key in {'obs', 'var'} else False)


def write_loom(filename, adata):
    from loompy import create
    row_attrs = adata.var.to_dict('list')
    row_attrs['var_names'] = adata.var_names
    col_attrs = adata.obs.to_dict('list')
    col_attrs['obs_names'] = adata.obs_names
    lc = create(
        filename,
        matrix=adata.X.T,
        row_attrs=row_attrs,
        col_attrs=col_attrs)
    lc.close()
