import os
import warnings
import pandas as pd
import numpy as np
from scipy.sparse import issparse


def write_csvs(dirname, adata, skip_data=True, sep=','):
    """See :class:`~anndata.AnnData.write_csvs`.
    """
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
            value = np.array(value)
            if np.ndim(value) == 0:
                value = value[None]
            try:
                df = pd.DataFrame(value)
            except:
                warnings.warn('Omitting to write \'{}\'.'.format(key))
                continue
        df.to_csv(filename, sep=sep,
                  header=True if key in {'obs', 'var', 'obsm', 'varm'} else False,
                  index=True if key in {'obs', 'var'} else False)


def write_loom(filename, adata):
    from loompy import create
    row_attrs = adata.var.to_dict('list')
    row_attrs['var_names'] = adata.var_names.values
    col_attrs = adata.obs.to_dict('list')
    col_attrs['obs_names'] = adata.obs_names.values
    lc = create(
        filename,
        matrix=adata.X.T,
        row_attrs=row_attrs,
        col_attrs=col_attrs)
    lc.close()
