from importlib.util import find_spec
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix
import anndata as ad


HERE = Path(__file__).parent


# -------------------------------------------------------------------------------
# Some test data
# -------------------------------------------------------------------------------

X_list = [    # data matrix of shape n_obs x n_vars
    [1, 0], [3, 0], [5, 6]]

obs_dict = {  # annotation of observations / rows
    'row_names': ['name1', 'name2', 'name3'],  # row annotation
    'oanno1': ['cat1', 'cat2', 'cat2'],        # categorical annotation
    'oanno1b': ['cat1', 'cat1', 'cat1'],       # categorical annotation with one category
    'oanno2': ['o1', 'o2', 'o3'],              # string annotation
    'oanno3': [2.1, 2.2, 2.3]}                 # float annotation

var_dict = {  # annotation of variables / columns
    'vanno1': [3.1, 3.2],
    'vanno2': ['cat1', 'cat1']}        # categorical annotation

uns_dict = {  # unstructured annotation
    'oanno1_colors': ['#000000', '#FFFFFF'],
    'uns2': ['some annotation']}


# -------------------------------------------------------------------------------
# The test functions
# -------------------------------------------------------------------------------


def test_readwrite_h5ad():
    for typ in [np.array, csr_matrix]:
        X = typ(X_list)
        adata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
        assert pd.api.types.is_string_dtype(adata.obs['oanno1'])
        adata.raw = adata
        adata.write('./test.h5ad')
        adata = ad.read('./test.h5ad')
        assert pd.api.types.is_categorical(adata.obs['oanno1'])
        assert pd.api.types.is_string_dtype(adata.obs['oanno2'])
        assert adata.obs.index.tolist() == ['name1', 'name2', 'name3']
        assert adata.obs['oanno1'].cat.categories.tolist() == ['cat1', 'cat2']
        assert pd.api.types.is_categorical(adata.raw.var['vanno2'])


def test_readwrite_dynamic():
    for typ in [np.array, csr_matrix]:
        X = typ(X_list)
        adata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
        adata.filename = './test.h5ad'  # change to backed mode
        adata.write()
        adata = ad.read('./test.h5ad')
        assert pd.api.types.is_categorical(adata.obs['oanno1'])
        assert pd.api.types.is_string_dtype(adata.obs['oanno2'])
        assert adata.obs.index.tolist() == ['name1', 'name2', 'name3']
        assert adata.obs['oanno1'].cat.categories.tolist() == ['cat1', 'cat2']


def test_readwrite_zarr():
    for typ in [np.array, csr_matrix]:
        X = typ(X_list)
        adata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
        assert pd.api.types.is_string_dtype(adata.obs['oanno1'])
        adata.write_zarr('./test_zarr_dir', chunks=True)
        adata = ad.read_zarr('./test_zarr_dir')
        assert pd.api.types.is_categorical(adata.obs['oanno1'])
        assert pd.api.types.is_string_dtype(adata.obs['oanno2'])
        assert adata.obs.index.tolist() == ['name1', 'name2', 'name3']
        assert adata.obs['oanno1'].cat.categories.tolist() == ['cat1', 'cat2']


@pytest.mark.skipif(not find_spec('loompy'), reason='Loompy is not installed (expected on Python 3.5)')
def test_readwrite_loom():
    for i, typ in enumerate([np.array, csr_matrix]):
        X = typ(X_list)
        adata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
        adata.write_loom('./test.loom')
        adata = ad.read_loom('./test.loom', sparse=(i == 1))
        if isinstance(X, np.ndarray):
            assert np.allclose(adata.X, X)
        else:
            # TODO: this should not be necessary
            assert np.allclose(adata.X.toarray(), X.toarray())


def test_read_csv():
    adata = ad.read_csv(HERE / 'adata.csv')
    assert adata.obs_names.tolist() == ['r1', 'r2', 'r3']
    assert adata.var_names.tolist() == ['c1', 'c2']
    assert adata.X.tolist() == X_list


def test_read_tsv_strpath():
    adata = ad.read_text(str(HERE / 'adata-comments.tsv'), '\t')
    assert adata.obs_names.tolist() == ['r1', 'r2', 'r3']
    assert adata.var_names.tolist() == ['c1', 'c2']
    assert adata.X.tolist() == X_list


def test_read_tsv_iter():
    with (HERE / 'adata-comments.tsv').open() as f:
        adata = ad.read_text(f, '\t')
        assert adata.obs_names.tolist() == ['r1', 'r2', 'r3']
        assert adata.var_names.tolist() == ['c1', 'c2']
        assert adata.X.tolist() == X_list


def test_write_csv():
    for typ in [np.array, csr_matrix]:
        X = typ(X_list)
        adata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
        adata.write_csvs('./test_csv_dir', skip_data=False)
