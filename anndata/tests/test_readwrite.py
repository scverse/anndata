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


@pytest.mark.parametrize('typ', [np.array, csr_matrix])
def test_readwrite_h5ad(typ):
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


@pytest.mark.parametrize('typ', [np.array, csr_matrix])
def test_readwrite_h5ad_one_dimensino(typ):
    X = typ(X_list)
    adata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    adata = adata[:, 0].copy()
    adata.write('./test.h5ad')
    adata = ad.read('./test.h5ad')


@pytest.mark.parametrize('typ', [np.array, csr_matrix])
def test_readwrite_dynamic(typ):
    X = typ(X_list)
    adata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    adata.filename = './test.h5ad'  # change to backed mode
    adata.write()
    adata = ad.read('./test.h5ad')
    assert pd.api.types.is_categorical(adata.obs['oanno1'])
    assert pd.api.types.is_string_dtype(adata.obs['oanno2'])
    assert adata.obs.index.tolist() == ['name1', 'name2', 'name3']
    assert adata.obs['oanno1'].cat.categories.tolist() == ['cat1', 'cat2']


@pytest.mark.skipif(not find_spec('zarr'), reason='Zarr is not installed')
@pytest.mark.parametrize('typ', [np.array, csr_matrix])
def test_readwrite_zarr(typ):
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
@pytest.mark.parametrize('typ', [np.array, csr_matrix])
def test_readwrite_loom(typ):
    X = typ(X_list)
    adata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    adata.write_loom('./test.loom')
    adata = ad.read_loom('./test.loom', sparse=typ is csr_matrix)
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


@pytest.mark.parametrize('typ', [np.array, csr_matrix])
def test_write_csv(typ):
    X = typ(X_list)
    adata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    adata.write_csvs('./test_csv_dir', skip_data=False)


@pytest.mark.parametrize(['read', 'write'], [
    (ad.read_h5ad, ad.AnnData.write_h5ad),
    # Loom can’t handle 0×0 matrices
    # (ad.read_loom, ad.AnnData.write_loom),
    # TODO: only run this if zarr is installed
    # (ad.read_zarr, lambda a, f: ad.readwrite.write_zarr(f, a)),
])
def test_readwrite_hdf5_empty(read, write):
    adata = ad.AnnData(uns=dict(empty=np.array([], dtype=float)))
    write(adata, HERE / 'test.empty.h5')
    ad_read = read(HERE / 'test.empty.h5')
    assert ad_read.uns['empty'].shape == (0,)
