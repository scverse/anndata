from importlib.util import find_spec
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd
from pandas.api.types import is_categorical
import pytest
from scipy.sparse import csr_matrix, issparse

import anndata as ad

from anndata.tests.helpers import gen_adata, asarray

HERE = Path(__file__).parent


# -------------------------------------------------------------------------------
# Some test data
# -------------------------------------------------------------------------------

X_sp = csr_matrix([
    [1, 0, 0],
    [3, 0, 0],
    [5, 6, 0],
    [0, 0, 0],
    [0, 0, 0]
])

X_list = [    # data matrix of shape n_obs x n_vars
    [1, 0],
    [3, 0],
    [5, 6],
]

obs_dict = dict(  # annotation of observations / rows
    row_names=['name1', 'name2', 'name3'],  # row annotation
    oanno1=['cat1', 'cat2', 'cat2'],        # categorical annotation
    oanno1b=['cat1', 'cat1', 'cat1'],       # categorical annotation with one category
    oanno1c=['cat1', 'cat1', np.nan],       # categorical annotation with a missing value
    oanno2=['o1', 'o2', 'o3'],              # string annotation
    oanno3=[2.1, 2.2, 2.3],                 # float annotation
)

var_dict = dict(  # annotation of variables / columns
    vanno1=[3.1, 3.2],
    vanno2=['cat1', 'cat1'],  # categorical annotation
)

uns_dict = dict(  # unstructured annotation
    oanno1_colors=['#000000', '#FFFFFF'],
    uns2=['some annotation'],
    uns3="another annotation",
    uns4={
        "a": 1,
        "b": [2, 3],
        "c": "4"
    }
)


@pytest.fixture(params=[{}, {"compression": "gzip"}])
def dataset_kwargs(request):
    return request.param

# -------------------------------------------------------------------------------
# The test functions
# -------------------------------------------------------------------------------


@pytest.mark.parametrize('typ', [np.array, csr_matrix])
def test_readwrite_h5ad(typ, dataset_kwargs, backing_h5ad):
    tmpdir = tempfile.TemporaryDirectory()
    tmpdirpth = Path(tmpdir.name)
    mid_pth = tmpdirpth / "mid.h5ad"

    X = typ(X_list)
    adata_src = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    assert not is_categorical(adata_src.obs['oanno1'])
    adata_src.raw = adata_src
    adata_src.write(backing_h5ad, **dataset_kwargs)

    adata_mid = ad.read(backing_h5ad)
    adata_mid.write(mid_pth, **dataset_kwargs)

    adata = ad.read_h5ad(mid_pth)
    assert is_categorical(adata.obs['oanno1'])
    assert not is_categorical(adata.obs['oanno2'])
    assert adata.obs.index.tolist() == ['name1', 'name2', 'name3']
    assert adata.obs['oanno1'].cat.categories.tolist() == ['cat1', 'cat2']
    assert is_categorical(adata.raw.var['vanno2'])
    assert np.all(adata.obs == adata_src.obs)
    assert np.all(adata.var == adata_src.var)
    assert np.all(adata.var.index == adata_src.var.index)
    assert adata.var.index.dtype == adata_src.var.index.dtype
    assert type(adata.raw.X) == type(adata_src.raw.X)
    assert np.allclose(asarray(adata.raw.X), asarray(adata_src.raw.X))
    assert np.all(adata.raw.var == adata_src.raw.var)


@pytest.mark.parametrize('typ', [np.array, csr_matrix])
def test_readwrite_maintain_X_dtype(typ, backing_h5ad):
    X = typ(X_list)
    adata_src = ad.AnnData(X, dtype='int8')
    adata_src.write(backing_h5ad)

    adata = ad.read(backing_h5ad)
    assert adata.X.dtype == adata_src.X.dtype

def test_read_write_maintain_obsmvarm_dtypes(backing_h5ad):
    M, N = 100, 101
    orig = gen_adata((M, N))

    orig.write(backing_h5ad)
    curr = ad.read(backing_h5ad)

    assert type(orig.obsm["array"]) is type(curr.obsm["array"])
    assert np.all(orig.obsm["array"] == curr.obsm["array"])
    assert np.all(orig.varm["array"] == curr.varm["array"])
    assert type(orig.obsm["sparse"]) is type(curr.obsm["sparse"])
    assert np.all((orig.obsm["sparse"] == curr.obsm["sparse"]).toarray())
    assert np.all((orig.varm["sparse"] == curr.varm["sparse"]).toarray())
    assert type(orig.obsm["df"]) is type(curr.obsm["df"])
    assert np.all(orig.obsm["df"] == curr.obsm["df"])
    assert np.all(orig.varm["df"] == curr.varm["df"])


def test_maintain_layers(backing_h5ad):
    M, N = 100, 101
    orig = gen_adata((M, N))
    orig.write(backing_h5ad)
    curr = ad.read(backing_h5ad)

    assert type(orig.layers["array"]) is type(curr.layers["array"])
    assert np.all(orig.layers["array"] == curr.layers["array"])
    assert type(orig.layers["sparse"]) is type(curr.layers["sparse"])
    assert np.all((orig.layers["sparse"] == curr.layers["sparse"]).toarray())
    assert type(orig.layers["df"]) is type(curr.layers["df"])
    assert np.all(orig.layers["df"] == curr.layers["df"])


def test_readwrite_sparse_as_dense(backing_h5ad):
    adata_src = ad.AnnData(X_sp)
    adata_src.write(backing_h5ad, force_dense=True)

    adata = ad.read(backing_h5ad, chunk_size=2)
    assert issparse(adata.X)
    assert np.allclose(X_sp.toarray(), adata.X.toarray())


@pytest.mark.parametrize('typ', [np.array, csr_matrix])
def test_readwrite_h5ad_one_dimension(typ, backing_h5ad):
    X = typ(X_list)
    adata_src = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    adata_one = adata_src[:, 0].copy()
    adata_one.write(backing_h5ad)
    adata = ad.read(backing_h5ad)
    assert adata.shape == (3, 1)


@pytest.mark.parametrize('typ', [np.array, csr_matrix])
def test_readwrite_dynamic(typ, backing_h5ad):
    X = typ(X_list)
    adata_src = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    adata_src.filename = backing_h5ad  # change to backed mode
    adata_src.write()

    adata = ad.read(backing_h5ad)
    assert is_categorical(adata.obs['oanno1'])
    assert not is_categorical(adata.obs['oanno2'])
    assert adata.obs.index.tolist() == ['name1', 'name2', 'name3']
    assert adata.obs['oanno1'].cat.categories.tolist() == ['cat1', 'cat2']


@pytest.mark.skipif(not find_spec('zarr'), reason='Zarr is not installed')
@pytest.mark.parametrize('typ', [np.array, csr_matrix])
def test_readwrite_zarr(typ, tmp_path):
    X = typ(X_list)
    adata_src = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    assert not is_categorical(adata_src.obs['oanno1'])
    adata_src.write_zarr(tmp_path / 'test_zarr_dir', chunks=True)

    adata = ad.read_zarr(tmp_path / 'test_zarr_dir')
    assert is_categorical(adata.obs['oanno1'])
    assert not is_categorical(adata.obs['oanno2'])
    assert adata.obs.index.tolist() == ['name1', 'name2', 'name3']
    assert adata.obs['oanno1'].cat.categories.tolist() == ['cat1', 'cat2']


@pytest.mark.skipif(not find_spec('loompy'), reason='Loompy is not installed (expected on Python 3.5)')
@pytest.mark.parametrize('typ', [np.array, csr_matrix])
def test_readwrite_loom(typ, tmp_path):
    X = typ(X_list)
    adata_src = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    adata_src.obsm['X_a'] = np.zeros((adata_src.n_obs, 2))
    adata_src.varm['X_b'] = np.zeros((adata_src.n_vars, 3))
    adata_src.write_loom(tmp_path / 'test.loom', write_obsm_varm=True)

    adata = ad.read_loom(tmp_path / 'test.loom', sparse=typ is csr_matrix, cleanup=True)
    if isinstance(X, np.ndarray):
        assert np.allclose(adata.X, X)
    else:
        # TODO: this should not be necessary
        assert np.allclose(adata.X.toarray(), X.toarray())
    assert 'X_a' in adata.obsm_keys() and adata.obsm['X_a'].shape[1] == 2
    assert 'X_b' in adata.varm_keys() and adata.varm['X_b'].shape[1] == 3
    # as we called with `cleanup=True`
    assert 'oanno1b' in adata.uns['loom-obs']
    assert 'vanno2' in adata.uns['loom-var']


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
def test_write_csv(typ, tmp_path):
    X = typ(X_list)
    adata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict)
    adata.write_csvs(tmp_path / 'test_csv_dir', skip_data=False)


@pytest.mark.parametrize(['read', 'write', 'name'], [
    pytest.param(ad.read_h5ad, ad.readwrite.write._write_h5ad, 'test_empty.h5ad'),
    pytest.param(ad.read_loom, ad.readwrite.write_loom, 'test_empty.loom', marks=pytest.mark.xfail(reason='Loom can’t handle 0×0 matrices')),
    pytest.param(ad.read_zarr, ad.readwrite.write_zarr, 'test_empty.zarr'),
    pytest.param(ad.read_zarr, ad.readwrite.write_zarr, 'test_empty.zip', marks=pytest.mark.xfail(reason='Zarr zip storage doesn’t seem to work…')),
])
def test_readwrite_hdf5_empty(read, write, name, tmp_path):
    if read is ad.read_zarr:
        pytest.importorskip('zarr')
    adata = ad.AnnData(uns=dict(empty=np.array([], dtype=float)))
    write(tmp_path / name, adata)
    ad_read = read(tmp_path / name)
    assert ad_read.uns['empty'].shape == (0,)


def test_read_excel():
    adata = ad.read_excel(HERE / 'data/excel.xlsx', 'Sheet1', dtype=int)
    assert adata.X.tolist() == X_list


def test_write_categorical(tmp_path):
    adata_pth = tmp_path / "adata.h5ad"
    orig = ad.AnnData(
        X=np.ones((5, 5)),
        obs=pd.DataFrame({
            "cat1": ["a", "a", "b", np.nan, np.nan],
            "cat2": pd.Categorical(["a", "a", "b", np.nan, np.nan])
        })
    )
    orig.write_h5ad(adata_pth)
    curr = ad.read_h5ad(adata_pth)
    assert np.all(orig.obs.notna() == curr.obs.notna())
    assert np.all(orig.obs.stack().dropna() == curr.obs.stack().dropna())
