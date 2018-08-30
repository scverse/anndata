import numpy as np
import anndata as ad
import os

X = np.array([[1, 2, 3],
              [4, 5, 6],
              [7, 8, 9]])
L = np.array([[10, 11, 12],
              [13, 14, 15],
              [16, 17, 18]])

def test_creation():
    adata = ad.AnnData(X=X, layers={'L':L.copy()})

    assert list(adata.layers.keys()) == ['L']
    assert (adata.layers['L'] == L).all()

def test_views():
    adata = ad.AnnData(X=X, layers={'L':L.copy()})
    adata_view = adata[1:, 1:]

    assert adata_view.layers.isview
    assert adata_view.layers._adata_ref == adata

    assert adata_view.layers.keys() == adata.layers.keys()
    assert (adata_view.layers['L'] == adata.layers['L'][1:, 1:]).all()

    adata.layers['S'] = X

    assert adata_view.layers.keys() == adata.layers.keys()
    assert (adata_view.layers['S'] == adata.layers['S'][1:, 1:]).all()

    adata_view.layers['T'] = X[1:, 1:]

    assert not adata_view.layers.isview
    assert not adata_view.isview

def test_readwrite():
    adata = ad.AnnData(X=X, layers={'L':L.copy()})
    adata.write('test.h5ad')
    adata_read = ad.read_h5ad('test.h5ad')

    assert adata.layers.keys() == adata_read.layers.keys()
    assert (adata.layers['L'] == adata_read.layers['L']).all()

    os.remove('test.h5ad')

def test_readwrite_loom():
    adata = ad.AnnData(X=X, layers={'L': L.copy()})
    adata.write_loom('test.loom')
    adata_read = ad.read('test.loom')

    assert adata.layers.keys() == adata_read.layers.keys()
    assert (adata.layers['L'] == adata_read.layers['L']).all()

    os.remove('test.loom')

def test_backed():
    #backed mode for layers isn't implemented, layers stay in memory
    pass
