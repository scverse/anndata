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
    adata = ad.AnnData(X=X, layers_X={'L':L.copy()})

    assert list(adata.layers_X.keys()) == ['L']
    assert (adata.layers_X['L'] == L).all()

def test_views():
    adata = ad.AnnData(X=X, layers_X={'L':L.copy()})
    adata_view = adata[1:, 1:]

    assert adata_view.layers_X.isview
    assert adata_view.layers_X._adata_ref == adata

    assert adata_view.layers_X.keys() == adata.layers_X.keys()
    assert (adata_view.layers_X['L'] == adata.layers_X['L'][1:, 1:]).all()

    adata.layers_X['S'] = X

    assert adata_view.layers_X.keys() == adata.layers_X.keys()
    assert (adata_view.layers_X['S'] == adata.layers_X['S'][1:, 1:]).all()

    adata_view.layers_X['T'] = X[1:, 1:]

    assert not adata_view.layers_X.isview
    assert not adata_view.isview

def test_readwrite():
    adata = ad.AnnData(X=X, layers_X={'L':L.copy()})
    adata.write('test.h5ad')

    adata_read = ad.read_h5ad('test.h5ad')

    assert adata.layers_X.keys() == adata_read.layers_X.keys()
    assert (adata.layers_X['L'] == adata_read.layers_X['L']).all()

    os.remove('test.h5ad')

def test_backed():
    #backed mode for layers isn't implemented, layers stay in memory
    pass
