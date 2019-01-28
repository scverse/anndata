from importlib.util import find_spec
from pathlib import Path

import pytest
import numpy as np
import anndata as ad


X = np.array([
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9],
])
L = np.array([
    [10, 11, 12],
    [13, 14, 15],
    [16, 17, 18],
])


def test_creation():
    adata = ad.AnnData(X=X, layers={'L': L.copy()})

    assert list(adata.layers.keys()) == ['L']
    assert (adata.layers['L'] == L).all()


def test_views():
    adata = ad.AnnData(X=X, layers={'L': L.copy()})
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


def test_readwrite(backing_h5ad):
    adata = ad.AnnData(X=X, layers={'L': L.copy()})
    adata.write(backing_h5ad)
    adata_read = ad.read_h5ad(backing_h5ad)

    assert adata.layers.keys() == adata_read.layers.keys()
    assert (adata.layers['L'] == adata_read.layers['L']).all()


@pytest.mark.skipif(find_spec('loompy') is None, reason="loompy not installed")
def test_readwrite_loom():
    loom_path = Path('test.loom')
    if loom_path.is_file():
        loom_path.unlink()
    adata = ad.AnnData(X=X, layers={'L': L.copy()})
    adata.write_loom(loom_path)
    adata_read = ad.read_loom(loom_path, X_name='')

    assert adata.layers.keys() == adata_read.layers.keys()
    assert (adata.layers['L'] == adata_read.layers['L']).all()


def test_backed():
    # backed mode for layers isn't implemented, layers stay in memory
    pass
