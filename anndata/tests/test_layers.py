from importlib.util import find_spec
from pathlib import Path

import pytest
import numpy as np
import anndata as ad


X = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
L = np.array([[10, 11, 12], [13, 14, 15], [16, 17, 18]])


def test_creation():
    adata = ad.AnnData(X=X, layers={'L': L.copy()})

    assert list(adata.layers.keys()) == ['L']
    assert "L" in adata.layers
    assert "X" not in adata.layers
    assert "some_other_thing" not in adata.layers
    assert (adata.layers['L'] == L).all()


def test_views():
    adata = ad.AnnData(X=X, layers={'L': L.copy()})
    adata_view = adata[1:, 1:]

    assert adata_view.layers.isview
    assert adata_view.layers.parent_mapping == adata.layers

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
def test_readwrite_loom(tmp_path):
    loom_path = Path(tmp_path / 'test.loom')
    adata = ad.AnnData(X=X, layers={'L': L.copy()})
    adata.write_loom(loom_path)
    adata_read = ad.read_loom(loom_path, X_name='')

    assert adata.layers.keys() == adata_read.layers.keys()
    assert (adata.layers['L'] == adata_read.layers['L']).all()


def test_backed():
    # backed mode for layers isn't implemented, layers stay in memory
    pass


def test_copy():
    adata = ad.AnnData(X=X, layers={'L': L.copy()})
    bdata = adata.copy()
    adata.layers['L'] += 10
    assert np.all(adata.layers['L'] != bdata.layers['L'])  # 201
