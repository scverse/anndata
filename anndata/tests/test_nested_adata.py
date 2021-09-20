"""Tests for storing adata in uns """

import tempfile
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csr_matrix, csc_matrix

import anndata as ad
from anndata.tests.helpers import gen_adata, assert_equal


@pytest.mark.parametrize("typ", [np.array, csr_matrix, csc_matrix])
def test_readwrite_nested_adata(typ):
    tmpdir = tempfile.TemporaryDirectory()
    tmpdirpth = Path(tmpdir.name)
    h5ad_pth = tmpdirpth / "adata.h5ad"
    M, N = 100, 101
    adata = gen_adata((M, N), X_type=typ)
    M_nested, N_nested = 10, 11
    nested_adata = gen_adata((M_nested, N_nested), X_type=typ)
    adata.uns['test'] = nested_adata

    adata.write_h5ad(h5ad_pth)
    from_h5ad = ad.read_h5ad(h5ad_pth)

    assert_equal(from_h5ad.uns['test'], nested_adata, exact=True)

@pytest.mark.parametrize("typ", [np.array, csr_matrix, csc_matrix])
def test_readwrite_nested_multilevel_adata(typ):
    tmpdir = tempfile.TemporaryDirectory()
    tmpdirpth = Path(tmpdir.name)
    h5ad_pth = tmpdirpth / "adata.h5ad"
    M, N = 100, 101
    adata = gen_adata((M, N), X_type=typ)

    nested_adata1 = gen_adata((10, 11), X_type=typ)
    nested_adata2 = gen_adata((40, 20), X_type=typ)
    nested_adata1.uns['test'] = nested_adata2
    adata.uns['test'] = nested_adata1

    adata.write_h5ad(h5ad_pth)
    from_h5ad = ad.read_h5ad(h5ad_pth)

    assert_equal(from_h5ad.uns['test'].uns['test'], nested_adata2, exact=True)
