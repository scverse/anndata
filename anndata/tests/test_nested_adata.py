"""Tests for storing adata in uns """

import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from anndata.tests.helpers import gen_adata, assert_equal
from scipy.sparse import csr_matrix, csc_matrix


@pytest.mark.parametrize("typ", [np.array, csr_matrix, csc_matrix])
@pytest.mark.parametrize("file_format", ['h5ad', 'zarr'])
def test_readwrite_nested_adata(typ, file_format):
    tmpdir = tempfile.TemporaryDirectory()
    tmpdirpth = Path(tmpdir.name)
    output_pth = tmpdirpth / "adata.{}".format(file_format)
    M, N = 100, 101
    adata = gen_adata((M, N), X_type=typ)
    M_nested, N_nested = 10, 11
    nested_adata = gen_adata((M_nested, N_nested), X_type=typ)
    adata.uns['test'] = nested_adata
    if file_format == 'h5ad':
        adata.write_h5ad(output_pth)
        from_disk = ad.read_h5ad(output_pth)
    elif file_format == 'zarr':
        adata.write_zarr(output_pth)
        from_disk = ad.read_zarr(output_pth)

    assert_equal(from_disk.uns['test'], nested_adata, exact=True)


@pytest.mark.parametrize("typ", [np.array, csr_matrix, csc_matrix])
@pytest.mark.parametrize("file_format", ['h5ad', 'zarr'])
def test_readwrite_nested_multilevel_adata(typ, file_format):
    tmpdir = tempfile.TemporaryDirectory()
    tmpdirpth = Path(tmpdir.name)
    output_pth = tmpdirpth / "adata.{}".format(file_format)
    M, N = 100, 101
    adata = gen_adata((M, N), X_type=typ)

    nested_adata1 = gen_adata((10, 11), X_type=typ)
    nested_adata2 = gen_adata((40, 20), X_type=typ)
    nested_adata1.uns['test'] = nested_adata2
    adata.uns['test'] = nested_adata1

    if file_format == 'h5ad':
        adata.write_h5ad(output_pth)
        from_disk = ad.read_h5ad(output_pth)
    elif file_format == 'zarr':
        adata.write_zarr(output_pth)
        from_disk = ad.read_zarr(output_pth)

    assert_equal(from_disk.uns['test'].uns['test'], nested_adata2, exact=True)
