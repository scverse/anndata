import h5py
import numpy as np
import pytest
from scipy import sparse

import anndata as ad
from anndata.core.sparsedataset import SparseDataset
from anndata.tests.helpers import assert_equal


@pytest.mark.parametrize(
    ['sparse_format', 'append_method'],
    [
        pytest.param(sparse.csr_matrix, sparse.vstack),
        pytest.param(sparse.csc_matrix, sparse.hstack),
    ],
)
def test_dataset_append_memory(tmp_path, sparse_format, append_method):
    h5_path = tmp_path / 'test.h5'
    a = sparse_format(sparse.random(100, 100))
    b = sparse_format(sparse.random(100, 100))

    with h5py.File(h5_path, "a") as f:
        ad.readwrite.h5ad.write_attribute(f, "mtx", a)
        diskmtx = SparseDataset(f["mtx"])

        diskmtx.append(b)
        fromdisk = diskmtx.tomemory()

    frommem = append_method([a, b])

    assert_equal(fromdisk, frommem)


@pytest.mark.parametrize(
    ['sparse_format', 'append_method'],
    [
        pytest.param(sparse.csr_matrix, sparse.vstack),
        pytest.param(sparse.csc_matrix, sparse.hstack),
    ],
)
def test_dataset_append_disk(tmp_path, sparse_format, append_method):
    h5_path = tmp_path / 'test.h5'
    a = sparse_format(sparse.random(10, 10))
    b = sparse_format(sparse.random(10, 10))

    with h5py.File(h5_path, "a") as f:
        ad.readwrite.h5ad.write_attribute(f, "a", a)
        ad.readwrite.h5ad.write_attribute(f, "b", b)
        a_disk = SparseDataset(f["a"])
        b_disk = SparseDataset(f["b"])

        a_disk.append(b_disk)
        fromdisk = a_disk.tomemory()

    frommem = append_method([a, b])

    assert_equal(fromdisk, frommem)


def test_wrong_formats(tmp_path):
    h5_path = tmp_path / "base.h5"
    base = sparse.random(100, 100, format="csr")

    with h5py.File(h5_path, "a") as f:
        ad.readwrite.h5ad.write_attribute(f, "base", base)
        disk_mtx = SparseDataset(f["base"])
        pre_checks = disk_mtx.tomemory()

        with pytest.raises(ValueError):
            disk_mtx.append(sparse.random(100, 100, format="csc"))
        with pytest.raises(ValueError):
            disk_mtx.append(sparse.random(100, 100, format="coo"))
        with pytest.raises(NotImplementedError):
            disk_mtx.append(np.random.random((100, 100)))
        disk_dense = f.create_dataset(
            "dense", data=np.random.random((100, 100))
        )
        with pytest.raises(NotImplementedError):
            disk_mtx.append(disk_dense)

        post_checks = disk_mtx.tomemory()

    # Check nothing changed
    assert not np.any((pre_checks != post_checks).toarray())
