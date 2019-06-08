# Had to rename the test module to h5py_ so that it wouldn't conflict with the
# h5py import upon testing.

import numpy as np
import scipy.sparse as ss

from anndata import h5py


def test_create_and_read_dataset(tmp_path):
    h5_path = tmp_path / 'test.h5'
    sparse_matrix = ss.csr_matrix([
        [0, 1, 0],
        [0, 0, 1],
        [0, 0, 0],
        [1, 1, 0],
    ], dtype=np.float64)
    with h5py.File(h5_path) as h5f:
        h5f.create_dataset('sparse/matrix', data=sparse_matrix)
    with h5py.File(h5_path) as h5f:
        assert (h5f['sparse']['matrix'][1:3] != sparse_matrix[1:3]).size == 0
        assert (h5f['sparse']['matrix'][2:] != sparse_matrix[2:]).size == 0
        assert (h5f['sparse']['matrix'][:2] != sparse_matrix[:2]).size == 0
        assert (h5f['sparse']['matrix'][-2:] != sparse_matrix[-2:]).size == 0
        assert (h5f['sparse']['matrix'][:-2] != sparse_matrix[:-2]).size == 0
        assert (h5f['sparse']['matrix'].value != sparse_matrix).size == 0


def test_create_dataset_from_dataset(tmp_path):
    from_h5_path = tmp_path / 'from.h5'
    to_h5_path = tmp_path / 'to.h5'
    sparse_matrix = ss.csr_matrix([
        [0, 1, 0],
        [0, 0, 1],
        [0, 0, 0],
        [1, 1, 0],
    ], dtype=np.float64)
    with h5py.File(from_h5_path) as from_h5f:
        from_dset = from_h5f.create_dataset('sparse/matrix', data=sparse_matrix)

        with h5py.File(to_h5_path) as to_h5f:
            to_h5f.create_dataset('sparse/matrix', data=from_dset)
            assert (to_h5f['sparse/matrix'].value != sparse_matrix).size == 0


def test_dataset_append(tmp_path):
    h5_path = tmp_path / 'test.h5'
    sparse_matrix = ss.csr_matrix([
        [0, 1, 0],
        [0, 0, 1],
        [0, 0, 0],
        [1, 1, 0],
    ], dtype=np.float64)
    to_append = ss.csr_matrix([
        [0, 1, 1],
        [1, 0, 0],
    ], dtype=np.float64)
    appended_matrix = ss.vstack((sparse_matrix, to_append))

    with h5py.File(h5_path) as h5f:
        h5f.create_dataset('matrix', data=sparse_matrix, chunks=(100000,))
        h5f['matrix'].append(to_append)
        assert (h5f['matrix'].value != appended_matrix).size == 0
