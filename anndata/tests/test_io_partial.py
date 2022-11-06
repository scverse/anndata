from anndata import AnnData
from anndata._io.specs.registry import read_elem_partial
from anndata._io import write_zarr
from scipy.sparse import csr_matrix
import numpy as np
import pytest
from pathlib import Path
import zarr

X = np.array([[1.0, 0.0, 3.0], [4.0, 0.0, 6.0], [0.0, 8.0, 0.0]], dtype="float32")
X_check = np.array([[4.0, 0.0], [0.0, 8.0]], dtype="float32")


@pytest.mark.parametrize("typ", [np.asarray, csr_matrix])
def test_zarr_read_partial(tmp_path, typ):

    adata = AnnData(X=typ(X))

    print(tmp_path)

    path = Path(tmp_path) / "test_zarr.zarr"

    write_zarr(path, adata)

    with zarr.open(path, mode="r") as store:
        X_part = read_elem_partial(store["X"], indices=([1, 2], [0, 1]))

    assert np.all(X_check == X_part)
