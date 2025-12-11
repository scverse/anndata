from __future__ import annotations

import warnings
from importlib.util import find_spec
from pathlib import Path

import h5py
import numpy as np
import pytest
import zarr
from scipy.sparse import csr_matrix

from anndata import AnnData, settings
from anndata._io.specs.registry import read_elem_partial
from anndata.io import read_elem, write_h5ad, write_zarr

X = np.array([[1.0, 0.0, 3.0], [4.0, 0.0, 6.0], [0.0, 8.0, 0.0]], dtype="float32")
X_check = np.array([[4.0, 0.0], [0.0, 8.0]], dtype="float32")

WRITER = dict(h5ad=write_h5ad, zarr=write_zarr)
READER = dict(h5ad=h5py.File, zarr=zarr.open)


@pytest.mark.parametrize("typ", [np.asarray, csr_matrix])
def test_read_partial_X(tmp_path, typ, diskfmt):
    adata = AnnData(X=typ(X))

    path = Path(tmp_path) / ("test_tp_X." + diskfmt)

    WRITER[diskfmt](path, adata)

    store = READER[diskfmt](path, mode="r")
    if diskfmt == "zarr":
        X_part = read_elem_partial(store["X"], indices=([1, 2], [0, 1]))
    else:
        # h5py doesn't allow fancy indexing across multiple dimensions
        X_part = read_elem_partial(store["X"], indices=([1, 2],))
        X_part = X_part[:, [0, 1]]
        store.close()

    assert np.all(X_check == X_part)


@pytest.mark.skipif(not find_spec("scanpy"), reason="Scanpy is not installed")
def test_read_partial_adata(tmp_path, diskfmt):
    import scanpy as sc

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", r"Moving element.*uns.*to.*obsp", FutureWarning
        )
        adata = sc.datasets.pbmc68k_reduced()

    path = Path(tmp_path) / ("test_rp." + diskfmt)

    # we’re not adding things to read_partial anymore, so it can only read non-nullable strings.
    # therefore force writing old format here
    with settings.override(allow_write_nullable_strings=False):
        WRITER[diskfmt](path, adata)

    storage = READER[diskfmt](path, mode="r")

    obs_idx = [1, 2]
    var_idx = [0, 3]
    adata_sbs = adata[obs_idx, var_idx]

    if diskfmt == "zarr":
        part = read_elem_partial(storage["X"], indices=(obs_idx, var_idx))
    else:
        # h5py doesn't allow fancy indexing across multiple dimensions
        part = read_elem_partial(storage["X"], indices=(obs_idx,))
        part = part[:, var_idx]
    assert np.all(part == adata_sbs.X)

    part = read_elem_partial(storage["obs"], indices=(obs_idx,))
    assert np.all(part.keys() == adata_sbs.obs.keys())
    assert np.all(part.index == adata_sbs.obs.index)

    part = read_elem_partial(storage["var"], indices=(var_idx,))
    assert np.all(part.keys() == adata_sbs.var.keys())
    assert np.all(part.index == adata_sbs.var.index)

    for key in storage["obsm"]:
        part = read_elem_partial(storage["obsm"][key], indices=(obs_idx,))
        assert np.all(part == adata_sbs.obsm[key])

    for key in storage["varm"]:
        part = read_elem_partial(storage["varm"][key], indices=(var_idx,))
        np.testing.assert_equal(part, adata_sbs.varm[key])

    for key in storage["obsp"]:
        part = read_elem_partial(storage["obsp"][key], indices=(obs_idx, obs_idx))
        part = part.toarray()
        assert np.all(part == adata_sbs.obsp[key])

    # check uns just in case
    np.testing.assert_equal(read_elem(storage["uns"]).keys(), adata.uns.keys())
