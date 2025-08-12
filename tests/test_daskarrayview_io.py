# tests/test_daskarrayview_io.py
from __future__ import annotations

import tempfile
from pathlib import Path

import dask.array as da
import numpy as np

import anndata as ad
from anndata._core.views import DaskArrayView


def test_dask_array_view_to_dask_array():
    """Test the to_dask_array method - covers views.py line."""
    X = np.random.random((100, 50))
    adata = ad.AnnData(X=da.from_array(X, chunks=(10, 10)))
    subset = adata[:10, :]

    converted = subset.X.to_dask_array()

    assert hasattr(converted, "compute")
    assert converted.shape == (10, 50)
    assert not isinstance(converted, DaskArrayView)


def test_dask_array_view_zarr_write():
    """Test DaskArrayView Zarr write - covers methods.py lines."""
    X = np.random.random((100, 50))
    adata = ad.AnnData(X=da.from_array(X, chunks=(10, 10)))
    subset = adata[:10, :]

    with tempfile.TemporaryDirectory() as tmpdir:
        zarr_path = Path(tmpdir) / "test.zarr"
        subset.write_zarr(zarr_path)

        loaded = ad.read_zarr(zarr_path)
        assert loaded.X.shape == (10, 50)


def test_dask_array_view_h5ad_write():
    """Test DaskArrayView H5AD write - covers methods.py lines."""
    X = np.random.random((100, 50))
    adata = ad.AnnData(X=da.from_array(X, chunks=(10, 10)))
    subset = adata[:10, :]

    with tempfile.TemporaryDirectory() as tmpdir:
        h5ad_path = Path(tmpdir) / "test.h5ad"
        subset.write_h5ad(h5ad_path)

        loaded = ad.read_h5ad(h5ad_path)
        assert loaded.X.shape == (10, 50)


def test_dask_array_view_layers_write():
    """Test DaskArrayView in layers - covers additional method paths."""
    X = np.random.random((50, 30))
    layer_data = np.random.random((50, 30))

    adata = ad.AnnData(
        X=da.from_array(X, chunks=(10, 10)),
        layers={"test": da.from_array(layer_data, chunks=(10, 10))},
    )
    subset = adata[:20, :]

    with tempfile.TemporaryDirectory() as tmpdir:
        zarr_path = Path(tmpdir) / "test_layers.zarr"
        subset.write_zarr(zarr_path)

        loaded = ad.read_zarr(zarr_path)
        assert loaded.X.shape == (20, 30)
        assert "test" in loaded.layers
