import numpy as np
import dask.array as da
import anndata as ad
import tempfile
import os
from anndata._core.views import DaskArrayView


def test_dask_array_view_to_dask_array():
    """Test the to_dask_array method - covers views.py line."""
    X = np.random.random((100, 50))
    adata = ad.AnnData(X=da.from_array(X, chunks=(10, 10)))
    subset = adata[:10, :]  # Creates DaskArrayView
    
    # Test the new method directly
    converted = subset.X.to_dask_array()
    
    assert hasattr(converted, 'compute')
    assert converted.shape == (10, 50)
    assert not isinstance(converted, DaskArrayView)


def test_dask_array_view_zarr_write():
    """Test DaskArrayView Zarr write - covers methods.py lines."""
    X = np.random.random((100, 50))
    adata = ad.AnnData(X=da.from_array(X, chunks=(10, 10)))
    subset = adata[:10, :]  # Creates DaskArrayView
    
    with tempfile.TemporaryDirectory() as tmpdir:
        zarr_path = os.path.join(tmpdir, "test.zarr")
        subset.write_zarr(zarr_path)
        
        # Verify it worked
        loaded = ad.read_zarr(zarr_path)
        assert loaded.X.shape == (10, 50)


def test_dask_array_view_h5ad_write():
    """Test DaskArrayView H5AD write - covers methods.py lines."""
    X = np.random.random((100, 50))
    adata = ad.AnnData(X=da.from_array(X, chunks=(10, 10)))
    subset = adata[:10, :]  # Creates DaskArrayView
    
    with tempfile.TemporaryDirectory() as tmpdir:
        h5ad_path = os.path.join(tmpdir, "test.h5ad")
        subset.write_h5ad(h5ad_path)
        
        # Verify it worked
        loaded = ad.read_h5ad(h5ad_path)
        assert loaded.X.shape == (10, 50)


def test_dask_array_view_layers_write():
    """Test DaskArrayView in layers - covers additional method paths."""
    X = np.random.random((50, 30))
    layer_data = np.random.random((50, 30))
    
    adata = ad.AnnData(
        X=da.from_array(X, chunks=(10, 10)),
        layers={"test": da.from_array(layer_data, chunks=(10, 10))}
    )
    subset = adata[:20, :]  # Creates DaskArrayView in both X and layers
    
    with tempfile.TemporaryDirectory() as tmpdir:
        zarr_path = os.path.join(tmpdir, "test_layers.zarr")
        subset.write_zarr(zarr_path)  # Tests layers with DaskArrayView
        
        loaded = ad.read_zarr(zarr_path)
        assert loaded.X.shape == (20, 30)
        assert "test" in loaded.layers