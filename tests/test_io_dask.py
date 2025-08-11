"""
Clean test demonstrating the DaskArrayView write fix for GitHub issue #2022.
"""

import numpy as np
import pandas as pd  # ← Add this import
import tempfile
from pathlib import Path
import anndata as ad
import dask.array as da


def test_daskarray_view_write_fix():
    """Test that DaskArrayView objects can now be written to storage formats."""
    print("Testing DaskArrayView write fix for GitHub issue #2022")
    print("=" * 60)
    
    # Create test data
    X = np.random.random((20, 10)).astype(np.float32)
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(20)])
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(10)])
    
    # Creating AnnData with DaskArray
    dask_array = da.from_array(X, chunks=(5, 5))
    adata = ad.AnnData(X=dask_array, obs=obs, var=var)
    
    print(f"Created AnnData with DaskArray: {type(adata.X)}")
    
    # Create subset - this creates DaskArrayView
    subset = adata[:8, :6]
    print(f"Created subset: {type(subset.X)}")
    
    # Verify it's actually a DaskArrayView
    from anndata._core.views import DaskArrayView
    assert isinstance(subset.X, DaskArrayView), "Subset should create DaskArrayView"
    print(f"Confirmed DaskArrayView creation")
    
    # Test the new to_dask_array() method
    converted = subset.X.to_dask_array()
    from anndata.compat import DaskArray
    assert isinstance(converted, DaskArray), "Conversion should produce DaskArray"
    assert not isinstance(converted, DaskArrayView), "Should not be a view anymore"
    print(f"to_dask_array() method works correctly")
    
    # Verify data integrity
    np.testing.assert_array_equal(converted.compute(), X[:8, :6])
    print(f"Data integrity preserved during conversion")
    
    # Test writing to Zarr (this was failing before the fix)
    with tempfile.TemporaryDirectory() as tmpdir:
        zarr_path = Path(tmpdir) / "test.zarr"
        print(f"Writing DaskArrayView to Zarr: {zarr_path}")
        
        try:
            subset.write_zarr(str(zarr_path))
            print("Zarr write successful!")
            
            # Verify by reading back
            loaded = ad.read_zarr(str(zarr_path))
            np.testing.assert_array_equal(loaded.X, X[:8, :6])
            assert loaded.shape == (8, 6)
            print("Zarr roundtrip verified")
            
        except Exception as e:
            print(f"Zarr write failed: {e}")
            return False
    
    # Test writing to HDF5 (this was also failing before the fix)
    with tempfile.TemporaryDirectory() as tmpdir:
        h5_path = Path(tmpdir) / "test.h5ad"
        print(f"Writing DaskArrayView to HDF5: {h5_path}")
        
        try:
            subset.write_h5ad(str(h5_path))
            print("HDF5 write successful!")
            
            # Verify by reading back
            loaded = ad.read_h5ad(str(h5_path))
            np.testing.assert_array_equal(loaded.X, X[:8, :6])
            assert loaded.shape == (8, 6)
            print("HDF5 roundtrip verified")
            
        except Exception as e:
            print(f"HDF5 write failed: {e}")
            return False
    
    return True


def test_daskarray_view_in_layers():
    """Test that DaskArrayView objects in layers can also be written."""
    print("\nTesting DaskArrayView in layers")
    print("-" * 40)
    
    # Create test data
    X = np.random.random((15, 8)).astype(np.float32)
    layer_data = np.random.random((15, 8)).astype(np.float32)
    
    # Create AnnData with DaskArray in layers
    dask_X = da.from_array(X, chunks=(5, 4))
    dask_layer = da.from_array(layer_data, chunks=(5, 4))
    
    adata = ad.AnnData(
        X=dask_X,
        layers={"test_layer": dask_layer}
    )
    
    # Create subset - this creates DaskArrayView in both X and layers
    subset = adata[:10, :5]
    
    from anndata._core.views import DaskArrayView
    assert isinstance(subset.X, DaskArrayView)
    assert isinstance(subset.layers["test_layer"], DaskArrayView)
    print("Created DaskArrayView in both X and layers")
    
    # Test writing
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "test_layers.zarr"
        
        try:
            subset.write_zarr(str(output_path))
            print("Write with DaskArrayView in layers successful!")
            
            # Verify
            loaded = ad.read_zarr(str(output_path))
            np.testing.assert_array_equal(loaded.X, X[:10, :5])
            np.testing.assert_array_equal(loaded.layers["test_layer"], layer_data[:10, :5])
            print("Layer data integrity verified")
            
            return True
            
        except Exception as e:
            print(f"Write with layers failed: {e}")
            return False


def main():
    """Run all tests to demonstrate the fix."""
    print("DaskArrayView Write Fix Demonstration")
    print("GitHub Issue #2022: No Write Method for DaskArrayView")
    print("=" * 70)
    
    try:
        # Test basic functionality
        result1 = test_daskarray_view_write_fix()
        
        # Test with layers
        result2 = test_daskarray_view_in_layers()
        
        print("\n" + "=" * 70)
        print("SUMMARY")
        print("=" * 70)
        
        if result1 and result2:
            print("ALL TESTS PASSED!")
            print("DaskArrayView objects can now be written to both Zarr and HDF5")
            print("Fix works for DaskArrayView in main data (X) and layers")
            print("Data integrity is preserved in all cases")
            print("GitHub issue #2022 is resolved")
        else:
            print("Some tests failed")
            
    except ImportError as e:
        print(f"Missing dependencies: {e}")
        print("Please install: pip install dask[array]")
    except Exception as e:
        print(f"Unexpected error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()