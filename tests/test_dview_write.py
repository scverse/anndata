import numpy as np
import anndata as ad
from anndata.compat import DaskArray

# Create test data
X = np.random.random((20, 10)).astype(np.float32)
adata = ad.AnnData(X=DaskArray.from_array(X, chunks=(5, 5)))

# Create subset (this creates DaskArrayView)  
subset = adata[:8, :6]
print(f"Subset X type: {type(subset.X)}")  # Should be DaskArrayView

# Try to write (this should work with the fix)
subset.write_zarr("test_output.zarr")
print("Write successful!")

# Read back and verify
loaded = ad.read_zarr("test_output.zarr")
print(f"Read successful! Shape: {loaded.shape}")