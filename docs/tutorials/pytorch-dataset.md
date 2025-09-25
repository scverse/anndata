# PyTorch Dataset for AnnData

The {class}`~anndata.experimental.pytorch.AnnDataset` class provides a PyTorch {class}`~torch.utils.data.Dataset` implementation for {class}`~anndata.AnnData` objects. This enables seamless integration with PyTorch DataLoaders for machine-learning workflows.

```{warning}
The `AnnDataset` class is part of the experimental API and may change in future versions.
```

## Key Features

- **Memory efficiency**: Works with both in-memory AnnData objects and streams from backed HDF5 data without loading into memory
- **Data transforms**: Applies custom transformations to loaded data
- **Multiprocessing safety**: Safe HDF5 handling following h5py best practices
- **Optimized batch loading**: Sorts indices for efficient sequential disk access

## Basic Usage

### Creating a Dataset

```python
import anndata as ad
from anndata.experimental.pytorch import AnnDataset
from torch.utils.data import DataLoader
import scanpy as sc

# Load example adata
adata = sc.datasets.pbmc3k_processed()

# Create a basic dataset
dataset = AnnDataset(adata)

print(f"Dataset length: {len(dataset)}")
print(f"Dataset shape: {dataset.shape}")

# Get a single observation
item = dataset[0]
print(f"observation keys: {list(item.keys())}")
print(f"Expression data shape: {item['X'].shape}")
```

### Using with DataLoader

#### Standard DataLoader (works for all data)

```python
# Create a DataLoader for training
dataloader = DataLoader(
    dataset,
    batch_size=32,
    shuffle=True,
    num_workers=4
)

# Iterate through batches
for batch in dataloader:
    X_batch = batch["X"]  # Shape: (batch_size, n_vars)
    # Your training code here
    break
```

#### Optimized DataLoader for Backed Data

For backed HDF5 data, use the optimized dataloader for better disk I/O performance:

```python
# Optimized DataLoader - sorts batch indices for sequential disk access
dataloader = dataset.get_optimized_dataloader(
    batch_size=32,
    shuffle=True,
    num_workers=4
)

# Alternative: Use standard DataLoader with optimized collate function
dataloader = DataLoader(
    dataset,
    batch_size=32,
    shuffle=True,
    num_workers=4,
    collate_fn=dataset.get_collate_fn()  # Enables sorted index loading
)
```

The optimized methods provide significant performance improvements for backed data by sorting batch indices for efficient sequential disk access, while maintaining the original batch order in returned data.

## Data Transformations

Transforms are applied to the expression data (`X`) after it's loaded from the AnnData object. The `transform` parameter accepts only Transform class instances that inherit from {class}`~anndata.experimental.pytorch.Transform`. This unified API ensures all transforms work seamlessly with multiprocessing by design.

### Creating Custom Transform Classes

Transform classes receive the loaded expression data as a PyTorch tensor and return the transformed tensor:

Example:

```python
from anndata.experimental.pytorch import Transform, Compose
import torch

class NormalizeRowSum(Transform):
    def __init__(self, target_sum=1e4):
        self.target_sum = target_sum

    def __call__(self, x):
        x = torch.clamp(x, min=0)  # Ensure positive values
        row_sum = torch.sum(x, dim=-1, keepdim=True) + 1e-8
        return x * (self.target_sum / row_sum)

    def __repr__(self):
        return f"NormalizeRowSum(target_sum={self.target_sum})"

class LogTransform(Transform):
    def __call__(self, x):
        return torch.log1p(x)


# Compose transforms
training_transform = Compose([
    NormalizeRowSum(target_sum=1e4),
    LogTransform(),
])

# Create dataset
dataset = AnnDataset(adata, transform=training_transform)

# Use optimized dataloader for backed data
dataloader = dataset.get_optimized_dataloader(batch_size=32, num_workers=4)
```


## Working with Subsets

You can create datasets that use only a subset of observations (rows) from your AnnData object using the `obs_subset` parameter.

```python
# Create dataset with specific observations
import numpy as np

# Select first 100 observations
subset_indices = np.arange(100)
subset_dataset = AnnDataset(adata, obs_subset=subset_indices)

print(f"Original dataset: {len(adata)} observations")
print(f"Subset dataset: {len(subset_dataset)} observations")

# Use with train/test splits
from sklearn.model_selection import train_test_split
train_idx, test_idx = train_test_split(np.arange(len(adata)), test_size=0.2)

train_dataset = AnnDataset(adata, obs_subset=train_idx, transform=transform)
test_dataset = AnnDataset(adata, obs_subset=test_idx, transform=transform)
```

## Loading from File Paths

AnnDataset can load directly from h5ad file paths, which is useful for large datasets that don't fit in memory:

```python
# Load from file path (recommended for large backed datasets)
dataset = AnnDataset("/path/to/large_dataset.h5ad")

# The dataset will automatically handle backed data efficiently
dataloader = dataset.get_optimized_dataloader(batch_size=64)
```

## Advanced Configuration

### Chunk Size Parameter

For sparse datasets stored as dense, you can control memory usage with the `chunk_size` parameter:

```python
# Adjust chunk size for memory management
dataset = AnnDataset(
    adata,
    chunk_size=1000,  # Process 1000 rows at a time (default: 6000)
    transform=transform
)
```

### Dataset Properties

The dataset provides useful properties for introspection:

```python
print(f"Dataset shape: {dataset.shape}")  # (n_obs, n_vars)
print(f"Number of observations: {len(dataset)}")
print(f"Number of variables: {dataset.shape[1]}")
```


### Accessing Metadata

The dataset includes observation metadata in the returned items:

```python
# Item contains both expression data and metadata
item = dataset[0]
print(f"Available keys: {list(item.keys())}")

# Expression data
X = item["X"]

# Metadata (if available) - prefixed with 'obs_'
if "obs_cell_type" in item:
    cell_type = item["obs_cell_type"]

# Different data types are handled automatically:
# - Numeric values: converted to torch.Tensor
# - String values: kept as strings
# - Other types: converted to strings
```

## Multiprocessing Safety

AnnDataset is designed to work safely with PyTorch's multiprocessing DataLoader:

```python
# Safe to use num_workers > 0
dataloader = DataLoader(
    dataset,
    batch_size=32,
    shuffle=True,
    num_workers=4  # Works correctly with backed data
)

# Or with optimized dataloader
dataloader = dataset.get_optimized_dataloader(
    batch_size=32,
    shuffle=True,
    num_workers=4  # Multiprocessing-safe
)
```

**Key multiprocessing features:**
- **Transform classes**: Required for multiprocessing compatibility (functions/lambdas won't work)
- **File reopening**: Each worker process opens backed files independently following h5py best practices
- **State management**: Dataset state is properly serialized for worker processes

## Error Handling and Validation

The dataset includes comprehensive input validation:

```python
# Invalid inputs will raise helpful errors
try:
    # Transform must be a Transform class instance
    dataset = AnnDataset(adata, transform=lambda x: x)  # Will raise TypeError
except TypeError as e:
    print(f"Error: {e}")

try:
    # obs_subset indices must be valid
    dataset = AnnDataset(adata, obs_subset=[99999])  # Will raise ValueError
except ValueError as e:
    print(f"Error: {e}")

try:
    # chunk_size must be positive
    dataset = AnnDataset(adata, chunk_size=0)  # Will raise ValueError
except ValueError as e:
    print(f"Error: {e}")
```
