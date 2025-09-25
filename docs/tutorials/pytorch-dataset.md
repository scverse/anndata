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
dataloader = DataLoader(dataset, batch_size=32, num_workers=4)
```


## Working with Subsets

You can create datasets that use only a subset of observations (rows) from your AnnData object using the `obs_subset` parameter.


### Accessing Metadata

The dataset includes observation metadata in the returned items:

```python
# Item contains both expression data and metadata
item = dataset[0]
print(f"Available keys: {list(item.keys())}")

# Expression data
X = item["X"]

# Metadata (if available)
if "obs_cell_type" in item:
    cell_type = item["obs_cell_type"]
```
