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

#### Loading from File Paths

AnnDataset can load directly from h5ad file paths, which is useful for large datasets that don't fit in memory:

```python
# Load from file path (recommended for large backed datasets)
dataset = AnnDataset("/path/to/large_dataset.h5ad")

# The dataset will automatically handle backed data efficiently
dataloader = dataset.get_optimized_dataloader(batch_size=64)
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

For backed HDF5 data, AnnDataset provides **two methods** to create optimized DataLoaders with better disk I/O performance.

The optimized data loader provide significant performance improvements for backed data by sorting batch indices for efficient sequential disk access, while maintaining the original batch order in returned data.


**Method 1: `get_optimized_dataloader()`** (Recommended)
```python
optimized_dataloader = dataset.get_optimized_dataloader(
    batch_size=32,
    shuffle=True,
    num_workers=4
)
```

**Method 2: Manually pass `get_collate_fn()` to the DataLoader init**
```python
# Use standard DataLoader with optimized collate function
from torch.utils.data import DataLoader

# Alternative: pass collate function
dataloader = DataLoader(
    dataset,
    batch_size=32,
    shuffle=True,
    num_workers=4,
    collate_fn=dataset.get_collate_fn()
)
```

## Accessing Metadata

The dataset includes observation metadata in the returned items:

```python
# Item contains both expression data and metadata
item = dataset[0]
print(f"Available keys: {list(item.keys())}")

# Expression data
X = item["X"]

# Metadata (if available) - prefixed with 'obs_'
cell_type = item["obs_cell_type"]

# Different data types are handled automatically:
# - Numeric values: converted to torch.Tensor
# - String values: kept as strings
# - Other types: converted to strings
```

## Data Transformations

Transforms are applied to the complete data dictionary (including `X` and observation metadata) after it's loaded from the AnnData object. The `transform` parameter accepts only Transform class instances that inherit from {class}`~anndata.experimental.pytorch.Transform`. This unified API ensures all transforms work seamlessly with multiprocessing by design.

### Creating Custom Transform Classes

Transform classes receive the complete data dictionary containing expression data (`X`) and observation metadata and return the transformed dictionary:

The data dictionary has the following structure:
- `"X"`: Expression data as a PyTorch tensor (shape: `[n_vars]`)
- `"obs_{column_name}"`: Observation metadata from `adata.obs` columns
  - Numeric columns (int, float, bool) → PyTorch tensors
  - String columns → Python strings
  - Other types → converted to strings

For example, if your `adata.obs` has columns `["cell_type", "batch", "n_genes"]`, the data dictionary will contain:
```python
{
    "X": torch.tensor([...]),           # Expression data
    "obs_cell_type": "T_cell",          # String metadata
    "obs_batch": torch.tensor(1),       # Numeric metadata
    "obs_n_genes": torch.tensor(2500)   # Numeric metadata
}
```

Example:

```python
from anndata.experimental.pytorch import Transform, Compose
import torch

class NormalizeRowSum(Transform):
    def __init__(self, target_sum=1e4):
        self.target_sum = target_sum

    def __call__(self, data_dict):
        X = data_dict["X"]
        X = torch.clamp(X, min=0)  # Ensure positive values
        row_sum = torch.sum(X, dim=-1, keepdim=True) + 1e-8
        data_dict["X"] = X * (self.target_sum / row_sum)
        return data_dict

    def __repr__(self):
        return f"NormalizeRowSum(target_sum={self.target_sum})"

class LogTransform(Transform):
    def __call__(self, data_dict):
        data_dict["X"] = torch.log1p(data_dict["X"])
        return data_dict

class EncodeLabels(Transform):
    """Encode string labels to integers in the data dictionary."""
    def __init__(self, obs_key='cell_type', encoded_key=None, label_encoder=None):
        self.obs_key = obs_key
        self.encoded_key = encoded_key or f"obs_{obs_key}_encoded"
        self.label_encoder = label_encoder

    def __call__(self, data_dict):
        # Access metadata using obs_{column_name} key pattern
        obs_key = f"obs_{self.obs_key}"

        if obs_key in data_dict and self.label_encoder is not None:
            string_label = data_dict[obs_key]  # Get string label from metadata
            encoded_label = self.label_encoder.transform([string_label])[0]
            data_dict[self.encoded_key] = torch.tensor(encoded_label, dtype=torch.long)
        return data_dict

# Set up label encoder
from sklearn.preprocessing import LabelEncoder
label_encoder = LabelEncoder()
label_encoder.fit(adata.obs['cell_type'].values)

# Compose transforms
training_transform = Compose([
    NormalizeRowSum(target_sum=1e4),
    LogTransform(),
    EncodeLabels(obs_key='cell_type', label_encoder=label_encoder)
])

# Create dataset
dataset = AnnDataset(adata, transform=training_transform)

# Use optimized dataloader for backed data
dataloader = dataset.get_optimized_dataloader(batch_size=32, num_workers=4)

# Access both processed data and encoded labels
for batch in dataloader:
    X = batch["X"]  # Normalized and log-transformed expression data
    y = batch["obs_cell_type_encoded"]  # Encoded labels as integers
    # Your training code here
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
