# PyTorch Dataset for AnnData

**Author:** AnnData Development Team

The {class}`~anndata.experimental.pytorch.AnnDataset` class provides a PyTorch {class}`~torch.utils.data.Dataset` implementation for {class}`~anndata.AnnData` objects, enabling seamless integration with PyTorch DataLoaders and training loops.

```{warning}
The `AnnDataset` class is part of the experimental API and may change in future versions.
```

## Key Features

- **Generic transform interface**: Accepts any callable transform function following PyTorch conventions
- **Memory efficiency**: Streams from backed HDF5 data without loading into memory
- **Multiprocessing safety**: Safe HDF5 handling following h5py best practices
- **Configurable chunk processing**: User-defined chunk sizes for memory management
- **Optimized batch loading**: Sorts indices for efficient sequential disk access

## Basic Usage

### Creating a Dataset

```python
import anndata as ad
from anndata.experimental.pytorch import AnnDataset
from torch.utils.data import DataLoader
import scanpy as sc

# Load real single-cell data
adata = sc.datasets.pbmc3k_processed()

# Create a basic dataset
dataset = AnnDataset(adata)

print(f"Dataset length: {len(dataset)}")
print(f"Dataset shape: {dataset.shape}")

# Get a single sample
sample = dataset[0]
print(f"Sample keys: {list(sample.keys())}")
print(f"Expression data shape: {sample['X'].shape}")
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

The `transform` parameter accepts any callable that takes a {class}`torch.Tensor` and returns a {class}`torch.Tensor`, following standard PyTorch Dataset conventions.

### Simple Normalization

```python
import torch

def normalize_transform(X):
    """Normalize to 10,000 counts per cell and log-transform."""
    # Ensure we have positive values and no NaNs
    X = torch.clamp(X, min=0)

    # Normalize to 10,000 counts per cell (sum across genes for each cell)
    row_sum = torch.sum(X, dim=-1, keepdim=True) + 1e-8  # Sum across genes (last dimension)
    X = X * (1e4 / row_sum)

    # Log1p transformation
    X = torch.log1p(X)

    # Ensure no NaNs or infinities
    X = torch.nan_to_num(X, nan=0.0, posinf=10.0, neginf=0.0)

    return X

# Create dataset with transform
dataset = AnnDataset(adata, transform=normalize_transform)
```

### Training with Augmentation

```python
def training_transform(X):
    """Complete preprocessing and augmentation pipeline."""
    # Ensure we have positive values and no NaNs
    X = torch.clamp(X, min=0)

    # Normalize to 10k counts per cell (sum across genes for each cell)
    row_sum = torch.sum(X, dim=-1, keepdim=True) + 1e-8  # Sum across genes (last dimension)
    X = X * (1e4 / row_sum)

    # Log transformation
    X = torch.log1p(X)

    # Add small amount of noise for augmentation
    if torch.rand(1).item() < 0.5:  # 50% chance
        noise = torch.normal(0, 0.01, size=X.shape)  # Reduced noise
        X = X + noise

    # Ensure no NaNs or infinities
    X = torch.nan_to_num(X, nan=0.0, posinf=10.0, neginf=0.0)

    return X

training_dataset = AnnDataset(
    adata,
    transform=training_transform,
    multiprocessing_safe=True,
    chunk_size=1000
)
```

## Working with Subsets

### Index-based Subsetting

```python
# Create dataset with specific cell indices
train_indices = [0, 1, 2, 100, 101, 102]  # Your training indices
train_dataset = AnnDataset(adata, obs_subset=train_indices)
```

### Large Dataset Configuration

For large datasets, configure chunk size and multiprocessing settings:

```python
large_dataset = AnnDataset(
    adata,  # Or pass file path: "large_data.h5ad"
    transform=normalize_transform,
    multiprocessing_safe=True,
    chunk_size=2000,  # Larger chunks for better performance
)

# Use with multiple workers
dataloader = DataLoader(
    large_dataset,
    batch_size=64,
    num_workers=8,  # More workers for better throughput
    pin_memory=True  # If using GPU
)
```

## Multiprocessing and Performance

The `AnnDataset` follows h5py best practices for multiprocessing:

- Each worker process opens HDF5 files independently
- No complex retry mechanisms needed
- Safe concurrent access to backed data

```python
# This is safe and efficient
dataloader = DataLoader(
    dataset,
    batch_size=128,
    num_workers=4,  # Multiple workers work safely
    shuffle=True
)
```

### Multiprocessing with Custom Transforms

```{note}
When using custom transform functions with multiprocessing (`num_workers > 0`), define your transforms in a separate Python module rather than inline. Functions defined in Jupyter notebooks or the main script cannot be pickled for multiprocessing.

**Production pattern:**
```python
# In transforms.py
def training_transform(X):
    row_sum = torch.sum(X, dim=-1, keepdim=True) + 1e-8
    X = X * (1e4 / row_sum)
    return torch.log1p(X)

# In your main script
from transforms import training_transform
dataloader = DataLoader(dataset, num_workers=4)  # Now multiprocessing works!
```
```

## Integration with Training Loops

### Basic Training Loop

```python
import torch.nn as nn
import torch.optim as optim

# Create model, loss, optimizer
model = nn.Linear(dataset.shape[1], 10)  # n_vars -> 10 classes
criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(model.parameters())

# Training loop
model.train()
for epoch in range(num_epochs):
    for batch in dataloader:
        X_batch = batch["X"]

        # Forward pass
        outputs = model(X_batch)
        # Your loss computation here

        # Backward pass
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
```

### Accessing Metadata

The dataset includes observation metadata in the returned samples:

```python
# Sample contains both expression data and metadata
sample = dataset[0]
print(f"Available keys: {list(sample.keys())}")

# Expression data
X = sample["X"]

# Metadata (if available)
if "obs_cell_type" in sample:
    cell_type = sample["obs_cell_type"]
```

## Performance Tips

1. **Use appropriate chunk sizes**: Larger chunks (1000-2000) for better I/O performance
2. **Enable multiprocessing**: Set `multiprocessing_safe=True` for backed data
3. **Use multiple workers**: DataLoader with `num_workers > 0` for better throughput
4. **Pin memory**: Use `pin_memory=True` when transferring to GPU

```python
# Optimized configuration
dataset = AnnDataset(
    adata,
    transform=your_transform,
    multiprocessing_safe=True,
    chunk_size=2000
)

dataloader = DataLoader(
    dataset,
    batch_size=128,
    num_workers=4,
    pin_memory=True,
    shuffle=True
)
```

This provides a clean, efficient interface for using AnnData objects in PyTorch workflows while maintaining compatibility with standard PyTorch patterns and best practices.
