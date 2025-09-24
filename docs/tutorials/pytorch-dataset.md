# PyTorch Dataset for AnnData

The {class}`~anndata.experimental.pytorch.AnnDataset` class provides a PyTorch {class}`~torch.utils.data.Dataset` implementation for {class}`~anndata.AnnData` objects. This enables seamless integration with PyTorch DataLoaders and training workflows for single-cell data analysis.

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

The `transform` parameter accepts only Transform class instances that inherit from {class}`~anndata.experimental.pytorch.Transform`. This unified API ensures all transforms work seamlessly with multiprocessing by design.

All transforms must inherit from the Transform base class and implement the `__call__` method.

### Creating Custom Transform Classes

Create your own Transform classes that are serializable and work seamlessly with multiprocessing:

```python
from anndata.experimental.pytorch import Transform, Compose
import torch

class NormalizeTransform(Transform):
    """Normalize counts to target sum per cell."""
    def __init__(self, target_sum=1e4):
        self.target_sum = target_sum

    def __call__(self, x):
        # Ensure positive values
        x = torch.clamp(x, min=0)
        # Normalize to target sum
        row_sum = torch.sum(x, dim=-1, keepdim=True) + 1e-8
        return x * (self.target_sum / row_sum)

    def __repr__(self):
        return f"NormalizeTransform(target_sum={self.target_sum})"

class LogTransform(Transform):
    """Apply log1p transformation."""
    def __call__(self, x):
        return torch.log1p(x)

class AddNoise(Transform):
    """Add Gaussian noise for augmentation."""
    def __init__(self, noise_std=0.01, probability=0.5):
        self.noise_std = noise_std
        self.probability = probability

    def __call__(self, x):
        if torch.rand(1).item() < self.probability:
            noise = torch.normal(0, self.noise_std, size=x.shape)
            x = x + noise
            x = torch.clamp(x, min=0)  # Keep non-negative
        return x

# Compose transforms
training_transform = Compose([
    NormalizeTransform(target_sum=1e4),
    LogTransform(),
    AddNoise(noise_std=0.01, probability=0.5)
])

# Create dataset - works perfectly with multiprocessing!
dataset = AnnDataset(adata, transform=training_transform)
dataloader = DataLoader(dataset, batch_size=32, num_workers=4)
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

```{note}
**Multiprocessing with Transforms**

For the best multiprocessing experience, create your own Transform classes:

```python
from anndata.experimental.pytorch import Transform, Compose

# Transform classes work seamlessly with multiprocessing
class MyTransform(Transform):
    def __init__(self, param=1.0):
        self.param = param
    def __call__(self, x):
        return torch.log1p(x * self.param)

transform = MyTransform(param=1e4)
dataloader = DataLoader(dataset, num_workers=4)  # No pickling issues!
```

All transforms must inherit from the Transform base class - this ensures they work seamlessly with multiprocessing without any additional configuration.
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
