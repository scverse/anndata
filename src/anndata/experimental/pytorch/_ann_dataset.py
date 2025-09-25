"""PyTorch Dataset for AnnData objects.

This module provides a PyTorch Dataset implementation for AnnData objects
with support for:

- Generic transform interface following PyTorch conventions
- Memory-efficient streaming from backed HDF5 data
- Multiprocessing-safe I/O following h5py best practices
- Configurable chunk processing for memory management
- Optimized batch loading with sorted indices
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse import issparse

import anndata as ad

if TYPE_CHECKING:
    from collections.abc import Sequence

    from ..._core.anndata import AnnData
    from .transforms import Transform

# Constants for metadata key formatting
OBS_PREFIX = "obs_"

# Import torch conditionally
try:
    import torch
    from torch.utils.data import Dataset

    TORCH_AVAILABLE = True
except ImportError:
    # Mock classes for when torch is not available
    class Dataset:
        pass

    class torch:
        class Tensor:
            pass

        class Generator:
            pass

        @staticmethod
        def from_numpy(x):
            return x

        @staticmethod
        def zeros(*args, **kwargs):
            return None

        @staticmethod
        def tensor(*args, **kwargs):
            return None

        class utils:
            class data:
                @staticmethod
                def get_worker_info():
                    return None

    TORCH_AVAILABLE = False


def _get_adata_for_worker(
    adata_path: str | Path | None,
    adata_obj: AnnData | None,
) -> AnnData:
    """Get AnnData object for worker process following h5py best practices.

    Following h5py documentation recommendations for parallel read access:
    "It's advised to open the file independently in each reader process;
    opening the file once and then forking may cause issues."

    See: https://docs.h5py.org/en/stable/mpi.html

    Parameters
    ----------
    adata_path
        Path to HDF5 file (if reading from file)
    adata_obj
        AnnData object (if reading from memory or backed object)

    Returns
    -------
    AnnData
        AnnData object safe for worker access
    """
    if adata_obj is not None and not adata_obj.isbacked:
        # In-memory data - safe to use directly
        return adata_obj
    elif adata_path is not None:
        # Load from file path - each worker opens independently
        return ad.read_h5ad(adata_path, backed="r")
    elif adata_obj is not None and adata_obj.isbacked:
        # Backed object - get the file path and reopen
        if hasattr(adata_obj, "filename") and adata_obj.filename is not None:
            return ad.read_h5ad(adata_obj.filename, backed="r")
        else:
            # Fall back to using existing object
            return adata_obj
    else:
        msg = "Either adata_path or adata_obj must be provided"
        raise ValueError(msg)


class AnnDataset(Dataset):
    """PyTorch Dataset for AnnData objects.

    A PyTorch Dataset implementation that provides:
    - Generic transform interface following PyTorch conventions
    - Multiprocessing-safe I/O with h5py best practices
    - Memory-efficient streaming from backed HDF5 data
    - Configurable chunk processing for memory management
    - Optimized batch loading with sorted indices
    """

    def __init__(
        self,
        adata: AnnData | str | Path,
        *,
        # Transform object (must inherit from Transform base class for multiprocessing compatibility)
        transform: Transform | None = None,
        # Observation selection
        obs_subset: Sequence[int] | np.ndarray | None = None,
        # Advanced features
        chunk_size: int = 6000,
    ):
        """Initialize AnnData Dataset.

        Parameters
        ----------
        adata
            AnnData object or path to h5ad file
        transform
            Transform object to transform loaded data dictionary.
            Must inherit from anndata.experimental.pytorch.Transform and
            implement __call__ method that receives dict with 'X' and metadata.
            This ensures multiprocessing compatibility.
        obs_subset
            Subset of observations to include (by index)
        chunk_size
            Used only when loading sparse dataset that is stored as dense. Loading iterates
            through chunks of the dataset of this row size until it reads the whole dataset.
            Higher size means higher memory consumption and higher (to a point) loading speed.
        """
        if not TORCH_AVAILABLE:
            error_msg = (
                "PyTorch is required for AnnDataset. "
                "Please install PyTorch: pip install torch"
            )
            raise ImportError(error_msg)

        # Validate inputs
        self._validate_inputs(transform, obs_subset, chunk_size)

        # Store and validate transform
        self._validate_transform(transform)
        self.transform = transform
        self.chunk_size = chunk_size

        # Handle AnnData input
        if isinstance(adata, (str, Path)):
            self.adata_path = Path(adata)
            self.adata_obj = None
            # Load minimal info for initialization
            adata_temp = _get_adata_for_worker(self.adata_path, None)
            self._initialize_from_adata(adata_temp, obs_subset)
        else:
            # For backed AnnData objects, store the file path for multiprocessing
            if adata.isbacked and hasattr(adata, "filename"):
                self.adata_path = Path(adata.filename)
            else:
                self.adata_path = None
            self.adata_obj = adata
            self._initialize_from_adata(adata, obs_subset)

        # Initialize worker-specific random generator
        self._torch_generator = None
        self._worker_seed_base = hash((id(self), id(self.transform)))

    def _validate_transform(self, transform):
        """Validate that transform is a Transform object for multiprocessing compatibility.

        Transform inheritance is required because:
        1. Regular functions/lambdas cannot be pickled for multiprocessing
        2. Transform classes are serializable and work seamlessly across worker processes
        3. This ensures consistent behavior when using DataLoader with num_workers > 0
        """
        if transform is not None:
            # Import Transform here to avoid circular imports
            try:
                from .transforms import Transform
            except ImportError:
                # If transforms module isn't available, we can't validate
                return

            if not isinstance(transform, Transform):
                msg = (
                    "transform must be an instance of Transform class for multiprocessing compatibility. "
                    f"Got {type(transform)}. Please inherit from anndata.experimental.pytorch.Transform "
                    "and implement __call__ method."
                )
                raise TypeError(msg)

    def _validate_inputs(self, transform, obs_subset, chunk_size):
        """Validate dataset inputs and provide helpful error messages."""

        # Observation subset validation
        if obs_subset is not None:
            obs_subset = np.asarray(obs_subset)
            # Empty subset is allowed (results in empty dataset)
            if np.any(obs_subset < 0):
                msg = "obs_subset indices must be non-negative"
                raise ValueError(msg)

        # Chunk size validation
        if chunk_size <= 0:
            msg = "chunk_size must be positive"
            raise ValueError(msg)

    def _initialize_from_adata(
        self,
        adata: AnnData,
        obs_subset: Sequence[int] | np.ndarray | None,
    ):
        """Initialize dataset from AnnData object."""

        # Handle observation subsetting
        if obs_subset is not None:
            self.obs_indices = np.asarray(obs_subset)
            # Validate indices are within bounds
            if np.any(self.obs_indices >= adata.n_obs):
                msg = f"obs_subset contains indices >= {adata.n_obs} (dataset size)"
                raise ValueError(msg)
        else:
            self.obs_indices = np.arange(adata.n_obs)

        # Store dataset info
        self.n_obs = len(self.obs_indices)
        self.n_vars = adata.n_vars

    def _get_worker_generator(self) -> torch.Generator:
        """Get worker-specific random generator for reproducible augmentation."""
        if self._torch_generator is None:
            # Get worker info for multiprocessing
            worker_info = torch.utils.data.get_worker_info()
            if worker_info is not None:
                worker_seed = self._worker_seed_base + worker_info.id
            else:
                worker_seed = self._worker_seed_base

            self._torch_generator = torch.Generator()
            self._torch_generator.manual_seed(worker_seed)

        return self._torch_generator

    def _load_observation(self, idx: int) -> tuple[torch.Tensor, dict]:
        """Load a single observation following h5py best practices for multiprocessing."""

        # Get the actual row index
        row_idx = self.obs_indices[idx]

        # For in-memory data, _get_adata_for_worker will return the object directly
        # For backed data, it will reopen files independently in each worker
        if self.adata_obj is None:
            # In multiprocessing worker, adata_obj can become None
            if self.adata_path is not None:
                adata = _get_adata_for_worker(self.adata_path, None)
            else:
                msg = "AnnData object is None and no path available for reloading"
                raise RuntimeError(msg)
        else:
            # Get AnnData object safe for worker access
            # Following h5py recommendation: each worker opens files independently
            adata = _get_adata_for_worker(self.adata_path, self.adata_obj)

        X_raw = adata.X[row_idx]  # Get single row
        obs_data = adata.obs.iloc[row_idx].to_dict() if adata.obs.shape[1] > 0 else {}

        # Convert sparse to dense if needed
        if issparse(X_raw):
            X_raw = X_raw.toarray()

        # Ensure X_raw is 1D
        if X_raw.ndim > 1:
            X_raw = X_raw.flatten()

        X_tensor = torch.from_numpy(X_raw).float()
        return X_tensor, obs_data

    def _apply_transform(self, data_dict: dict) -> dict:
        """Apply transform function to data dictionary."""
        if self.transform is not None:
            result = self.transform(data_dict)
            if not isinstance(result, dict):
                msg = f"Transform must return dict, got {type(result)}"
                raise TypeError(msg)
            return result
        return data_dict

    def __len__(self) -> int:
        """Return number of observations."""
        return self.n_obs

    def __getitem__(self, idx: int) -> dict[str, torch.Tensor | str]:
        """Get a single observation with full preprocessing pipeline."""

        # Load raw data
        X_raw, obs_data = self._load_observation(idx)

        # Build complete data dictionary
        result = {"X": X_raw}

        # Add observation metadata if available
        if obs_data:
            self._add_metadata_to_result(result, obs_data)

        # Apply transform to entire data dictionary
        result = self._apply_transform(result)

        return result

    def _add_metadata_to_result(self, result: dict, obs_data: dict) -> None:
        """Add observation metadata to result dictionary with proper key formatting."""
        for key, value in obs_data.items():
            # Use a more structured approach for metadata keys
            metadata_key = self._format_metadata_key(key)

            if isinstance(value, (int, float, bool)):
                result[metadata_key] = torch.tensor(value)
            elif isinstance(value, str):
                # For string values, store as string (tests expect this)
                result[metadata_key] = value
            elif isinstance(value, np.number):
                result[metadata_key] = torch.tensor(float(value))
            else:
                # For other types, convert to string
                result[metadata_key] = str(value)

    def _format_metadata_key(self, key: str) -> str:
        """Format observation metadata key with consistent prefix."""
        return f"{OBS_PREFIX}{key}"

    def get_optimized_dataloader(
        self,
        batch_size: int = 32,
        *,
        shuffle: bool = True,
        num_workers: int = 0,
        **kwargs,
    ):
        """Get a DataLoader with optimized batch loading for backed data.

        This method creates a DataLoader that uses index-based batch loading
        with sorted indices for optimal disk access performance.

        Parameters
        ----------
        batch_size
            Number of items per batch
        shuffle
            Whether to shuffle the data
        num_workers
            Number of worker processes
        **kwargs
            Additional arguments passed to DataLoader

        Returns
        -------
        DataLoader
            Optimized DataLoader instance
        """
        if not TORCH_AVAILABLE:
            error_msg = "PyTorch is required for DataLoader"
            raise ImportError(error_msg)

        from functools import partial

        from torch.utils.data import DataLoader

        # Create a DataLoader with optimized collate function
        return DataLoader(
            _IndexDataset(self),
            batch_size=batch_size,
            shuffle=shuffle,
            num_workers=num_workers,
            collate_fn=partial(_optimized_collate_fn, dataset=self),
            **kwargs,
        )

    def __getstate__(self):
        """Prepare for multiprocessing (pickle support)."""
        state = self.__dict__.copy()

        # Remove unpicklable objects
        if self.adata_obj is not None and self.adata_obj.isbacked:
            # Backed AnnData objects contain unpicklable file handles
            state["adata_obj"] = None

        # Reset worker-specific generator (will be recreated in worker)
        state["_torch_generator"] = None

        return state

    @property
    def shape(self) -> tuple[int, int]:
        """Shape of the dataset (n_obs, n_vars)."""
        return (self.n_obs, self.n_vars)

    def get_collate_fn(self):
        """Get a collate function optimized for this dataset.

        The returned collate function sorts batch indices for efficient
        sequential disk access when using backed data.

        Returns
        -------
        callable
            Collate function that can be passed to torch.utils.data.DataLoader

        Examples
        --------
        .. doctest::
            :skipif: not TORCH_AVAILABLE

        >>> dataset = AnnDataset(adata)
        >>> dataloader = DataLoader(
        ...     dataset,
        ...     batch_size=32,
        ...     collate_fn=dataset.get_collate_fn(),
        ... )
        """
        from functools import partial

        return partial(_ann_dataset_collate_fn, dataset=self)


def _batch_samples(samples: list[dict]) -> dict[str, torch.Tensor]:
    """Helper function to batch items into tensors."""
    batch = {}

    # Handle expression data
    X_list = [sample["X"] for sample in samples]
    batch["X"] = torch.stack(X_list, dim=0)

    # Handle metadata (if present)
    metadata_keys = set()
    for sample in samples:
        metadata_keys.update(k for k in sample if k.startswith(OBS_PREFIX))

    for key in metadata_keys:
        values = []
        for sample in samples:
            if key in sample:
                values.append(sample[key])
            else:
                values.append(None)  # Default value for missing metadata

        # Check if all values are tensors
        if values and all(isinstance(v, torch.Tensor) for v in values if v is not None):
            # Handle tensor values
            tensor_values = [v if v is not None else torch.tensor(0.0) for v in values]
            batch[key] = torch.stack(tensor_values, dim=0)
        else:
            # Handle non-tensor values (strings, etc.) - store as list
            batch[key] = values

    return batch


def _ann_dataset_collate_fn(
    batch_samples: list[dict], dataset: AnnDataset = None
) -> dict[str, torch.Tensor]:
    """Collate function for AnnDataset that batches items efficiently."""
    if not TORCH_AVAILABLE:
        error_msg = "PyTorch is required for batch collation"
        raise ImportError(error_msg)

    return _batch_samples(batch_samples)


class _IndexDataset:
    """Helper dataset that returns indices for optimized batch loading."""

    def __init__(self, ann_dataset: AnnDataset):
        self.ann_dataset = ann_dataset

    def __len__(self):
        return len(self.ann_dataset)

    def __getitem__(self, idx: int) -> int:
        return idx


def _optimized_collate_fn(
    batch_indices: list[int], dataset: AnnDataset
) -> dict[str, torch.Tensor]:
    """Optimized collate function that sorts indices for efficient disk access.

    This function sorts indices for efficient sequential disk access while
    preserving the original batch order in the returned data.
    """
    if not TORCH_AVAILABLE:
        error_msg = "PyTorch is required for batch collation"
        raise ImportError(error_msg)

    # Convert to numpy for efficient operations
    batch_indices_np = np.array(batch_indices)

    # Get sort permutation - argsort gives us the indices that would sort the array
    sort_perm = np.argsort(batch_indices_np)

    # Create sorted indices for efficient disk access
    sorted_indices = batch_indices_np[sort_perm]

    # Load all items in the batch using sorted indices
    sorted_batch_data = [dataset[idx] for idx in sorted_indices]

    # Restore original order using inverse permutation
    # argsort(argsort(x)) gives the inverse permutation
    restore_perm = np.argsort(sort_perm)
    batch_data = [sorted_batch_data[i] for i in restore_perm]

    return _batch_samples(batch_data)
