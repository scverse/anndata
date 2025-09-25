"""Test multiprocessing functionality of AnnDataset.

This module tests that AnnDataset works correctly with DataLoader
when using num_workers > 0, ensuring transforms are properly pickled.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

import anndata as ad
from anndata.experimental.pytorch import AnnDataset

try:
    import torch
except ImportError:
    torch = None


# Module-level transform classes for multiprocessing compatibility
try:
    from anndata.experimental.pytorch import Transform

    class LogTransform(Transform):
        """Simple log transform."""

        def __call__(self, x):
            return torch.log1p(x)

    class TrainingTransform(Transform):
        """Training transform with normalization and augmentation."""

        def __call__(self, x):
            # Normalize to 10k counts per cell
            row_sum = torch.sum(x, dim=-1, keepdim=True) + 1e-8
            x = x * (1e4 / row_sum)
            # Log transform
            x = torch.log1p(x)
            # Add small amount of noise
            x = x + 0.01 * torch.randn_like(x)
            return x

    class NormalizeTransform(Transform):
        """Simple normalization transform."""

        def __call__(self, x):
            return torch.log1p(x)

    class CustomTransform(Transform):
        """Simple custom transform."""

        def __call__(self, x):
            return x * 1.5

    class ComposedTransform(Transform):
        """Composed transform with normalization and noise."""

        def __call__(self, x):
            # Normalize
            row_sum = torch.sum(x, dim=-1, keepdim=True) + 1e-8
            x = x * (1e4 / row_sum)
            # Add noise
            x = x + 0.01 * torch.randn_like(x)
            return x

    # Create instances for use in tests
    log_transform = LogTransform()
    training_transform = TrainingTransform()
    normalize_transform = NormalizeTransform()
    custom_transform = CustomTransform()
    composed_transform = ComposedTransform()

except ImportError:
    # Fallback for when Transform is not available
    log_transform = None
    training_transform = None
    normalize_transform = None
    custom_transform = None
    composed_transform = None


class TestAnnDatasetMultiprocessing:
    """Test multiprocessing functionality of AnnDataset."""

    @pytest.fixture
    def sample_adata(self):
        """Create sample AnnData for testing."""
        np.random.seed(42)
        n_obs, n_vars = 100, 50
        X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)

        obs_data = {
            "cell_type": np.random.choice(["A", "B", "C"], n_obs),
            "batch": np.random.choice(["batch_1", "batch_2"], n_obs),
        }

        return ad.AnnData(X=X, obs=obs_data)

    @pytest.fixture
    def backed_adata(self, sample_adata):
        """Create backed AnnData for testing."""
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as f:
            temp_path = f.name

        # Write data to file
        sample_adata.write(temp_path)

        # Load as backed
        backed_adata = ad.read_h5ad(temp_path, backed="r")

        yield backed_adata

        # Cleanup
        Path(temp_path).unlink()

    def test_multiprocessing_with_in_memory_data(self, sample_adata):
        """Test multiprocessing with in-memory data."""
        try:
            from torch.utils.data import DataLoader
        except ImportError:
            pytest.skip("PyTorch not available")

        # Create dataset with simple transform function
        dataset = AnnDataset(sample_adata, transform=normalize_transform)

        # Test DataLoader with multiprocessing
        dataloader = DataLoader(dataset, batch_size=16, shuffle=True, num_workers=2)

        # Test that we can iterate through the dataloader
        batches = []
        for batch in dataloader:
            batches.append(batch)
            if len(batches) >= 3:  # Only test first few batches
                break

        assert len(batches) > 0
        assert "X" in batches[0]
        assert batches[0]["X"].shape[0] == 16  # batch_size
        assert batches[0]["X"].shape[1] == 50  # n_vars

    def test_multiprocessing_with_backed_data(self, backed_adata):
        """Test multiprocessing with backed data."""
        try:
            from torch.utils.data import DataLoader
        except ImportError:
            pytest.skip("PyTorch not available")

        # Create dataset with backed data using transform function
        dataset = AnnDataset(
            backed_adata,
            transform=training_transform,
            chunk_size=20,
        )

        # Test DataLoader with multiprocessing
        dataloader = DataLoader(
            dataset,
            batch_size=10,
            shuffle=False,  # Easier to test deterministically
            num_workers=2,
        )

        # Test that we can iterate through the dataloader
        batches = []
        for batch in dataloader:
            batches.append(batch)
            if len(batches) >= 2:  # Only test first few batches
                break

        assert len(batches) > 0
        assert "X" in batches[0]
        assert batches[0]["X"].shape[0] == 10  # batch_size
        assert batches[0]["X"].shape[1] == 50  # n_vars

    def test_multiprocessing_with_file_path(self, sample_adata):
        """Test multiprocessing with file path input."""
        try:
            from torch.utils.data import DataLoader
        except ImportError:
            pytest.skip("PyTorch not available")

        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as f:
            temp_path = f.name

        try:
            # Write data to file
            sample_adata.write(temp_path)

            # Create dataset from file path
            dataset = AnnDataset(temp_path, transform=log_transform)

            # Test DataLoader with multiprocessing
            dataloader = DataLoader(dataset, batch_size=8, shuffle=False, num_workers=2)

            # Test that we can iterate through the dataloader
            batches = []
            for batch in dataloader:
                batches.append(batch)
                if len(batches) >= 2:
                    break

            assert len(batches) > 0
            assert "X" in batches[0]
            assert batches[0]["X"].shape[0] == 8  # batch_size
            assert batches[0]["X"].shape[1] == 50  # n_vars

        finally:
            # Cleanup
            Path(temp_path).unlink()

    def test_multiprocessing_with_optimized_dataloader(self, sample_adata):
        """Test multiprocessing with optimized DataLoader."""
        import importlib.util

        if importlib.util.find_spec("torch") is None:
            pytest.skip("PyTorch not available")

        # Create dataset
        dataset = AnnDataset(sample_adata, transform=training_transform)

        # Test optimized DataLoader with multiprocessing
        optimized_loader = dataset.get_optimized_dataloader(
            batch_size=12, shuffle=False, num_workers=2
        )

        # Test that we can iterate through the optimized loader
        batches = []
        for batch in optimized_loader:
            batches.append(batch)
            if len(batches) >= 2:
                break

        assert len(batches) > 0
        assert "X" in batches[0]
        assert batches[0]["X"].shape[0] == 12  # batch_size
        assert batches[0]["X"].shape[1] == 50  # n_vars

    def test_multiprocessing_with_obs_subset(self, sample_adata):
        """Test multiprocessing with observation subsetting."""
        try:
            from torch.utils.data import DataLoader
        except ImportError:
            pytest.skip("PyTorch not available")

        # Create subset indices
        subset_indices = np.arange(0, 50)  # First 50 observations

        # Create dataset with subset
        dataset = AnnDataset(
            sample_adata,
            obs_subset=subset_indices,
            transform=log_transform,
        )

        # Test DataLoader with multiprocessing
        dataloader = DataLoader(dataset, batch_size=10, shuffle=False, num_workers=2)

        # Test that we can iterate through the dataloader
        batches = []
        for batch in dataloader:
            batches.append(batch)
            if len(batches) >= 2:
                break

        assert len(batches) > 0
        assert "X" in batches[0]
        assert batches[0]["X"].shape[0] == 10  # batch_size
        assert batches[0]["X"].shape[1] == 50  # n_vars
        assert len(dataset) == 50  # Should have 50 observations

    def test_transform_consistency_across_workers(self, sample_adata):
        """Test that transforms are applied consistently across workers."""
        try:
            from torch.utils.data import DataLoader
        except ImportError:
            pytest.skip("PyTorch not available")

        # Create dataset with deterministic transform
        dataset = AnnDataset(
            sample_adata,
            transform=log_transform,
        )

        # Test with single worker first
        single_worker_loader = DataLoader(
            dataset, batch_size=16, shuffle=False, num_workers=0
        )

        # Test with multiple workers
        multi_worker_loader = DataLoader(
            dataset, batch_size=16, shuffle=False, num_workers=2
        )

        # Get first batch from each
        single_batch = next(iter(single_worker_loader))
        multi_batch = next(iter(multi_worker_loader))

        # Both should have same shape and reasonable values
        assert single_batch["X"].shape == multi_batch["X"].shape
        assert torch.allclose(single_batch["X"], multi_batch["X"], atol=1e-6)

    def test_multiprocessing_error_handling(self, sample_adata):
        """Test error handling in multiprocessing scenarios."""
        try:
            import torch  # noqa: F401
        except ImportError:
            pytest.skip("PyTorch not available")

        # This should raise an error during initialization
        with pytest.raises(ValueError, match="chunk_size must be positive"):
            AnnDataset(
                sample_adata,
                transform=log_transform,
                chunk_size=0,  # Invalid chunk size
            )

    def test_no_multiprocessing_fallback(self, sample_adata):
        """Test that dataset works without multiprocessing."""
        try:
            from torch.utils.data import DataLoader
        except ImportError:
            pytest.skip("PyTorch not available")

        # Create dataset without multiprocessing
        dataset = AnnDataset(
            sample_adata,
            transform=log_transform,
        )

        # Test DataLoader without workers
        dataloader = DataLoader(dataset, batch_size=16, shuffle=False, num_workers=0)

        # Test that we can iterate through the dataloader
        batch = next(iter(dataloader))
        assert "X" in batch
        assert batch["X"].shape[0] == 16
        assert batch["X"].shape[1] == 50

    def test_custom_transform_class_multiprocessing(self, sample_adata):
        """Test multiprocessing with custom transform classes."""
        try:
            from torch.utils.data import DataLoader
        except ImportError:
            pytest.skip("PyTorch not available")

        # Create dataset with custom transform
        dataset = AnnDataset(
            sample_adata,
            transform=custom_transform,
        )

        # Test DataLoader with multiprocessing
        dataloader = DataLoader(dataset, batch_size=8, shuffle=False, num_workers=2)

        # Test that we can iterate through the dataloader
        batch = next(iter(dataloader))
        assert "X" in batch
        assert batch["X"].shape[0] == 8
        assert batch["X"].shape[1] == 50

    def test_composed_transforms_multiprocessing(self, sample_adata):
        """Test multiprocessing with composed transforms."""
        try:
            from torch.utils.data import DataLoader
        except ImportError:
            pytest.skip("PyTorch not available")

        # Create dataset with composed transform
        dataset = AnnDataset(
            sample_adata,
            transform=composed_transform,
        )

        # Test DataLoader with multiprocessing
        dataloader = DataLoader(dataset, batch_size=16, shuffle=False, num_workers=2)

        # Test that we can iterate through the dataloader
        batch = next(iter(dataloader))
        assert "X" in batch
        assert batch["X"].shape[0] == 16
        assert batch["X"].shape[1] == 50
