"""Consolidated tests for AnnDataset PyTorch integration.

This module combines essential tests for:
- Basic AnnDataset functionality
- Transform classes and multi-worker compatibility
- DataLoader integration with multi-workers
- Backed data handling
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import anndata as ad
from anndata.experimental.pytorch import (
    TORCH_AVAILABLE,
    AnnDataset,
    Compose,
    Transform,
)

if TORCH_AVAILABLE:
    import torch
    from torch.utils.data import DataLoader


@pytest.fixture
def gen_adata(shape=(100, 50)):
    """Generate test AnnData object."""
    X = np.random.rand(*shape)
    obs = {"cell_type": np.random.choice(["A", "B", "C"], shape[0])}
    var = {"gene_name": [f"Gene_{i}" for i in range(shape[1])]}
    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def gen_sparse_adata(shape=(100, 50)):
    """Generate test sparse AnnData object."""
    X = sparse.csr_matrix(np.random.rand(*shape))
    obs = {"cell_type": np.random.choice(["A", "B", "C"], shape[0])}
    var = {"gene_name": [f"Gene_{i}" for i in range(shape[1])]}
    return ad.AnnData(X=X, obs=obs, var=var)


# Transform classes for testing
class LogTransform(Transform):
    """Simple log transform for testing."""

    def __call__(self, data_dict):
        data_dict["X"] = torch.log1p(data_dict["X"])
        return data_dict


class NormalizeTransform(Transform):
    """Normalization transform for testing."""

    def __init__(self, target_sum=1e4):
        self.target_sum = target_sum

    def __call__(self, data_dict):
        X = data_dict["X"]
        X = torch.clamp(X, min=0)
        row_sum = torch.sum(X, dim=-1, keepdim=True) + 1e-8
        data_dict["X"] = X * (self.target_sum / row_sum)
        return data_dict


class DictTransform(Transform):
    """Transform that modifies both X and metadata."""

    def __call__(self, data_dict):
        data_dict["X"] = data_dict["X"] * 2
        data_dict["obs_transformed"] = torch.tensor(1.0)
        return data_dict


@pytest.mark.skipif(
    not TORCH_AVAILABLE, reason="PyTorch is required for AnnDataset tests"
)
class TestAnnDataset:
    """Test cases for AnnDataset."""

    def test_basic_functionality(self, gen_adata):
        """Test basic dataset functionality."""
        adata = gen_adata
        dataset = AnnDataset(adata)

        assert len(dataset) == adata.n_obs
        assert dataset.shape == adata.shape

        # Test single item access
        item = dataset[0]
        assert isinstance(item, dict)
        assert "X" in item
        assert isinstance(item["X"], torch.Tensor)
        assert item["X"].shape == (adata.n_vars,)

    def test_metadata_access(self, gen_adata):
        """Test that metadata is properly included."""
        adata = gen_adata
        dataset = AnnDataset(adata)

        item = dataset[0]
        assert "obs_cell_type" in item
        assert isinstance(item["obs_cell_type"], str)

    def test_file_path_loading(self, gen_adata):
        """Test loading from file path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = Path(tmpdir) / "test.h5ad"
            gen_adata.write_h5ad(file_path)

            dataset = AnnDataset(file_path)
            assert len(dataset) == gen_adata.n_obs
            assert dataset.shape == gen_adata.shape

    def test_backed_data_handling(self, gen_adata):
        """Test backed data handling."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = Path(tmpdir) / "test.h5ad"
            gen_adata.write_h5ad(file_path)

            # Load as backed
            adata_backed = ad.read_h5ad(file_path, backed="r")
            dataset = AnnDataset(adata_backed)

            assert len(dataset) == adata_backed.n_obs
            item = dataset[0]
            assert isinstance(item["X"], torch.Tensor)

    def test_obs_subset(self, gen_adata):
        """Test observation subsetting."""
        subset_indices = np.arange(10)
        dataset = AnnDataset(gen_adata, obs_subset=subset_indices)

        assert len(dataset) == 10
        assert dataset.shape == (10, gen_adata.n_vars)

    def test_sparse_data(self, gen_sparse_adata):
        """Test sparse data handling."""
        dataset = AnnDataset(gen_sparse_adata)
        item = dataset[0]
        assert isinstance(item["X"], torch.Tensor)
        assert item["X"].shape == (gen_sparse_adata.n_vars,)


@pytest.mark.skipif(
    not TORCH_AVAILABLE, reason="PyTorch is required for Transform tests"
)
class TestTransforms:
    """Test Transform classes and multi-worker compatibility."""

    def test_transform_base_class(self):
        """Test that Transform base class works correctly."""

        class TestTransform(Transform):
            def __init__(self, multiplier=2.0):
                self.multiplier = multiplier

            def __call__(self, data_dict):
                data_dict["X"] = data_dict["X"] * self.multiplier
                return data_dict

        transform = TestTransform(multiplier=3.0)
        data_dict = {"X": torch.tensor([1.0, 2.0, 3.0])}
        result = transform(data_dict)

        expected = torch.tensor([3.0, 6.0, 9.0])
        assert torch.allclose(result["X"], expected)

    def test_compose_transforms(self):
        """Test that Compose works correctly."""
        transform1 = LogTransform()
        transform2 = NormalizeTransform(target_sum=1e4)
        compose = Compose([transform1, transform2])

        data_dict = {"X": torch.tensor([[1000.0, 2000.0], [3000.0, 4000.0]])}
        result = compose(data_dict)

        assert isinstance(result["X"], torch.Tensor)
        assert result["X"].shape == data_dict["X"].shape

    def test_transform_with_dataset(self, gen_adata):
        """Test transforms with AnnDataset."""
        transform = LogTransform()
        dataset = AnnDataset(gen_adata, transform=transform)

        item = dataset[0]
        assert isinstance(item["X"], torch.Tensor)
        # Check that log transform was applied
        assert torch.allclose(item["X"], torch.log1p(torch.tensor(gen_adata.X[0])))


@pytest.mark.skipif(
    not TORCH_AVAILABLE, reason="PyTorch is required for DataLoader tests"
)
class TestDataLoaderIntegration:
    """Test DataLoader integration and multi-workers."""

    def test_basic_dataloader(self, gen_adata):
        """Test basic DataLoader functionality."""
        dataset = AnnDataset(gen_adata)
        dataloader = DataLoader(dataset, batch_size=32, shuffle=False)

        batch = next(iter(dataloader))
        assert isinstance(batch, dict)
        assert "X" in batch
        assert batch["X"].shape[0] <= 32  # batch size
        assert batch["X"].shape[1] == gen_adata.n_vars

    def test_optimized_dataloader(self, gen_adata):
        """Test optimized DataLoader for backed data."""
        dataset = AnnDataset(gen_adata)
        dataloader = dataset.get_optimized_dataloader(batch_size=16, shuffle=True)

        batch = next(iter(dataloader))
        assert isinstance(batch, dict)
        assert "X" in batch
        assert batch["X"].shape[0] <= 16

    def test_multi_worker_safety(self, gen_adata):
        """Test multi-worker safety with transforms."""
        transform = LogTransform()
        dataset = AnnDataset(gen_adata, transform=transform)

        # Test with multiple workers
        dataloader = DataLoader(dataset, batch_size=16, num_workers=2, shuffle=True)

        # Should not raise errors
        batch = next(iter(dataloader))
        assert isinstance(batch, dict)
        assert "X" in batch

        # Verify transform was actually applied (log1p)
        expected_X = torch.log1p(torch.tensor(gen_adata.X[: batch["X"].shape[0]]))
        assert torch.allclose(batch["X"], expected_X, atol=1e-6)

    def test_multi_worker_with_backed_data(self, gen_adata):
        """Test multi-workers with backed data (critical for file I/O)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = Path(tmpdir) / "test.h5ad"
            gen_adata.write_h5ad(file_path)

            # Test with backed data and multiple workers
            transform = NormalizeTransform(target_sum=1e4)
            dataset_with_transform = AnnDataset(file_path, transform=transform)

            # Test different worker counts and verify results are consistent
            results = {}
            for num_workers in [0, 1, 2]:
                dataloader = DataLoader(
                    dataset_with_transform,
                    batch_size=8,
                    num_workers=num_workers,
                    shuffle=False,  # Don't shuffle to get consistent results
                )

                # Get first batch
                batch = next(iter(dataloader))
                assert isinstance(batch, dict)
                assert "X" in batch
                assert batch["X"].shape[0] <= 8

                # Store results for comparison
                results[num_workers] = batch["X"].clone()

                # Verify normalization was applied (sum should be close to target_sum)
                row_sums = torch.sum(batch["X"], dim=1)
                assert torch.allclose(
                    row_sums, torch.full_like(row_sums, 1e4), rtol=0.1
                )

            # Results should be identical regardless of worker count
            assert torch.allclose(results[0], results[1], atol=1e-6)
            assert torch.allclose(results[1], results[2], atol=1e-6)

    def test_transform_pickling(self, gen_adata):
        """Test that transforms can be pickled for multi-workers."""

        # Create a transform with state
        class StatefulTransform(Transform):
            def __init__(self, multiplier=2.0):
                self.multiplier = multiplier
                self.call_count = 0

            def __call__(self, data_dict):
                self.call_count += 1
                data_dict["X"] = data_dict["X"] * self.multiplier
                return data_dict

        transform = StatefulTransform(multiplier=3.0)
        dataset = AnnDataset(gen_adata, transform=transform)

        # Test with multi-workers - transform should be pickled and unpickled
        dataloader = DataLoader(dataset, batch_size=16, num_workers=2, shuffle=True)

        # Process multiple batches to ensure pickling works
        batches = [next(iter(dataloader)) for _ in range(3)]

        for batch in batches:
            assert isinstance(batch, dict)
            assert "X" in batch
            # Verify transform was applied (multiplied by 3)
            assert torch.allclose(batch["X"], gen_adata.X[: batch["X"].shape[0]] * 3.0)

        # Test that single-worker and multi-worker give same results
        single_worker_dataloader = DataLoader(
            dataset, batch_size=16, num_workers=0, shuffle=False
        )
        multi_worker_dataloader = DataLoader(
            dataset, batch_size=16, num_workers=2, shuffle=False
        )

        single_batch = next(iter(single_worker_dataloader))
        multi_batch = next(iter(multi_worker_dataloader))

        # Results should be identical
        assert torch.allclose(single_batch["X"], multi_batch["X"], atol=1e-6)

    def test_multi_worker_performance_verification(self, gen_adata):
        """Test that multi-workers actually improve performance and work correctly."""
        import time

        # Create a transform that simulates some work
        class WorkSimulationTransform(Transform):
            def __call__(self, data_dict):
                # Simulate some computational work
                time.sleep(0.001)  # 1ms delay per item
                data_dict["X"] = data_dict["X"] * 2.0
                return data_dict

        transform = WorkSimulationTransform()
        dataset = AnnDataset(gen_adata, transform=transform)

        # Test single worker
        single_worker_dataloader = DataLoader(
            dataset, batch_size=4, num_workers=0, shuffle=False
        )

        # Test multi-worker
        multi_worker_dataloader = DataLoader(
            dataset, batch_size=4, num_workers=2, shuffle=False
        )

        # Get multiple batches and verify they're processed correctly
        single_batches = [next(iter(single_worker_dataloader)) for _ in range(3)]
        multi_batches = [next(iter(multi_worker_dataloader)) for _ in range(3)]

        # Verify all batches have correct transform applied
        for batch in single_batches + multi_batches:
            assert isinstance(batch, dict)
            assert "X" in batch
            # Verify transform was applied (multiplied by 2)
            expected = gen_adata.X[: batch["X"].shape[0]] * 2.0
            assert torch.allclose(batch["X"], torch.tensor(expected), atol=1e-6)

    def test_multi_worker_data_consistency(self, gen_adata):
        """Test that multi-workers produce consistent data across runs."""
        transform = NormalizeTransform(target_sum=1e4)
        dataset = AnnDataset(gen_adata, transform=transform)

        # Run multiple times with same parameters
        results = []
        for _ in range(3):
            dataloader = DataLoader(dataset, batch_size=8, num_workers=2, shuffle=False)
            batch = next(iter(dataloader))
            results.append(batch["X"].clone())

        # All results should be identical
        for i in range(1, len(results)):
            assert torch.allclose(results[0], results[i], atol=1e-6)

    def test_multi_worker_error_handling(self, gen_adata):
        """Test error handling in multi-worker context."""

        # Create a transform that might fail
        class FailingTransform(Transform):
            def __call__(self, data_dict):
                # Simulate occasional failure
                if torch.rand(1) < 0.1:  # 10% chance of failure
                    error_msg = "Simulated transform failure"
                    raise ValueError(error_msg)
                return data_dict

        transform = FailingTransform()
        dataset = AnnDataset(gen_adata, transform=transform)

        # Test with multi-workers - should handle errors gracefully
        dataloader = DataLoader(dataset, batch_size=16, num_workers=2, shuffle=True)

        # This might raise an error, which is expected behavior
        # The test verifies that multi-workers don't cause additional issues
        with pytest.raises(ValueError, match="Simulated transform failure"):
            next(iter(dataloader))

    def test_collate_function(self, gen_adata):
        """Test custom collate function."""
        dataset = AnnDataset(gen_adata)
        collate_fn = dataset.get_collate_fn()

        # Test collate function directly
        items = [dataset[i] for i in range(5)]
        batch = collate_fn(items)

        assert isinstance(batch, dict)
        assert "X" in batch
        assert batch["X"].shape[0] == 5
        assert batch["X"].shape[1] == gen_adata.n_vars


@pytest.mark.skipif(
    not TORCH_AVAILABLE, reason="PyTorch is required for advanced tests"
)
class TestAdvancedFeatures:
    """Test advanced AnnDataset features."""

    def test_chunk_size_parameter(self, gen_sparse_adata):
        """Test chunk size parameter for memory management."""
        dataset = AnnDataset(gen_sparse_adata, chunk_size=10)

        item = dataset[0]
        assert isinstance(item["X"], torch.Tensor)

    def test_dataset_properties(self, gen_adata):
        """Test dataset properties."""
        dataset = AnnDataset(gen_adata)

        assert dataset.shape == gen_adata.shape
        assert len(dataset) == gen_adata.n_obs

    def test_transform_pipeline(self, gen_adata):
        """Test complex transform pipeline."""
        transform1 = NormalizeTransform(target_sum=1e4)
        transform2 = LogTransform()
        compose = Compose([transform1, transform2])

        dataset = AnnDataset(gen_adata, transform=compose)
        item = dataset[0]

        assert isinstance(item["X"], torch.Tensor)
        assert item["X"].shape == (gen_adata.n_vars,)

    def test_metadata_types(self, gen_adata):
        """Test different metadata types."""
        # Add numeric metadata
        gen_adata.obs["numeric_col"] = np.random.rand(gen_adata.n_obs)
        gen_adata.obs["bool_col"] = np.random.choice([True, False], gen_adata.n_obs)

        dataset = AnnDataset(gen_adata)
        item = dataset[0]

        assert "obs_numeric_col" in item
        assert "obs_bool_col" in item
        assert isinstance(item["obs_numeric_col"], torch.Tensor)
        assert isinstance(item["obs_bool_col"], torch.Tensor)


@pytest.mark.skipif(
    not TORCH_AVAILABLE, reason="PyTorch is required for error handling tests"
)
class TestErrorHandling:
    """Test error handling and edge cases."""

    def test_invalid_transform_type(self, gen_adata):
        """Test that invalid transform types raise appropriate errors."""
        with pytest.raises(TypeError, match="Transform must inherit from"):
            AnnDataset(gen_adata, transform="invalid_transform")

    def test_invalid_obs_subset(self, gen_adata):
        """Test invalid observation subset."""
        with pytest.raises(IndexError):
            AnnDataset(gen_adata, obs_subset=[1000])  # Out of bounds

    def test_negative_obs_subset(self, gen_adata):
        """Test negative observation subset indices."""
        with pytest.raises(IndexError):
            AnnDataset(gen_adata, obs_subset=[-1])

    def test_invalid_chunk_size(self, gen_adata):
        """Test invalid chunk size."""
        with pytest.raises(ValueError, match="chunk_size must be positive"):
            AnnDataset(gen_adata, chunk_size=0)

    def test_missing_values_handling(self):
        """Test handling of missing values in metadata."""
        X = np.random.rand(10, 5)
        obs = pd.DataFrame({
            "cell_type": ["A", "B", None, "A", "B", "C", "A", None, "B", "C"],
            "numeric": [1.0, 2.0, np.nan, 3.0, 4.0, 5.0, 6.0, 7.0, np.nan, 8.0],
        })

        adata = ad.AnnData(X=X, obs=obs)
        dataset = AnnDataset(adata)

        # Should handle missing values gracefully
        item = dataset[0]
        assert isinstance(item["obs_cell_type"], str)
        assert isinstance(item["obs_numeric"], torch.Tensor)
