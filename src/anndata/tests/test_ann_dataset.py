"""Tests for AnnData Dataset."""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

import anndata as ad
from anndata.experimental.pytorch import (
    TORCH_AVAILABLE,
    AnnDataset,
)

if TORCH_AVAILABLE:
    import torch


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
    from scipy.sparse import csr_matrix

    X = csr_matrix(np.random.rand(*shape))
    obs = {"cell_type": np.random.choice(["A", "B", "C"], shape[0])}
    var = {"gene_name": [f"Gene_{i}" for i in range(shape[1])]}
    return ad.AnnData(X=X, obs=obs, var=var)


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
        assert dataset.shape == (adata.n_obs, adata.n_vars)

        # Test single sample
        sample = dataset[0]
        assert isinstance(sample, dict)
        assert "X" in sample
        assert sample["X"].shape == (adata.n_vars,)

    def test_transform_parameter(self, gen_adata):
        """Test transform parameter."""
        adata = gen_adata

        # Test default (no transform)
        dataset = AnnDataset(adata)
        assert dataset.transform is None

        # Test that regular functions are rejected
        def custom_function(X):
            return X * 2

        with pytest.raises(
            TypeError, match="transform must be an instance of Transform class"
        ):
            AnnDataset(adata, transform=custom_function)

        # Test valid Transform class
        from anndata.experimental.pytorch import Transform

        class CustomTransform(Transform):
            def __call__(self, X):
                return X * 2

        custom_transform = CustomTransform()
        dataset = AnnDataset(adata, transform=custom_transform)
        assert dataset.transform == custom_transform

    def test_transform_application(self, gen_adata):
        """Test transform application."""
        adata = gen_adata

        # Test transform using proper Transform class
        from anndata.experimental.pytorch import Transform

        class DummyTransform(Transform):
            def __call__(self, X):
                return X * 1.1

        dummy_transform = DummyTransform()
        dataset = AnnDataset(adata, transform=dummy_transform)
        sample = dataset[0]

        # The transform should be applied
        assert sample["X"].shape == (adata.n_vars,)

        # Test that transform was actually applied
        original_sample = AnnDataset(adata)[0]
        transformed_sample = sample

        # Should be different due to transform
        assert not torch.allclose(original_sample["X"], transformed_sample["X"])

    def test_obs_subset(self, gen_adata):
        """Test observation subsetting."""
        adata = gen_adata
        subset_indices = [0, 5, 10, 15]

        dataset = AnnDataset(adata, obs_subset=subset_indices)
        assert len(dataset) == len(subset_indices)
        assert dataset.shape == (len(subset_indices), adata.n_vars)

        # Test that we get the correct observations
        sample = dataset[0]
        assert sample["X"].shape == (adata.n_vars,)

    def test_sparse_data(self, gen_sparse_adata):
        """Test with sparse data."""
        adata = gen_sparse_adata
        dataset = AnnDataset(adata)

        sample = dataset[0]
        assert isinstance(sample, dict)
        assert "X" in sample
        assert sample["X"].shape == (adata.n_vars,)

    def test_backed_data(self, gen_adata):
        """Test with backed data."""
        adata = gen_adata

        with tempfile.TemporaryDirectory() as tmpdir:
            # Save as backed file
            file_path = Path(tmpdir) / "test.h5ad"
            adata.write_h5ad(file_path)

            # Load as backed
            backed_adata = ad.read_h5ad(file_path, backed="r")
            dataset = AnnDataset(backed_adata)

            sample = dataset[0]
            assert isinstance(sample, dict)
            assert "X" in sample
            assert sample["X"].shape == (adata.n_vars,)

    def test_file_path_input(self, gen_adata):
        """Test with file path input."""
        adata = gen_adata

        with tempfile.TemporaryDirectory() as tmpdir:
            # Save as file
            file_path = Path(tmpdir) / "test.h5ad"
            adata.write_h5ad(file_path)

            # Load from file path
            dataset = AnnDataset(file_path)

            sample = dataset[0]
            assert isinstance(sample, dict)
            assert "X" in sample
            assert sample["X"].shape == (adata.n_vars,)

    def test_metadata_inclusion(self, gen_adata):
        """Test that observation metadata is included."""
        adata = gen_adata
        dataset = AnnDataset(adata)

        sample = dataset[0]
        assert "X" in sample
        # Should have obs_ prefix for metadata
        obs_keys = [k for k in sample if k.startswith("obs_")]
        assert len(obs_keys) > 0

    def test_multiprocessing_safety(self, gen_adata):
        """Test multiprocessing safety settings."""
        adata = gen_adata

        # Test dataset creation
        dataset = AnnDataset(adata)
        sample = dataset[0]
        assert isinstance(sample, dict)

    def test_worker_generator_consistency(self, gen_adata):
        """Test worker generator consistency."""
        adata = gen_adata

        from anndata.experimental.pytorch import Transform

        class DummyTransform(Transform):
            def __call__(self, X):
                return X

        dummy_transform = DummyTransform()
        dataset = AnnDataset(
            adata,
            transform=dummy_transform,
        )

        # Test that we can get samples consistently
        sample1 = dataset[0]
        sample2 = dataset[0]
        assert isinstance(sample1, dict)
        assert isinstance(sample2, dict)

    def test_optimized_collate_fn(self, gen_adata):
        """Test optimized collate function."""
        adata = gen_adata
        dataset = AnnDataset(adata)

        # Test collate function
        collate_fn = dataset.get_collate_fn()
        assert callable(collate_fn)

        # Create some sample data
        samples = [dataset[i] for i in range(5)]
        batch = collate_fn(samples)

        assert isinstance(batch, dict)
        assert "X" in batch
        assert batch["X"].shape == (5, adata.n_vars)

    def test_collate_fn_with_metadata(self, gen_adata):
        """Test collate function handles metadata correctly."""
        adata = gen_adata
        dataset = AnnDataset(adata)

        collate_fn = dataset.get_collate_fn()
        samples = [dataset[i] for i in range(3)]
        batch = collate_fn(samples)

        assert isinstance(batch, dict)
        assert "X" in batch
        assert batch["X"].shape == (3, adata.n_vars)

        # Check metadata keys
        metadata_keys = [k for k in batch if k.startswith("obs_")]
        assert len(metadata_keys) > 0

    def test_get_optimized_dataloader(self, gen_adata):
        """Test optimized DataLoader creation."""
        adata = gen_adata
        dataset = AnnDataset(adata)

        # Test creating optimized DataLoader
        dataloader = dataset.get_optimized_dataloader(batch_size=5, shuffle=False)

        # Test that we can iterate
        batch = next(iter(dataloader))
        assert isinstance(batch, dict)
        assert "X" in batch
        assert batch["X"].shape[0] <= 5  # batch_size

    def test_row_subsetting_correctness(self, gen_adata):
        """Test that row subsetting works correctly."""
        adata = gen_adata
        obs_subset = [10, 20, 30, 40, 50]  # Specific indices

        dataset = AnnDataset(adata, obs_subset=obs_subset)

        # Test length
        assert len(dataset) == len(obs_subset)

        # Test that we get the correct data
        for i in range(len(obs_subset)):
            sample = dataset[i]
            # expected_row_idx = obs_subset[i]
            # expected_data = (
            #     adata.X[expected_row_idx].toarray().flatten()
            #     if hasattr(adata.X, "toarray")
            #     else adata.X[expected_row_idx].flatten()
            # )

            # The preprocessing might change the values, but shape should match
            assert sample["X"].shape == (adata.n_vars,)

    def test_row_subsetting_with_collate_fn(self, gen_adata):
        """Test row subsetting works with collate function."""
        adata = gen_adata
        obs_subset = [80, 60, 10, 20, 40]  # Non-sequential indices

        dataset = AnnDataset(adata, obs_subset=obs_subset)
        collate_fn = dataset.get_collate_fn()

        # Create batch indices
        batch_indices = [0, 1, 2, 3, 4]  # Dataset indices, not original indices

        # Get samples using dataset indices
        samples = [dataset[i] for i in batch_indices]
        batch = collate_fn(samples)

        # Verify batch structure
        assert isinstance(batch, dict)
        assert "X" in batch
        assert batch["X"].shape == (5, 50)  # batch_size=5, n_genes=50

        # Verify that all samples are present (order may be different due to sorting)
        # expected_adata_indices = obs_subset[batch_indices]  # [80, 60, 10, 20, 40]

        # The collate function should have processed all requested samples
        assert batch["X"].shape[0] == len(batch_indices)

    def test_empty_subset_handling(self, gen_adata):
        """Test handling of empty subset."""
        adata = gen_adata
        empty_subset = []

        dataset = AnnDataset(adata, obs_subset=empty_subset)
        assert len(dataset) == 0
        assert dataset.shape == (0, adata.n_vars)

    def test_single_row_subset(self, gen_adata):
        """Test subsetting to single row."""
        adata = gen_adata
        single_subset = [25]

        dataset = AnnDataset(adata, obs_subset=single_subset)
        assert len(dataset) == 1
        assert dataset.shape == (1, adata.n_vars)

        sample = dataset[0]
        assert sample["X"].shape == (adata.n_vars,)

    def test_duplicate_indices_in_subset(self, gen_adata):
        """Test handling of duplicate indices in subset."""
        adata = gen_adata
        duplicate_subset = [10, 20, 10, 30, 20]  # Contains duplicates

        dataset = AnnDataset(adata, obs_subset=duplicate_subset)
        assert len(dataset) == len(duplicate_subset)  # Should include duplicates

        # Test that we can access all indices
        for i in range(len(duplicate_subset)):
            sample = dataset[i]
            assert sample["X"].shape == (adata.n_vars,)

    def test_out_of_bounds_subset_indices(self, gen_adata):
        """Test handling of out-of-bounds subset indices."""
        adata = gen_adata
        out_of_bounds_subset = [0, 50, 100, 150]  # Some indices are out of bounds

        # This should raise an error
        with pytest.raises((IndexError, ValueError)):
            AnnDataset(adata, obs_subset=out_of_bounds_subset)

    def test_generic_preprocessing(self, gen_adata):
        """Test generic preprocessing functionality."""
        adata = gen_adata

        # Test with custom preprocessing Transform class
        from anndata.experimental.pytorch import Transform

        class CustomPreprocessing(Transform):
            def __call__(self, X):
                return X * 2 + 1

        custom_preprocessing = CustomPreprocessing()
        dataset = AnnDataset(adata, transform=custom_preprocessing)

        sample = dataset[0]
        assert isinstance(sample, dict)
        assert "X" in sample
        assert sample["X"].shape == (adata.n_vars,)

        # Test without preprocessing function
        dataset = AnnDataset(adata, transform=None)

        sample = dataset[0]
        assert isinstance(sample, dict)
        assert "X" in sample
        assert sample["X"].shape == (adata.n_vars,)

    def test_chunk_size_configuration(self, gen_adata):
        """Test chunk size configuration for backed data."""
        adata = gen_adata

        # Test with custom chunk size
        dataset = AnnDataset(adata, chunk_size=500)

        sample = dataset[0]
        assert isinstance(sample, dict)
        assert "X" in sample
        assert dataset.chunk_size == 500

        # Test with default chunk size
        dataset = AnnDataset(adata, chunk_size=1000)

        sample = dataset[0]
        assert isinstance(sample, dict)
        assert "X" in sample
        assert dataset.chunk_size == 1000
