"""Tests for AnnDataset Transform classes and multiprocessing compatibility."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest


class TestTransforms:
    """Test Transform classes and multiprocessing compatibility."""

    def test_transform_base_class(self, gen_adata):
        """Test that Transform base class works correctly."""
        try:
            import torch

            from anndata.experimental.pytorch import Transform
        except ImportError:
            pytest.skip("PyTorch not available")

        class TestTransform(Transform):
            def __init__(self, multiplier=2.0):
                self.multiplier = multiplier

            def __call__(self, x):
                return x * self.multiplier

            def __repr__(self):
                return f"TestTransform(multiplier={self.multiplier})"

        # Test transform creation and usage
        transform = TestTransform(multiplier=3.0)
        x = torch.tensor([1.0, 2.0, 3.0])
        result = transform(x)

        expected = torch.tensor([3.0, 6.0, 9.0])
        assert torch.allclose(result, expected)
        assert str(transform) == "TestTransform(multiplier=3.0)"

    def test_compose_transforms(self, gen_adata):
        """Test that Compose works correctly."""
        try:
            import torch

            from anndata.experimental.pytorch import Compose, Transform
        except ImportError:
            pytest.skip("PyTorch not available")

        class MultiplyTransform(Transform):
            def __init__(self, multiplier=2.0):
                self.multiplier = multiplier

            def __call__(self, x):
                return x * self.multiplier

            def __repr__(self):
                return f"MultiplyTransform(multiplier={self.multiplier})"

        class AddTransform(Transform):
            def __init__(self, value=1.0):
                self.value = value

            def __call__(self, x):
                return x + self.value

            def __repr__(self):
                return f"AddTransform(value={self.value})"

        # Test composition
        compose = Compose([MultiplyTransform(multiplier=2.0), AddTransform(value=1.0)])

        x = torch.tensor([1.0, 2.0, 3.0])
        result = compose(x)

        # Should first multiply by 2, then add 1: [1,2,3] -> [2,4,6] -> [3,5,7]
        expected = torch.tensor([3.0, 5.0, 7.0])
        assert torch.allclose(result, expected)

        # Test repr
        repr_str = str(compose)
        assert "Compose(" in repr_str
        assert "MultiplyTransform(multiplier=2.0)" in repr_str
        assert "AddTransform(value=1.0)" in repr_str

    def test_compose_empty_transforms(self):
        """Test that Compose raises error with empty transforms."""
        try:
            from anndata.experimental.pytorch import Compose
        except ImportError:
            pytest.skip("PyTorch not available")

        with pytest.raises(ValueError, match="At least one transform must be provided"):
            Compose([])

    def test_transforms_with_anndataset(self, gen_adata):
        """Test that Transform classes work with AnnDataset."""
        try:
            import torch

            from anndata.experimental.pytorch import AnnDataset, Compose, Transform
        except ImportError:
            pytest.skip("PyTorch not available")

        class NormalizeTransform(Transform):
            def __init__(self, target_sum=1e4):
                self.target_sum = target_sum

            def __call__(self, x):
                x = torch.clamp(x, min=0)
                row_sum = torch.sum(x, dim=-1, keepdim=True) + 1e-8
                return x * (self.target_sum / row_sum)

        class LogTransform(Transform):
            def __call__(self, x):
                return torch.log1p(x)

        # Create composed transform
        transform = Compose([NormalizeTransform(target_sum=1e4), LogTransform()])

        # Create dataset with transform
        adata = gen_adata
        dataset = AnnDataset(adata, transform=transform)

        # Test that we can get samples
        sample = dataset[0]
        assert isinstance(sample, dict)
        assert "X" in sample
        assert isinstance(sample["X"], torch.Tensor)

        # Test that transform was applied (values should be log-transformed)
        assert sample["X"].min() >= 0  # log1p ensures non-negative values

    def test_multiprocessing_with_custom_transforms(self, gen_adata):
        """Test that custom Transform classes work with multiprocessing."""
        try:
            import torch
            from torch.utils.data import DataLoader

            from anndata.experimental.pytorch import AnnDataset, Transform
        except ImportError:
            pytest.skip("PyTorch not available")

        class SerializableTransform(Transform):
            """A serializable transform for testing multiprocessing."""

            def __init__(self, target_sum=1e4):
                self.target_sum = target_sum

            def __call__(self, x):
                # Simple normalization
                x = torch.clamp(x, min=0)
                row_sum = torch.sum(x, dim=-1, keepdim=True) + 1e-8
                x = x * (self.target_sum / row_sum)
                return torch.log1p(x)

            def __repr__(self):
                return f"SerializableTransform(target_sum={self.target_sum})"

        # Create dataset with custom transform
        adata = gen_adata
        transform = SerializableTransform(target_sum=1e4)
        dataset = AnnDataset(adata, transform=transform)

        # Test with multiprocessing (num_workers=2)
        # This should not raise any pickling errors
        dataloader = DataLoader(
            dataset,
            batch_size=5,
            shuffle=False,
            num_workers=2,  # This tests multiprocessing
        )

        # Test that we can iterate through the dataloader
        batches = []
        for batch in dataloader:
            batches.append(batch)
            if len(batches) >= 2:  # Just test first couple of batches
                break

        assert len(batches) > 0
        assert "X" in batches[0]
        assert isinstance(batches[0]["X"], torch.Tensor)
        assert batches[0]["X"].shape[0] <= 5  # batch_size

    def test_multiprocessing_with_backed_data_and_transforms(self, backed_adata):
        """Test multiprocessing with backed data and custom transforms."""
        try:
            import torch
            from torch.utils.data import DataLoader

            from anndata.experimental.pytorch import AnnDataset, Transform
        except ImportError:
            pytest.skip("PyTorch not available")

        class BackedTransform(Transform):
            """Transform for testing with backed data."""

            def __init__(self, scale=1e4):
                self.scale = scale

            def __call__(self, x):
                x = torch.clamp(x, min=0)
                # Simple scaling and log transform
                x = x * (self.scale / (torch.sum(x, dim=-1, keepdim=True) + 1e-8))
                return torch.log1p(x)

        # Create dataset with backed data and transform
        transform = BackedTransform(scale=1e4)
        dataset = AnnDataset(backed_adata, transform=transform, chunk_size=50)

        # Test multiprocessing with backed data
        dataloader = DataLoader(
            dataset,
            batch_size=10,
            shuffle=False,
            num_workers=2,  # Test multiprocessing
        )

        # Verify it works
        batch = next(iter(dataloader))
        assert isinstance(batch, dict)
        assert "X" in batch
        assert batch["X"].shape[0] <= 10

    def test_transform_with_file_path(self, gen_adata):
        """Test that transforms work when AnnDataset is created from file path."""
        try:
            import torch

            from anndata.experimental.pytorch import AnnDataset, Transform
        except ImportError:
            pytest.skip("PyTorch not available")

        class FileTransform(Transform):
            def __call__(self, x):
                return torch.log1p(torch.clamp(x, min=0))

        # Save data to temporary file
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
            temp_path = tmp.name

        try:
            # Write data to file
            gen_adata.write(temp_path)

            # Create dataset from file path with transform
            transform = FileTransform()
            dataset = AnnDataset(temp_path, transform=transform)

            # Test that it works
            sample = dataset[0]
            assert isinstance(sample, dict)
            assert "X" in sample
            assert isinstance(sample["X"], torch.Tensor)

        finally:
            # Cleanup
            Path(temp_path).unlink()

    def test_transform_pickling(self, gen_adata):
        """Test that Transform classes can be pickled and unpickled."""
        try:
            import pickle

            import torch

            from anndata.experimental.pytorch import Compose, Transform
        except ImportError:
            pytest.skip("PyTorch not available")

        class PicklableTransform(Transform):
            def __init__(self, param=1.0):
                self.param = param

            def __call__(self, x):
                return x * self.param

        # Test pickling individual transform
        transform = PicklableTransform(param=2.5)
        pickled = pickle.dumps(transform)
        unpickled = pickle.loads(pickled)

        assert unpickled.param == 2.5

        # Test that unpickled transform works
        x = torch.tensor([1.0, 2.0])
        result1 = transform(x)
        result2 = unpickled(x)
        assert torch.allclose(result1, result2)

        # Test pickling Compose
        compose = Compose([
            PicklableTransform(param=2.0),
            PicklableTransform(param=3.0),
        ])

        pickled_compose = pickle.dumps(compose)
        unpickled_compose = pickle.loads(pickled_compose)

        # Test that composed transforms work the same
        result1 = compose(x)
        result2 = unpickled_compose(x)
        assert torch.allclose(result1, result2)
