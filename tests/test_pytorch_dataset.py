"""Tests for PyTorch Dataset integration."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import anndata as ad

# Import PyTorch components conditionally
try:
    import torch
    from torch.utils.data import DataLoader

    from anndata.experimental.pytorch import AnnDataset, Compose, Transform

    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False

pytestmark = pytest.mark.skipif(not TORCH_AVAILABLE, reason="PyTorch not available")


class SimpleTransform(Transform):
    """Simple transform for testing that just passes data through."""

    def __call__(self, data):
        return data


class DictTransform(Transform):
    """Dict-based transform that modifies both X and metadata."""

    def __call__(self, data_dict):
        if not isinstance(data_dict, dict):
            msg = "This transform requires a dict"
            raise TypeError(msg)

        # Modify X
        data_dict["X"] = data_dict["X"] * 2

        # Add new metadata field
        data_dict["obs_transformed"] = torch.tensor(1.0)

        return data_dict


class LabelEncodingTransform(Transform):
    """Transform that encodes string labels to integers."""

    def __init__(self, obs_key="cell_type", mapping=None):
        self.obs_key = obs_key
        self.mapping = mapping or {"A": 0, "B": 1, "C": 2}

    def __call__(self, data_dict):
        # Access metadata using obs_{column_name} key pattern
        obs_key = f"obs_{self.obs_key}"
        obs_encoded_key = f"{obs_key}_encoded"

        if obs_key in data_dict and isinstance(data_dict[obs_key], str):
            encoded = self.mapping.get(data_dict[obs_key], -1)
            data_dict[obs_encoded_key] = torch.tensor(encoded, dtype=torch.long)
        return data_dict


@pytest.fixture
def simple_adata():
    """Create simple AnnData for testing."""
    X = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], dtype=np.float32)
    obs = pd.DataFrame({
        "cell_type": ["A", "B", "A"],
        "n_genes": [100, 200, 150],
        "batch": [1, 2, 1],
        "quality": [0.8, 0.9, 0.7],
    })
    var = pd.DataFrame(index=["gene1", "gene2"])
    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def sparse_adata():
    """Create sparse AnnData for testing."""
    X = sparse.csr_matrix([[1.0, 0.0], [0.0, 4.0], [5.0, 0.0]], dtype=np.float32)
    obs = pd.DataFrame({"cell_type": ["A", "B", "C"], "n_genes": [50, 100, 75]})
    return ad.AnnData(X=X, obs=obs)


class TestAnnDatasetBasic:
    """Test basic AnnDataset functionality."""

    def test_creation(self, simple_adata):
        """Test dataset creation."""
        dataset = AnnDataset(simple_adata)
        assert len(dataset) == 3
        assert dataset.shape == (3, 2)

    def test_getitem_no_transform(self, simple_adata):
        """Test __getitem__ without transforms."""
        dataset = AnnDataset(simple_adata)
        item = dataset[0]

        assert isinstance(item, dict)
        assert "X" in item
        assert torch.allclose(item["X"], torch.tensor([1.0, 2.0]))

        # Check metadata
        assert item["obs_cell_type"] == "A"
        assert item["obs_n_genes"] == torch.tensor(100)
        assert item["obs_batch"] == torch.tensor(1)
        assert torch.allclose(item["obs_quality"], torch.tensor(0.8))

    def test_obs_subset(self, simple_adata):
        """Test observation subsetting."""
        dataset = AnnDataset(simple_adata, obs_subset=[0, 2])
        assert len(dataset) == 2

        item0 = dataset[0]
        item1 = dataset[1]

        assert item0["obs_cell_type"] == "A"
        assert item1["obs_cell_type"] == "A"
        assert torch.allclose(item0["X"], torch.tensor([1.0, 2.0]))
        assert torch.allclose(item1["X"], torch.tensor([5.0, 6.0]))

    def test_sparse_data(self, sparse_adata):
        """Test with sparse data."""
        dataset = AnnDataset(sparse_adata)
        item = dataset[0]

        assert torch.allclose(item["X"], torch.tensor([1.0, 0.0]))
        assert item["obs_cell_type"] == "A"


class TestTransforms:
    """Test transform functionality."""

    def test_no_transform(self, simple_adata):
        """Test dataset without transforms."""
        dataset = AnnDataset(simple_adata, transform=None)
        item = dataset[0]
        assert torch.allclose(item["X"], torch.tensor([1.0, 2.0]))

    def test_simple_transform(self, simple_adata):
        """Test simple pass-through transform."""
        dataset = AnnDataset(simple_adata, transform=SimpleTransform())
        item = dataset[0]
        assert torch.allclose(item["X"], torch.tensor([1.0, 2.0]))

    def test_dict_transform(self, simple_adata):
        """Test dict-based transform."""
        dataset = AnnDataset(simple_adata, transform=DictTransform())
        item = dataset[0]

        # X should be doubled
        assert torch.allclose(item["X"], torch.tensor([2.0, 4.0]))

        # Should have new metadata field
        assert "obs_transformed" in item
        assert item["obs_transformed"] == torch.tensor(1.0)

        # Original metadata should be preserved
        assert item["obs_cell_type"] == "A"

    def test_label_encoding_transform(self, simple_adata):
        """Test label encoding within transforms."""
        transform = LabelEncodingTransform(obs_key="cell_type")
        dataset = AnnDataset(simple_adata, transform=transform)

        item = dataset[0]
        assert item["obs_cell_type"] == "A"
        assert item["obs_cell_type_encoded"] == torch.tensor(0)

        item = dataset[1]
        assert item["obs_cell_type"] == "B"
        assert item["obs_cell_type_encoded"] == torch.tensor(1)

    def test_compose_transforms(self, simple_adata):
        """Test composing multiple transforms."""
        transform = Compose([
            DictTransform(),  # Doubles X, adds obs_transformed
            LabelEncodingTransform(obs_key="cell_type"),  # Adds encoded labels
        ])
        dataset = AnnDataset(simple_adata, transform=transform)

        item = dataset[0]

        # X should be doubled
        assert torch.allclose(item["X"], torch.tensor([2.0, 4.0]))

        # Should have both new fields
        assert item["obs_transformed"] == torch.tensor(1.0)
        assert item["obs_cell_type_encoded"] == torch.tensor(0)

        # Original metadata preserved
        assert item["obs_cell_type"] == "A"


class TestDataTypes:
    """Test handling of different data types."""

    def test_various_obs_types(self):
        """Test various observation data types."""
        X = np.array([[1.0, 2.0]], dtype=np.float32)
        obs = pd.DataFrame({
            "int_col": [42],
            "float_col": [3.14],
            "bool_col": [True],
            "str_col": ["test"],
            "category_col": pd.Categorical(["cat1"]),
            "nullable_int": pd.array([10], dtype="Int64"),
        })
        adata = ad.AnnData(X=X, obs=obs)
        dataset = AnnDataset(adata)
        item = dataset[0]

        # Check type conversions
        assert item["obs_int_col"] == torch.tensor(42)
        assert torch.allclose(item["obs_float_col"], torch.tensor(3.14))
        assert item["obs_bool_col"] == torch.tensor(True)  # noqa: FBT003
        assert item["obs_str_col"] == "test"
        assert item["obs_category_col"] == "cat1"  # Should be converted to string
        assert item["obs_nullable_int"] == torch.tensor(10.0)  # np.number -> float

    def test_missing_values(self):
        """Test handling of missing values."""
        X = np.array([[1.0, 2.0]], dtype=np.float32)
        obs = pd.DataFrame({
            "with_nan": [np.nan],
            "with_none": [None],
        })
        adata = ad.AnnData(X=X, obs=obs)
        dataset = AnnDataset(adata)
        item = dataset[0]

        # NaN should be converted to string
        assert item["obs_with_nan"] == "nan"
        assert item["obs_with_none"] == "None"


class TestDataLoader:
    """Test DataLoader integration."""

    def test_standard_dataloader(self, simple_adata):
        """Test standard DataLoader."""
        dataset = AnnDataset(simple_adata)
        dataloader = DataLoader(dataset, batch_size=2, shuffle=False)

        batch = next(iter(dataloader))

        assert "X" in batch
        assert batch["X"].shape == (2, 2)  # batch_size=2, n_vars=2

        # Check metadata batching
        assert "obs_cell_type" in batch
        assert batch["obs_cell_type"] == ["A", "B"]

        assert "obs_n_genes" in batch
        assert torch.allclose(batch["obs_n_genes"], torch.tensor([100, 200]))

    def test_optimized_dataloader(self, simple_adata):
        """Test optimized DataLoader."""
        dataset = AnnDataset(simple_adata)
        dataloader = dataset.get_optimized_dataloader(batch_size=2, shuffle=False)

        batch = next(iter(dataloader))

        assert "X" in batch
        assert batch["X"].shape == (2, 2)

        # Should have same structure as standard dataloader
        assert "obs_cell_type" in batch
        assert batch["obs_cell_type"] == ["A", "B"]

    def test_dataloader_with_transforms(self, simple_adata):
        """Test DataLoader with transforms."""
        transform = DictTransform()
        dataset = AnnDataset(simple_adata, transform=transform)
        dataloader = DataLoader(dataset, batch_size=2, shuffle=False)

        batch = next(iter(dataloader))

        # X should be doubled by transform
        expected_X = torch.tensor([[2.0, 4.0], [6.0, 8.0]])
        assert torch.allclose(batch["X"], expected_X)

        # Should have transform-added metadata
        assert "obs_transformed" in batch
        assert torch.allclose(batch["obs_transformed"], torch.tensor([1.0, 1.0]))


class TestErrorHandling:
    """Test error handling and edge cases."""

    def test_invalid_transform_type(self, simple_adata):
        """Test error when transform is not Transform instance."""
        with pytest.raises(
            TypeError, match="transform must be an instance of Transform"
        ):
            AnnDataset(simple_adata, transform=lambda x: x)

    def test_invalid_obs_subset(self, simple_adata):
        """Test error with invalid observation subset."""
        with pytest.raises(ValueError, match="obs_subset contains indices"):
            AnnDataset(simple_adata, obs_subset=[0, 1, 10])  # Index 10 doesn't exist

    def test_negative_obs_subset(self, simple_adata):
        """Test error with negative indices."""
        with pytest.raises(ValueError, match="obs_subset indices must be non-negative"):
            AnnDataset(simple_adata, obs_subset=[-1, 0])

    def test_invalid_chunk_size(self, simple_adata):
        """Test error with invalid chunk size."""
        with pytest.raises(ValueError, match="chunk_size must be positive"):
            AnnDataset(simple_adata, chunk_size=0)


class TestFileLoading:
    """Test loading from file paths."""

    def test_load_from_path(self, simple_adata, tmp_path):
        """Test loading dataset from file path."""
        file_path = tmp_path / "test.h5ad"
        simple_adata.write(file_path)

        dataset = AnnDataset(str(file_path))
        assert len(dataset) == 3

        item = dataset[0]
        assert torch.allclose(item["X"], torch.tensor([1.0, 2.0]))
        assert item["obs_cell_type"] == "A"

    def test_load_from_path_with_transform(self, simple_adata, tmp_path):
        """Test loading from path with transforms."""
        file_path = tmp_path / "test.h5ad"
        simple_adata.write(file_path)

        dataset = AnnDataset(str(file_path), transform=DictTransform())
        item = dataset[0]

        # Transform should work on file-loaded data
        assert torch.allclose(item["X"], torch.tensor([2.0, 4.0]))
        assert item["obs_transformed"] == torch.tensor(1.0)


if __name__ == "__main__":
    pytest.main([__file__])
