"""Tests for AnnDataset batch collation fixes.

This module tests the fixes for:
- Custom transform keys being properly handled in batch collation
- _batch_samples function handling all keys, not just obs_ prefixed ones
- EncodeLabels transform working with proper obs key prefixes
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest
from sklearn.preprocessing import LabelEncoder

import anndata as ad
from anndata.experimental.pytorch import (
    TORCH_AVAILABLE,
    AnnDataset,
)

if TORCH_AVAILABLE:
    import torch
    from torch.utils.data import DataLoader

    from anndata.experimental.pytorch import (
        Compose,
        Transform,
    )
    from anndata.experimental.pytorch._ann_dataset import _batch_samples

    # Test transforms
    class EncodeLabels(Transform):
        """Encode string labels to integers using proper obs key prefix."""

        def __init__(self, obs_key: str, encoded_key: str, label_encoder: LabelEncoder):
            self.obs_key = obs_key
            self.encoded_key = encoded_key
            self.label_encoder = label_encoder

        def __call__(self, data_dict):
            # Get string label from metadata and encode it
            # Note: AnnDataset prefixes obs keys with "obs_"
            prefixed_obs_key = f"obs_{self.obs_key}"
            string_label = data_dict[prefixed_obs_key]
            encoded_label = self.label_encoder.transform([string_label])[0]
            data_dict[self.encoded_key] = torch.tensor(encoded_label, dtype=torch.long)
            return data_dict

        def __repr__(self):
            return f"EncodeLabels(obs_key='{self.obs_key}', encoded_key='{self.encoded_key}')"

    class AddCustomKey(Transform):
        """Transform that adds custom keys without obs_ prefix."""

        def __init__(self, custom_key: str, value: float):
            self.custom_key = custom_key
            self.value = value

        def __call__(self, data_dict):
            data_dict[self.custom_key] = torch.tensor(self.value)
            return data_dict

    class MultipleCustomKeys(Transform):
        """Transform that adds multiple custom keys."""

        def __call__(self, data_dict):
            data_dict["custom_tensor"] = torch.tensor([1.0, 2.0, 3.0])
            data_dict["custom_scalar"] = torch.tensor(42.0)
            data_dict["custom_string"] = "test_value"
            data_dict["custom_int"] = torch.tensor(123, dtype=torch.long)
            return data_dict


@pytest.fixture
def gen_adata_with_labels(shape=(100, 50)):
    """Generate test AnnData object with categorical labels."""
    X = np.random.rand(*shape)
    cell_types = np.random.choice(["B cells", "T cells", "NK cells"], shape[0])
    obs = {
        "cell_type": cell_types,
        "numeric_meta": np.random.rand(shape[0]),
        "bool_meta": np.random.choice([True, False], shape[0]),
    }
    var = {"gene_name": [f"Gene_{i}" for i in range(shape[1])]}
    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.mark.skipif(
    not TORCH_AVAILABLE, reason="PyTorch is required for batch collation tests"
)
class TestBatchCollationFixes:
    """Test fixes for batch collation handling custom transform keys."""

    def test_batch_samples_handles_all_keys(self):
        """Test that _batch_samples handles all keys, not just obs_ prefixed ones."""
        # Create sample data with mixed key types
        samples = [
            {
                "X": torch.tensor([1.0, 2.0]),
                "obs_cell_type": "A",
                "obs_numeric": torch.tensor(1.5),
                "custom_key": torch.tensor(10.0),
                "encoded_label": torch.tensor(0, dtype=torch.long),
                "string_key": "test",
            },
            {
                "X": torch.tensor([3.0, 4.0]),
                "obs_cell_type": "B",
                "obs_numeric": torch.tensor(2.5),
                "custom_key": torch.tensor(20.0),
                "encoded_label": torch.tensor(1, dtype=torch.long),
                "string_key": "test2",
            },
        ]

        batch = _batch_samples(samples)

        # Check that all keys are present
        expected_keys = {
            "X",
            "obs_cell_type",
            "obs_numeric",
            "custom_key",
            "encoded_label",
            "string_key",
        }
        assert set(batch.keys()) == expected_keys

        # Check tensor keys are properly batched
        assert torch.allclose(batch["X"], torch.tensor([[1.0, 2.0], [3.0, 4.0]]))
        assert torch.allclose(batch["obs_numeric"], torch.tensor([1.5, 2.5]))
        assert torch.allclose(batch["custom_key"], torch.tensor([10.0, 20.0]))
        assert torch.allclose(batch["encoded_label"], torch.tensor([0, 1]))

        # Check string keys are stored as lists
        assert batch["obs_cell_type"] == ["A", "B"]
        assert batch["string_key"] == ["test", "test2"]

    def test_batch_samples_handles_missing_keys(self):
        """Test that _batch_samples handles samples with missing keys."""
        samples = [
            {
                "X": torch.tensor([1.0, 2.0]),
                "obs_cell_type": "A",
                "custom_key": torch.tensor(10.0),
            },
            {
                "X": torch.tensor([3.0, 4.0]),
                "obs_cell_type": "B",
                # Missing custom_key
            },
            {
                "X": torch.tensor([5.0, 6.0]),
                "obs_cell_type": "C",
                "custom_key": torch.tensor(30.0),
                "extra_key": torch.tensor(100.0),  # Only in this sample
            },
        ]

        batch = _batch_samples(samples)

        # Check that all keys from all samples are present
        expected_keys = {"X", "obs_cell_type", "custom_key", "extra_key"}
        assert set(batch.keys()) == expected_keys

        # Check X is properly batched
        expected_X = torch.tensor([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
        assert torch.allclose(batch["X"], expected_X)

        # Check string keys
        assert batch["obs_cell_type"] == ["A", "B", "C"]

        # Check keys with missing values - should use default tensor for missing values
        expected_custom = torch.tensor([10.0, 0.0, 30.0])  # 0.0 is default for missing
        assert torch.allclose(batch["custom_key"], expected_custom)

        expected_extra = torch.tensor([0.0, 0.0, 100.0])  # 0.0 is default for missing
        assert torch.allclose(batch["extra_key"], expected_extra)

    def test_encode_labels_transform_with_batch_loading(self, gen_adata_with_labels):
        """Test EncodeLabels transform works correctly with batch loading."""
        adata = gen_adata_with_labels

        # Create and fit label encoder
        label_encoder = LabelEncoder()
        label_encoder.fit(adata.obs["cell_type"].values)

        # Create transform
        encode_transform = EncodeLabels(
            obs_key="cell_type",
            encoded_key="cell_type_encoded",
            label_encoder=label_encoder,
        )

        # Create dataset with transform
        dataset = AnnDataset(adata, transform=encode_transform)

        # Test single item access
        item = dataset[0]
        assert "cell_type_encoded" in item
        assert isinstance(item["cell_type_encoded"], torch.Tensor)
        assert item["cell_type_encoded"].dtype == torch.long

        # Test batch loading
        dataloader = DataLoader(dataset, batch_size=8, shuffle=False)
        batch = next(iter(dataloader))

        # Check that encoded labels are in the batch
        assert "cell_type_encoded" in batch
        assert isinstance(batch["cell_type_encoded"], torch.Tensor)
        assert batch["cell_type_encoded"].shape == (8,)
        assert batch["cell_type_encoded"].dtype == torch.long

        # Verify encoding is correct
        original_labels = [adata.obs["cell_type"].iloc[i] for i in range(8)]
        expected_encoded = torch.tensor(label_encoder.transform(original_labels))
        assert torch.equal(batch["cell_type_encoded"], expected_encoded)

    def test_multiple_custom_transforms_batch_loading(self, gen_adata_with_labels):
        """Test multiple custom transforms with batch loading."""
        adata = gen_adata_with_labels

        # Create label encoder
        label_encoder = LabelEncoder()
        label_encoder.fit(adata.obs["cell_type"].values)

        # Create multiple transforms
        transforms = Compose([
            EncodeLabels(
                obs_key="cell_type",
                encoded_key="cell_type_encoded",
                label_encoder=label_encoder,
            ),
            AddCustomKey("custom_value", 3.14),
            MultipleCustomKeys(),
        ])

        dataset = AnnDataset(adata, transform=transforms)

        # Test single item
        item = dataset[0]
        expected_keys = {
            "X",
            "obs_cell_type",
            "obs_numeric_meta",
            "obs_bool_meta",
            "cell_type_encoded",
            "custom_value",
            "custom_tensor",
            "custom_scalar",
            "custom_string",
            "custom_int",
        }
        assert set(item.keys()) == expected_keys

        # Test batch loading
        dataloader = DataLoader(dataset, batch_size=4, shuffle=False)
        batch = next(iter(dataloader))

        # All keys should be present in batch
        assert set(batch.keys()) == expected_keys

        # Check tensor keys are properly batched
        assert batch["cell_type_encoded"].shape == (4,)
        assert torch.allclose(batch["custom_value"], torch.tensor([3.14] * 4))
        assert batch["custom_tensor"].shape == (4, 3)
        assert torch.allclose(batch["custom_scalar"], torch.tensor([42.0] * 4))
        assert torch.allclose(batch["custom_int"], torch.tensor([123] * 4))

        # Check string keys are lists
        assert len(batch["custom_string"]) == 4
        assert all(val == "test_value" for val in batch["custom_string"])

    def test_optimized_dataloader_with_custom_keys(self, gen_adata_with_labels):
        """Test optimized dataloader works with custom transform keys."""
        adata = gen_adata_with_labels

        # Save to file for backed loading
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = Path(tmpdir) / "test.h5ad"
            adata.write_h5ad(file_path)

            # Create label encoder
            label_encoder = LabelEncoder()
            label_encoder.fit(adata.obs["cell_type"].values)

            # Create transform
            encode_transform = EncodeLabels(
                obs_key="cell_type",
                encoded_key="encoded_labels",
                label_encoder=label_encoder,
            )

            # Create dataset with file path (backed loading)
            dataset = AnnDataset(file_path, transform=encode_transform)

            # Test optimized dataloader
            dataloader = dataset.get_optimized_dataloader(
                batch_size=6, shuffle=True, num_workers=0
            )

            batch = next(iter(dataloader))

            # Check custom key is present and properly batched
            assert "encoded_labels" in batch
            assert isinstance(batch["encoded_labels"], torch.Tensor)
            assert batch["encoded_labels"].shape == (6,)
            assert batch["encoded_labels"].dtype == torch.long

    def test_mixed_data_types_in_batch(self, gen_adata_with_labels):
        """Test batch collation with mixed data types."""
        adata = gen_adata_with_labels

        class MixedDataTransform(Transform):
            def __call__(self, data_dict):
                data_dict["tensor_float"] = torch.tensor(1.5)
                data_dict["tensor_int"] = torch.tensor(42, dtype=torch.long)
                data_dict["tensor_bool"] = torch.tensor(data=True, dtype=torch.bool)
                data_dict["string_val"] = "test_string"
                data_dict["list_val"] = [1, 2, 3]  # Non-tensor value
                return data_dict

        dataset = AnnDataset(adata, transform=MixedDataTransform())
        dataloader = DataLoader(dataset, batch_size=3, shuffle=False)
        batch = next(iter(dataloader))

        # Tensor values should be batched
        assert torch.allclose(batch["tensor_float"], torch.tensor([1.5, 1.5, 1.5]))
        assert torch.equal(batch["tensor_int"], torch.tensor([42, 42, 42]))
        assert torch.equal(batch["tensor_bool"], torch.tensor([True, True, True]))

        # Non-tensor values should be lists
        assert batch["string_val"] == ["test_string", "test_string", "test_string"]
        # Lists get transposed by PyTorch's default collation behavior when using DataLoader
        assert isinstance(batch["list_val"], list)
        expected_transposed = [
            torch.tensor([1, 1, 1]),
            torch.tensor([2, 2, 2]),
            torch.tensor([3, 3, 3]),
        ]
        assert len(batch["list_val"]) == 3
        for i, tensor in enumerate(batch["list_val"]):
            assert torch.equal(tensor, expected_transposed[i])

    def test_multiworker_with_custom_keys(self, gen_adata_with_labels):
        """Test multiworker DataLoader with custom transform keys."""
        adata = gen_adata_with_labels

        # Create label encoder
        label_encoder = LabelEncoder()
        label_encoder.fit(adata.obs["cell_type"].values)

        encode_transform = EncodeLabels(
            obs_key="cell_type", encoded_key="labels", label_encoder=label_encoder
        )

        dataset = AnnDataset(adata, transform=encode_transform)

        # Test with multiple workers
        dataloader = DataLoader(dataset, batch_size=4, shuffle=False, num_workers=2)

        batch = next(iter(dataloader))

        # Custom key should be present and properly handled
        assert "labels" in batch
        assert isinstance(batch["labels"], torch.Tensor)
        assert batch["labels"].shape == (4,)

    def test_edge_case_empty_batch(self):
        """Test edge case with empty batch."""
        samples = []

        # Should handle empty samples gracefully
        batch = _batch_samples(samples)

        # Should return empty dict for empty samples
        assert isinstance(batch, dict)
        assert len(batch) == 0

    def test_edge_case_single_sample_batch(self, gen_adata_with_labels):
        """Test batch with single sample."""
        adata = gen_adata_with_labels

        transform = AddCustomKey("single_key", 99.9)
        dataset = AnnDataset(adata, transform=transform)

        # Create batch with single sample
        dataloader = DataLoader(dataset, batch_size=1, shuffle=False)
        batch = next(iter(dataloader))

        # Should handle single sample correctly
        assert "single_key" in batch
        assert torch.allclose(batch["single_key"], torch.tensor([99.9]))
        assert batch["X"].shape == (1, adata.n_vars)

    def test_backward_compatibility_obs_keys_only(self, gen_adata_with_labels):
        """Test that original obs_ prefixed keys still work correctly."""
        adata = gen_adata_with_labels

        # Transform that only uses obs_ keys (original behavior)
        class ObsOnlyTransform(Transform):
            def __call__(self, data_dict):
                # Just modify existing obs keys, don't add custom ones
                if "obs_numeric_meta" in data_dict:
                    data_dict["obs_numeric_meta"] = data_dict["obs_numeric_meta"] * 2
                return data_dict

        dataset = AnnDataset(adata, transform=ObsOnlyTransform())
        dataloader = DataLoader(dataset, batch_size=5, shuffle=False)
        batch = next(iter(dataloader))

        # Should still work as before
        assert "obs_cell_type" in batch
        assert "obs_numeric_meta" in batch
        assert "obs_bool_meta" in batch

        # Check that transformation was applied
        original_values = torch.tensor(
            adata.obs["numeric_meta"].iloc[:5].values, dtype=torch.float32
        )
        expected_values = original_values * 2
        assert torch.allclose(batch["obs_numeric_meta"], expected_values)
