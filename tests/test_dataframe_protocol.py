"""Tests for the DataFrameLike protocol."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import pytest

from anndata._types import DataFrameLike, DataFrameLikeIlocIndexer

if TYPE_CHECKING:
    from typing import Any, Literal, Self


class MockIlocIndexer:
    """Mock iloc indexer for testing."""

    def __init__(self, data: pd.DataFrame):
        self._data = data

    def __getitem__(self, idx: Any) -> MockDataFrame:
        result = self._data.iloc[idx]
        if isinstance(result, pd.DataFrame):
            return MockDataFrame(result)
        # For single row selection, wrap in DataFrame
        return MockDataFrame(pd.DataFrame([result]))


class MockDataFrame:
    """A minimal DataFrame-like class for testing the protocol."""

    def __init__(self, data: pd.DataFrame):
        self._data = data

    @property
    def index(self) -> pd.Index:
        return self._data.index

    @property
    def columns(self) -> pd.Index:
        return self._data.columns

    @property
    def shape(self) -> tuple[int, int]:
        return self._data.shape

    @property
    def iloc(self) -> MockIlocIndexer:
        return MockIlocIndexer(self._data)

    def reindex(
        self,
        index: pd.Index | None = None,
        axis: Literal[0] = 0,
        fill_value: Any = np.nan,
    ) -> Self:
        # axis=0 is the default; don't pass it when index is specified
        # since pandas doesn't allow both keyword arguments together
        return MockDataFrame(self._data.reindex(index=index, fill_value=fill_value))


class TestDataFrameLikeProtocol:
    """Test the DataFrameLike protocol with different implementations."""

    @pytest.fixture
    def sample_df(self) -> pd.DataFrame:
        """Create a sample pandas DataFrame for testing."""
        return pd.DataFrame(
            {"a": [1, 2, 3], "b": [4.0, 5.0, 6.0], "c": ["x", "y", "z"]},
            index=["row1", "row2", "row3"],
        )

    def test_pandas_dataframe_is_dataframe_like(self, sample_df: pd.DataFrame):
        """pd.DataFrame should satisfy the DataFrameLike protocol."""
        assert isinstance(sample_df, DataFrameLike)

    def test_pandas_iloc_is_iloc_indexer(self, sample_df: pd.DataFrame):
        """pd.DataFrame.iloc should satisfy the DataFrameLikeIlocIndexer protocol."""
        assert isinstance(sample_df.iloc, DataFrameLikeIlocIndexer)

    def test_mock_dataframe_is_dataframe_like(self, sample_df: pd.DataFrame):
        """MockDataFrame should satisfy the DataFrameLike protocol."""
        mock_df = MockDataFrame(sample_df)
        assert isinstance(mock_df, DataFrameLike)

    def test_mock_iloc_is_iloc_indexer(self, sample_df: pd.DataFrame):
        """MockDataFrame.iloc should satisfy the DataFrameLikeIlocIndexer protocol."""
        mock_df = MockDataFrame(sample_df)
        assert isinstance(mock_df.iloc, DataFrameLikeIlocIndexer)

    def test_dataframe_like_has_required_properties(self, sample_df: pd.DataFrame):
        """Verify DataFrameLike objects have the required properties."""
        for df in [sample_df, MockDataFrame(sample_df)]:
            assert hasattr(df, "index")
            assert hasattr(df, "columns")
            assert hasattr(df, "shape")
            assert hasattr(df, "iloc")
            assert hasattr(df, "reindex")

    def test_dataframe_like_index(self, sample_df: pd.DataFrame):
        """Verify index property returns a pd.Index."""
        mock_df = MockDataFrame(sample_df)
        assert isinstance(mock_df.index, pd.Index)
        pd.testing.assert_index_equal(mock_df.index, sample_df.index)

    def test_dataframe_like_columns(self, sample_df: pd.DataFrame):
        """Verify columns property returns a pd.Index."""
        mock_df = MockDataFrame(sample_df)
        assert isinstance(mock_df.columns, pd.Index)
        pd.testing.assert_index_equal(mock_df.columns, sample_df.columns)

    def test_dataframe_like_shape(self, sample_df: pd.DataFrame):
        """Verify shape property returns correct tuple."""
        mock_df = MockDataFrame(sample_df)
        assert mock_df.shape == sample_df.shape
        assert mock_df.shape == (3, 3)

    def test_dataframe_like_iloc(self, sample_df: pd.DataFrame):
        """Verify iloc indexer works correctly."""
        mock_df = MockDataFrame(sample_df)

        # Test single row selection
        result = mock_df.iloc[0]
        assert isinstance(result, DataFrameLike)

        # Test slice selection
        result = mock_df.iloc[0:2]
        assert isinstance(result, DataFrameLike)
        assert result.shape[0] == 2

    def test_dataframe_like_reindex(self, sample_df: pd.DataFrame):
        """Verify reindex method works correctly."""
        mock_df = MockDataFrame(sample_df)
        new_index = pd.Index(["row1", "row2", "row4"])

        result = mock_df.reindex(index=new_index, fill_value=-1)
        assert isinstance(result, DataFrameLike)
        pd.testing.assert_index_equal(result.index, new_index)

    def test_non_dataframe_is_not_dataframe_like(self):
        """Objects that don't implement the protocol should not match."""
        assert not isinstance([], DataFrameLike)
        assert not isinstance({}, DataFrameLike)
        assert not isinstance("string", DataFrameLike)
        assert not isinstance(42, DataFrameLike)
        assert not isinstance(np.array([1, 2, 3]), DataFrameLike)


@pytest.mark.usefixtures("xr_available")
class TestDataset2DIsDataFrameLike:
    """Test that Dataset2D satisfies the DataFrameLike protocol."""

    @pytest.fixture
    def xr_available(self):
        """Skip tests if xarray is not available."""
        pytest.importorskip("xarray")

    @pytest.fixture
    def sample_dataset2d(self):
        """Create a sample Dataset2D for testing."""
        from anndata._core.xarray import Dataset2D
        from anndata.compat import XDataset

        ds = XDataset(
            {
                "a": (["idx"], [1, 2, 3]),
                "b": (["idx"], [4.0, 5.0, 6.0]),
            },
            coords={"idx": ["row1", "row2", "row3"]},
        )
        return Dataset2D(ds)

    def test_dataset2d_is_dataframe_like(self, sample_dataset2d):
        """Dataset2D should satisfy the DataFrameLike protocol."""
        assert isinstance(sample_dataset2d, DataFrameLike)

    def test_dataset2d_iloc_is_iloc_indexer(self, sample_dataset2d):
        """Dataset2D.iloc should satisfy the DataFrameLikeIlocIndexer protocol."""
        assert isinstance(sample_dataset2d.iloc, DataFrameLikeIlocIndexer)

    def test_dataset2d_has_required_properties(self, sample_dataset2d):
        """Verify Dataset2D has the required properties."""
        assert hasattr(sample_dataset2d, "index")
        assert hasattr(sample_dataset2d, "columns")
        assert hasattr(sample_dataset2d, "shape")
        assert hasattr(sample_dataset2d, "iloc")
        assert hasattr(sample_dataset2d, "reindex")

    def test_dataset2d_properties_return_correct_types(self, sample_dataset2d):
        """Verify Dataset2D properties return correct types."""
        assert isinstance(sample_dataset2d.index, pd.Index)
        assert isinstance(sample_dataset2d.columns, pd.Index)
        assert isinstance(sample_dataset2d.shape, tuple)
        assert len(sample_dataset2d.shape) == 2


class TestDataFrameLikeWithAnnData:
    """Test that DataFrameLike protocol works correctly with AnnData objects."""

    @pytest.fixture
    def simple_adata(self):
        """Create a simple AnnData object for testing."""
        import anndata as ad

        return ad.AnnData(
            X=np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
            obs=pd.DataFrame(
                {"cell_type": ["A", "B", "C"]},
                index=["cell1", "cell2", "cell3"],
            ),
            var=pd.DataFrame(
                {"gene_name": ["g1", "g2", "g3"]},
                index=["gene1", "gene2", "gene3"],
            ),
        )

    def test_adata_obs_is_dataframe_like(self, simple_adata):
        """AnnData.obs should satisfy the DataFrameLike protocol."""
        assert isinstance(simple_adata.obs, DataFrameLike)

    def test_adata_var_is_dataframe_like(self, simple_adata):
        """AnnData.var should satisfy the DataFrameLike protocol."""
        assert isinstance(simple_adata.var, DataFrameLike)

    def test_adata_obs_has_required_properties(self, simple_adata):
        """Verify AnnData.obs has all required DataFrameLike properties."""
        obs = simple_adata.obs
        assert hasattr(obs, "index")
        assert hasattr(obs, "columns")
        assert hasattr(obs, "shape")
        assert hasattr(obs, "iloc")
        assert hasattr(obs, "reindex")

    def test_adata_obs_iloc_subsetting(self, simple_adata):
        """Verify iloc subsetting works on AnnData.obs."""
        obs = simple_adata.obs
        subset = obs.iloc[0:2]
        assert isinstance(subset, DataFrameLike)
        assert subset.shape[0] == 2

    def test_adata_subset_preserves_dataframe_like(self, simple_adata):
        """Verify subsetting AnnData preserves DataFrameLike for obs/var."""
        adata_subset = simple_adata[0:2, 0:2]
        assert isinstance(adata_subset.obs, DataFrameLike)
        assert isinstance(adata_subset.var, DataFrameLike)
        assert adata_subset.obs.shape[0] == 2
        assert adata_subset.var.shape[0] == 2

    def test_adata_copy_preserves_dataframe_like(self, simple_adata):
        """Verify copying AnnData preserves DataFrameLike for obs/var."""
        adata_copy = simple_adata.copy()
        assert isinstance(adata_copy.obs, DataFrameLike)
        assert isinstance(adata_copy.var, DataFrameLike)

    def test_set_obs_with_dataframe(self, simple_adata):
        """Verify setting obs with pd.DataFrame works."""
        new_obs = pd.DataFrame(
            {"new_col": [1, 2, 3]},
            index=["cell1", "cell2", "cell3"],
        )
        simple_adata.obs = new_obs
        assert isinstance(simple_adata.obs, DataFrameLike)
        assert "new_col" in simple_adata.obs.columns

    def test_set_var_with_dataframe(self, simple_adata):
        """Verify setting var with pd.DataFrame works."""
        new_var = pd.DataFrame(
            {"new_col": [1, 2, 3]},
            index=["gene1", "gene2", "gene3"],
        )
        simple_adata.var = new_var
        assert isinstance(simple_adata.var, DataFrameLike)
        assert "new_col" in simple_adata.var.columns


class TestCustomDataFrameLikeWithAnnData:
    """Test that custom DataFrameLike implementations work with AnnData."""

    def test_init_adata_with_custom_dataframe_like_obs(self):
        """Verify AnnData can be initialized with a custom DataFrameLike obs."""
        import anndata as ad

        # Create a custom DataFrameLike object
        obs_df = pd.DataFrame(
            {"cell_type": ["A", "B", "C"]},
            index=["cell1", "cell2", "cell3"],
        )
        mock_obs = MockDataFrame(obs_df)

        # Verify MockDataFrame satisfies the protocol
        assert isinstance(mock_obs, DataFrameLike)

        # Create AnnData with custom DataFrameLike
        adata = ad.AnnData(
            X=np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
            obs=mock_obs,
            var=pd.DataFrame(
                {"gene_name": ["g1", "g2", "g3"]},
                index=["gene1", "gene2", "gene3"],
            ),
        )

        # Verify obs is the MockDataFrame (unchanged)
        assert isinstance(adata.obs, DataFrameLike)
        assert adata.obs.shape[0] == 3

    def test_init_adata_with_custom_dataframe_like_var(self):
        """Verify AnnData can be initialized with a custom DataFrameLike var."""
        import anndata as ad

        # Create a custom DataFrameLike object
        var_df = pd.DataFrame(
            {"gene_name": ["g1", "g2", "g3"]},
            index=["gene1", "gene2", "gene3"],
        )
        mock_var = MockDataFrame(var_df)

        # Verify MockDataFrame satisfies the protocol
        assert isinstance(mock_var, DataFrameLike)

        # Create AnnData with custom DataFrameLike var
        adata = ad.AnnData(
            X=np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
            obs=pd.DataFrame(
                {"cell_type": ["A", "B", "C"]},
                index=["cell1", "cell2", "cell3"],
            ),
            var=mock_var,
        )

        # Verify var is the MockDataFrame (unchanged)
        assert isinstance(adata.var, DataFrameLike)
        assert adata.var.shape[0] == 3

    def test_custom_dataframe_like_length_validation(self):
        """Verify length validation works for custom DataFrameLike."""
        import anndata as ad

        # Create a custom DataFrameLike with wrong length
        obs_df = pd.DataFrame(
            {"cell_type": ["A", "B"]},  # Only 2 rows, but X has 3
            index=["cell1", "cell2"],
        )
        mock_obs = MockDataFrame(obs_df)

        # Should raise ValueError due to length mismatch
        with pytest.raises(ValueError, match="must have as many rows"):
            ad.AnnData(
                X=np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
                obs=mock_obs,
                var=pd.DataFrame(
                    {"gene_name": ["g1", "g2", "g3"]},
                    index=["gene1", "gene2", "gene3"],
                ),
            )


@pytest.mark.usefixtures("xr_available")
class TestDataset2DWithAnnData:
    """Test that Dataset2D works correctly as obs/var in AnnData."""

    @pytest.fixture
    def xr_available(self):
        """Skip tests if xarray is not available."""
        pytest.importorskip("xarray")

    @pytest.fixture
    def adata_with_dataset2d_obs(self):
        """Create an AnnData with Dataset2D obs."""
        import anndata as ad
        from anndata._core.xarray import Dataset2D
        from anndata.compat import XDataset

        X = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

        # Create Dataset2D for obs
        obs_ds = XDataset(
            {"cell_type": (["idx"], ["A", "B", "C"])},
            coords={"idx": ["cell1", "cell2", "cell3"]},
        )
        obs = Dataset2D(obs_ds)

        var = pd.DataFrame(
            {"gene_name": ["g1", "g2", "g3"]},
            index=["gene1", "gene2", "gene3"],
        )

        return ad.AnnData(X=X, obs=obs, var=var)

    def test_adata_with_dataset2d_obs_is_dataframe_like(self, adata_with_dataset2d_obs):
        """AnnData with Dataset2D obs should satisfy DataFrameLike."""
        assert isinstance(adata_with_dataset2d_obs.obs, DataFrameLike)

    def test_adata_with_dataset2d_subset(self, adata_with_dataset2d_obs):
        """Subsetting AnnData with Dataset2D obs should work."""
        adata_subset = adata_with_dataset2d_obs[0:2]
        assert isinstance(adata_subset.obs, DataFrameLike)
        assert adata_subset.obs.shape[0] == 2

    def test_adata_with_dataset2d_obs_index(self, adata_with_dataset2d_obs):
        """Dataset2D obs should have correct index."""
        obs = adata_with_dataset2d_obs.obs
        pd.testing.assert_index_equal(
            obs.index, pd.Index(["cell1", "cell2", "cell3"]), check_names=False
        )
