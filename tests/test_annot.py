"""Test handling of values in `obs`/ `var`"""

from __future__ import annotations

from contextlib import nullcontext

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
from natsort import natsorted

import anndata as ad
from anndata._warnings import ImplicitModificationWarning
from anndata.tests.helpers import assert_equal, get_multiindex_columns_df


@pytest.mark.parametrize("dtype", [object, "string"])
def test_str_to_categorical(dtype):
    obs = pd.DataFrame(
        {"str": ["a", "a", None, "b", "b"]}, index=[f"cell-{i}" for i in range(5)]
    )
    obs["str"] = obs["str"].astype(dtype)
    a = ad.AnnData(obs=obs.copy())

    a.strings_to_categoricals()
    expected = obs["str"].astype(pd.CategoricalDtype(pd.Index(["a", "b"], dtype=dtype)))
    pd.testing.assert_series_equal(expected, a.obs["str"])


@pytest.mark.parametrize("dtype", [object, "string"])
def test_to_categorical_ordering(dtype):
    obs = pd.DataFrame(
        {"str": ["10", "11", "3", "9", "10", "10"]},
        index=[f"cell-{i}" for i in range(6)],
    )
    obs["str"] = obs["str"].astype(dtype)
    a = ad.AnnData(obs=obs.copy())

    a.strings_to_categoricals()

    expected = obs["str"].astype(
        pd.CategoricalDtype(categories=natsorted(obs["str"].unique()))
    )
    pd.testing.assert_series_equal(expected, a.obs["str"])


def test_non_str_to_not_categorical():
    # Test case based on https://github.com/scverse/anndata/issues/141#issuecomment-802105259
    obs = pd.DataFrame(index=[f"cell-{i}" for i in range(5)]).assign(
        str_with_nan=["foo", "bar", None, np.nan, "foo"],
        boolean_with_nan_and_none=[True, False, np.nan, None, True],
        boolean_with_nan=[True, False, np.nan, np.nan, True],
        boolean_with_none=[True, False, None, None, True],
    )
    adata = ad.AnnData(obs=obs.copy())

    orig_dtypes = {k: v.name for k, v in obs.dtypes.items()}
    expected_dtypes = orig_dtypes.copy()
    expected_dtypes["str_with_nan"] = "category"

    adata.strings_to_categoricals()
    result_dtypes = {k: v.name for k, v in adata.obs.dtypes.items()}

    assert expected_dtypes == result_dtypes

    expected_non_transformed = obs.drop(columns=["str_with_nan"])
    result_non_transformed = adata.obs.drop(columns=["str_with_nan"])

    pd.testing.assert_frame_equal(expected_non_transformed, result_non_transformed)


def test_error_col_multiindex():
    adata = ad.AnnData(np.random.rand(100, 10))
    df = get_multiindex_columns_df((adata.shape[0], 20))
    with pytest.raises(ValueError, match=r"MultiIndex columns are not supported"):
        adata.obs = df


@pytest.mark.parametrize(
    "restrict_index_types", [True, False], ids=["error", "no_error"]
)
def test_row_multiindex(restrict_index_types):
    ad.settings.restrict_index_types = restrict_index_types
    df = pd.DataFrame(
        {"x": [1, 2, 3]},
        index=pd.MultiIndex.from_tuples([("a", 1), ("b", 2), ("c", 3)]),
    )
    with (
        pytest.raises(
            ValueError, match=r"pandas.MultiIndex not supported as index for obs or var"
        )
        if restrict_index_types
        else nullcontext()
    ):
        adata = ad.AnnData(df)
        assert_equal(adata[0], adata[pd.MultiIndex.from_tuples([("a", 1)])])


@pytest.mark.parametrize(
    "restrict_index_types", [True, False], ids=["error", "no_error"]
)
def test_row_multiindex_setter(restrict_index_types):
    ad.settings.restrict_index_types = restrict_index_types
    df = pd.DataFrame(
        {"x": [1, 2, 3]},
        index=pd.MultiIndex.from_tuples([("a", 1), ("b", 2), ("c", 3)]),
    )
    adata = ad.AnnData(np.random.rand(3, 10))
    with (
        pytest.raises(
            ValueError, match=r"pandas.MultiIndex not supported as index for obs or var"
        )
        if restrict_index_types
        else nullcontext()
    ):
        adata.obs = df
        assert_equal(adata[0], adata[pd.MultiIndex.from_tuples([("a", 1)])])


@pytest.mark.parametrize(
    "restrict_index_types", [True, False], ids=["restrict", "unrestricted"]
)
def test_arrow_index(restrict_index_types):
    ad.settings.restrict_index_types = restrict_index_types
    df = pd.DataFrame(
        {"x": [1, 2, 3]},
        index=pd.Index(
            pd.arrays.ArrowExtensionArray(
                pa.array([
                    {"x": 1, "y": True},
                    {"x": 2, "y": False},
                    {"x": 3, "y": False},
                ])
            )
        ),
    )
    with (
        pytest.warns(ImplicitModificationWarning, match=r"Transforming to str index.")
        if restrict_index_types
        else nullcontext()
    ):
        adata = ad.AnnData(np.random.rand(3, 10), obs=df)
    # If we restricted the type, the index is now a string
    with (
        pytest.raises(TypeError, match=r"unhashable type: 'dict'")
        if restrict_index_types
        else nullcontext()
    ):
        assert_equal(
            adata[pd.arrays.ArrowExtensionArray(pa.array([{"x": 1, "y": True}]))],
            adata[0],
        )
