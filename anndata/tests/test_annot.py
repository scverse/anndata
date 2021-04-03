"""Test handling of values in `obs`/ `var`"""
import numpy as np
import pandas as pd

import anndata as ad

import pytest


@pytest.mark.parametrize("dtype", [object, "string"])
def test_str_to_categorical(dtype):
    obs = pd.DataFrame(
        {"str": ["a", "a", None, "b", "b"]}, index=[f"cell-{i}" for i in range(5)]
    )
    obs["str"] = obs["str"].astype(dtype)
    a = ad.AnnData(obs=obs.copy())

    a.strings_to_categoricals()
    expected = obs["str"].astype("category")
    pd.testing.assert_series_equal(expected, a.obs["str"])


def test_non_str_to_not_categorical():
    # Test case based on https://github.com/theislab/anndata/issues/141#issuecomment-802105259
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
