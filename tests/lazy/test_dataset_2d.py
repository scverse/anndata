from __future__ import annotations

from importlib.util import find_spec

import numpy as np
import pandas as pd
import pytest

from anndata._core.xarray import Dataset2D
from anndata.compat import XDataArray

pytestmark = pytest.mark.skipif(not find_spec("xarray"), reason="xarray not installed")


@pytest.fixture
def dataset_2d():
    return Dataset2D(
        {"foo": ("obs_names", pd.array(["a", "b", "c"], dtype="category"))},
        coords={"obs_names": [1, 2, 3]},
    )


def test_dataset_2d_set_array_different_dims(dataset_2d):
    da = XDataArray(
        np.arange(3),
        coords=[np.arange(3)],
        dims="not_obs_names",
    )
    with pytest.raises(ValueError, match="dimension obs_names, found not_obs_names"):
        dataset_2d["foo"] = da


def test_dataset_2d_set_array_multiple_dims(dataset_2d):
    da = XDataArray(
        np.arange(9).reshape(3, 3),
        coords={"obs_names": np.arange(3), "not_obs_names": np.arange(3)},
        dims=("obs_names", "not_obs_names"),
    )
    with pytest.raises(ValueError, match="should have only one dimension"):
        dataset_2d["foo"] = da


def test_dataset_2d_set_dataarray(dataset_2d):
    da = XDataArray(
        np.arange(3), coords={"obs_names": [1, 2, 3]}, dims=("obs_names"), name="bar"
    )
    dataset_2d["bar"] = da
    assert dataset_2d["bar"].dims == ("obs_names",)
    assert dataset_2d["bar"].equals(da)


def test_dataset_2d_set_dataarray_bad_coord_name(dataset_2d):
    da = XDataArray(
        np.arange(3),
        coords={"foo": ("obs_names", np.arange(3))},
        dims=("obs_names"),
        name="bar",
    )
    with pytest.raises(ValueError, match="DataArray should have coordinate obs_names"):
        dataset_2d["bar"] = da


def test_dataset_2d_set_dataarray_bad_name(dataset_2d):
    da = XDataArray(
        np.arange(3),
        coords={"obs_names": np.arange(3)},
        dims=("obs_names"),
        name="not_bar",
    )
    with pytest.raises(
        ValueError, match="DataArray should have name bar, found not_bar"
    ):
        dataset_2d["bar"] = da


def test_dataset_2d_set_extension_array(dataset_2d):
    array = pd.array(["e", "f", "g"], dtype="category")
    dataset_2d["bar"] = array
    assert dataset_2d["bar"].dims == ("obs_names",)
    assert dataset_2d["bar"].data is array


@pytest.mark.parametrize(
    "data", [np.arange(3), XDataArray(np.arange(3), dims="obs_names", name="obs_names")]
)
def test_dataset_2d_set_index(data, dataset_2d):
    with pytest.raises(
        KeyError,
        match="Cannot set obs_names as a variable",
    ):
        dataset_2d["obs_names"] = data
