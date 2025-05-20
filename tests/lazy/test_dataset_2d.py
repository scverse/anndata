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


def test_dataset_2d_set_dataarray(dataset_2d):
    da = XDataArray(
        np.arange(3), coords={"obs_names": [1, 2, 3]}, dims=("obs_names"), name="bar"
    )
    dataset_2d["bar"] = da
    assert dataset_2d["bar"].dims == ("obs_names",)
    assert dataset_2d["bar"].equals(da)


def test_dataset_2d_set_extension_array(dataset_2d):
    array = pd.array(["e", "f", "g"], dtype="category")
    dataset_2d["bar"] = array
    assert dataset_2d["bar"].dims == ("obs_names",)
    assert dataset_2d["bar"].data is array


@pytest.mark.parametrize(
    ("da", "msg"),
    [
        pytest.param(
            XDataArray(
                np.arange(3),
                coords={"foo": ("obs_names", np.arange(3))},
                dims="obs_names",
                name="bar",
            ),
            "DataArray should have coordinate obs_names",
            id="coord_name",
        ),
        pytest.param(
            XDataArray(
                np.arange(3),
                coords={"obs_names": np.arange(3)},
                dims=("obs_names",),
                name="not_bar",
            ),
            "DataArray should have name bar, found not_bar",
            id="dataarray_name",
        ),
        pytest.param(
            XDataArray(
                np.arange(9).reshape(3, 3),
                coords={"obs_names": np.arange(3), "not_obs_names": np.arange(3)},
                dims=("obs_names", "not_obs_names"),
            ),
            "should have only one dimension",
            id="multiple_dims",
        ),
        pytest.param(
            XDataArray(
                np.arange(3),
                coords=[np.arange(3)],
                dims="not_obs_names",
            ),
            "dimension obs_names, found not_obs_names",
            id="name_conflict",
        ),
    ],
)
def test_dataset_2d_set_with_bad_dataarray(da, msg, dataset_2d):
    with pytest.raises(ValueError, match=msg):
        dataset_2d["bar"] = da


@pytest.mark.parametrize(
    "data", [np.arange(3), XDataArray(np.arange(3), dims="obs_names", name="obs_names")]
)
def test_dataset_2d_set_index(data, dataset_2d):
    with pytest.raises(
        KeyError,
        match="Cannot set obs_names as a variable",
    ):
        dataset_2d["obs_names"] = data


def test_dataset_2d_set_tuple(dataset_2d):
    with pytest.raises(
        TypeError,
        match="Setting with a tuple is not permitted",
    ):
        dataset_2d["foo"] = ("obs_names", [1, 2, 3])
