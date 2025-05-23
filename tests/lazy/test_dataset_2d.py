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


@pytest.mark.parametrize(
    "setter",
    [
        pd.array(["e", "f", "g"], dtype="category"),
        ("obs_names", pd.array(["e", "f", "g"], dtype="category")),
    ],
    ids=["array", "tuple_with_array"],
)
def test_dataset_2d_set_extension_array(dataset_2d, setter):
    dataset_2d["bar"] = setter
    assert dataset_2d["bar"].dims == ("obs_names",)
    assert dataset_2d["bar"].data is setter[1] if isinstance(setter, tuple) else setter


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
        pytest.param(
            ("not_obs_names", [1, 2, 3]),
            "Setting value tuple should have first entry",
            id="tuple_bad_dim",
        ),
        pytest.param(
            (("not_obs_names",), [1, 2, 3]),
            "Dimension tuple should have only",
            id="nested_tuple_bad_dim",
        ),
        pytest.param(
            (("obs_names", "bar"), [1, 2, 3]),
            "Dimension tuple is too long",
            id="nested_tuple_too_long",
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
