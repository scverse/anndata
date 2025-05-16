from __future__ import annotations

from importlib.util import find_spec

import numpy as np
import pandas as pd
import pytest

from anndata._core.xarray import Dataset2D

pytestmark = pytest.mark.skipif(not find_spec("xarray"), reason="xarray not installed")


def test_dataset_2d_set_array_different_dims():
    import xarray as xr

    ds = Dataset2D(
        {"foo": ("obs_names", pd.array(["a", "b", "c"], dtype="category"))},
        coords={"obs_names": [1, 2, 3]},
    )
    da = xr.DataArray(
        np.arange(3),
        coords=[np.arange(3)],
        dims="not_obs_names",
    )
    with pytest.raises(ValueError, match="dimension obs_names, found not_obs_names"):
        ds["foo"] = da


def test_dataset_2d_set_extension_array():
    ds = Dataset2D(
        {"foo": ("obs_names", pd.array(["a", "b", "c"], dtype="category"))},
        coords={"obs_names": [1, 2, 3]},
    )
    array = pd.array(["e", "f", "g"], dtype="category")
    ds["foo"] = array
    assert ds["foo"].dims == ("obs_names",)
    assert ds["foo"].data is array
