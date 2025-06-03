from __future__ import annotations

import string

import numpy as np
import pandas as pd
import pytest

from anndata._core.xarray import Dataset2D
from anndata.compat import XDataArray, XDataset, XVariable
from anndata.tests.helpers import gen_typed_df

pytestmark = pytest.importorskip("xarray")


@pytest.fixture
def df():
    return gen_typed_df(10)


@pytest.fixture
def dataset2d(df):
    return Dataset2D.from_dataframe(df)


def test_shape(df, dataset2d):
    assert dataset2d.shape == df.shape


def test_columns(df, dataset2d):
    assert np.all(dataset2d.columns.sort_values() == df.columns.sort_values())


def test_to_memory(df, dataset2d):
    memory_df = dataset2d.to_memory()
    assert np.all(df == memory_df)
    assert np.all(df.index == memory_df.index)
    assert np.all(df.columns.sort_values() == memory_df.columns.sort_values())


def test_getitem(df, dataset2d):
    col = df.columns[0]
    assert np.all(dataset2d[col] == df[col])

    empty_dset = dataset2d[[]]
    assert empty_dset.shape == (df.shape[0], 0)
    assert np.all(empty_dset.index == dataset2d.index)


def test_backed_property(dataset2d):
    assert not dataset2d.is_backed

    dataset2d.is_backed = True
    assert dataset2d.is_backed

    dataset2d.is_backed = False
    assert not dataset2d.is_backed


def test_index_dim(dataset2d):
    assert dataset2d.index_dim == "index"
    assert dataset2d.true_index_dim == dataset2d.index_dim

    col = next(iter(dataset2d.keys()))
    dataset2d.true_index_dim = col
    assert dataset2d.index_dim == "index"
    assert dataset2d.true_index_dim == col

    with pytest.raises(ValueError, match=r"Unknown variable `test`\."):
        dataset2d.true_index_dim = "test"

    dataset2d.true_index_dim = None
    assert dataset2d.true_index_dim == dataset2d.index_dim


def test_index(dataset2d):
    alphabet = np.asarray(
        list(string.ascii_letters + string.digits + string.punctuation)
    )
    new_idx = pd.Index(
        [
            "".join(np.random.choice(alphabet, size=10))
            for _ in range(dataset2d.shape[0])
        ],
        name="test_index",
    )

    col = next(iter(dataset2d.keys()))
    dataset2d.true_index_dim = col

    dataset2d.index = new_idx
    assert np.all(dataset2d.index == new_idx)
    assert dataset2d.true_index_dim == dataset2d.index_dim == new_idx.name
    assert list(dataset2d.coords.keys()) == [new_idx.name]


@pytest.fixture
def dataset_2d_one_column():
    return Dataset2D(
        {"foo": ("obs_names", pd.array(["a", "b", "c"], dtype="category"))},
        coords={"obs_names": [1, 2, 3]},
    )


def test_dataset_2d_set_dataarray(dataset_2d_one_column):
    da = XDataArray(
        np.arange(3), coords={"obs_names": [1, 2, 3]}, dims=("obs_names"), name="bar"
    )
    dataset_2d_one_column["bar"] = da
    assert dataset_2d_one_column["bar"].dims == ("obs_names",)
    assert dataset_2d_one_column["bar"].equals(da)


def test_dataset_2d_set_dataset(dataset_2d_one_column):
    ds = XDataset(
        data_vars={
            "foo": ("obs_names", np.arange(3)),
            "bar": ("obs_names", np.arange(3) + 3),
        },
        coords={"obs_names": [1, 2, 3]},
    )
    key = ["foo", "bar"]
    dataset_2d_one_column[key] = ds
    assert tuple(dataset_2d_one_column[key].sizes.keys()) == ("obs_names",)
    assert dataset_2d_one_column[key].equals(ds)


@pytest.mark.parametrize(
    "setter",
    [
        pd.array(["e", "f", "g"], dtype="category"),
        ("obs_names", pd.array(["e", "f", "g"], dtype="category")),
    ],
    ids=["array", "tuple_with_array"],
)
def test_dataset_2d_set_extension_array(dataset_2d_one_column, setter):
    dataset_2d_one_column["bar"] = setter
    assert dataset_2d_one_column["bar"].dims == ("obs_names",)
    assert (
        dataset_2d_one_column["bar"].data is setter[1]
        if isinstance(setter, tuple)
        else setter
    )


@pytest.mark.parametrize(
    ("da", "pattern"),
    [
        pytest.param(
            XDataset(
                data_vars={"bar": ("obs_names", np.arange(3))},
                coords={"foo": ("obs_names", np.arange(3))},
            ),
            r"Dataset should have coordinate obs_names",
            id="coord_name_dataset",
        ),
        pytest.param(
            XDataArray(
                np.arange(3),
                coords={"foo": ("obs_names", np.arange(3))},
                dims="obs_names",
                name="bar",
            ),
            r"DataArray should have coordinate obs_names",
            id="coord_name",
        ),
        pytest.param(
            XDataArray(
                np.arange(3),
                coords={"obs_names": np.arange(3)},
                dims=("obs_names",),
                name="not_bar",
            ),
            r"DataArray should have name bar, found not_bar",
            id="dataarray_name",
        ),
        pytest.param(
            XDataset(
                data_vars={
                    "foo": (["obs_names", "not_obs_names"], np.arange(9).reshape(3, 3))
                },
                coords={"obs_names": np.arange(3), "not_obs_names": np.arange(3)},
            ),
            r"Dataset should have only one dimension",
            id="multiple_dims_dataset",
        ),
        pytest.param(
            XDataArray(
                np.arange(9).reshape(3, 3),
                coords={"obs_names": np.arange(3), "not_obs_names": np.arange(3)},
                dims=("obs_names", "not_obs_names"),
            ),
            r"DataArray should have only one dimension",
            id="multiple_dims_dataarray",
        ),
        pytest.param(
            XVariable(
                data=np.arange(9).reshape(3, 3),
                dims=("obs_names", "not_obs_names"),
            ),
            r"Variable should have only one dimension",
            id="multiple_dims_variable",
        ),
        pytest.param(
            XDataset(
                data_vars={"foo": ("other", np.arange(3))},
                coords={"obs_names": ("other", np.arange(3))},
            ),
            r"Dataset should have dimension obs_names",
            id="name_conflict_dataset",
        ),
        pytest.param(
            XVariable(
                data=np.arange(3),
                dims="not_obs_names",
            ),
            r"Variable should have dimension obs_names, found not_obs_names",
            id="name_conflict_variable",
        ),
        pytest.param(
            XDataArray(
                np.arange(3),
                coords=[np.arange(3)],
                dims="not_obs_names",
            ),
            r"DataArray should have dimension obs_names, found not_obs_names",
            id="name_conflict_dataarray",
        ),
        pytest.param(
            ("not_obs_names", [1, 2, 3]),
            r"Setting value tuple should have first entry",
            id="tuple_bad_dim",
        ),
        pytest.param(
            (("not_obs_names",), [1, 2, 3]),
            r"Dimension tuple should have only",
            id="nested_tuple_bad_dim",
        ),
        pytest.param(
            (("obs_names", "bar"), [1, 2, 3]),
            r"Dimension tuple is too long",
            id="nested_tuple_too_long",
        ),
    ],
)
def test_dataset_2d_set_with_bad_obj(da, pattern, dataset_2d_one_column):
    with pytest.raises(ValueError, match=pattern):
        dataset_2d_one_column["bar"] = da


@pytest.mark.parametrize(
    "data", [np.arange(3), XDataArray(np.arange(3), dims="obs_names", name="obs_names")]
)
def test_dataset_2d_set_index(data, dataset_2d_one_column):
    with pytest.raises(
        KeyError,
        match="Cannot set obs_names as a variable",
    ):
        dataset_2d_one_column["obs_names"] = data
