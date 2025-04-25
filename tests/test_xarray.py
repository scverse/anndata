from __future__ import annotations

import string

import numpy as np
import pandas as pd
import pytest

from anndata._core.xarray import Dataset2D
from anndata.tests.helpers import gen_typed_df


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
