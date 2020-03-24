import pandas as pd

from anndata.utils import make_index_unique


def test_make_index_unique():
    index = pd.Index(["val", "val", "val-1", "val-1"])
    result = make_index_unique(index)
    expected = pd.Index(["val", "val-2", "val-1", "val-1-1"])
    assert all(left == right for left, right in zip(result, expected))
