import pandas as pd

from anndata.utils import make_index_unique


def test_make_index_unique():
    index = pd.Index(["SNORD113", "SNORD113", "SNORD113-1"])
    result = make_index_unique(index)
    expected = pd.Index(["SNORD113", "SNORD113-2", "SNORD113-1"])
    assert all(l == r for l, r in zip(result, expected))
