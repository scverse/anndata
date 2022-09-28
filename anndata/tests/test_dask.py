"""
For tests using dask
"""
import anndata as ad
import pandas as pd

import pytest

from anndata.tests.helpers import assert_equal
from anndata.compat import DaskDataFrame


pytest.importorskip("dask.array")


def test_dask_X_view():
    import dask.array as da

    M, N = 50, 30
    adata = ad.AnnData(
        obs=pd.DataFrame(index=[f"cell{i:02}" for i in range(M)]),
        var=pd.DataFrame(index=[f"gene{i:02}" for i in range(N)]),
    )
    adata.X = da.ones((M, N))
    view = adata[:30]
    view.copy()


def test_assert_equal_dask_dataframes():
    import dask.dataframe as dd

    data_1 = {"col1": [1, 2, 3, 4], "col2": [5, 6, 7, 8]}
    data_2 = {"col1": [1, 2, 3, 4], "col2": [5, 6, 7, 8]}

    a = dd.from_pandas(pd.DataFrame(data=data_1), npartitions=1)

    b = dd.from_pandas(pd.DataFrame(data=data_2), npartitions=1)

    assert_equal(a, b)
