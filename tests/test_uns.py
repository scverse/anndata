from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from anndata import AnnData
from anndata.tests.helpers import assert_equal


def test_uns_color_subset():
    # Tests for https://github.com/scverse/anndata/issues/257
    obs = pd.DataFrame(
        {
            "cat1": pd.Categorical(list("aabcd")),
            "cat2": pd.Categorical(list("aabbb")),
        },
        index=[f"cell{i}" for i in range(5)],
    )

    # If number of categories does not match number of colors, they should be reset
    wrong_color_length_adata = AnnData(
        np.ones((5, 5)),
        obs=obs,
        uns={
            "cat1_colors": ["red", "green", "blue"],
            "cat2_colors": ["red", "green", "blue"],
        },
    )
    v = wrong_color_length_adata[:, [0, 1]]

    assert "cat1_colors" not in v.uns
    assert "cat2_colors" not in v.uns

    # Otherwise the colors should still match after resetting
    cat1_colors = np.array(["red", "green", "blue", "yellow"], dtype=object)
    adata = AnnData(np.ones((5, 5)), obs=obs, uns={"cat1_colors": cat1_colors.copy()})

    for color, idx in [("red", [0, 1]), ("green", [2]), ("blue", [3]), ("yellow", [4])]:
        v = adata[idx, :]
        assert len(v.uns["cat1_colors"]) == 1
        assert v.uns["cat1_colors"][0] == color

        c = v.copy()
        assert_equal(v.uns, c.uns, elem_name="uns")
        with pytest.raises(AssertionError):
            assert_equal(adata.uns, c.uns, elem_name="uns")

    # But original object should not change
    assert list(adata.uns["cat1_colors"]) == list(cat1_colors)
