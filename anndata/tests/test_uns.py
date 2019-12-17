import numpy as np
import pandas as pd

from anndata import AnnData


def test_uns_color_subset():
    # Tests for https://github.com/theislab/anndata/issues/257
    obs = pd.DataFrame(index=[f"cell{i}" for i in range(5)])
    obs["cat1"] = pd.Series(list("aabcd"), index=obs.index, dtype="category")
    obs["cat2"] = pd.Series(list("aabbb"), index=obs.index, dtype="category")
    uns = dict(
        cat1_colors=["red", "green", "blue"], cat2_colors=["red", "green", "blue"],
    )

    adata = AnnData(np.ones((5, 5)), obs=obs, uns=uns)

    # If number of categories does not match number of colors,
    # they should be reset
    v = adata[:, [0, 1]]
    assert "cat1_colors" not in v.uns
    assert "cat2_colors" not in v.uns

    # Otherwise the colors should still match after reseting
    adata.uns["cat1_colors"] = ["red", "green", "blue", "yellow"]
    v = adata[[0, 1], :]
    assert len(v.uns["cat1_colors"]) == 1
    assert v.uns["cat1_colors"][0] == "red"
