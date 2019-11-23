import numpy as np
import pandas as pd

from anndata import AnnData


def test_uns_color_subset():
    # Tests for https://github.com/theislab/anndata/issues/257
    # Tests when the length of the color vector doesn't match with the number of categories
    obs = pd.DataFrame(index=[f"cell{i}" for i in range(5)])
    obs["cat"] = pd.Series(list("aabcd"), index=obs.index, dtype="category")
    uns = {"cat_colors": ["red", "green", "blue"]}

    adata = AnnData(np.ones((5, 5)), obs=obs, uns=uns)

    v = adata[:, [0, 1]]
    assert "cat_colors" not in v.uns

    adata.uns["cat_colors"] = ["red", "green", "blue", "yellow"]
    v = adata[[0, 1], :]
    assert len(v.uns["cat_colors"]) == 1
    assert v.uns["cat_colors"][0] == "red"
