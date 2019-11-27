import numpy as np
import pytest

import anndata as ad
from anndata.core.anndata import ImplicitModificationWarning


# -------------------------------------------------------------------------------
# Some test data
# -------------------------------------------------------------------------------

X_list = [
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9],
]  # data matrix of shape n_obs x n_vars

obs_dict = dict(  # annotation of observations / rows
    row_names=['name1', 'name2', 'name3'],  # row annotation
    oanno1=['cat1', 'cat2', 'cat2'],  # categorical annotation
    oanno2=['o1', 'o2', 'o3'],  # string annotation
    oanno3=[2.1, 2.2, 2.3],  # float annotation
)

var_dict = dict(  # annotation of variables / columns
    col_names=['var1', 'var2', 'var3'], vanno1=[3.1, 3.2, 3.3]
)

uns_dict = dict(  # unstructured annotation
    oanno1_colors=['#000000', '#FFFFFF'], uns2=['some annotation']
)


# -------------------------------------------------------------------------------
# The test functions
# -------------------------------------------------------------------------------


def test_raw(backing_h5ad):
    X = np.array(X_list)
    adata = ad.AnnData(
        X, obs=obs_dict, var=var_dict, uns=uns_dict, dtype='int32'
    )

    # init raw
    adata.raw = adata

    assert adata.raw[:, 0].X.tolist() == [[1], [4], [7]]

    adata = adata[:, [0, 1]]

    assert adata.var_names.tolist() == ['var1', 'var2']
    assert adata.raw.var_names.tolist() == ['var1', 'var2', 'var3']

    # read write
    with pytest.warns(
        ImplicitModificationWarning, match="Initializing view as actual"
    ):  # TODO: don’t modify adata just to write it
        adata.write(backing_h5ad)
    adata = ad.read(backing_h5ad)

    assert adata.raw[:, 0].X.tolist() == [[1], [4], [7]]
    assert adata.raw.var_names.tolist() == ['var1', 'var2', 'var3']
    assert adata.var_names.tolist() == ['var1', 'var2']

    del adata.raw
    assert adata.raw is None
