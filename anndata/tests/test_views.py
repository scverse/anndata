import numpy as np
import anndata as ad


# -------------------------------------------------------------------------------
# Some test data
# -------------------------------------------------------------------------------

X_list = [    # data matrix of shape n_obs x n_vars
    [1, 2, 3], [4, 5, 6], [7, 8, 9]]

obs_dict = {  # annotation of observations / rows
    'row_names': ['name1', 'name2', 'name3'],  # row annotation
    'oanno1': ['cat1', 'cat2', 'cat2'],        # categorical annotation
    'oanno2': ['o1', 'o2', 'o3'],              # string annotation
    'oanno3': [2.1, 2.2, 2.3]}                 # float annotation

var_dict = {  # annotation of variables / columns
    'vanno1': [3.1, 3.2, 3.3]}

uns_dict = {  # unstructured annotation
    'oanno1_colors': ['#000000', '#FFFFFF'],
    'uns2': ['some annotation']}


# -------------------------------------------------------------------------------
# The test functions
# -------------------------------------------------------------------------------


def test_views():
    X = np.array(X_list)
    adata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict, dtype='int32')

    assert adata[:, 0].isview
    assert adata[:, 0].X.tolist() == [1, 4, 7]

    adata[:2, 0].X = [0, 0]

    assert adata[:, 0].X.tolist() == [0, 0, 7]

    adata_subset = adata[:2, [0, 1]]

    assert adata_subset.isview
    # now transition to actual object
    adata_subset.obs['foo'] = range(2)
    assert not adata_subset.isview

    assert adata_subset.obs['foo'].tolist() == list(range(2))


def test_slice_copy():
    adata = ad.AnnData(np.empty((100, 100)))
    adata.obsm['o'] = np.empty((100, 50))

    adata = adata[:50]
    adata.obsm['o'] = np.ones((50, 20))
