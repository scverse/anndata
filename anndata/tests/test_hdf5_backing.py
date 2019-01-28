from pathlib import Path

import pytest
import numpy as np

import anndata as ad


# -------------------------------------------------------------------------------
# Some test data
# -------------------------------------------------------------------------------

X_list = [    # data matrix of shape n_obs x n_vars
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9],
]

obs_dict = dict(  # annotation of observations / rows
    row_names=['name1', 'name2', 'name3'],  # row annotation
    oanno1=['cat1', 'cat2', 'cat2'],        # categorical annotation
    oanno2=['o1', 'o2', 'o3'],              # string annotation
    oanno3=[2.1, 2.2, 2.3],                 # float annotation
)

var_dict = dict(  # annotation of variables / columns
    vanno1=[3.1, 3.2, 3.3],
)

uns_dict = dict(  # unstructured annotation
    oanno1_colors=['#000000', '#FFFFFF'],
    uns2=['some annotation'],
)


# -------------------------------------------------------------------------------
# The test functions
# -------------------------------------------------------------------------------

# this is very similar to the views test
def test_backing(backing_h5ad):
    X = np.array(X_list)
    adata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict, dtype='int32')
    assert not adata.isbacked

    adata.filename = backing_h5ad
    adata.write()
    assert not adata.file.isopen
    assert adata.isbacked
    assert adata[:, 0].isview
    assert adata[:, 0].X.tolist() == [1, 4, 7]
    # this might give us a trouble as the user might not
    # know that the file is open again....
    assert adata.file.isopen

    adata[:2, 0].X = [0, 0]
    assert adata[:, 0].X.tolist() == [0, 0, 7]

    adata_subset = adata[:2, [0, 1]]
    assert adata_subset.isview
    with pytest.raises(ValueError):
        # cannot set view in backing mode...
        adata_subset.obs['foo'] = range(2)

    # need to copy first
    copy = Path('./test.subset.h5ad')
    adata_subset = adata_subset.copy(copy)
    copy.unlink()
    # now transition to actual object
    adata_subset.obs['foo'] = range(2)
    assert not adata_subset.isview
    assert adata_subset.obs['foo'].tolist() == list(range(2))

    # save
    adata_subset.write()


def test_double_index(backing_h5ad):
    X = np.array(X_list)
    adata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict, dtype='int32')
    adata.filename = backing_h5ad
    with pytest.raises(ValueError):
        # no view of view of backed object currently
        adata[:2][:, 0]

    # close backing file
    adata.write()


def test_return_to_memory_mode(backing_h5ad):
    X = np.array(X_list)
    adata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict, dtype='int32')
    adata.filename = backing_h5ad
    assert adata.isbacked

    adata.filename = None
    assert not adata.isbacked

    bdata = ad.AnnData(X, obs=obs_dict, var=var_dict, uns=uns_dict, dtype='int32')
    # make sure the previous file had been properly closed
    # when setting `adata.filename = None`
    # if it hadn't the following line would throw an error
    bdata.filename = backing_h5ad
    # close the file
    bdata.filename = None
