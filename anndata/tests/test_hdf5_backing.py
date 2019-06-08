import joblib
import pytest
import numpy as np
from scipy import sparse

import anndata as ad


# -------------------------------------------------------------------------------
# Some test data
# -------------------------------------------------------------------------------

@pytest.fixture
def adata():
    X_list = [    # data matrix of shape n_obs x n_vars
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
    ]
    X = np.array(X_list)
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
    return ad.AnnData(
        X,
        obs=obs_dict,
        var=var_dict,
        uns=uns_dict,
        obsm={"o1": np.zeros((X.shape[0], 10))},
        varm={"v1": np.ones((X.shape[1], 20))},
        layers={
            "float": X.astype(float),
            "sparse": sparse.csr_matrix(X)
        },
        dtype='int32'
    )
# -------------------------------------------------------------------------------
# The test functions
# -------------------------------------------------------------------------------


# this is very similar to the views test
def test_backing(adata, tmp_path, backing_h5ad):
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
    subset_hash = joblib.hash(adata_subset)

    # cannot set view in backing mode...
    with pytest.raises(ValueError):
        adata_subset.obs['foo'] = range(2)
    with pytest.raises(ValueError):
        adata_subset.var['bar'] = -12
    with pytest.raises(ValueError):
        adata_subset.obsm['o2'] = np.ones((2, 2))
    with pytest.raises(ValueError):
        adata_subset.varm['v2'] = np.zeros((2, 2))
    with pytest.raises(ValueError):
        adata_subset.layers['float2'] = adata_subset.layers['float'].copy()

    # Things should stay the same after failed operations
    assert subset_hash == joblib.hash(adata_subset)
    assert adata_subset.isview

    # need to copy first
    adata_subset = adata_subset.copy(tmp_path / 'test.subset.h5ad')
    # now transition to actual object
    adata_subset.obs['foo'] = range(2)
    assert not adata_subset.isview
    assert adata_subset.obs['foo'].tolist() == list(range(2))

    # save
    adata_subset.write()


def test_double_index(adata, backing_h5ad):
    adata.filename = backing_h5ad
    with pytest.raises(ValueError):
        # no view of view of backed object currently
        adata[:2][:, 0]

    # close backing file
    adata.write()


def test_return_to_memory_mode(adata, backing_h5ad):
    bdata = adata.copy()
    adata.filename = backing_h5ad
    assert adata.isbacked

    adata.filename = None
    assert not adata.isbacked

    # make sure the previous file had been properly closed
    # when setting `adata.filename = None`
    # if it hadn't the following line would throw an error
    bdata.filename = backing_h5ad
    # close the file
    bdata.filename = None
