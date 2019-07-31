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


@pytest.fixture(params=[sparse.csr_matrix, sparse.csc_matrix])
def sparse_format(request):
    return request.param

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
    assert adata[:, 0].X.tolist() == np.reshape([1, 4, 7], (3, 1)).tolist()
    # this might give us a trouble as the user might not
    # know that the file is open again....
    assert adata.file.isopen

    adata[:2, 0].X = [0, 0]
    assert adata[:, 0].X.tolist() == np.reshape([0, 0, 7], (3, 1)).tolist()

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


def test_backed_modification(adata, backing_h5ad):
    adata.X[:, 1] = 0  # Make it a little sparse
    adata.X = sparse.csr_matrix(adata.X)
    assert not adata.isbacked

    # While this currently makes the file backed, it doesn't write it as sparse
    adata.filename = backing_h5ad
    adata.write()
    assert not adata.file.isopen
    assert adata.isbacked

    adata.X[0, [0, 2]] = 10
    adata.X[1, [0, 2]] = [11, 12]
    adata.X[2, 1] = 13  # If it were written as sparse, this should fail

    assert adata.isbacked

    assert np.all(adata.X[0, :] == np.array([10, 0, 10]))
    assert np.all(adata.X[1, :] == np.array([11, 0, 12]))
    assert np.all(adata.X[2, :] == np.array([7, 13, 9]))


def test_backed_modification_sparse(adata, backing_h5ad, sparse_format):
    adata.X[:, 1] = 0  # Make it a little sparse
    adata.X = sparse_format(adata.X)
    assert not adata.isbacked

    adata.write(backing_h5ad)
    adata = ad.read_h5ad(backing_h5ad, backed="r+")

    assert adata.filename == backing_h5ad
    assert adata.isbacked

    adata.X[0, [0, 2]] = 10
    adata.X[1, [0, 2]] = [11, 12]
    with pytest.raises(ValueError):
        adata.X[2, 1] = 13

    assert adata.isbacked

    assert np.all(adata.X[0, :] == np.array([10, 0, 10]))
    assert np.all(adata.X[1, :] == np.array([11, 0, 12]))
    assert np.all(adata.X[2, :] == np.array([7, 0, 9]))


# TODO: Work around h5py not supporting this
# def test_backed_view_modification(adata, backing_h5ad):
#     adata.write(backing_h5ad)
#     backed_adata = ad.read_h5ad(backing_h5ad, backed=True)

#     backed_view = backed_adata[[1, 2], :]
#     backed_view.X = 0

#     assert np.all(backed_adata.X[:3, :] == 0)


# TODO: Implement
# def test_backed_view_modification_sparse(adata, backing_h5ad, sparse_format):
#     adata[:, 1] = 0  # Make it a little sparse
#     adata.X = sparse_format(adata.X)
#     adata.write(backing_h5ad)
#     backed_adata = ad.read_h5ad(backing_h5ad, backed=True)

#     backed_view = backed_adata[[1,2], :]
#     backed_view.X = 0
#     assert np.all(backed_adata.X[[1,2], :] == 0)
