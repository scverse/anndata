import numpy as np
from numpy import ma
import pandas as pd
import pickle
from scipy import sparse as sp

from anndata import AnnData


def test_creation():
    AnnData(np.array([[1, 2], [3, 4]]))
    AnnData(np.array([[1, 2], [3, 4]]), {}, {})
    AnnData(ma.array([[1, 2], [3, 4]]), uns={'mask': [0, 1, 1, 0]})
    AnnData(sp.eye(2))
    X = np.array([[1, 2, 3], [4, 5, 6]])
    adata = AnnData(
        X=X,
        obs=dict(Obs=['A', 'B']),
        var=dict(Feat=['a', 'b', 'c']),
        obsm=dict(X_pca=np.array([[1, 2], [3, 4]])),
        raw=dict(X=X, var={'var_names': ['a', 'b', 'c']}))

    assert adata.raw.X.tolist() == X.tolist()
    assert adata.raw.var_names.tolist() == ['a', 'b', 'c']

    assert AnnData(np.array([1, 2])).X.shape == (2,)

    from pytest import raises
    raises(ValueError, AnnData,
           np.array([[1, 2], [3, 4]]),
           dict(TooLong=[1, 2, 3, 4]))

    # init with empty data matrix
    shape = (3, 5)
    adata = AnnData(None, uns={'test': np.array((3, 3))}, shape=shape)
    assert adata.X is None
    assert adata.shape == shape
    assert 'test' in adata.uns


def test_names():
    adata = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(obs_names=['A', 'B']),
        dict(var_names=['a', 'b', 'c']))

    assert adata.obs_names.tolist() == 'A B'.split()
    assert adata.var_names.tolist() == 'a b c'.split()

    adata = AnnData(np.array([[1, 2], [3, 4], [5, 6]]),
                    var={'var_names': ['a', 'b']})
    assert adata.var_names.tolist() == ['a', 'b']


def test_indices_dtypes():
    adata = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(obs_names=['A', 'B']),
        dict(var_names=['a', 'b', 'c']))
    adata.obs_names = ['ö', 'a']
    assert adata.obs_names.tolist() == ['ö', 'a']


def test_creation_from_vector():
    adata = AnnData(np.array([1, 2, 3]))
    assert adata.X.shape == (3,)
    adata = AnnData(np.array([[1], [2], [3]]))
    assert adata.X.shape == (3,)


def test_slicing():
    adata = AnnData(np.array([[1, 2, 3],
                              [4, 5, 6]]))

    assert np.all(adata[:, 0].X == adata.X[:, 0])

    assert adata[0, 0].X.tolist() == 1
    assert adata[0, :].X.tolist() == [1, 2, 3]
    assert adata[:, 0].X.tolist() == [1, 4]

    assert adata[:, [0, 1]].X.tolist() == [[1, 2], [4, 5]]
    assert adata[:, np.array([0, 2])].X.tolist() == [[1, 3], [4, 6]]
    assert adata[:, np.array([False, True, True])].X.tolist() == [[2, 3], [5, 6]]
    assert adata[:, 1:3].X.tolist() == [[2, 3], [5, 6]]


def test_slicing_strings():
    adata = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(obs_names=['A', 'B']),
        dict(var_names=['a', 'b', 'c']))

    assert adata['A', 'a'].X.tolist() == 1
    assert adata['A', :].X.tolist() == [1, 2, 3]
    assert adata[:, 'a'].X.tolist() == [1, 4]
    assert adata[:, ['a', 'b']].X.tolist() == [[1, 2], [4, 5]]
    assert adata[:, np.array(['a', 'c'])].X.tolist() == [[1, 3], [4, 6]]
    assert adata[:, 'b':'c'].X.tolist() == [[2, 3], [5, 6]]

    from pytest import raises
    with raises(IndexError): _ = adata[:, 'X']
    with raises(IndexError): _ = adata['X', :]
    with raises(IndexError): _ = adata['A':'X', :]
    with raises(IndexError): _ = adata[:, 'a':'X']


def test_slicing_graphs():
    adata = AnnData(
        np.array([[1, 2], [3, 4], [5, 6]]),
        uns={'neighbors':
             {'connectivities': sp.csr_matrix(np.ones((3, 3)))}})
    adata_sub = adata[[0, 1], :]
    assert adata_sub.uns['neighbors']['connectivities'].shape[0] == 2
    assert adata.uns['neighbors']['connectivities'].shape[0] == 3
    assert adata_sub.copy().uns['neighbors']['connectivities'].shape[0] == 2


def test_slicing_series():
    adata = AnnData(
        np.array([[1, 2], [3, 4], [5, 6]]),
        dict(obs_names=['A', 'B', 'C']),
        dict(var_names=['a', 'b']))
    df = pd.DataFrame({'a': ['1', '2', '2']})
    df1 = pd.DataFrame({'b': ['1', '2']})
    assert (adata[df['a'].values == '2'].X.tolist()
            == adata[df['a'] == '2'].X.tolist())
    assert (adata[:, df1['b'].values == '2'].X.tolist()
            == adata[:, df1['b'] == '2'].X.tolist())


def test_slicing_remove_unused_categories():
    adata = AnnData(
        np.array([[1, 2], [3, 4], [5, 6], [7, 8]]),
        dict(k=['a', 'a', 'b', 'b']))
    adata._sanitize()
    assert adata[3:5].obs['k'].cat.categories.tolist() == ['b']


def test_slicing_integer_index():
    adata = AnnData(
        np.array([[0, 1, 2], [3, 4, 5]]),
        var=dict(var_names=[10, 11, 12]))
    sliced = adata[:, adata.X.sum(0) > 3]  # This used to fail
    assert sliced.shape == (2, 2)


def test_get_subset_annotation():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
                    dict(S=['A', 'B']),
                    dict(F=['a', 'b', 'c']))

    assert adata[0, 0].obs['S'].tolist() == ['A']
    assert adata[0, 0].var['F'].tolist() == ['a']


def test_transpose():
    adata = AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(obs_names=['A', 'B']),
        dict(var_names=['a', 'b', 'c']))

    adata1 = adata.T

    # make sure to not modify the original!
    assert adata.obs_names.tolist() == ['A', 'B']
    assert adata.var_names.tolist() == ['a', 'b', 'c']

    assert adata1.obs_names.tolist() == ['a', 'b', 'c']
    assert adata1.var_names.tolist() == ['A', 'B']
    assert adata1.X.shape == adata.X.T.shape

    adata2 = adata.transpose()
    assert np.array_equal(adata1.X, adata2.X)
    assert np.array_equal(adata1.obs, adata2.obs)
    assert np.array_equal(adata1.var, adata2.var)


def test_append_col():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))

    adata.obs['new'] = [1, 2]
    # this worked in the initial AnnData, but not with a dataframe
    # adata.obs[['new2', 'new3']] = [['A', 'B'], ['c', 'd']]

    from pytest import raises
    with raises(ValueError):
        adata.obs['new4'] = 'far too long'.split()


def test_delete_col():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]), dict(o1=[1, 2], o2=[3, 4]))
    assert ['o1', 'o2'] == adata.obs_keys()

    del adata.obs['o1']
    assert ['o2'] == adata.obs_keys()
    assert [3, 4] == adata.obs['o2'].tolist()


def test_set_obs():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))

    adata.obs = pd.DataFrame({'a': [3, 4]})
    assert adata.obs_names.tolist() == [0, 1]

    from pytest import raises
    with raises(ValueError):
        adata.obs = pd.DataFrame({'a': [3, 4, 5]})
        adata.obs = {'a': [1, 2]}


def test_multicol():
    adata = AnnData(np.array([[1, 2, 3], [4, 5, 6]]))
    # 'c' keeps the columns as should be
    adata.obsm['c'] = np.array([[0., 1.], [2, 3]])
    assert adata.obsm_keys() == ['c']
    assert adata.obsm['c'].tolist() == [[0., 1.], [2, 3]]


def test_n_obs():
    adata = AnnData(np.array([[1, 2], [3, 4], [5, 6]]))
    assert adata.n_obs == 3
    adata1 = adata[:2, ]
    assert adata1.n_obs == 2


def test_concatenate():
    # dense data
    adata1 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
                     {'obs_names': ['s1', 's2'],
                      'anno1': ['c1', 'c2']},
                     {'var_names': ['a', 'b', 'c'],
                      'annoA': [0, 1, 2]})
    adata2 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
                     {'obs_names': ['s3', 's4'],
                      'anno1': ['c3', 'c4']},
                     {'var_names': ['d', 'c', 'b'],
                      'annoA': [0, 1, 2]})
    adata3 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
                     {'obs_names': ['s1', 's2'],
                      'anno2': ['d3', 'd4']},
                     {'var_names': ['d', 'c', 'b'],
                      'annoB': [0, 1, 2]})

    # inner join
    adata = adata1.concatenate(adata2, adata3)
    assert adata.X.astype(int).tolist() == [[2, 3], [5, 6], [3, 2], [6, 5], [3, 2], [6, 5]]
    assert adata.obs_keys() == ['anno1', 'anno2', 'batch']
    assert adata.var_keys() == ['annoA-0', 'annoA-1', 'annoB-2']
    assert adata.var.values.tolist() == [[1, 2, 2], [2, 1, 1]]
    adata = adata1.concatenate(adata2, adata3, batch_key='batch1')
    assert adata.obs_keys() == ['anno1', 'anno2', 'batch1']
    adata = adata1.concatenate(adata2, adata3, batch_categories=['a1', 'a2', 'a3'])
    assert adata.obs['batch'].cat.categories.tolist() == ['a1', 'a2', 'a3']
    assert adata.var_names.tolist() == ['b', 'c']

    # outer join
    adata = adata1.concatenate(adata2, adata3, join='outer')
    from numpy import ma
    Xma = ma.masked_invalid(adata.X)
    Xma_ref = ma.masked_invalid(np.array([
        [1.0, 2.0, 3.0, np.nan],
        [4.0, 5.0, 6.0, np.nan],
        [np.nan, 3.0, 2.0, 1.0],
        [np.nan, 6.0, 5.0, 4.0],
        [np.nan, 3.0, 2.0, 1.0],
        [np.nan, 6.0, 5.0, 4.0]]))
    assert np.array_equal(Xma.mask, Xma_ref.mask)
    assert np.allclose(Xma.compressed(), Xma_ref.compressed())
    var_ma = ma.masked_invalid(adata.var.values.tolist())
    var_ma_ref = ma.masked_invalid(np.array(
        [[0.0, np.nan, np.nan], [1.0, 2.0, 2.0], [2.0, 1.0, 1.0], [np.nan, 0.0, 0.0]]))
    assert np.array_equal(var_ma.mask, var_ma_ref.mask)
    assert np.allclose(var_ma.compressed(), var_ma_ref.compressed())

    # sparse data
    from scipy.sparse import csr_matrix
    adata1 = AnnData(csr_matrix([[0, 2, 3], [0, 5, 6]]),
                     {'obs_names': ['s1', 's2'],
                      'anno1': ['c1', 'c2']},
                     {'var_names': ['a', 'b', 'c']})
    adata2 = AnnData(csr_matrix([[0, 2, 3], [0, 5, 6]]),
                     {'obs_names': ['s3', 's4'],
                      'anno1': ['c3', 'c4']},
                     {'var_names': ['d', 'c', 'b']})
    adata3 = AnnData(csr_matrix([[1, 2, 0], [0, 5, 6]]),
                     {'obs_names': ['s5', 's6'],
                      'anno2': ['d3', 'd4']},
                     {'var_names': ['d', 'c', 'b']})

    # inner join
    adata = adata1.concatenate(adata2, adata3)
    assert adata.X.toarray().astype(int).tolist() == [[2, 3], [5, 6], [3, 2], [6, 5], [0, 2], [6, 5]]

    # outer join
    adata = adata1.concatenate(adata2, adata3, join='outer')
    assert adata.X.toarray().tolist() == [
        [0.0, 2.0, 3.0, 0.0],
        [0.0, 5.0, 6.0, 0.0],
        [0.0, 3.0, 2.0, 0.0],
        [0.0, 6.0, 5.0, 0.0],
        [0.0, 0.0, 2.0, 1.0],
        [0.0, 6.0, 5.0, 0.0]]


def test_rename_categories():
    X = np.ones((6, 3))
    obs = pd.DataFrame(
        {'cat_anno': pd.Categorical(['a', 'a', 'a', 'a', 'b', 'a'])})
    adata = AnnData(X=X, obs=obs)
    adata.uns['tool'] = {}
    adata.uns['tool']['cat_array'] = np.rec.fromarrays(
        [np.ones(2) for cat in adata.obs['cat_anno'].cat.categories],
        dtype=[(cat, 'float32') for cat in adata.obs['cat_anno'].cat.categories])
    adata.uns['tool']['params'] = {'groupby': 'cat_anno'}

    new_categories = ['c', 'd']
    adata.rename_categories('cat_anno', new_categories)

    assert list(adata.obs['cat_anno'].cat.categories) == new_categories
    assert list(adata.uns['tool']['cat_array'].dtype.names) == new_categories

# def test_pickle():
#     adata = AnnData()
#     adata2 = pickle.loads(pickle.dumps(adata))
#     assert adata2.obsm._parent == adata2
