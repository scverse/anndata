import numpy as np
from numpy import ma
import pandas as pd
from scipy import sparse as sp

from anndata import AnnData


def test_creation():
    AnnData(np.array([[1, 2], [3, 4]]))
    AnnData(np.array([[1, 2], [3, 4]]), {}, {})
    AnnData(ma.array([[1, 2], [3, 4]]), uns={'mask': [0, 1, 1, 0]})
    AnnData(sp.eye(2))
    AnnData(
        np.array([[1, 2, 3], [4, 5, 6]]),
        dict(Obs=['A', 'B']),
        dict(Feat=['a', 'b', 'c']))

    assert AnnData(np.array([1, 2])).X.shape == (2,)

    from pytest import raises
    raises(ValueError, AnnData,
           np.array([[1, 2], [3, 4]]),
           dict(TooLong=[1, 2, 3, 4]))


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
    adata1 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
                     {'obs_names': ['s1', 's2'],
                      'anno1': ['c1', 'c2']},
                     {'var_names': ['a', 'b', 'c']})
    adata2 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
                     {'obs_names': ['s3', 's4'],
                      'anno1': ['c3', 'c4']},
                     {'var_names': ['d', 'c', 'b']})
    adata3 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
                     {'obs_names': ['s5', 's6'],
                      'anno2': ['d3', 'd4']},
                     {'var_names': ['d', 'c', 'b']})
    adata = adata1.concatenate(adata2, adata3)
    assert adata.n_vars == 2
    assert adata.obs_keys() == ['anno1', 'anno2', 'batch']
    adata = adata1.concatenate(adata2, adata3, batch_key='batch1')
    assert adata.obs_keys() == ['anno1', 'anno2', 'batch1']
    adata = adata1.concatenate(adata2, adata3, batch_categories=['a1', 'a2', 'a3'])
    assert adata.obs['batch'].cat.categories.tolist() == ['a1', 'a2', 'a3']
    assert adata.var_names.tolist() == ['b', 'c']
    

def test_concatenate_sparse():
    from scipy.sparse import csr_matrix
    adata1 = AnnData(csr_matrix([[0, 2, 3], [0, 5, 6]]),
                     {'obs_names': ['s1', 's2'],
                      'anno1': ['c1', 'c2']},
                     {'var_names': ['a', 'b', 'c']})
    adata2 = AnnData(csr_matrix([[0, 2, 3], [0, 5, 6]]),
                     {'obs_names': ['s3', 's4'],
                      'anno1': ['c3', 'c4']},
                     {'var_names': ['b', 'c', 'd']})
    adata3 = AnnData(csr_matrix([[1, 2, 0], [0, 5, 6]]),
                     {'obs_names': ['s5', 's6'],
                      'anno2': ['d3', 'd4']},
                     {'var_names': ['b', 'c', 'd']})
    adata = adata1.concatenate(adata2, adata3)
    assert adata.n_vars == 2


def test_concatenate_outer():
    adata1 = AnnData(np.array([[1, 2, 3], [4, 5, 6]]),
                     {'obs_names': ['s1', 's2'],
                      'anno1': ['c1', 'c2']},
                     {'var_names': ['a', 'b', 'c']})
    adata2 = AnnData(np.array([[1, 2, 3], [4, 5, 6], [7,8,9]]),
                     {'obs_names': ['s3', 's4', 's5'],
                      'anno2': ['c3', 'c4', 'c5']},
                     {'var_names': ['b', 'c', 'd']})
    adata = adata1.concatenate(adata2, join='outer')
    assert adata.n_vars == 4
    assert adata.obs_keys() == ['anno1', 'anno2', 'batch']

# TODO: remove logging and actually test values
# from scanpy import logging as logg

# def test_profile_memory():
#     import gc
#     dim = 10  # increase this when profiling
#     print()
#     logg.print_memory_usage('start profiling')
#     X = np.random.rand(dim, dim).astype('float32')
#     logg.print_memory_usage('allocated X')
#     var_filter = np.array([0, 1])
#     X = X[:, var_filter]
#     logg.print_memory_usage('sliced X')
#     X = np.random.rand(dim, dim).astype('float32')
#     logg.print_memory_usage('allocated X')
#     adata = AnnData(X)
#     logg.print_memory_usage('init adata with reference to X')
#     adata.var['multi'] = np.random.rand(dim, 3)
#     logg.print_memory_usage('added some annotation')
#     # ------------------------------------------------
#     # compare adata.__getitem__ with adata.filter_var
#     # ------------------------------------------------
#     # here, it doesn't make a difference in other scenarios
#     # (e.g. sc.preprocess.weinreb16), filter_var seems to invoke earlier garbage
#     # collection than slicing
#     # adata.filter_var(var_filter)  # inplace
#     adata = adata[:, var_filter]  # with copy
#     logg.print_memory_usage('sliced adata')
#     gc.collect()
#     logg.print_memory_usage('after calling gc.collect()')
#     return adata


# def test_profile_memory_2():
#     adata = test_profile_memory()
#     logg.print_memory_usage('after leaving function')
