import logging as logg
import pandas as pd


def make_index_unique(index, join='-'):
    """Makes the index unique by appending '1', '2', etc.

    The first occurance of a non-unique value is ignored.

    Parameters
    ----------
    join : `str`, optional (default: '')
         The connecting string between name and integer.

    Examples
    --------
    >>> adata1 = sc.AnnData(np.ones((3, 2)), {'obs_names': ['a', 'b', 'c']})
    >>> adata2 = sc.AnnData(np.zeros((3, 2)), {'obs_names': ['d', 'b', 'b']})
    >>> adata = adata1.concatenate(adata2)
    >>> adata.obs_names
    Index(['a', 'b', 'c', 'd', 'b', 'b'], dtype='object')
    >>> adata.obs_names_make_unique()
    >>> adata.obs_names
    Index(['a', 'b', 'c', 'd', 'b-1', 'b-2'], dtype='object')
    """
    if index.is_unique:
        return index
    from collections import defaultdict
    values = index.values
    indices_dup = index.duplicated(keep='first')
    values_dup = values[indices_dup]
    counter = defaultdict(lambda: 0)
    for i, v in enumerate(values_dup):
        counter[v] += 1
        values_dup[i] += join + str(counter[v])
    values[indices_dup] = values_dup
    index = pd.Index(values)
    return index


def warn_names_duplicates(string, df):
    names = 'Observation' if string == 'obs' else 'Variable'
    logg.info(
        '{} names are not unique. '
        'To make them unique, call `.{}_names_make_unique`.'
        .format(names, string))
