import logging as logg
import warnings
from functools import wraps
from typing import Mapping, Any, Sequence

import pandas as pd
import numpy as np


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


def convert_dictionary_to_structured_array(source: Mapping[str, Sequence[Any]]):

    names = list(source.keys())
    try:  # transform to byte-strings
        cols = [np.asarray(col) if np.array(col[0]).dtype.char not in {'U', 'S'}
                else np.asarray(col).astype('U') for col in source.values()]
    except UnicodeEncodeError:
        raise ValueError(
            'Currently only support ascii strings. Don\'t use "ö" etc. for sample annotation.')

    # if old_index_key not in source:
    #     names.append(new_index_key)
    #     cols.append(np.arange(len(cols[0]) if cols else n_row).astype('U'))
    # else:
    #     names[names.index(old_index_key)] = new_index_key
    #     cols[names.index(old_index_key)] = cols[names.index(old_index_key)].astype('U')
    dtype_list = list(zip(names, [str(c.dtype) for c in cols], [(c.shape[1],) for c in cols]))
    # might be unnecessary
    dtype = np.dtype(dtype_list)

    arr = np.zeros((len(cols[0]),), dtype)
    # here, we do not want to call BoundStructArray.__getitem__
    # but np.ndarray.__getitem__, therefore we avoid the following line
    # arr = np.ndarray.__new__(cls, (len(cols[0]),), dtype)
    for i, name in enumerate(dtype.names):
        arr[name] = np.array(cols[i], dtype=dtype_list[i][1])

    return arr


def deprecated(new_name):
    """
    This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used.
    """
    def decorator(func):
        @wraps(func)
        def new_func(*args, **kwargs):
            warnings.simplefilter('always', DeprecationWarning)  # turn off filter
            warnings.warn(
                'Use {0} instead of {1}, {1} will be removed in the future.'
                .format(new_name, func.__name__),
                category=DeprecationWarning,
                stacklevel=2,
            )
            warnings.simplefilter('default', DeprecationWarning)  # reset filter
            return func(*args, **kwargs)
        setattr(new_func, '__deprecated', True)
        return new_func
    return decorator


class DeprecationMixinMeta(type):
    """
    Use this as superclass so deprecated methods and properties
    do not appear in vars(MyClass)/dir(MyClass)
    """
    def __dir__(cls):
        def is_deprecated(attr):
            if isinstance(attr, property):
                attr = attr.fget
            return getattr(attr, '__deprecated', False)
        
        return [
            item for item in type.__dir__(cls)
            if not is_deprecated(getattr(cls, item, None))
        ]
