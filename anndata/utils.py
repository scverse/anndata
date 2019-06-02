import warnings
from functools import wraps, singledispatch
from typing import Mapping, Any, Sequence, Union, Tuple

import pandas as pd
import numpy as np
from scipy.sparse import spmatrix

from .logging import get_logger

logger = get_logger(__name__)


@singledispatch
def convert_to_dict(obj) -> dict:
    return dict(obj)

@convert_to_dict.register(dict)
def convert_to_dict_dict(obj: dict):
    return obj

@convert_to_dict.register(np.ndarray)
def convert_to_dict_ndarray(obj: np.ndarray):
    if obj.dtype.fields is None:
        raise TypeError(
            "Can only convert np.ndarray with compound dtypes to dict, passed "
            "array had '{}'.".format(obj.dtype)
        )
    return {k: obj[k] for k in obj.dtype.fields.keys()}

@convert_to_dict.register(type(None))
def convert_to_dict_nonetype(obj: None):
    return dict()


def make_index_unique(index: pd.Index, join: str = '-'):
    """Makes the index unique by appending '1', '2', etc.

    The first occurance of a non-unique value is ignored.

    Parameters
    ----------
    join
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


def warn_names_duplicates(attr: str):
    names = 'Observation' if attr == 'obs' else 'Variable'
    logger.info(
        '{} names are not unique. '
        'To make them unique, call `.{}_names_make_unique`.'
        .format(names, attr))


def warn_no_string_index(names: Sequence[Any]):
    if not isinstance(names[0], str):
        logger.warning(
            'AnnData expects string indices for some functionality, but your first two indices are: {}. '
            .format(names[:2]))


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


def deprecated(new_name: str):
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


Index = Union[slice, int, np.int64, np.ndarray, spmatrix]


def get_n_items_idx(idx: Index, l: int):
    if isinstance(idx, np.ndarray) and idx.dtype == bool:
        return idx.sum()
    elif isinstance(idx, slice):
        start = 0 if idx.start is None else idx.start
        stop = l if idx.stop is None else idx.stop
        step = 1 if idx.step is None else idx.step
        return (stop - start) // step
    elif isinstance(idx, (int, np.int_, np.int64, np.int32)):
        return 1
    else:
        return len(idx)


def unpack_index(index: Union[Index, Tuple[Index, Index]]) -> Tuple[Index, Index]:
    # handle indexing with boolean matrices
    if (
        isinstance(index, (spmatrix, np.ndarray))
        and index.ndim == 2
        and index.dtype.kind == 'b'
    ): return index.nonzero()

    if not isinstance(index, tuple):
        return index, slice(None)
    elif len(index) == 2:
        return index
    elif len(index) == 1:
        return index[0], slice(None)
    else:
        raise IndexError('invalid number of indices')
