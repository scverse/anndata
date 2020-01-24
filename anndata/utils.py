import warnings
from functools import wraps, singledispatch
from typing import Mapping, Any, Sequence

import pandas as pd
import numpy as np

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
            "Can only convert np.ndarray with compound dtypes to dict, "
            f"passed array had “{obj.dtype}”."
        )
    return {k: obj[k] for k in obj.dtype.fields.keys()}


@convert_to_dict.register(type(None))
def convert_to_dict_nonetype(obj: None):
    return dict()


def make_index_unique(index: pd.Index, join: str = "-"):
    """\
    Makes the index unique by appending '1', '2', etc.

    The first occurance of a non-unique value is ignored.

    Parameters
    ----------
    join
         The connecting string between name and integer.

    Examples
    --------
    >>> from anndata import AnnData
    >>> adata1 = AnnData(np.ones((3, 2)), dict(obs_names=['a', 'b', 'c']))
    >>> adata2 = AnnData(np.zeros((3, 2)), dict(obs_names=['d', 'b', 'b']))
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
    indices_dup = index.duplicated(keep="first")
    values_dup = values[indices_dup]
    counter = defaultdict(lambda: 0)
    for i, v in enumerate(values_dup):
        counter[v] += 1
        values_dup[i] += join + str(counter[v])
    values[indices_dup] = values_dup
    index = pd.Index(values)
    return index


def warn_names_duplicates(attr: str):
    names = "Observation" if attr == "obs" else "Variable"
    logger.info(
        f"{names} names are not unique. "
        f"To make them unique, call `.{attr}_names_make_unique`."
    )


def ensure_df_homogeneous(df: pd.DataFrame, name: str) -> np.ndarray:
    arr = df.to_numpy()
    if df.dtypes.nunique() != 1:
        warnings.warn(f"{name} converted to numpy array with dtype {arr.dtype}")
    return arr


def convert_dictionary_to_structured_array(source: Mapping[str, Sequence[Any]]):
    names = list(source.keys())
    try:  # transform to byte-strings
        cols = [
            np.asarray(col)
            if np.array(col[0]).dtype.char not in {"U", "S"}
            else np.asarray(col).astype("U")
            for col in source.values()
        ]
    except UnicodeEncodeError:
        raise ValueError(
            "Currently only support ascii strings. "
            "Don’t use “ö” etc. for sample annotation."
        )

    # if old_index_key not in source:
    #     names.append(new_index_key)
    #     cols.append(np.arange(len(cols[0]) if cols else n_row).astype("U"))
    # else:
    #     names[names.index(old_index_key)] = new_index_key
    #     cols[names.index(old_index_key)] = cols[names.index(old_index_key)].astype("U")
    dtype_list = list(
        zip(names, [str(c.dtype) for c in cols], [(c.shape[1],) for c in cols])
    )
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
    """\
    This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used.
    """

    def decorator(func):
        @wraps(func)
        def new_func(*args, **kwargs):
            # turn off filter
            warnings.simplefilter("always", DeprecationWarning)
            warnings.warn(
                f"Use {new_name} instead of {func.__name__}, "
                f"{func.__name__} will be removed in the future.",
                category=DeprecationWarning,
                stacklevel=2,
            )
            warnings.simplefilter("default", DeprecationWarning)  # reset filter
            return func(*args, **kwargs)

        setattr(new_func, "__deprecated", True)
        return new_func

    return decorator


class DeprecationMixinMeta(type):
    """\
    Use this as superclass so deprecated methods and properties
    do not appear in vars(MyClass)/dir(MyClass)
    """

    def __dir__(cls):
        def is_deprecated(attr):
            if isinstance(attr, property):
                attr = attr.fget
            return getattr(attr, "__deprecated", False)

        return [
            item
            for item in type.__dir__(cls)
            if not is_deprecated(getattr(cls, item, None))
        ]
