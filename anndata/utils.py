import warnings
from functools import wraps, singledispatch
from typing import Mapping, Any, Sequence, Union, Callable

import h5py
import pandas as pd
import numpy as np
from scipy import sparse

from .logging import get_logger
from ._core.sparse_dataset import SparseDataset

logger = get_logger(__name__)


@singledispatch
def asarray(x):
    """Convert x to a numpy array"""
    return np.asarray(x)


@asarray.register(sparse.spmatrix)
def asarray_sparse(x):
    return x.toarray()


@asarray.register(SparseDataset)
def asarray_sparse_dataset(x):
    return asarray(x.value)


@asarray.register(h5py.Dataset)
def asarray_h5py_dataset(x):
    return x[...]


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


@singledispatch
def dim_len(x, axis):
    """\
    Return the size of an array in dimension `axis`.

    Returns None if `x` is an awkward array with variable length in the requested dimension.
    """
    return x.shape[axis]


try:
    from .compat import awkward as ak

    def _size_at_depth(layout, depth, lateral_context, **kwargs):
        """Callback function for dim_len_awkward, resolving the dim_len for a given level"""
        if layout.is_numpy:
            # if it's an embedded rectilinear array, we have to deal with its shape
            # which might not be 1-dimensional
            if layout.is_unknown:
                shape = (0,)
            else:
                shape = layout.shape
            numpy_axis = lateral_context["axis"] - depth + 1
            if not (1 <= numpy_axis < len(shape)):
                raise TypeError(f"axis={lateral_context['axis']} is too deep")
            lateral_context["out"] = shape[numpy_axis]
            return ak.contents.EmptyArray()

        elif layout.is_list and depth == lateral_context["axis"]:
            if layout.parameter("__array__") in ("string", "bytestring"):
                # Strings are implemented like an array of lists of uint8 (ListType(NumpyType(...)))
                # which results in an extra hierarchy-level that shouldn't show up in dim_len
                # See https://github.com/scikit-hep/awkward/discussions/1654#discussioncomment-3736747
                raise TypeError(f"axis={lateral_context['axis']} is too deep")

            if layout.is_regular:
                # if it's a regular list, you want the size
                lateral_context["out"] = layout.size
            else:
                # if it's an irregular list, you want a null token
                lateral_context["out"] = -1
            return ak.contents.EmptyArray()

        elif layout.is_record and depth == lateral_context["axis"]:
            lateral_context["out"] = len(layout.fields)
            return ak.contents.EmptyArray()

        elif layout.is_record:
            # currently, we don't recurse into records
            # in theory we could, just not sure how to do it at the moment
            # Would need to consider cases like: scalars, unevenly sized values
            raise TypeError(
                f"Cannot recurse into record type found at axis={lateral_context['axis']}"
            )

        elif layout.is_union:
            # if it's a union, you could get the result of each union branch
            # separately and see if they're all the same; if not, it's an error
            result = None
            for content in layout.contents:
                context = {"axis": lateral_context["axis"]}
                ak.transform(
                    _size_at_depth,
                    content,
                    lateral_context=context,
                )
                if result is None:
                    result = context["out"]
                elif result != context["out"]:
                    # Union branches have different lengths -> return null token
                    lateral_context["out"] = -1
                    return ak.contents.EmptyArray()
            lateral_context["out"] = result
            return ak.contents.EmptyArray()

    @dim_len.register(ak.Array)
    def dim_len_awkward(array, axis):
        """Get the length of an awkward array in a given dimension

        Returns None if the dimension is of variable length.

        Code adapted from @jpivarski's solution in https://github.com/scikit-hep/awkward/discussions/1654#discussioncomment-3521574
        """
        if axis < 0:  # negative axis is another can of worms... maybe later
            raise NotImplementedError("Does not support negative axis")
        elif axis == 0:
            return len(array)
        else:
            # communicate with the recursive function using a context (lateral)
            context = {"axis": axis}

            # "transform" but we don't care what kind of array it returns
            ak.transform(
                _size_at_depth,
                array,
                lateral_context=context,
            )

            # Use `None` as null token.
            return None if context["out"] == -1 else context["out"]

    @asarray.register(ak.Array)
    def asarray_awkward(x):
        return x

except ImportError:
    pass


def make_index_unique(index: pd.Index, join: str = "-"):
    """
    Makes the index unique by appending a number string to each duplicate index element:
    '1', '2', etc.

    If a tentative name created by the algorithm already exists in the index, it tries
    the next integer in the sequence.

    The first occurrence of a non-unique value is ignored.

    Parameters
    ----------
    join
         The connecting string between name and integer.

    Examples
    --------
    >>> from anndata import AnnData
    >>> adata = AnnData(np.ones((2, 3)), var=pd.DataFrame(index=["a", "a", "b"]))
    >>> adata.var_names
    Index(['a', 'a', 'b'], dtype='object')
    >>> adata.var_names_make_unique()
    >>> adata.var_names
    Index(['a', 'a-1', 'b'], dtype='object')
    """
    if index.is_unique:
        return index
    from collections import Counter

    values = index.values.copy()
    indices_dup = index.duplicated(keep="first")
    values_dup = values[indices_dup]
    values_set = set(values)
    counter = Counter()
    issue_interpretation_warning = False
    example_colliding_values = []
    for i, v in enumerate(values_dup):
        while True:
            counter[v] += 1
            tentative_new_name = v + join + str(counter[v])
            if tentative_new_name not in values_set:
                values_set.add(tentative_new_name)
                values_dup[i] = tentative_new_name
                break
            issue_interpretation_warning = True
            if len(example_colliding_values) < 5:
                example_colliding_values.append(tentative_new_name)

    if issue_interpretation_warning:
        warnings.warn(
            f"Suffix used ({join}[0-9]+) to deduplicate index values may make index "
            + "values difficult to interpret. There values with a similar suffixes in "
            + "the index. Consider using a different delimiter by passing "
            + "`join={delimiter}`"
            + "Example key collisions generated by the make_index_unique algorithm: "
            + str(example_colliding_values)
        )
    values[indices_dup] = values_dup
    index = pd.Index(values, name=index.name)
    return index


def warn_names_duplicates(attr: str):
    names = "Observation" if attr == "obs" else "Variable"
    warnings.warn(
        f"{names} names are not unique. "
        f"To make them unique, call `.{attr}_names_make_unique`.",
        UserWarning,
        stacklevel=2,
    )


def ensure_df_homogeneous(
    df: pd.DataFrame, name: str
) -> Union[np.ndarray, sparse.csr_matrix]:
    # TODO: rename this function, I would not expect this to return a non-dataframe
    if all(isinstance(dt, pd.SparseDtype) for dt in df.dtypes):
        arr = df.sparse.to_coo().tocsr()
    else:
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


def import_function(module: str, name: str) -> Callable:
    """\
    Try to import function from module. If the module is not installed or
    function is not part of the module, it returns a dummy function that raises
    the respective import error once the function is called. This could be a
    ModuleNotFoundError if the module is missing or an AttributeError if the
    module is installed but the function is not exported by it.

    Params
    -------
    module
        Module to import from. Can be nested, e.g. "sklearn.utils".
    name
        Name of function to import from module.
    """
    from importlib import import_module

    try:
        module = import_module(module)
        func = getattr(module, name)
    except (ImportError, AttributeError) as e:
        error = e

        def func(*_, **__):
            raise error

    return func
