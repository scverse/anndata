from __future__ import annotations

import re
import sys
import warnings
from functools import partial, singledispatch
from types import FunctionType, UnionType
from typing import TYPE_CHECKING, Literal, TypeAliasType, get_args, get_origin

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

import anndata

from ._core.sparse_dataset import BaseCompressedSparseDataset
from ._warnings import warn
from .compat import CSArray, CupyArray, CupySparseMatrix, DaskArray
from .logging import get_logger

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Mapping, Sequence
    from typing import Any, LiteralString

logger = get_logger(__name__)


def import_name(full_name: str) -> Any:
    from importlib import import_module

    parts = full_name.split(".")
    obj = import_module(parts[0])
    for _i, name in enumerate(parts[1:]):
        i = _i
        try:
            obj = import_module(f"{obj.__name__}.{name}")
        except ModuleNotFoundError:
            break
    else:
        i = len(parts)
    for name in parts[i + 1 :]:
        try:
            obj = getattr(obj, name)
        except AttributeError as e:
            msg = f"{parts[:i]}, {parts[i + 1 :]}, {obj} {name}"
            raise RuntimeError(msg) from e
    return obj


@singledispatch
def asarray(x):
    """Convert x to a numpy array"""
    return np.asarray(x)


@asarray.register(CSArray)
@asarray.register(sparse.spmatrix)
def asarray_sparse(x):
    return x.toarray()


@asarray.register(BaseCompressedSparseDataset)
def asarray_sparse_dataset(x):
    return asarray(x.to_memory())


@asarray.register(h5py.Dataset)
def asarray_h5py_dataset(x):
    return x[...]


@asarray.register(CupyArray)
def asarray_cupy(x):
    return x.get()


@asarray.register(CupySparseMatrix)
def asarray_cupy_sparse(x):
    return x.toarray().get()


@asarray.register(DaskArray)
def asarray_dask(x):
    return asarray(x.compute())


@singledispatch
def convert_to_dict(obj) -> dict:
    return dict(obj)


@convert_to_dict.register(dict)
def convert_to_dict_dict(obj: dict):
    return obj


@convert_to_dict.register(np.ndarray)
def convert_to_dict_ndarray(obj: np.ndarray):
    if obj.dtype.fields is None:
        msg = (
            "Can only convert np.ndarray with compound dtypes to dict, "
            f"passed array had “{obj.dtype}”."
        )
        raise TypeError(msg)
    return {k: obj[k] for k in obj.dtype.fields}


@convert_to_dict.register(type(None))
def convert_to_dict_nonetype(obj: None):
    return dict()


@singledispatch
def axis_len(x, axis: Literal[0, 1]) -> int | None:
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
            shape = (0,) if layout.is_unknown else layout.shape
            numpy_axis = lateral_context["axis"] - depth + 1
            if not (1 <= numpy_axis < len(shape)):
                msg = f"axis={lateral_context['axis']} is too deep"
                raise TypeError(msg)
            lateral_context["out"] = shape[numpy_axis]
            return ak.contents.EmptyArray()

        elif layout.is_list and depth == lateral_context["axis"]:
            if layout.parameter("__array__") in {"string", "bytestring"}:
                # Strings are implemented like an array of lists of uint8 (ListType(NumpyType(...)))
                # which results in an extra hierarchy-level that shouldn't show up in dim_len
                # See https://github.com/scikit-hep/awkward/discussions/1654#discussioncomment-3736747
                msg = f"axis={lateral_context['axis']} is too deep"
                raise TypeError(msg)

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
            msg = f"Cannot recurse into record type found at axis={lateral_context['axis']}"
            raise TypeError(msg)

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

    @axis_len.register(ak.Array)
    def axis_len_awkward(array, axis: Literal[0, 1]) -> int | None:
        """Get the length of an awkward array in a given axis

        Returns None if the axis is of variable length.

        Code adapted from @jpivarski's solution in https://github.com/scikit-hep/awkward/discussions/1654#discussioncomment-3521574
        """
        if axis < 0:  # negative axis is another can of worms... maybe later
            msg = "Does not support negative axis"
            raise NotImplementedError(msg)
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


def make_index_unique(index: pd.Index[str], join: str = "-") -> pd.Index[str]:
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
    >>> adata.var_names.astype("string")
    Index(['a', 'a', 'b'], dtype='string')
    >>> adata.var_names_make_unique()
    >>> adata.var_names.astype("string")
    Index(['a', 'a-1', 'b'], dtype='string')
    """
    if index.is_unique:
        return index
    from collections import Counter

    values = index.array.copy()
    indices_dup = index.duplicated(keep="first") & ~index.isna()
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
        msg = (
            f"Suffix used ({join}[0-9]+) to deduplicate index values may make index values difficult to interpret. "
            "There values with a similar suffixes in the index. "
            "Consider using a different delimiter by passing `join={delimiter}`. "
            "Example key collisions generated by the make_index_unique algorithm: "
            f"{example_colliding_values}"
        )
        warn(msg, UserWarning)
    values[indices_dup] = values_dup
    index = pd.Index(values, name=index.name)
    return index


def join_english(words: Iterable[str], conjunction: str = "or") -> str:
    words = list(words)  # no need to be efficient
    if len(words) == 0:
        return ""
    if len(words) == 1:
        return words[0]
    if len(words) == 2:
        return f"{words[0]} {conjunction} {words[1]}"
    return ", ".join(words[:-1]) + f", {conjunction} {words[-1]}"


def warn_names_duplicates(attr: str):
    names = "Observation" if attr == "obs" else "Variable"
    msg = (
        f"{names} names are not unique. "
        f"To make them unique, call `.{attr}_names_make_unique`."
    )
    warn(msg, UserWarning)


def ensure_df_homogeneous(
    df: pd.DataFrame, name: str
) -> np.ndarray | sparse.csr_matrix:
    # TODO: rename this function, I would not expect this to return a non-dataframe
    if all(isinstance(dt, pd.SparseDtype) for dt in df.dtypes):
        arr = df.sparse.to_coo().tocsr()
    else:
        arr = df.to_numpy()
    if df.dtypes.nunique() != 1:
        msg = f"{name} converted to numpy array with dtype {arr.dtype}"
        warn(msg, UserWarning)
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
    except UnicodeEncodeError as e:
        msg = (
            "Currently only support ascii strings. "
            "Don’t use “ö” etc. for sample annotation."
        )
        raise ValueError(msg) from e

    # if old_index_key not in source:
    #     names.append(new_index_key)
    #     cols.append(np.arange(len(cols[0]) if cols else n_row).astype("U"))
    # else:
    #     names[names.index(old_index_key)] = new_index_key
    #     cols[names.index(old_index_key)] = cols[names.index(old_index_key)].astype("U")
    dtype_list = list(
        zip(
            names,
            [str(c.dtype) for c in cols],
            [(c.shape[1],) for c in cols],
            strict=True,
        )
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


def warn_once(msg: str, category: type[Warning]) -> None:
    warn(msg, category)
    # Prevent from showing up every time an awkward array is used
    # You'd think `'once'` works, but it doesn't at the repl and in notebooks
    warnings.filterwarnings("ignore", category=category, message=re.escape(msg))


if TYPE_CHECKING:
    from warnings import deprecated
else:
    if sys.version_info >= (3, 13):
        from warnings import deprecated as _deprecated
    else:
        from typing_extensions import deprecated as _deprecated
    deprecated = partial(_deprecated, category=FutureWarning)


def deprecation_msg(
    name: LiteralString, new_name: LiteralString, add_msg: LiteralString | None = None
) -> LiteralString:
    msg = (
        f"Use {new_name} instead of {name}, "
        f"{name} is deprecated and will be removed in the future."
    )
    if add_msg is not None:
        msg += f" {add_msg}"
    return msg


class DeprecationMixinMeta(type):
    """\
    Use this as superclass so deprecated methods and properties
    do not appear in vars(MyClass)/dir(MyClass)
    """

    def __dir__(cls):
        dont_hide = getattr(cls, "_DONT_HIDE_DEPRECATED", set())

        def is_hidden(attr: object) -> bool:
            if isinstance(attr, property):
                attr = attr.fget
            is_deprecated = bool(getattr(attr, "__deprecated__", None))
            return is_deprecated and getattr(attr, "__name__", None) not in dont_hide

        return [
            item
            for item in type.__dir__(cls)
            if not is_hidden(getattr(cls, item, None))
        ]


def set_module[C: FunctionType | type](name: str, /) -> Callable[[C], C]:
    def decorator(f: C) -> C:
        f.__module__ = name
        return f

    return decorator


def get_union_members[T](
    typ: TypeAliasType | UnionType | T,
    can_get_args: Callable[[object], bool] = lambda _: False,
) -> tuple[T, ...]:
    """Get args of a nested union of types."""
    while isinstance(typ, TypeAliasType):
        typ = typ.__value__
    if not isinstance(typ, UnionType) and not can_get_args(typ):
        return (typ,)
    return tuple(
        t for sub in get_args(typ) for t in get_union_members(sub, can_get_args)
    )


get_literal_members = partial(
    get_union_members, can_get_args=lambda x: get_origin(x) is Literal
)
"""Get args of a nested union of Literal types."""


def raise_value_error_if_multiindex_columns(df: pd.DataFrame, attr: str):
    if isinstance(df.columns, pd.MultiIndex):
        msg = (
            "MultiIndex columns are not supported in AnnData. "
            f"Please use a single-level index for {attr}."
        )
        raise ValueError(msg)


def module_get_attr_redirect(
    attr_name: str,
    deprecated_mapping: Mapping[str, str],
    old_module_path: str | None = None,
) -> Any:
    full_old_module_path = (
        f"anndata{'.' + old_module_path if old_module_path is not None else ''}"
    )
    if new_path := deprecated_mapping.get(attr_name):
        msg = (
            f"Importing {attr_name} from `{full_old_module_path}` is deprecated. "
            f"Import anndata.{new_path} instead."
        )
        warn(msg, FutureWarning)
        # hacky import_object_by_name, but we test all these
        mod = anndata
        while "." in new_path:
            mod_name, new_path = new_path.split(".", 1)
            mod = getattr(mod, mod_name)
        return getattr(mod, new_path)
    msg = f"module {full_old_module_path} has no attribute {attr_name!r}"
    raise AttributeError(msg)
