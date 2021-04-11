from copy import deepcopy
from functools import reduce, wraps
from inspect import signature, Parameter
from typing import Collection, Union, Mapping, MutableMapping, Optional
from warnings import warn

import h5py
from scipy.sparse import spmatrix
import numpy as np
import pandas as pd

from ._overloaded_dict import _overloaded_uns, OverloadedDict
from .._core.index import _subset

# try importing zarr, dask, and zappy
from packaging import version

try:
    from zarr.core import Array as ZarrArray
except ImportError:

    class ZarrArray:
        @staticmethod
        def __repr__():
            return "mock zarr.core.Array"


try:
    from zappy.base import ZappyArray
except ImportError:

    class ZappyArray:
        @staticmethod
        def __repr__():
            return "mock zappy.base.ZappyArray"


try:
    from dask.array import Array as DaskArray
except ImportError:

    class DaskArray:
        @staticmethod
        def __repr__():
            return "mock dask.array.core.Array"


try:
    from typing import Literal
except ImportError:
    try:
        from typing_extensions import Literal
    except ImportError:

        class LiteralMeta(type):
            def __getitem__(cls, values):
                if not isinstance(values, tuple):
                    values = (values,)
                return type("Literal_", (Literal,), dict(__args__=values))

        class Literal(metaclass=LiteralMeta):
            pass


def _from_fixed_length_strings(value):
    """\
    Convert from fixed length strings to unicode.

    For backwards compatability with older h5ad and zarr files.
    """
    new_dtype = []
    for dt in value.dtype.descr:
        dt_list = list(dt)
        dt_type = dt[1]
        # could probably match better
        is_annotated = isinstance(dt_type, tuple)
        if is_annotated:
            dt_type = dt_type[0]
        # Fixing issue introduced with h5py v2.10.0, see:
        # https://github.com/h5py/h5py/issues/1307
        if issubclass(np.dtype(dt_type).type, np.string_):
            dt_list[1] = f"U{int(dt_type[2:])}"
        elif is_annotated or np.issubdtype(np.dtype(dt_type), np.str_):
            dt_list[1] = "O"  # Assumption that it’s a vlen str
        new_dtype.append(tuple(dt_list))
    return value.astype(new_dtype)


def _decode_structured_array(
    arr: np.ndarray, dtype: Optional[np.dtype] = None, copy: bool = False
) -> np.ndarray:
    """
    h5py 3.0 now reads all strings as bytes. There is a helper method which can convert these to strings,
    but there isn't anything for fields of structured dtypes.

    Params
    ------
    arr
        An array with structured dtype
    dtype
        dtype of the array. This is checked for h5py string data types.
        Passing this is allowed for cases where array may have been processed by another function before hand.
    """
    if copy:
        arr = arr.copy()
    if dtype is None:
        dtype = arr.dtype
    # codecs.decode is 2x slower than this lambda, go figure
    decode = np.frompyfunc(lambda x: x.decode("utf-8"), 1, 1)
    for k, (dt, _) in dtype.fields.items():
        check = h5py.check_string_dtype(dt)
        if check is not None and check.encoding == "utf-8":
            decode(arr[k], out=arr[k])
    return arr


def _to_fixed_length_strings(value: np.ndarray) -> np.ndarray:
    """\
    Convert variable length strings to fixed length.

    Currently a workaround for
    https://github.com/zarr-developers/zarr-python/pull/422
    """
    new_dtype = []
    for dt_name, (dt_type, dt_offset) in value.dtype.fields.items():
        if dt_type.kind == "O":
            #  Assuming the objects are str
            size = max(len(x.encode()) for x in value.getfield("O", dt_offset))
            new_dtype.append((dt_name, ("U", size)))
        else:
            new_dtype.append((dt_name, dt_type))
    return value.astype(new_dtype)


#############################
# Dealing with uns
#############################


def _clean_uns(d: Mapping[str, MutableMapping[str, Union[pd.Series, str, int]]]):
    """
    Compat function for when categorical keys were stored in uns.
    This used to be buggy because when storing categorical columns in obs and var with
    the same column name, only one `<colname>_categories` is retained.
    """
    k_to_delete = set()
    for cats_name, cats in d.get("uns", {}).items():
        if not cats_name.endswith("_categories"):
            continue
        name = cats_name.replace("_categories", "")
        # fix categories with a single category
        if isinstance(cats, (str, int)):
            cats = [cats]
        for ann in ["obs", "var"]:
            if name not in d[ann]:
                continue
            codes: np.ndarray = d[ann][name].values
            # hack to maybe find the axis the categories were for
            if not np.all(codes < len(cats)):
                continue
            d[ann][name] = pd.Categorical.from_codes(codes, cats)
            k_to_delete.add(cats_name)
    for cats_name in k_to_delete:
        del d["uns"][cats_name]


def _move_adj_mtx(d):
    """
    Read-time fix for moving adjacency matrices from uns to obsp
    """
    n = d.get("uns", {}).get("neighbors", {})
    obsp = d.setdefault("obsp", {})

    for k in ("distances", "connectivities"):
        if (
            (k in n)
            and isinstance(n[k], (spmatrix, np.ndarray))
            and len(n[k].shape) == 2
        ):
            warn(
                f"Moving element from .uns['neighbors']['{k}'] to .obsp['{k}'].\n\n"
                "This is where adjacency matrices should go now.",
                FutureWarning,
            )
            obsp[k] = n.pop(k)


def _find_sparse_matrices(d: Mapping, n: int, keys: tuple, paths: list):
    """Find paths to sparse matrices with shape (n, n)."""
    for k, v in d.items():
        if isinstance(v, Mapping):
            _find_sparse_matrices(v, n, (*keys, k), paths)
        elif isinstance(v, spmatrix) and v.shape == (n, n):
            paths.append((*keys, k))
    return paths


def _slice_uns_sparse_matrices(uns: MutableMapping, oidx: "Index1d", orig_n_obs: int):
    """slice sparse spatrices of n_obs × n_obs in self.uns"""
    if isinstance(oidx, slice) and len(range(*oidx.indices(orig_n_obs))) == orig_n_obs:
        return uns  # slice of entire dimension is a no-op

    paths = _find_sparse_matrices(uns, orig_n_obs, (), [])

    if not paths:
        return uns

    uns = deepcopy(uns)
    for path in paths:
        str_path = "".join(f"['{key}']" for key in path)
        warn(
            f"During AnnData slicing, found matrix at .uns{str_path} that happens"
            f" to be dimensioned at n_obs×n_obs ({orig_n_obs}×{orig_n_obs}).\n\n"
            "These matrices should now be stored in the .obsp attribute.\n"
            "This slicing behavior will be removed in anndata 0.8.",
            FutureWarning,
        )
        d = reduce(lambda d, k: d[k], path[:-1], uns)
        d[path[-1]] = _subset(d[path[-1]], (oidx, oidx))
    return uns


# This function was adapted from scikit-learn
# github.com/scikit-learn/scikit-learn/blob/master/sklearn/utils/validation.py
def _deprecate_positional_args(func=None, *, version: str = "1.0 (renaming of 0.25)"):
    """Decorator for methods that issues warnings for positional arguments.
    Using the keyword-only argument syntax in pep 3102, arguments after the
    * will issue a warning when passed as a positional argument.

    Parameters
    ----------
    func
        Function to check arguments on.
    version
        The version when positional arguments will result in error.
    """

    def _inner_deprecate_positional_args(f):
        sig = signature(f)
        kwonly_args = []
        all_args = []

        for name, param in sig.parameters.items():
            if param.kind == Parameter.POSITIONAL_OR_KEYWORD:
                all_args.append(name)
            elif param.kind == Parameter.KEYWORD_ONLY:
                kwonly_args.append(name)

        @wraps(f)
        def inner_f(*args, **kwargs):
            extra_args = len(args) - len(all_args)
            if extra_args <= 0:
                return f(*args, **kwargs)

            # extra_args > 0
            args_msg = [
                "{}={}".format(name, arg)
                for name, arg in zip(kwonly_args[:extra_args], args[-extra_args:])
            ]
            args_msg = ", ".join(args_msg)
            warn(
                f"Pass {args_msg} as keyword args. From version {version} passing "
                "these as positional arguments will result in an error",
                FutureWarning,
            )
            kwargs.update(zip(sig.parameters, args))
            return f(**kwargs)

        return inner_f

    if func is not None:
        return _inner_deprecate_positional_args(func)

    return _inner_deprecate_positional_args
