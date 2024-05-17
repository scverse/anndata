from __future__ import annotations

import os
import sys
from codecs import decode
from collections.abc import Mapping
from contextlib import AbstractContextManager
from dataclasses import dataclass, field
from functools import singledispatch, wraps
from inspect import Parameter, signature
from pathlib import Path
from typing import Any, Union
from warnings import warn

import h5py
import numpy as np
import pandas as pd
from packaging.version import Version
from scipy.sparse import issparse, spmatrix

from .exceptiongroups import add_note  # noqa: F401


class Empty:
    pass


Index1D = Union[slice, int, str, np.int64, np.ndarray]
Index = Union[Index1D, tuple[Index1D, Index1D], spmatrix]
H5Group = h5py.Group
H5Array = h5py.Dataset


#############################
# stdlib
#############################


if sys.version_info >= (3, 11):
    from contextlib import chdir
else:

    @dataclass
    class chdir(AbstractContextManager):
        path: Path
        _old_cwd: list[Path] = field(default_factory=list)

        def __enter__(self) -> None:
            self._old_cwd.append(Path())
            os.chdir(self.path)

        def __exit__(self, *_exc_info) -> None:
            os.chdir(self._old_cwd.pop())


if sys.version_info >= (3, 10):
    from itertools import pairwise
else:

    def pairwise(iterable):
        from itertools import tee

        a, b = tee(iterable)
        next(b, None)
        return zip(a, b)


#############################
# Optional deps
#############################

try:
    from zarr.core import Array as ZarrArray
    from zarr.hierarchy import Group as ZarrGroup
except ImportError:

    class ZarrArray:
        @staticmethod
        def __repr__():
            return "mock zarr.core.Array"

    class ZarrGroup:
        @staticmethod
        def __repr__():
            return "mock zarr.core.Group"


try:
    import awkward

    AwkArray = awkward.Array

except ImportError:

    class AwkArray:
        @staticmethod
        def __repr__():
            return "mock awkward.highlevel.Array"


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
    from cupy import ndarray as CupyArray
    from cupyx.scipy.sparse import (
        csc_matrix as CupyCSCMatrix,
    )
    from cupyx.scipy.sparse import (
        csr_matrix as CupyCSRMatrix,
    )
    from cupyx.scipy.sparse import (
        spmatrix as CupySparseMatrix,
    )
except ImportError:

    class CupySparseMatrix:
        @staticmethod
        def __repr__():
            return "mock cupyx.scipy.sparse.spmatrix"

    class CupyCSRMatrix:
        @staticmethod
        def __repr__():
            return "mock cupyx.scipy.sparse.csr_matrix"

    class CupyCSCMatrix:
        @staticmethod
        def __repr__():
            return "mock cupyx.scipy.sparse.csc_matrix"

    class CupyArray:
        @staticmethod
        def __repr__():
            return "mock cupy.ndarray"


#############################
# IO helpers
#############################


@singledispatch
def _read_attr(attrs: Mapping, name: str, default: Any | None = Empty):
    if default is Empty:
        return attrs[name]
    else:
        return attrs.get(name, default=default)


@_read_attr.register(h5py.AttributeManager)
def _read_attr_hdf5(
    attrs: h5py.AttributeManager, name: str, default: Any | None = Empty
):
    """
    Read an HDF5 attribute and perform all necessary conversions.

    At the moment, this only implements conversions for string attributes, other types
    are passed through. String conversion is needed compatibility with other languages.
    For example Julia's HDF5.jl writes string attributes as fixed-size strings, which
    are read as bytes by h5py.
    """
    if name not in attrs and default is not Empty:
        return default
    attr = attrs[name]
    attr_id = attrs.get_id(name)
    dtype = h5py.check_string_dtype(attr_id.dtype)
    if dtype is None:
        return attr
    else:
        if dtype.length is None:  # variable-length string, no problem
            return attr
        elif len(attr_id.shape) == 0:  # Python bytestring
            return attr.decode("utf-8")
        else:  # NumPy array
            return [decode(s, "utf-8") for s in attr]


def _from_fixed_length_strings(value):
    """\
    Convert from fixed length strings to unicode.

    For backwards compatibility with older h5ad and zarr files.
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
        if issubclass(np.dtype(dt_type).type, np.bytes_):
            dt_list[1] = f"U{int(dt_type[2:])}"
        elif is_annotated or np.issubdtype(np.dtype(dt_type), np.str_):
            dt_list[1] = "O"  # Assumption that it’s a vlen str
        new_dtype.append(tuple(dt_list))
    return value.astype(new_dtype)


def _decode_structured_array(
    arr: np.ndarray, dtype: np.dtype | None = None, copy: bool = False
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


def _clean_uns(adata: AnnData):  # noqa: F821
    """
    Compat function for when categorical keys were stored in uns.
    This used to be buggy because when storing categorical columns in obs and var with
    the same column name, only one `<colname>_categories` is retained.
    """
    k_to_delete = set()
    for cats_name, cats in adata.uns.items():
        if not cats_name.endswith("_categories"):
            continue
        name = cats_name.replace("_categories", "")
        # fix categories with a single category
        if isinstance(cats, (str, int)):
            cats = [cats]
        for ann in [adata.obs, adata.var]:
            if name not in ann:
                continue
            codes: np.ndarray = ann[name].values
            # hack to maybe find the axis the categories were for
            if not np.all(codes < len(cats)):
                continue
            ann[name] = pd.Categorical.from_codes(codes, cats)
            k_to_delete.add(cats_name)
    for cats_name in k_to_delete:
        del adata.uns[cats_name]


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
                f"{name}={arg}"
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


def _transpose_by_block(dask_array: DaskArray) -> DaskArray:
    import dask.array as da

    b = dask_array.blocks
    b_raveled = b.ravel()
    block_layout = np.zeros(b.shape, dtype=object)

    for i in range(block_layout.size):
        block_layout.flat[i] = b_raveled[i].map_blocks(
            lambda x: x.T, chunks=b_raveled[i].chunks[::-1]
        )

    return da.block(block_layout.T.tolist())


def _safe_transpose(x):
    """Safely transpose x

    This is a workaround for: https://github.com/scipy/scipy/issues/19161
    """

    if isinstance(x, DaskArray) and issparse(x._meta):
        return _transpose_by_block(x)
    else:
        return x.T


def _map_cat_to_str(cat: pd.Categorical) -> pd.Categorical:
    if Version(pd.__version__) >= Version("2.1"):
        # Argument added in pandas 2.1
        return cat.map(str, na_action="ignore")
    else:
        return cat.map(str)
