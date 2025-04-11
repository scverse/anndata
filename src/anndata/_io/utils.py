from __future__ import annotations

from functools import WRAPPER_ASSIGNMENTS, wraps
from itertools import pairwise
from typing import TYPE_CHECKING, cast
from warnings import warn

import h5py
from packaging.version import Version

from .._core.sparse_dataset import BaseCompressedSparseDataset

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping
    from typing import Any, Literal

    from .._types import StorageType, _WriteInternal
    from ..compat import H5Group, ZarrGroup
    from ..typing import RWAble
    from .specs.registry import Writer

    Storage = StorageType | BaseCompressedSparseDataset

# For allowing h5py v3
# https://github.com/scverse/anndata/issues/442
H5PY_V3 = Version(h5py.__version__).major >= 3

# -------------------------------------------------------------------------------
# Type conversion
# -------------------------------------------------------------------------------


# Could be numba’d if it returned tuples instead of slices
def idx_chunks_along_axis(shape: tuple, axis: int, chunk_size: int):
    """\
    Gives indexer tuples chunked along an axis.

    Params
    ------
    shape
        Shape of array to be chunked
    axis
        Axis to chunk along
    chunk_size
        Size of chunk along axis

    Returns
    -------
    An iterator of tuples for indexing into an array of passed shape.
    """
    total = shape[axis]
    cur = 0
    mutable_idx = [slice(None) for i in range(len(shape))]
    while cur + chunk_size < total:
        mutable_idx[axis] = slice(cur, cur + chunk_size)
        yield tuple(mutable_idx)
        cur += chunk_size
    mutable_idx[axis] = slice(cur, None)
    yield tuple(mutable_idx)


def is_float(string):
    """\
    Check whether string is float.

    See also
    --------
    http://stackoverflow.com/questions/736043/checking-if-a-string-can-be-converted-to-float-in-python
    """
    try:
        float(string)
        return True
    except ValueError:
        return False


def is_int(string):
    """Check whether string is integer."""
    try:
        int(string)
        return True
    except ValueError:
        return False


def convert_bool(string):
    """Check whether string is boolean."""
    if string == "True":
        return True, True
    elif string == "False":
        return True, False
    else:
        return False, False


def convert_string(string):
    """Convert string to int, float or bool."""
    if is_int(string):
        return int(string)
    elif is_float(string):
        return float(string)
    elif convert_bool(string)[0]:
        return convert_bool(string)[1]
    elif string == "None":
        return None
    else:
        return string


def check_key(key):
    """Checks that passed value is a valid h5py key.

    Should convert it if there is an obvious conversion path, error otherwise.
    """
    typ = type(key)
    if issubclass(typ, str):
        return str(key)
    # TODO: Should I try to decode bytes? It's what h5py would do,
    # but it will be read out as a str.
    # elif issubclass(typ, bytes):
    # return key
    else:
        msg = f"{key} of type {typ} is an invalid key. Should be str."
        raise TypeError(msg)


# -------------------------------------------------------------------------------
# Generic functions
# -------------------------------------------------------------------------------


def read_attribute(*args, **kwargs):
    from .specs import read_elem

    msg = "This internal function has been deprecated, please use read_elem instead"
    warn(msg, FutureWarning, stacklevel=2)
    return read_elem(*args, **kwargs)


def write_attribute(*args, **kwargs):
    from .specs import write_elem

    msg = "This internal function has been deprecated, please use write_elem instead"
    warn(msg, FutureWarning, stacklevel=2)
    return write_elem(*args, **kwargs)


# -------------------------------------------------------------------------------
# Errors handling
# -------------------------------------------------------------------------------
# TODO: Is there a consistent way to do this which just modifies the previously
# thrown error? Could do a warning?


class AnnDataReadError(OSError):
    """Error caused while trying to read in AnnData."""


def _get_display_path(store: Storage) -> str:
    """Return an absolute path of an element (always starts with “/”)."""
    if isinstance(store, BaseCompressedSparseDataset):
        store = store.group
    path = store.name or "??"  # can be None
    return f"/{path.removeprefix('/')}"


def add_key_note(
    e: BaseException, store: Storage, path: str, key: str, op: Literal["read", "writ"]
) -> None:
    if any(
        f"Error raised while {op}ing key" in note
        for note in getattr(e, "__notes__", [])
    ):
        return

    dir = "to" if op == "writ" else "from"
    msg = f"Error raised while {op}ing key {key!r} of {type(store)} {dir} {path}"
    e.add_note(msg)


def report_read_key_on_error(func):
    """\
    A decorator for hdf5/zarr element reading which makes keys involved in errors get reported.

    Example
    -------
    >>> import zarr
    >>> import numpy as np
    >>> @report_read_key_on_error
    ... def read_arr(group):
    ...     raise NotImplementedError()
    >>> z = zarr.open("tmp.zarr", mode="w")
    >>> z["X"] = np.array([1, 2, 3])
    >>> read_arr(z["X"])  # doctest: +SKIP
    """

    @wraps(func)
    def func_wrapper(*args, **kwargs):
        from anndata._io.specs import Reader

        # Figure out signature (method vs function) by going through args
        for arg in args:
            if not isinstance(arg, Reader):
                store = cast("Storage", arg)
                break
        else:
            msg = "No element found in args."
            raise ValueError(msg)
        try:
            return func(*args, **kwargs)
        except Exception as e:
            path, key = _get_display_path(store).rsplit("/", 1)
            add_key_note(e, store, path or "/", key, "read")
            raise

    return func_wrapper


def report_write_key_on_error(func):
    """\
    A decorator for hdf5/zarr element writing which makes keys involved in errors get reported.

    Example
    -------
    >>> import zarr
    >>> @report_write_key_on_error
    ... def write_arr(group, key, val):
    ...     raise NotImplementedError()
    >>> z = zarr.open("tmp.zarr", mode="w")
    >>> X = [1, 2, 3]
    >>> write_arr(z, "X", X)  # doctest: +SKIP
    """

    @wraps(func)
    def func_wrapper(*args, **kwargs):
        from anndata._io.specs import Writer

        # Figure out signature (method vs function) by going through args
        for arg, _key in pairwise(args):
            key = _key
            if not isinstance(arg, Writer):
                store = cast("Storage", arg)
                break
        else:
            msg = "No element found in args."
            raise ValueError(msg)
        try:
            return func(*args, **kwargs)
        except Exception as e:
            path = _get_display_path(store)
            add_key_note(e, store, path, key, "writ")
            raise

    return func_wrapper


# -------------------------------------------------------------------------------
# Common h5ad/zarr stuff
# -------------------------------------------------------------------------------


def _read_legacy_raw(
    f: ZarrGroup | H5Group,
    modern_raw,  # TODO: type
    read_df: Callable,
    read_attr: Callable,
    *,
    attrs=("X", "var", "varm"),
) -> dict:
    """\
    Backwards compat for reading legacy raw.
    Makes sure that no modern raw group coexists with legacy raw.* groups.
    """
    if modern_raw:
        if any(k.startswith("raw.") for k in f):
            what = f"File {f.filename}" if hasattr(f, "filename") else "Store"
            msg = f"{what} has both legacy and current raw formats."
            raise ValueError(msg)
        return modern_raw

    raw = {}
    if "X" in attrs and "raw.X" in f:
        raw["X"] = read_attr(f["raw.X"])
    if "var" in attrs and "raw.var" in f:
        raw["var"] = read_df(f["raw.var"])  # Backwards compat
    if "varm" in attrs and "raw.varm" in f:
        raw["varm"] = read_attr(f["raw.varm"])
    return raw


def zero_dim_array_as_scalar(func: _WriteInternal):
    """\
    A decorator for write_elem implementations of arrays where zero-dimensional arrays need special handling.
    """

    @wraps(func, assigned=(*WRAPPER_ASSIGNMENTS, "__defaults__", "__kwdefaults__"))
    def func_wrapper(
        f: StorageType,
        k: str,
        elem: RWAble,
        *,
        _writer: Writer,
        dataset_kwargs: Mapping[str, Any],
    ):
        if elem.shape == ():
            _writer.write_elem(f, k, elem[()], dataset_kwargs=dataset_kwargs)
        else:
            func(f, k, elem, _writer=_writer, dataset_kwargs=dataset_kwargs)

    return func_wrapper


def no_write_dataset_2d(write):
    def raise_error_if_dataset_2d_present(store, adata, *args, **kwargs):
        from anndata.experimental.backed._compat import has_dataset_2d

        if has_dataset_2d(adata):
            msg = (
                "Writing AnnData objects with a Dataset2D not supported yet. "
                "Please use `ds.to_memory` to bring the dataset into memory. "
                "Note that if you have generated this object by concatenating several `AnnData` objects"
                "the original types may be lost."
            )
            raise NotImplementedError(msg)
        return write(store, adata, *args, **kwargs)

    return raise_error_if_dataset_2d_present
