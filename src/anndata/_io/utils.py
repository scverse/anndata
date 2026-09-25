from __future__ import annotations

import re
import sys
from collections.abc import Callable, Mapping
from functools import WRAPPER_ASSIGNMENTS, cache, wraps
from itertools import pairwise
from typing import TYPE_CHECKING, Literal, cast

import numpy as np
import pandas as pd

from .._core.sparse_dataset import BaseCompressedSparseDataset
from .._write_compat import WriteCompat
from ..utils import warn

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, Literal

    import h5py
    import zarr
    from pandas.core.dtypes.dtypes import BaseMaskedDtype

    from .._types import StorageType, _ArrayStorageType, _WriteInternal
    from .specs.registry import Writer

    Storage = StorageType | BaseCompressedSparseDataset


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


@cache
def pandas_nullable_dtype(dtype: np.dtype) -> BaseMaskedDtype:
    """Infer nullable dtype from numpy dtype.

    There is no public pandas API for this, so this is the cleanest way.
    See <https://github.com/pandas-dev/pandas/issues/63608>
    """
    try:
        from pandas.core.dtypes.dtypes import BaseMaskedDtype
    except ImportError:
        pass
    else:
        if hasattr(BaseMaskedDtype, "from_numpy_dtype"):
            return BaseMaskedDtype.from_numpy_dtype(dtype)

    array_type: type[pd.arrays.BooleanArray | pd.arrays.IntegerArray]
    match dtype.kind:
        case "b":
            array_type = pd.arrays.BooleanArray
        case "i" | "u":
            array_type = pd.arrays.IntegerArray
        case _:
            raise NotImplementedError
    return array_type(np.ones(1, dtype), np.ones(1, bool)).dtype


# -------------------------------------------------------------------------------
# Generic functions
# -------------------------------------------------------------------------------


def read_attribute(*args, **kwargs):
    from .specs import read_elem

    msg = "This internal function has been deprecated, please use read_elem instead"
    warn(msg, FutureWarning)
    return read_elem(*args, **kwargs)


def write_attribute(*args, **kwargs):
    from .specs import write_elem

    msg = "This internal function has been deprecated, please use write_elem instead"
    warn(msg, FutureWarning)
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
    group = store.group if isinstance(store, BaseCompressedSparseDataset) else store
    path = group.name or "??"  # can be None
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
        __tracebackhide__ = True
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


def _check_has_no_slash_key(attr: str, elem: object, *, compat: WriteCompat) -> None:
    """Only attempt to write slash keys where people rely on it for backwards compatibility."""

    if attr in {"obs", "var", "uns", "raw"}:
        return  # `write_elem` checks these against `settings.disallow_forward_slash_in_h5ad`
    if compat >= WriteCompat.V0_14:
        return  # keys get escaped, see `escape_key`
    assert isinstance(elem, Mapping)
    if any("/" in k for k in elem if k not in {"/", None}):
        msg = f"Forward slashes are not allowed in keys in {attr}"
        raise ValueError(msg)


# The characters Windows forbids in file names, per “Naming Files, Paths, and Namespaces”:
# https://learn.microsoft.com/en-us/windows/win32/fileio/naming-a-file#naming-conventions
# That is a superset of what POSIX (`/` and NUL) and macOS (those plus `:`, the classic
# HFS separator) forbid, so escaping these makes a key safe as a path segment on any of
# them – which matters for zarr stores that map keys to real files.
# `%` is in there because it introduces an escape sequence.
UNSAFE_KEY_CHARS = frozenset('%/\\:*?"<>|') | frozenset(map(chr, range(32)))
_ESCAPE_RE = re.compile("%([0-9A-F]{2})")

if sys.version_info >= (3, 13):
    from ntpath import isreserved as _is_reserved_win
else:  # pragma: no cover
    _WIN_DEVICE_NAMES = frozenset(
        {"CON", "PRN", "AUX", "NUL", "CONIN$", "CONOUT$"}
        | {f"{dev}{n}" for dev in ("COM", "LPT") for n in "123456789¹²³"}
    )

    def _is_reserved_win(name: str) -> bool:
        """Backport of `ntpath.isreserved`, minus its check for the characters we escape."""
        if name[-1:] in (".", " "):  # trailing dots and spaces are reserved
            return name not in (".", "..")
        return name.partition(".")[0].rstrip(" ").upper() in _WIN_DEVICE_NAMES


def _is_reserved_key(key: str) -> bool:
    """Whether `key` is unusable as a child key even though its characters are safe.

    The first four are reserved by the `zarr` v3 spec
    (https://zarr-specs.readthedocs.io/en/latest/v3/core/index.html#node-names),
    the last one catches Windows device names like `CON` and trailing dots/spaces.
    """
    return (
        not key
        or set(key) == {"."}  # `.`, `..`, `...`, …
        or key.startswith("__")
        or key == "zarr.json"
        or _is_reserved_win(key)
    )


def key_needs_escaping(key: str) -> bool:
    """Whether `key` is unusable as a child key as-is."""
    return not UNSAFE_KEY_CHARS.isdisjoint(key) or _is_reserved_key(key)


def escape_key(key: str) -> str:
    """Make `key` usable as a child key.

    Characters a file system may choke on are percent-escaped,
    and a key that is reserved even without them is wrapped in a `%` on either end.
    That `%` cannot be mistaken for an escape sequence,
    as those are always a `%` followed by two hexadecimal digits.
    Wrapping unconditionally would be just as readable back,
    but we only do it where needed, so that keys stay legible on disk.

    >>> escape_key("foo/bar")
    'foo%2Fbar'
    >>> escape_key("CON")  # a Windows device name
    '%CON%'
    >>> escape_key("NUL.txt")  # `NUL` stays the device whatever follows the period
    '%NUL.txt%'
    >>> [unescape_key(escape_key(k)) for k in ["100%/day", "CON", "..", "__x", ""]]
    ['100%/day', 'CON', '..', '__x', '']
    """
    escaped = "".join(f"%{ord(c):02X}" if c in UNSAFE_KEY_CHARS else c for c in key)
    # escaping can neither create nor remove a reserved key, so this is stable.
    # Wrapping breaks every rule at once: the result is non-empty, is not all periods,
    # starts with neither `__` nor a device name, ends in neither a period nor a space,
    # and is not `zarr.json` – so it is never reserved itself, whatever `key` was.
    return f"%{escaped}%" if _is_reserved_key(escaped) else escaped


def unescape_key(key: str) -> str:
    """Invert :func:`escape_key`."""
    if key.endswith("%"):  # a `%` can never be the 2nd or 3rd character of an escape
        key = key[1:-1]
    return _ESCAPE_RE.sub(lambda m: chr(int(m[1], 16)), key)


# -------------------------------------------------------------------------------
# Common h5ad/zarr stuff
# -------------------------------------------------------------------------------


def _read_legacy_raw(
    f: zarr.Group | h5py.Group,
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


def zero_dim_array_as_scalar[S: StorageType, T: np.ndarray | _ArrayStorageType](
    func: _WriteInternal[S, T],
) -> _WriteInternal[S, T]:
    """\
    A decorator for write_elem implementations of arrays where zero-dimensional arrays need special handling.
    """

    @wraps(func, assigned=(*WRAPPER_ASSIGNMENTS, "__defaults__", "__kwdefaults__"))
    def func_wrapper(
        f: S,
        k: str,
        elem: T,
        *,
        _writer: Writer,
        dataset_kwargs: Mapping[str, Any],
    ) -> None:
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
