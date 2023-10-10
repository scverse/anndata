from __future__ import annotations

from functools import wraps
from typing import Callable, Literal
from warnings import warn

from packaging import version
import h5py

from .._core.sparse_dataset import BaseCompressedSparseDataset
from anndata.compat import H5Group, ZarrGroup, add_note

# For allowing h5py v3
# https://github.com/scverse/anndata/issues/442
H5PY_V3 = version.parse(h5py.__version__).major >= 3

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
        raise TypeError(f"{key} of type {typ} is an invalid key. Should be str.")


# -------------------------------------------------------------------------------
# Generic functions
# -------------------------------------------------------------------------------


def read_attribute(*args, **kwargs):
    from .specs import read_elem

    warn(
        "This internal function has been deprecated, please use read_elem instead",
        DeprecationWarning,
    )
    return read_elem(*args, **kwargs)


def write_attribute(*args, **kwargs):
    from .specs import write_elem

    warn(
        "This internal function has been deprecated, please use write_elem instead",
        DeprecationWarning,
    )
    return write_elem(*args, **kwargs)


# -------------------------------------------------------------------------------
# Errors handling
# -------------------------------------------------------------------------------
# TODO: Is there a consistent way to do this which just modifies the previously
# thrown error? Could do a warning?


class AnnDataReadError(OSError):
    """Error caused while trying to read in AnnData."""

    pass


def _get_parent(elem):
    try:
        import zarr
    except ImportError:
        zarr = None
    if zarr and isinstance(elem, (zarr.Group, zarr.Array)):
        parent = elem.store  # Not sure how to always get a name out of this
    elif isinstance(elem, BaseCompressedSparseDataset):
        parent = elem.group.file.name
    else:
        parent = elem.file.name
    return parent


def re_raise_error(e, elem, key, op=Literal["read", "writ"]):
    if any(
        f"Error raised while {op}ing key" in note
        for note in getattr(e, "__notes__", [])
    ):
        raise
    else:
        parent = _get_parent(elem)
        add_note(
            e,
            f"Error raised while {op}ing key {key!r} of {type(elem)} to " f"{parent}",
        )
        raise e


def report_read_key_on_error(func):
    """\
    A decorator for zarr element reading which makes keys involved in errors get reported.

    Example
    -------
    >>> import zarr
    >>> @report_read_key_on_error
    ... def read_arr(group):
    ...     raise NotImplementedError()
    >>> z = zarr.open("tmp.zarr")
    >>> z["X"] = [1, 2, 3]
    >>> read_arr(z["X"])  # doctest: +SKIP
    """

    @wraps(func)
    def func_wrapper(*args, **kwargs):
        from anndata._io.specs import Reader

        # Figure out signature (method vs function) by going through args
        for elem in args:
            if not isinstance(elem, Reader):
                break
        try:
            return func(*args, **kwargs)
        except Exception as e:
            re_raise_error(e, elem, elem.name, "read")

    return func_wrapper


def report_write_key_on_error(func):
    """\
    A decorator for zarr element reading which makes keys involved in errors get reported.

    Example
    -------
    >>> import zarr
    >>> @report_write_key_on_error
    ... def write_arr(group, key, val):
    ...     raise NotImplementedError()
    >>> z = zarr.open("tmp.zarr")
    >>> X = [1, 2, 3]
    >>> write_arr(z, "X", X)  # doctest: +SKIP
    """

    @wraps(func)
    def func_wrapper(*args, **kwargs):
        from anndata._io.specs import Writer

        # Figure out signature (method vs function) by going through args
        for i in range(len(args)):
            elem = args[i]
            key = args[i + 1]
            if not isinstance(elem, Writer):
                break
        try:
            return func(*args, **kwargs)
        except Exception as e:
            re_raise_error(e, elem, key, "writ")

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
            raise ValueError(f"{what} has both legacy and current raw formats.")
        return modern_raw

    raw = {}
    if "X" in attrs and "raw.X" in f:
        raw["X"] = read_attr(f["raw.X"])
    if "var" in attrs and "raw.var" in f:
        raw["var"] = read_df(f["raw.var"])  # Backwards compat
    if "varm" in attrs and "raw.varm" in f:
        raw["varm"] = read_attr(f["raw.varm"])
    return raw
