from functools import wraps, singledispatch
from typing import Mapping, Any

from .. import AnnData, Raw
from ..core.sparsedataset import SparseDataset

# -------------------------------------------------------------------------------
# Type conversion
# -------------------------------------------------------------------------------


# Could be numba'd if it returned tuples instead of slices
def idx_chunks_along_axis(shape: tuple, axis: int, chunk_size: int):
    """
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
    """Check whether string is float.

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
    """Check whether string is integer.
    """
    try:
        int(string)
        return True
    except ValueError:
        return False


def convert_bool(string):
    """Check whether string is boolean.
    """
    if string == 'True':
        return True, True
    elif string == 'False':
        return True, False
    else:
        return False, False


def convert_string(string):
    """Convert string to int, float or bool.
    """
    if is_int(string):
        return int(string)
    elif is_float(string):
        return float(string)
    elif convert_bool(string)[0]:
        return convert_bool(string)[1]
    elif string == 'None':
        return None
    else:
        return string


# -------------------------------------------------------------------------------
# Generic functions
# -------------------------------------------------------------------------------


@singledispatch
def write_attribute(*args, **kwargs):
    raise NotImplementedError(
        "Unrecognized argument types for `write_attribute`."
    )


@singledispatch
def read_attribute(*args, **kwargs):
    raise NotImplementedError(
        "Unrecognized argument types for `read_attribute`."
    )


@read_attribute.register(type(None))
def read_attribute_none(value) -> None:
    return None


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
    elif isinstance(elem, SparseDataset):
        parent = elem.group.file.name
    else:
        parent = elem.file.name
    return parent


def report_read_key_on_error(func):
    """
    A decorator for zarr element reading which makes keys involved in errors get reported.

    Example
    -------
    >>> import zarr
    >>> @report_read_key_on_error
    ... def read_arr(group):
    ...     raise NotImplementedError()
    >>> z = zarr.open("tmp.zarr")
    >>> z["X"] = [1, 2, 3]
    >>> read_arr(z["X"])
    """

    @wraps(func)
    def func_wrapper(elem, *args, **kwargs):
        try:
            return func(elem, *args, **kwargs)
        except Exception as e:
            if isinstance(e, AnnDataReadError):
                raise e
            else:
                parent = _get_parent(elem)
                raise AnnDataReadError(
                    f"Above error raised while reading key '{elem.name}' of "
                    f"type {type(elem)} from {parent}."
                )

    return func_wrapper


def report_write_key_on_error(func):
    """
    A decorator for zarr element reading which makes keys involved in errors get reported.

    Example
    -------
    >>> import zarr
    >>> @report_write_key_on_error
    ... def write_arr(group, key, val):
    ...     raise NotImplementedError()
    >>> z = zarr.open("tmp.zarr")
    >>> X = [1, 2, 3]
    >>> write_arr(z, "X", X)
    """

    @wraps(func)
    def func_wrapper(elem, key, val, *args, **kwargs):
        try:
            return func(elem, key, val, *args, **kwargs)
        except Exception as e:
            parent = _get_parent(elem)
            raise type(e)(
                f"{e}\n\n"
                f"Above error raised while writing key '{key}' of {type(elem)}"
                f" from {parent}."
            ) from e

    return func_wrapper


def _read_legacy_raw_into(f, raw, read_df, read_attr, filename=None):
    what = f'File {filename}' if filename else 'Store'
    # Backwards compat for reading legacy raw
    assert not (
        raw and any(k.startswith("raw.") for k in f)
    ), f"{what} has both legacy and current raw formats."
    raw = {}
    if "raw.var" in f:
        raw["var"] = read_df(f["raw.var"])  # Backwards compat
    if "raw.varm" in f:
        raw["varm"] = read_attr(f["raw.varm"])
    if "raw.X" in f:
        raw["X"] = read_attr(f["raw.X"])
