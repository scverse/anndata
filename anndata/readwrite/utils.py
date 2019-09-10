from functools import wraps
from .. import h5py as patched_h5py

# -------------------------------------------------------------------------------
# Type conversion
# -------------------------------------------------------------------------------


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


class AnnDataReadError(OSError):
    """Error caused while trying to read in AnnData."""

    pass


def report_key_on_error(func):
    """
    A decorator for zarr element reading which makes keys involved in errors get reported.

    Example
    -------
    >>> import zarr
    >>> @report_key_on_error
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
                try:
                    import zarr
                except ImportError:
                    zarr = None
                if zarr and isinstance(elem, (zarr.Group, zarr.Array)):
                    parent = (
                        elem.store
                    )  # Not sure how to always get a name out of this
                elif isinstance(elem, patched_h5py.Group):
                    parent = elem.h5py_group.file.name
                else:
                    parent = elem.file.name
                raise AnnDataReadError(
                    f"Above error raised while reading key '{elem.name}' of type {type(elem)} from {parent}."
                )

    return func_wrapper
