from functools import wraps
from typing import Optional, TypeVar


def format_msg(elem_name: Optional[str]) -> str:
    if elem_name is not None:
        return f"Error raised from element {elem_name!r}."
    else:
        return ""


# TODO: it would be better to modify the other exception
def report_name(func):
    """Report name of element being tested if test fails."""

    @wraps(func)
    def func_wrapper(*args, _elem_name: Optional[str] = None, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            if _elem_name is not None and not hasattr(e, "_name_attached"):
                msg = format_msg(_elem_name)
                args = list(e.args)
                if len(args) == 0:
                    args = [msg]
                else:
                    args[0] = f"{args[0]}\n\n{msg}"
                e.args = tuple(args)
                e._name_attached = True
            raise e

    return func_wrapper
