import sys


def add_note(err: BaseException, msg: str) -> BaseException:
    """
    Adds a note to an exception inplace and returns it.
    """
    if sys.version_info < (3, 11):
        err.__notes__ = getattr(err, "__notes__", []) + [msg]
    else:
        err.add_note(msg)
    return err
