from abc import ABC
from pathlib import PurePath
from typing import Union

try:
    from os import PathLike, fspath
except ImportError:
    class PathLike(ABC):
        @classmethod
        def __subclasshook__(cls, subclass):
            return issubclass(subclass, PurePath) or hasattr(subclass, '__fspath__')

    def fspath(path: Union[PathLike, str, bytes]) -> Union[str, bytes]:
        """Copied from https://github.com/python/cpython/blob/3.6/Lib/os.py"""
        if isinstance(path, (str, bytes)):
            return path
        if isinstance(path, PurePath):
            return str(path)

        path_type = type(path)
        try:
            path_repr = path_type.__fspath__(path)
        except AttributeError:
            if hasattr(path_type, '__fspath__'):
                raise
            else:
                raise TypeError("expected str, bytes or os.PathLike object, not " + path_type.__name__)
        if isinstance(path_repr, (str, bytes)):
            return path_repr
        else:
            raise TypeError(
                "expected {}.__fspath__() to return str or bytes, not {}"
                .format(path_type.__name__, type(path_repr).__name__)
            )
