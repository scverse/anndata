from collections.abc import Callable as Callable
from collections.abc import Generator, Iterable
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Literal, TypeVar

_T = TypeVar("_T")

@dataclass
class SettingsManager:
    __doc_tmpl__: str = ...
    def describe(
        self,
        option: str | Iterable[str] | None = None,
        *,
        should_print_description: bool = True,
        as_rst: bool = False,
    ) -> str: ...
    def deprecate(
        self, option: str, removal_version: str, message: str | None = None
    ) -> None: ...
    def register(
        self,
        option: str,
        *,
        default_value: _T,
        description: str,
        validate: Callable[[_T], None],
        option_type: object | None = None,
        get_from_env: Callable[[str, _T], _T] = ...,
    ) -> None: ...
    def __setattr__(self, option: str, val: object) -> None: ...
    def __getattr__(self, option: str) -> object: ...
    def __dir__(self) -> Iterable[str]: ...
    def reset(self, option: Iterable[str] | str) -> None: ...
    @contextmanager
    def override(self, **overrides) -> Generator[None]: ...
    @property
    def __doc__(self): ...

class _AnnDataSettingsManager(SettingsManager):
    remove_unused_categories: bool = True
    check_uniqueness: bool = True
    allow_write_nullable_strings: bool = False
    zarr_write_format: Literal[2, 3] = 2
    use_sparse_array_on_read: bool = False
    min_rows_for_chunked_h5_copy: int = 1000

settings: _AnnDataSettingsManager
