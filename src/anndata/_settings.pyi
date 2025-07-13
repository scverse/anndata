from collections.abc import Callable as Callable
from collections.abc import Generator, Iterable
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Literal

@dataclass
class SettingsManager:
    __doc_tmpl__: str = ...
    remove_unused_categories: bool = True
    check_uniqueness: bool = True
    allow_write_nullable_strings: bool = False
    zarr_write_format: Literal[2, 3] = 2
    use_sparse_array_on_read: bool = False
    min_rows_for_chunked_h5_copy: int = 1000

    def __setattr__(self, option: str, val: object) -> None: ...
    def __getattr__(self, option: str) -> object: ...
    def reset(self, option: Iterable[str] | str) -> None: ...
    @contextmanager
    def override(self, **overrides) -> Generator[None]: ...

settings: SettingsManager
