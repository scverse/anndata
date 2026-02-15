from collections.abc import Callable as Callable
from collections.abc import Generator, Iterable
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Literal, Self

@dataclass
class SettingsManager[T]:
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
        default_value: T,
        description: str,
        validate: Callable[[T, Self], None],
        option_type: object | None = None,
        get_from_env: Callable[[str, T], T] = ...,
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
    copy_on_write_X: bool = False
    allow_write_nullable_strings: bool | None = None
    zarr_write_format: Literal[2, 3] = 2
    use_sparse_array_on_read: bool = False
    min_rows_for_chunked_h5_copy: int = 1000
    disallow_forward_slash_in_h5ad: bool = False
    write_csr_csc_indices_with_min_possible_dtype: bool = False
    auto_shard_zarr_v3: bool = False
    repr_html_enabled: bool = True
    repr_html_fold_threshold: int = 5
    repr_html_max_depth: int = 3
    repr_html_max_items: int = 200
    repr_html_max_categories: int = 20
    repr_html_unique_limit: int = 1_000_000
    repr_html_dataframe_expand: bool = False
    repr_html_max_field_width: int = 400
    repr_html_type_width: int = 220
    repr_html_max_lazy_categories: int = 100
    repr_html_max_readme_size: int = 100_000

settings: _AnnDataSettingsManager
