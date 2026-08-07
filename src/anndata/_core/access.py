from __future__ import annotations

from functools import reduce
from typing import TYPE_CHECKING, NamedTuple

if TYPE_CHECKING:
    from anndata import AnnData


class ElementRef[K: (str, str | None)](NamedTuple):
    parent: AnnData
    attrname: str
    keys: tuple[K, ...] = ()

    def __str__(self) -> str:
        return f".{self.attrname}" + "".join(f"[{x!r}]" for x in self.keys)

    @property
    def _container(self):
        """The collection holding `keys[-1]`, i.e. one level up from the value."""
        return reduce(
            lambda d, k: d[k], self.keys[:-1], getattr(self.parent, self.attrname)
        )

    @property
    def unsubset_element(self):
        """Get referenced value as it exists before `self.parent`’s subsetting.

        `None` if `self.parent` is not a view, and so has no unsubset value.
        """
        if (adata_ref := self.parent._adata_ref) is None:
            return None
        return self._get_in(adata_ref)

    def _get_in(self, adata: AnnData):
        return reduce(lambda d, k: d[k], self.keys, getattr(adata, self.attrname))

    def get(self):
        """Get referenced value in self.parent."""
        return self._get_in(self.parent)

    def set(self, val):
        """Set referenced value in self.parent."""
        self._container[self.keys[-1]] = val

    def delete(self):
        del self._container[self.keys[-1]]
