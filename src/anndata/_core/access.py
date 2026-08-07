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
    def _parent_el(self):
        return reduce(
            lambda d, k: d[k], self.keys[:-1], getattr(self.parent, self.attrname)
        )

    def _get_in(self, adata: AnnData):
        return reduce(lambda d, k: d[k], self.keys, getattr(adata, self.attrname))

    def get(self):
        """Get referenced value in self.parent."""
        return self._get_in(self.parent)

    def get_unsubset(self):
        """Get referenced value in the AnnData `self.parent` is a view of.

        That is, the element before the view’s subsetting was applied.
        `None` if `self.parent` is not a view.
        """
        if (adata_ref := self.parent._adata_ref) is None:
            return None
        return self._get_in(adata_ref)

    def set(self, val):
        """Set referenced value in self.parent."""
        self._parent_el[self.keys[-1]] = val

    def delete(self):
        del self._parent_el[self.keys[-1]]
