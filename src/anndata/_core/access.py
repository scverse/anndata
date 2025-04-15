from __future__ import annotations

from functools import reduce
from typing import TYPE_CHECKING, NamedTuple

if TYPE_CHECKING:
    from anndata import AnnData


class ElementRef(NamedTuple):
    parent: AnnData
    attrname: str
    keys: tuple[str, ...] = ()

    def __str__(self) -> str:
        return f".{self.attrname}" + "".join(f"[{x!r}]" for x in self.keys)

    @property
    def _parent_el(self):
        return reduce(
            lambda d, k: d[k], self.keys[:-1], getattr(self.parent, self.attrname)
        )

    def get(self):
        """Get referenced value in self.parent."""
        return reduce(lambda d, k: d[k], self.keys, getattr(self.parent, self.attrname))

    def set(self, val):
        """Set referenced value in self.parent."""
        self._parent_el[self.keys[-1]] = val

    def delete(self):
        del self._parent_el[self.keys[-1]]
