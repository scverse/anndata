from functools import reduce
from typing import NamedTuple, Tuple

from . import anndata


class ElementRef(NamedTuple):
    parent: "anndata.AnnData"
    attrname: str
    keys: Tuple[str, ...] = ()

    def __str__(self) -> str:
        return f".{self.attrname}" + "".join(map(lambda x: f"['{x}']", self.keys))

    def get(self):
        """Get referenced value in self.parent."""
        return reduce(lambda d, k: d[k], self.keys, getattr(self.parent, self.attrname))

    def set(self, val):
        """Set referenced value in self.parent."""
        m = reduce(
            lambda d, k: d[k], self.keys[:-1], getattr(self.parent, self.attrname)
        )
        m[self.keys[-1]] = val

    def delete(self):
        m = reduce(
            lambda d, k: d[k], self.keys[:-1], getattr(self.parent, self.attrname)
        )
        del m[self.keys[-1]]
