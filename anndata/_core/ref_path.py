from enum import Enum
from typing import Iterable, Union, Sequence, Dict, Tuple, Optional, TYPE_CHECKING

import pandas as pd
import numpy as np

from .index import get_x_vector
from .raw import Raw

if TYPE_CHECKING:
    from .anndata import AnnData


# TODO: allow sequences: ("obsm", "X_pca", [0, 1])


class Attr(Enum):
    layers = ((str, ("obs", "var"), str), dict(l=(), X=("X",)))
    obs = ((str,), dict(o=()))
    obsm = ((str, (int, str)), dict(om=()))
    obsp = ((str, str, (0, 1)), dict(op=()))
    var = ((str,), dict(v=()))
    varm = ((str, (int, str)), dict(vm=()))
    varp = ((str, str, (0, 1)), dict(vp=()))
    raw = ((("X", "obs"), ...), dict(rX=("X",), ro=("obs",)))

    def __repr__(self):
        return f"{self.__class__.__name__}.{self.name}"

    @property
    def validation(self):
        return self.value[0]

    @property
    def short_codes(self):
        return self.value[1]

    def validate(self, path: Sequence[Union[str, int]]):
        if self is not Attr.raw and not len(path) == len(self.validation):
            raise ValueError(
                f"Path length mismatch: required ({len(self.validation)}) ≠ "
                f"got ({len(path)}) in path {path!r} for attr {self.name}."
            )
        for i, (elem, check) in enumerate(zip(path, self.validation)):
            err_prefix = f"Element path[{i}]={path[i]!r} in path {path!r}"
            if isinstance(check, tuple):
                if isinstance(check[0], type):
                    if not any(isinstance(elem, c) for c in check):
                        raise ValueError(
                            f"{err_prefix} is not of one of the types {check!r}"
                        )
                elif not any(elem == c for c in check):
                    raise ValueError(f"{err_prefix} is not of one of {check!r}")
            elif isinstance(check, type):
                if not isinstance(elem, check):
                    raise ValueError(f"{err_prefix} is non of type {check}")
            else:
                assert check is Ellipsis
                assert self is Attr.raw
                sub_path = RefPath.parse(path)
                # layers or obs
                sub_path.attr.validate(sub_path.path)

    @staticmethod
    def prefix(prefix: str) -> Tuple[str, Tuple[str, ...]]:
        if prefix in Attr.__members__:  # name
            return prefix, ()
        for attr in Attr:
            path_prefix = attr.short_codes.get(prefix)
            if path_prefix:
                return attr.name, path_prefix


class RefPath:
    attr: Attr
    path: Sequence[Union[str, int]]

    def __init__(self, attr: str, *path: Union[str, int]):
        self.attr = Attr[attr]
        self.path = path
        self.attr.validate(path)

    @staticmethod
    def parse(full_path: Union["RefPath", str, Sequence[Union[str, int]]]) -> "RefPath":
        if isinstance(full_path, RefPath):
            return full_path
        if isinstance(full_path, str):
            attr_or_code, *path = full_path.split("/")
        else:
            attr_or_code, *path = full_path
        # path can be shorthand: "obs/Foo" and "o/Foo"
        attr, path_prefix = Attr.prefix(attr_or_code)
        return RefPath(attr, *path_prefix, *path)

    def __repr__(self):
        return f"RefPath({self.attr.name!r}, {', '.join(map(repr, self.path))})"

    def make_name(self, length: int = 1) -> str:
        path = (self.attr.name, *self.path)
        if length == 1 and self.attr in {Attr.obsp, Attr.varp}:
            return path[-2]  # just key
        if length in {1, 2} and self.attr.name in {Attr.obsm, Attr.varm}:
            if isinstance(path[-1], int):
                return f"{path[-2]}{path[-1] + 1}"  # X_pca1
            else:  # normal
                return path[-1] if length == 1 else f"{path[-2]}-{path[-1]}"
        if length <= 2:
            return "-".join(path[-length:])
        return f"{path[:-2]}{self.make_name(length=2)}"

    def get_vector(self, adata: Union["AnnData", Raw], alias_col: Optional[str] = None):
        attr = getattr(adata, self.attr.name, None)
        if attr is None:
            raise ValueError(
                f"Unknown AnnData attribute {self.attr.name!r} (path={self.path!r})"
            )
        if attr is adata.layers:  # X is here after normalizing
            layer_name, dim, key = self.path
            layer_name = None if layer_name == "X" else layer_name
            if alias_col is not None:
                idx = getattr(adata, dim)[alias_col]
                key = idx.index[idx == key]
            return get_x_vector(adata, dim, key, layer_name)
        if attr is adata.obs or attr is adata.var:
            (col,) = self.path
            return attr[col]
        assert not isinstance(adata, Raw)
        if attr is adata.obsm or attr is adata.varm:
            m_name, col = self.path
            m = attr[m_name]
            return m[col] if isinstance(m, pd.DataFrame) else m[:, col]
        if attr is adata.obsp or attr is adata.varp:
            p_name, obs_name, orient = self.path
            p = attr[p_name]
            # TODO: take from index
            return np.ravel(p[:, obs_name] if orient == 1 else p[obs_name, :])
        if attr is adata.raw:
            raw_attr, *raw_path = self.path
            return RefPath(raw_attr, *raw_path).get_vector(adata.raw, alias_col)
        else:
            assert False, f"Unhandled attr {self.attr.name!r}"


def get_df(
    adata: "AnnData",
    paths: Iterable[Union[str, RefPath]] = (),
    *,
    gene_symbols: RefPath = None,
    # no use_raw, no layer.
) -> pd.DataFrame:
    paths = list(map(RefPath.parse, paths))
    names = paths_to_names(paths)
    columns = {n: p.get_vector(adata, gene_symbols) for n, p in names.items()}
    return pd.DataFrame(columns, adata.obs_names)


def paths_to_names(paths: Sequence[RefPath], length: int = 1) -> Dict[str, RefPath]:
    names = {}
    dupes = {}
    for path, name in zip(paths, (p.make_name(length) for p in paths)):
        dupes.setdefault(name, []).append(path)
    for name, paths_dup in dupes.items():
        if len(paths_dup) == 1:
            names[name] = paths_dup[0]
        elif any(len(p) > length for p in paths_dup):
            names.update(paths_to_names(paths_dup, length + 1))
        else:
            raise ValueError(f"Not sure how {name} can be extended for {paths_dup}")
    return names
