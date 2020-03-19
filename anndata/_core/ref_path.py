from enum import Enum
from typing import Iterable, Union, Sequence, Dict, Tuple, Optional, TYPE_CHECKING

import pandas as pd
import numpy as np
import scipy.sparse as ssp

from .index import get_x_vector
from .raw import Raw
from ..compat import Literal

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
    raw = ((("X", "var", "varm"), ...), dict(rX=("X",), rv=("var",), rvm=("varm",)))

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
        elif self is Attr.raw:
            sub_path = RefPath.parse(path)
            # layers or obs
            sub_path.attr.validate(sub_path.path)
            return

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
                assert False, f"Unhandled check {check!r} for {self}"

    @staticmethod
    def prefix(prefix: str) -> Tuple[str, Tuple[str, ...]]:
        if prefix in Attr.__members__:  # name
            return prefix, ()
        for attr in Attr:
            path_prefix = attr.short_codes.get(prefix)
            if path_prefix is not None:
                return attr.name, path_prefix
        raise ValueError(f"Unknown attr name or short code {prefix!r}.")


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
        from_str = isinstance(full_path, str)
        if from_str:
            # TODO: don’t split off gene names with slashes
            full_path = full_path.split("/")

        if not full_path:
            raise ValueError(f"No path specified.")
        attr_or_code, *path = full_path
        # path can be shorthand: "obs/Foo" and "o/Foo"
        attr, path_prefix = Attr.prefix(attr_or_code)
        if from_str and attr in {"obsp", "varp"}:
            if path[-1] not in "01":
                raise ValueError(f"Invalid last segment of {attr} path: {path[-1]!r}.")
            path[-1] = int(path[-1])
        return RefPath(attr, *path_prefix, *path)

    @property
    def dim(self) -> Literal["obs", "var"]:
        if self.attr in {Attr.obs, Attr.obsm, Attr.obsp}:
            return "obs"
        if self.attr in {Attr.var, Attr.varm, Attr.varp}:
            return "var"
        if self.attr is Attr.layers:
            idx_dim = self.path[1]
            return "obs" if idx_dim == "var" else "var"
        if self.attr is Attr.raw:
            return self.raw_subpath.dim
        assert False, f"Unimplemented attr {self.attr}"

    @property
    def raw_subpath(self):
        try:
            assert self.attr is Attr.raw
            raw_attr, *raw_path = self.path
            return RefPath.parse((raw_attr, *raw_path))
        except (AssertionError, ValueError):
            raise AttributeError(
                f"{self.__class__.__name__} with attr={self.attr.name} "
                "has no `raw_subpath`"
            )

    def __repr__(self):
        return f"RefPath({self.attr.name!r}, {', '.join(map(repr, self.path))})"

    def __eq__(self, other: "RefPath"):
        if not isinstance(other, RefPath):
            try:
                return self == RefPath.parse(other)
            except ValueError:
                return False
        return self.attr == other.attr and self.path == other.path

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

    def get_vector(self, adata: Union["AnnData", Raw]):
        attr = getattr(adata, self.attr.name)
        if self.attr is Attr.layers:  # X is here after normalizing
            layer_name, dim, key = self.path
            layer_name = None if layer_name == "X" else layer_name
            return get_x_vector(adata, dim, key, layer_name)
        if self.attr is Attr.obs or self.attr is Attr.var:
            (col,) = self.path
            return _to_vector(attr[col])
        assert not isinstance(adata, Raw)
        if self.attr is Attr.obsm or self.attr is Attr.varm:
            m_name, col = self.path
            m = attr[m_name]
            return _to_vector(m[col] if isinstance(m, pd.DataFrame) else m[:, col])
        if self.attr is Attr.obsp or self.attr is Attr.varp:
            p_name, key, orient = self.path
            p = attr[p_name]
            idx = getattr(adata, f"{attr.dim}_names") == key
            return _to_vector(p[:, idx] if orient == 1 else p[idx, :])
        if self.attr is Attr.raw:
            return self.raw_subpath.get_vector(adata.raw)
        else:
            assert False, f"Unhandled attr {self.attr.name!r}"


def _to_vector(v: Union[ssp.spmatrix, np.ndarray, pd.Series]) -> np.ndarray:
    v = v.toarray() if ssp.issparse(v) else v
    v = v.values if isinstance(v, pd.Series) else v
    return np.ravel(v)


# AnnData methods


def resolve_path(
    adata: "AnnData",
    *path: Union[str, RefPath, int],
    alias_col: Optional[str],
    dim: Optional[Literal["obs", "var"]] = None,
    use_raw: bool = False,
    layer: Optional[str] = None,
) -> RefPath:
    """\
    """
    try:
        return RefPath.parse(*path)
    except ValueError:
        pass

    # TODO
    raise NotImplementedError()

    if alias_col is not None:
        idx = getattr(adata, dim)[alias_col]
        key = idx.index[idx == key]


def get_vector(
    adata: "AnnData",
    *path: Union[str, RefPath, int],
    alias_col: Optional[str],
    dim: Optional[Literal["obs", "var"]] = None,
    use_raw: bool = False,
    layer: Optional[str] = None,
) -> np.ndarray:
    """\
    For the syntax see :meth:`AnnData.resolve_path`
    """
    return resolve_path(**locals()).get_vector(adata)


def get_df(
    adata: "AnnData",
    paths: Iterable[Union[str, Tuple[Union[str, int], ...], RefPath]] = (),
    *,
    alias_col: Optional[str] = None,
    dim: Optional[Literal["obs", "var"]] = None,
    use_raw: bool = False,
    layer: Optional[str] = None,
) -> pd.DataFrame:
    """\
    """
    kwargs = locals()
    del kwargs["paths"], kwargs["adata"]
    # TODO: split paths in sequences before resolving: ("X", "var", ["GeneA", "GeneB"])
    paths = [resolve_path(adata, p, **kwargs) for p in paths]
    names = paths_to_names(paths)
    columns = {n: p.get_vector(adata) for n, p in names.items()}
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
