from enum import Enum
from itertools import product
from typing import Union, Optional  # Special
from typing import Iterable, Sequence, Generator  # ABCs
from typing import Dict, Tuple  # Classes

import pandas as pd
import numpy as np
import scipy.sparse as ssp

from ..compat import Literal
from . import anndata
from .index import get_x_vector
from .raw import Raw


# TODO: allow sequences: ("obsm", "X_pca", [0, 1])
from ..utils import _doc_params


class AttrInfo(Enum):
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
        if self is not AttrInfo.raw and not len(path) == len(self.validation):
            raise ValueError(
                f"Path length mismatch: required ({len(self.validation)}) ≠ "
                f"got ({len(path)}) in path {path!r} for attr {self.name}."
            )
        elif self is AttrInfo.raw:
            sub_path = RefPath.parse(path)
            # layers or obs
            sub_path._attr.validate(sub_path.path)
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
        if prefix in AttrInfo.__members__:  # name
            return prefix, ()
        for attr in AttrInfo:
            path_prefix = attr.short_codes.get(prefix)
            if path_prefix is not None:
                return attr.name, path_prefix
        raise ValueError(f"Unknown attr name or short code {prefix!r}.")


SHORTCUTS = ", ".join(
    f"`{attr.name}` ({', '.join(f'`{c}`' for c in attr.short_codes)})"
    for attr in AttrInfo
)


class RefPath:
    """\
    A fully resolved path referring to a vector in an :attr:`~anndata.AnnData` object.
    The vector will have the length of the `AnnData`’s :attr:`dim`,
    :attr:`~anndata.AnnData.n_obs` or :attr:`~anndata.AnnData.n_vars`.

    Depending on :attr:`attr`, :attr:`path` can be:

    layers
        `(name: str, dim: "obs"|"var", <dim>_name: str)`

        e.g. `RefPath("X", "var", "Actb")`
    obs, var
        `(column: str,)`
    obsm, varm
        `(name: str, column: int|str)`
    obsp, varp
        `(name: str, <dim>_name: str, axis: 0|1)`
    raw
        `(attr: "X"|"var"|"varm", ...)`

        :attr:`path` is a subpath into :attr:`~anndata.AnnData.raw`.
        See also :attr:`raw_subpath`.
    """

    _attr: AttrInfo
    _path: Tuple[Union[str, int], ...]

    def __init__(self, attr: str, *path: Union[str, int]):
        self._attr = AttrInfo[attr]
        self._path = path
        self._attr.validate(path)

    @staticmethod
    @_doc_params(shortcuts=SHORTCUTS)
    def parse(full_path: Union["RefPath", str, Sequence[Union[str, int]]]) -> "RefPath":
        """\
        Converts tuples or strings with path specifications to `RefPath`\\ s.

        Parameters
        ----------
        full_path
            Attribute, followed by path (see :class:`~anndata.RefPath`).
            Either a string containing a `'/'`-delimited path, or a tuple.
            Both can contain shortcuts for the attribute:

            {shortcuts}
        """
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
        attr, path_prefix = AttrInfo.prefix(attr_or_code)
        if from_str and attr in {"obsp", "varp"}:
            if path[-1] not in "01":
                raise ValueError(f"Invalid last segment of {attr} path: {path[-1]!r}.")
            path[-1] = int(path[-1])
        return RefPath(attr, *path_prefix, *path)

    @property
    def attr(self) -> str:
        """Name of the referred to attribute in an :class:`~anndata.AnnData` object."""
        return self._attr.name

    @property
    def path(self) -> Tuple[Union[str, int], ...]:
        """Path to a vector in :attr:`attr`. See :class:`~anndata.RefPath`."""
        return self._path

    @property
    def dim(self) -> Literal["obs", "var"]:
        """\
        Dimension this path refers to (`obs` or `var`).

        Returns e.g. `var` for an :attr:`attr` of `var`, `varm`, `varp`,
        and `obs` for `RefPath("layers", "var", "GeneName")`,
        as a vector for a certain `var` has :attr:`~anndata.AnnData.n_obs` entries.

        .. note::
           Returns `var` for `RefPath("raw", "var", "ColName")`,
           but note that :attr:`~anndata.AnnData.raw` can have a higher
           :attr:`~anndata.AnnData.n_vars` than its parent :class:`~anndata.AnnData`.
        """
        if self._attr in {AttrInfo.obs, AttrInfo.obsm, AttrInfo.obsp}:
            return "obs"
        if self._attr in {AttrInfo.var, AttrInfo.varm, AttrInfo.varp}:
            return "var"
        if self._attr is AttrInfo.layers:
            idx_dim = self.path[1]
            return "obs" if idx_dim == "var" else "var"
        if self._attr is AttrInfo.raw:
            return self.raw_subpath.dim
        assert False, f"Unimplemented attr {self._attr}"

    @property
    def raw_subpath(self) -> "RefPath":
        """Returns a sub-path that resolves within a :attr:`~anndata.AnnData.raw`."""
        try:
            assert self._attr is AttrInfo.raw
            raw_attr, *raw_path = self.path
            return RefPath.parse((raw_attr, *raw_path))
        except (AssertionError, ValueError):
            raise AttributeError(
                f"{self.__class__.__name__} with attr={self.attr} "
                "has no `raw_subpath`"
            )

    def __repr__(self):
        return f"RefPath({self.attr!r}, {', '.join(map(repr, self.path))})"

    def __eq__(self, other: "RefPath"):
        if not isinstance(other, RefPath):
            try:
                return self == RefPath.parse(other)
            except ValueError:
                return False
        return self._attr == other._attr and self.path == other.path

    # TODO: make public?
    def _make_name(self, length: int = 1) -> str:
        path = (self.attr, *self.path)
        if length == 1 and self._attr in {AttrInfo.obsp, AttrInfo.varp}:
            return path[-2]  # just key
        if length in {1, 2} and self.attr in {AttrInfo.obsm, AttrInfo.varm}:
            if isinstance(path[-1], int):
                return f"{path[-2]}{path[-1] + 1}"  # X_pca1
            else:  # normal
                return path[-1] if length == 1 else f"{path[-2]}-{path[-1]}"
        if length <= 2:
            return "-".join(path[-length:])
        return f"{path[:-2]}{self._make_name(length=2)}"

    def get_vector(self, adata: Union["anndata.AnnData", Raw]):
        """Returns the referred-to vector from an :class:`~anndata.AnnData` object."""
        attr = getattr(adata, self.attr)
        if self._attr is AttrInfo.layers:  # X is here after normalizing
            layer_name, dim, key = self.path
            layer_name = None if layer_name == "X" else layer_name
            return get_x_vector(adata, dim, key, layer_name)
        if self._attr is AttrInfo.obs or self._attr is AttrInfo.var:
            (col,) = self.path
            return _to_vector(attr[col])
        assert not isinstance(adata, Raw)
        if self._attr is AttrInfo.obsm or self._attr is AttrInfo.varm:
            m_name, col = self.path
            m = attr[m_name]
            return _to_vector(m[col] if isinstance(m, pd.DataFrame) else m[:, col])
        if self._attr is AttrInfo.obsp or self._attr is AttrInfo.varp:
            p_name, key, orient = self.path
            p = attr[p_name]
            idx = getattr(adata, f"{attr.dim}_names") == key
            return _to_vector(p[:, idx] if orient == 1 else p[idx, :])
        if self._attr is AttrInfo.raw:
            return self.raw_subpath.get_vector(adata.raw)
        else:
            assert False, f"Unhandled attr {self.attr!r}"


def _to_vector(v: Union[ssp.spmatrix, np.ndarray, pd.Series]) -> np.ndarray:
    v = v.toarray() if ssp.issparse(v) else v
    v = v.values if isinstance(v, pd.Series) else v
    return np.ravel(v)


# AnnData methods


RefPathLike = Union[str, Tuple[Union[str, int], ...], RefPath]
PARAMS_RESOLVE = """\
adata
    Annotated data object.
path
    This supports subpaths of the :class:`~anndata.RefPath` syntax.
    `str` keys or tuple subpaths (like `'GeneA'` or `('X_pca', 0)`) are resolved
    according to `dim`, `use_raw`, and `alias_col`.
    As `RefPath`\\ s are always unique, they get passed through.
dim
    Dimension to resolve paths in.
    If `dim=None`, both dimensions are tried and an error is thrown for duplicates.
    If e.g. `dim='obs'`, it would find an unique name in
    :attr:`~anndata.AnnData.obs_names` or :attr:`~anndata.AnnData.obs`\\ `.columns`.
use_raw
    Resolve partial paths for `X`, `var`, or `varm` in `adata.raw`
    instead of `adata`?
alias_col
    A column in `adata.<dim>` with gene names to use instead of `adata.<dim>_names`
    (autodetected if `dim=None`)
layer
    The layer to get the vector from if the path resolves to a `<dim>_name`.\
"""


@_doc_params(params_resolve=PARAMS_RESOLVE)
def resolve_path(
    adata: "anndata.AnnData",
    *path: Union[str, RefPath, int],
    dim: Optional[Literal["obs", "var"]] = None,
    use_raw: bool = False,
    alias_col: Optional[str],
    layer: Optional[str] = None,
) -> RefPath:
    """\
    Resolves a :class:`~anndata.RefPath`-like :class:`tuple` or :class:`str` key.

    Parameters
    ----------
    {params_resolve}
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


@_doc_params(params_resolve=PARAMS_RESOLVE)
def get_vector(
    adata: "anndata.AnnData",
    *path: Union[str, RefPath, int],
    dim: Optional[Literal["obs", "var"]] = None,
    use_raw: bool = False,
    alias_col: Optional[str],
    layer: Optional[str] = None,
) -> np.ndarray:
    """\
    Get a single 1D vector using the `path`.

    Parameters
    ----------
    {params_resolve}
    """
    return resolve_path(**locals()).get_vector(adata)


@_doc_params(params_resolve=PARAMS_RESOLVE)
def get_df(
    adata: "anndata.AnnData",
    paths: Iterable[RefPathLike],
    *,
    dim: Optional[Literal["obs", "var"]] = None,
    use_raw: bool = False,
    alias_col: Optional[str] = None,
    layer: Optional[str] = None,
) -> pd.DataFrame:
    """\
    Resolves multiple paths, gets vectors via :meth:`~anndata.AnnData.get_vector` and
    joins them to a :class:`pandas.DataFrame`.

    So becomes `("obs", ["A", "B"])` the paths `("obs", "A")` and `("obs", "B")`.
    The data frame column names are unique and as short as possible.

    Parameters
    ----------
    {params_resolve}
    """
    kwargs = locals()
    del kwargs["paths"], kwargs["adata"]
    paths = [resolve_path(adata, p, **kwargs) for p in split_paths(paths)]
    names = paths_to_names(paths)
    columns = {n: p.get_vector(adata) for n, p in names.items()}
    return pd.DataFrame(columns, adata.obs_names)


def split_paths(
    multipath: Union[RefPathLike, Iterable[RefPathLike]]
) -> Generator[RefPathLike, None, None]:
    if isinstance(multipath, RefPath):
        yield multipath  # validated, so no inner sequence!
    elif isinstance(multipath, str):
        # TODO: globs and stuff. probably needs resolving info
        yield multipath
    elif isinstance(multipath, tuple):
        yield from product(
            *([elem] if isinstance(elem, (str, int)) else elem for elem in multipath)
        )
    else:  # iterable
        for mp in multipath:
            yield from split_paths(mp)


def paths_to_names(paths: Sequence[RefPath], length: int = 1) -> Dict[str, RefPath]:
    names = {}
    dupes = {}
    for path, name in zip(paths, (p._make_name(length) for p in paths)):
        dupes.setdefault(name, []).append(path)
    for name, paths_dup in dupes.items():
        if len(paths_dup) == 1:
            names[name] = paths_dup[0]
        elif any(len(p) > length for p in paths_dup):
            names.update(paths_to_names(paths_dup, length + 1))
        else:
            raise ValueError(f"Not sure how {name} can be extended for {paths_dup}")
    return names
