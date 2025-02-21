from __future__ import annotations

import inspect

# from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Generic, TypeVar
from warnings import warn

if TYPE_CHECKING:
    from collections.abc import Callable
from anndata import AnnData

# Based off of the extension framework in Polars
# https://github.com/pola-rs/polars/blob/main/py-polars/polars/api.py

__all__ = ["register_anndata_namespace"]


def find_stacklevel() -> int:
    """
    Find the first place in the stack that is not inside AnnData.

    Taken from:
    https://github.com/pandas-dev/pandas/blob/ab89c53f48df67709a533b6a95ce3d911871a0a8/pandas/util/_exceptions.py#L30-L51
    """
    import anndata as ad

    pkg_dir = str(Path(ad.__file__).parent)

    # https://stackoverflow.com/questions/17407119/python-inspect-stack-is-slow
    frame = inspect.currentframe()
    n = 0
    try:
        while frame:
            fname = inspect.getfile(frame)
            if fname.startswith(pkg_dir) or (
                (qualname := getattr(frame.f_code, "co_qualname", None))
                # ignore @singledispatch wrappers
                and qualname.startswith("singledispatch.")
            ):
                frame = frame.f_back
                n += 1
            else:
                break
    finally:
        # https://docs.python.org/3/library/inspect.html
        # > Though the cycle detector will catch these, destruction of the frames
        # > (and local variables) can be made deterministic by removing the cycle
        # > in a finally clause.
        del frame
    return n


NS = TypeVar("NS")


# Currently, there are no reserved namespaces internally, but if there ever are,
# we can use this to not allow them to be overridden.
_reserved_namespaces: set[str] = set.union(*(cls._accessors for cls in (AnnData,)))


class AccessorNameSpace(Generic[NS]):
    """Establish property-like namespace object for user-defined functionality."""

    def __init__(self, name: str, namespace: type[NS]) -> None:
        self._accessor = name
        self._ns = namespace

    def __get__(self, instance: NS | None, cls: type[NS]) -> NS | type[NS]:
        if instance is None:
            return self._ns

        ns_instance = self._ns(instance)  # type: ignore[call-arg]
        setattr(instance, self._accessor, ns_instance)
        return ns_instance


def _create_namespace(name: str, cls: type[AnnData]) -> Callable[[type[NS]], type[NS]]:
    """Register custom namespace against the underlying AnnData class."""

    def namespace(ns_class: type[NS]) -> type[NS]:
        if name in _reserved_namespaces:
            msg = f"cannot override reserved namespace {name!r}"
            raise AttributeError(msg)

        elif hasattr(cls, name):
            warn(
                f"Overriding existing custom namespace {name!r} (on {cls.__name__!r})",
                UserWarning,
                stacklevel=find_stacklevel(),
            )

        setattr(cls, name, AccessorNameSpace(name, ns_class))
        cls._accessors.add(name)
        return ns_class

    return namespace


def register_anndata_namespace(name: str) -> Callable[[type[NS]], type[NS]]:
    """Decorator for registering custom functionality with an :class:`~anndata.AnnData` object.

    Parameters
    ----------
    name
        Name under which the accessor should be registered. A warning is issued
        if this name conflicts with a preexisting attribute.

    Examples
    --------
    >>> import anndata as ad
    >>> from scipy.sparse import csr_matrix
    >>> import numpy as np
    >>>
    >>>
    >>> @ad.extensions.register_anndata_namespace("transforms")
    ... class TransformX:
    ...     def __init__(self, adata: ad.AnnData):
    ...         self._adata = adata
    >>>     def arcsinh_cofactor(
    ...         self, shift: float, scale: float, layer: str, inplace: bool = False
    ...     ) -> ad.AnnData:
    ...         self._adata.layers[layer] = np.arcsinh(self._adata.X.toarray() / scale) + shift
    ...         return None if inplace else self._adata
    >>>
    >>> rng = default_rng(42)
    >>> adata = ad.AnnData(
    ...     X=csr_matrix(rng.poisson(1, size=(100, 2000)), dtype=np.float32),
    ... )
    >>> adata.transforms.arcsinh_cofactor(1, 1, "arcsinh", inplace=True)
    >>> adata
    AnnData object with n_obs × n_vars = 100 × 2000
         layers: 'arcsinh'
    """
    return _create_namespace(name, AnnData)
