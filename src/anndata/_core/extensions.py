from __future__ import annotations

from typing import TYPE_CHECKING

from scverse_misc import make_register_namespace_decorator

from ..types import ExtensionNamespace
from .anndata import AnnData

if TYPE_CHECKING:
    from collections.abc import Callable

__all__ = ["register_anndata_namespace"]


def register_anndata_namespace[NameSpT: ExtensionNamespace](
    name: str,
) -> Callable[[type[NameSpT]], type[NameSpT]]:
    """Decorator for registering custom functionality with an :class:`~anndata.AnnData` object.

    This decorator allows you to extend AnnData objects with custom methods and properties
    organized under a namespace. The namespace becomes accessible as an attribute on AnnData
    instances, providing a clean way to you to add domain-specific functionality without modifying
    the AnnData class itself, or extending the class with additional methods as you see fit in your workflow.

    Parameters
    ----------
    name
        Name under which the accessor should be registered. This will be the attribute name
        used to access your namespace's functionality on AnnData objects (e.g., `adata.{name}`).
        Cannot conflict with existing AnnData attributes like `obs`, `var`, `X`, etc. The list of reserved
        attributes includes everything outputted by `dir(AnnData)`.

    Returns
    -------
    A decorator that registers the decorated class as a custom namespace.

    Notes
    -----
    Implementation requirements:

    1. The decorated class must have an `__init__` method that accepts exactly one parameter
       (besides `self`) named `adata` and annotated with type :class:`~anndata.AnnData`.
    2. The namespace will be initialized with the AnnData object on first access and then
       cached on the instance.
    3. If the namespace name conflicts with an existing namespace, a warning is issued.
    4. If the namespace name conflicts with a built-in AnnData attribute, an AttributeError is raised.

    Examples
    --------
    Simple transformation namespace with two methods:

    >>> import anndata as ad
    >>> import numpy as np
    >>>
    >>> @ad.register_anndata_namespace("transform")
    ... class TransformX:
    ...     def __init__(self, adata: ad.AnnData):
    ...         self._adata = adata
    ...
    ...     def log1p(
    ...         self, layer: str = None, inplace: bool = False
    ...     ) -> ad.AnnData | None:
    ...         '''Log1p transform the data.'''
    ...         data = self._adata.layers[layer] if layer else self._adata.X
    ...         log1p_data = np.log1p(data)
    ...
    ...         if layer:
    ...             layer_name = f"{layer}_log1p" if not inplace else layer
    ...         else:
    ...             layer_name = "log1p"
    ...
    ...         self._adata.layers[layer_name] = log1p_data
    ...
    ...         if not inplace:
    ...             return self._adata
    ...
    ...     def arcsinh(
    ...         self, layer: str = None, scale: float = 1.0, inplace: bool = False
    ...     ) -> ad.AnnData | None:
    ...         '''Arcsinh transform the data with optional scaling.'''
    ...         data = self._adata.layers[layer] if layer else self._adata.X
    ...         asinh_data = np.arcsinh(data / scale)
    ...
    ...         if layer:
    ...             layer_name = f"{layer}_arcsinh" if not inplace else layer
    ...         else:
    ...             layer_name = "arcsinh"
    ...
    ...         self._adata.layers[layer_name] = asinh_data
    ...
    ...         if not inplace:
    ...             return self._adata
    >>>
    >>> # Create an AnnData object
    >>> rng = np.random.default_rng(42)
    >>> adata = ad.AnnData(X=rng.poisson(1, size=(100, 2000)))
    >>>
    >>> # Use the registered namespace
    >>> adata.transform.log1p()  # Transforms X and returns the AnnData object
    AnnData object with n_obs × n_vars = 100 × 2000
        layers: 'log1p'
    >>> adata.transform.arcsinh()  # Transforms X and returns the AnnData object
    AnnData object with n_obs × n_vars = 100 × 2000
        layers: 'log1p', 'arcsinh'
    """
    return make_register_namespace_decorator(AnnData, "adata", name, "numpy")
