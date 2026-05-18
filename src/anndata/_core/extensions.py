from __future__ import annotations

from scverse_misc import make_register_namespace_decorator

from .anndata import AnnData

__all__ = ["register_anndata_namespace"]


register_anndata_namespace = make_register_namespace_decorator(
    AnnData, "adata", "register_anndata_namespace", "numpy"
)
