from __future__ import annotations

from typing import TYPE_CHECKING

from ..._core.xarray import Dataset2D

if TYPE_CHECKING:
    from anndata import AnnData


def has_dataset_2d(adata: AnnData) -> bool:
    if any(isinstance(annot_df, Dataset2D) for annot_df in [adata.obs, adata.var]):
        return True
    for annot_m_key in ["varm", "obsm"]:
        annot_m = getattr(adata, annot_m_key)
        if any(isinstance(maybe_df, Dataset2D) for maybe_df in annot_m.values()):
            return True
    return False
