from typing import TYPE_CHECKING, MutableMapping, Optional
from collections import OrderedDict
from copy import copy

from anndata.compat import _overloaded_uns, OverloadedDict, _slice_uns_sparse_matrices
from typing_extensions import Literal
import pandas as pd
from anndata._core.views import DataFrameView, DictView

if TYPE_CHECKING:
    from anndata._core.anndata import AnnData


class DataFrameDescriptor:
    """\
    One-dimensional annotation of observations or features.

    Parameters
    ----------
    attr:
        Named of the attribute.
    """

    def __init__(self, attr: Literal["obs", "var"]):
        self._attr = attr

    def __get__(self, adata: "AnnData", objtype: Optional[type] = None) -> pd.DataFrame:
        if not adata.is_view:
            return adata._obs if self.attr == "obs" else adata._var

        adata_ref = adata._adata_ref
        if self.attr == "obs":
            container_full = adata_ref.obs
            container = container_full.iloc[adata._oidx]
        elif self.attr == "var":
            container_full = adata_ref.var
            container = adata_ref.var.iloc[adata._vidx]
        else:
            raise NotImplementedError(self.attr)

        # we don't care about uns at this point
        adata._remove_unused_categories(container_full, container, uns={})
        return DataFrameView(container, view_args=(adata, self.attr))

    def __set__(self, adata: "AnnData", value: pd.DataFrame) -> None:
        adata._set_dim_df(value, self.attr)

    def __delete__(self, adata: "AnnData") -> None:
        setattr(
            adata, self.attr, pd.DataFrame(index=getattr(adata, f"{self.attr}_names"))
        )

    @property
    def attr(self) -> Literal["obs", "var"]:
        return self._attr


class UnsDescriptor:
    """Unstructured annotation (ordered dictionary)."""

    def __get__(self, adata: "AnnData", objtype: Optional[type] = None):
        if not adata.is_view:
            return _overloaded_uns(adata, adata._uns)

        adata_ref = adata._adata_ref
        # special case for old neighbors, backwards compat. Remove in anndata 0.8.
        # fix categories
        uns = _slice_uns_sparse_matrices(
            copy(adata_ref._uns), adata._oidx, adata_ref.n_obs
        )

        adata._remove_unused_categories(adata_ref.obs, adata.obs, uns, update_sub=False)
        adata._remove_unused_categories(adata_ref.var, adata.var, uns, update_sub=False)
        uns = DictView(uns, view_args=(adata, "_uns"))

        return _overloaded_uns(adata, uns)

    def __set__(self, adata: "AnnData", value: MutableMapping) -> None:
        if not isinstance(value, MutableMapping):
            raise ValueError(
                "Only mutable mapping types (e.g. dict) are allowed for `.uns`."
            )
        if isinstance(value, (OverloadedDict, DictView)):
            value = value.copy()
        if adata.is_view:
            adata._init_as_actual(adata.copy())
        adata._uns = value

    def __delete__(self, adata: "AnnData") -> None:
        adata._uns = OrderedDict()
