from __future__ import annotations

from typing import TYPE_CHECKING, cast, overload

import numpy as np
import pandas as pd

from ..compat import CupyArray, CupySparseMatrix
from ..types import SupportsArrayApiBase
from ..utils import asarray
from .aligned_df import _gen_dataframe
from .aligned_mapping import (
    AlignedMappingProperty,
    AxisArrays,
    _copy_array,
    _on_disk_x,
)
from .index import (
    _as_numpy_idx,
    _get_vector_ambiguous,
    _normalize_index,
    _subset,
    unpack_index,
)

if TYPE_CHECKING:
    from collections.abc import Mapping, MutableMapping, Sequence
    from typing import ClassVar

    from ..abc import CSCDataset, CSRDataset
    from ..acc import AdRef
    from ..typing import Index, InMemoryArray, _Index1DNorm
    from .aligned_mapping import AxisArraysView, Value
    from .anndata import AnnData
    from .xarray import Dataset2D


# TODO: Implement views for Raw
class Raw:
    # `Raw` parents its `varm` like a non-view `AnnData` does, and never is a view
    is_view: ClassVar = False
    _adata_ref: ClassVar[None] = None
    _oidx: ClassVar[None] = None
    _vidx: ClassVar[None] = None

    _varm: MutableMapping[str, Value]
    """Backing store for the `varm` `AlignedMappingProperty`."""

    _X: InMemoryArray | None
    """In-memory `.X`; `None` if the parent `AnnData` is backed."""

    def __init__(
        self,
        adata: AnnData,
        X: InMemoryArray | None = None,
        var: pd.DataFrame | Dataset2D | Mapping[str, Sequence] | None = None,
        varm: Mapping[str, Value] | None = None,
    ) -> None:
        if X is not None and X.shape[0] != adata.n_obs:
            msg = f"X has {X.shape[0]} rows, but n_obs is {adata.n_obs}"
            raise ValueError(msg)

        self._adata = adata
        self._n_obs = adata.n_obs
        # construct manually
        if adata.isbacked == (X is None):
            # Move from GPU to CPU since it's large and not always used
            if isinstance(X, CupyArray | CupySparseMatrix):
                self._X = X.get()
            else:
                self._X = X
            n_var = None if self._X is None else self._X.shape[1]
            self._var = _gen_dataframe(
                var, ["var_names"], source="X", attr="var", length=n_var
            )
            self.varm = varm
        elif X is None:  # construct from adata
            # Move from GPU to CPU since it's large and not always used
            if adata.X is None:
                self._X = None
            elif isinstance(adata.X, CupyArray | CupySparseMatrix):
                self._X = adata.X.get()
            else:
                self._X = _copy_array(adata.X)
            if not isinstance(var_df := adata.var, pd.DataFrame):
                msg = "Cannot create `.raw` from an AnnData with a lazy `var`"
                raise NotImplementedError(msg)
            self._var = var_df.copy()
            self.varm = adata.varm.copy()
        elif adata.isbacked:
            msg = "Cannot specify X if adata is backed"
            raise ValueError(msg)

    def _get_X(self, layer=None):
        if layer is not None:
            raise ValueError()
        return self.X

    @property
    def X(self) -> InMemoryArray | CSRDataset | CSCDataset | None:
        # TODO: Handle unsorted array of integer indices for h5py.Datasets
        if not self._adata.isbacked:
            return self._X
        if not self._adata.file.is_open:
            self._adata.file.open()
        # Handle legacy file formats:
        if "raw/X" in self._adata.file:
            X = _on_disk_x(self._adata.file["raw/X"])
        elif "raw.X" in self._adata.file:
            X = _on_disk_x(self._adata.file["raw.X"])  # Backwards compat
        else:
            msg = (
                f"Could not find dataset for raw X in file: "
                f"{self._adata.file.filename}."
            )
            raise AttributeError(msg)
        # Check if we need to subset
        if (oidx := self._adata._oidx) is not None:  # i.e. `self._adata.is_view`
            # TODO: As noted above, implement views of raw
            #       so we can know if we need to subset by var
            return _subset(X, (oidx, slice(None)))
        return X

    @property
    def shape(self) -> tuple[int, int]:
        return self.n_obs, self.n_vars

    @property
    def var(self) -> pd.DataFrame:
        return self._var

    @property
    def n_vars(self) -> int:
        return self._var.shape[0]

    @property
    def n_obs(self) -> int:
        return self._n_obs

    varm: AlignedMappingProperty[AxisArrays, AxisArraysView, str] = (
        AlignedMappingProperty(AxisArrays, 1)
    )

    @property
    def var_names(self) -> pd.Index[str]:
        return self.var.index

    @property
    def obs_names(self) -> pd.Index[str]:
        return self._adata.obs_names

    @overload
    def __getitem__(self, index: AdRef) -> InMemoryArray: ...
    @overload
    def __getitem__(self, index: Index | tuple[AnnData, Index]) -> Raw: ...
    def __getitem__(
        self, index: Index | tuple[AnnData, Index] | AdRef
    ) -> Raw | InMemoryArray:
        from ..acc import AdRef
        from .anndata import AnnData

        if isinstance(index, AdRef):
            # `RefAcc.get` has no official `Raw` support, but works for the accessors it has
            return index.acc.get(self, index.idx)

        if (
            isinstance(index, tuple)
            and len(index) == 2
            and isinstance(parent := index[0], AnnData)
        ):
            adata = parent
            oidx, vidx = self._normalize_indices(index[1])
        else:
            oidx, vidx = self._normalize_indices(index)
            adata = self._adata[oidx]

        # To preserve two dimensional shape
        if isinstance(vidx, int | np.integer):
            vidx = slice(vidx, vidx + 1, 1)
        if isinstance(oidx, int | np.integer):
            oidx = slice(oidx, oidx + 1, 1)

        X = (
            None
            if self._adata.isbacked or self._X is None
            else _subset(self._X, (oidx, vidx))
        )

        var = self._var.iloc[_as_numpy_idx(vidx)]
        new = Raw(adata, X=X, var=var)
        if self.varm is not None:
            # Since there is no view of raws
            hack = cast("AnnData", _RawViewHack(self, vidx))
            new.varm = self.varm._view(hack, (vidx,)).copy()
        return new

    def __str__(self) -> str:
        descr = f"Raw AnnData with n_obs × n_vars = {self.n_obs} × {self.n_vars}"
        for attr in ["var", "varm"]:
            keys = getattr(self, attr).keys()
            if len(keys) > 0:
                descr += f"\n    {attr}: {str(list(keys))[1:-1]}"
        return descr

    def copy(self) -> Raw:
        return Raw(
            self._adata,
            # a backed `.X` lives in the file and cannot be passed to `Raw`
            X=None if self._adata.isbacked or self._X is None else _copy_array(self._X),
            var=self.var.copy(),
            varm=dict(self._varm),
        )

    def to_adata(self) -> AnnData:
        """Create full AnnData object."""
        from anndata import AnnData

        return AnnData(
            X=None if (X := self.X) is None else _copy_array(X),
            var=self.var.copy(),
            varm=dict(self._varm),
            obs=self._adata.obs.copy(),
            obsm=self._adata.obsm.copy(),
            obsp=self._adata.obsp.copy(),
            uns=dict(self._adata.uns),
        )

    def _normalize_indices(
        self, packed_index: Index
    ) -> tuple[_Index1DNorm | int | np.integer, _Index1DNorm | int | np.integer]:
        # deal with slicing with pd.Series
        if isinstance(packed_index, pd.Series):
            packed_index = packed_index.values
        if isinstance(packed_index, tuple):
            if len(packed_index) != 2:
                raise IndexDimError(len(packed_index))
            if isinstance(packed_index[1], pd.Series):
                packed_index = packed_index[0], packed_index[1].values
            if isinstance(packed_index[0], pd.Series):
                packed_index = packed_index[0].values, packed_index[1]
        obs, var = unpack_index(packed_index)
        obs = _normalize_index(obs, self._adata.obs_names)
        var = _normalize_index(var, self.var_names)
        return obs, var

    def var_vector(self, k: str, /) -> InMemoryArray:
        # TODO decorator to copy AnnData.var_vector docstring
        return _get_vector_ambiguous(self, k, "var")

    def obs_vector(self, k: str, /) -> InMemoryArray:
        # TODO decorator to copy AnnData.obs_vector docstring
        if (X := self.X) is None:
            msg = f"Cannot get vector {k!r} from a `Raw` without `X`."
            raise ValueError(msg)
        # non-unique names resolve to a slice or a mask rather than a position
        idx = self.var_names.get_loc(k)
        # a scalar position would drop the dimension we want to ravel
        a = _subset(X, (slice(None), np.array([idx]) if isinstance(idx, int) else idx))
        if isinstance(a, np.ndarray):
            return np.ravel(a)
        if isinstance(a, SupportsArrayApiBase):
            # the array API standard has no `ravel`, and we stay in the namespace
            return a.__array_namespace__().reshape(a, (a.size,))
        return np.ravel(asarray(a))


# This exists to accommodate AlignedMappings,
# until we implement a proper RawView or get rid of Raw in favor of modes.
class _RawViewHack:
    def __init__(self, raw: Raw, vidx: _Index1DNorm):
        self.parent_raw = raw
        self.vidx = vidx

    @property
    def shape(self) -> tuple[int, int]:
        return self.parent_raw.n_obs, len(self.var_names)

    @property
    def obs_names(self) -> pd.Index:
        return self.parent_raw.obs_names

    @property
    def var_names(self) -> pd.Index:
        return self.parent_raw.var_names[_as_numpy_idx(self.vidx)]


class IndexDimError(IndexError):
    MSG = (
        "You tried to slice an AnnData(View) object with an"
        "{}-dimensional index, but only 2 dimensions exist in such an object."
    )
    MSG_1D = (
        "\nIf you tried to slice cells using adata[cells, ], "
        "note that Python (unlike R) uses adata[cells, :] as slicing syntax."
    )

    def __init__(self, n_dims: int):
        msg = self.MSG.format(n_dims)
        if n_dims == 1:
            msg += self.MSG_1D
        super().__init__(msg)
