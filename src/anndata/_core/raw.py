from __future__ import annotations

from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import issparse

from ..compat import CupyArray, CupySparseMatrix
from .aligned_df import _gen_dataframe
from .aligned_mapping import AlignedMappingProperty, AxisArrays
from .index import _normalize_index, _subset, get_vector, unpack_index
from .sparse_dataset import sparse_dataset

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence
    from typing import ClassVar

    from ..compat import CSMatrix, Index, Index1DNorm
    from .aligned_mapping import AxisArraysView
    from .anndata import AnnData
    from .sparse_dataset import BaseCompressedSparseDataset


# TODO: Implement views for Raw
class Raw:
    is_view: ClassVar = False

    def __init__(
        self,
        adata: AnnData,
        X: np.ndarray | CSMatrix | None = None,
        var: pd.DataFrame | Mapping[str, Sequence] | None = None,
        varm: AxisArrays | Mapping[str, np.ndarray] | None = None,
    ):
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
            if isinstance(adata.X, CupyArray | CupySparseMatrix):
                self._X = adata.X.get()
            else:
                self._X = adata.X.copy()
            self._var = adata.var.copy()
            self.varm = adata.varm.copy()
        elif adata.isbacked:
            msg = "Cannot specify X if adata is backed"
            raise ValueError(msg)

    def _get_X(self, layer=None):
        if layer is not None:
            raise ValueError()
        return self.X

    @property
    def X(self) -> BaseCompressedSparseDataset | np.ndarray | CSMatrix:
        # TODO: Handle unsorted array of integer indices for h5py.Datasets
        if not self._adata.isbacked:
            return self._X
        if not self._adata.file.is_open:
            self._adata.file.open()
        # Handle legacy file formats:
        if "raw/X" in self._adata.file:
            X = self._adata.file["raw/X"]
        elif "raw.X" in self._adata.file:
            X = self._adata.file["raw.X"]  # Backwards compat
        else:
            msg = (
                f"Could not find dataset for raw X in file: "
                f"{self._adata.file.filename}."
            )
            raise AttributeError(msg)
        if isinstance(X, h5py.Group):
            X = sparse_dataset(X)
        # Check if we need to subset
        if self._adata.is_view:
            # TODO: As noted above, implement views of raw
            #       so we can know if we need to subset by var
            return _subset(X, (self._adata._oidx, slice(None)))
        else:
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

    varm: AlignedMappingProperty[AxisArrays | AxisArraysView] = AlignedMappingProperty(
        "varm", AxisArrays, 1
    )

    @property
    def var_names(self) -> pd.Index[str]:
        return self.var.index

    @property
    def obs_names(self) -> pd.Index[str]:
        return self._adata.obs_names

    def __getitem__(self, index: Index) -> Raw:
        oidx, vidx = self._normalize_indices(index)

        # To preserve two dimensional shape
        if isinstance(vidx, int | np.integer):
            vidx = slice(vidx, vidx + 1, 1)
        if isinstance(oidx, int | np.integer):
            oidx = slice(oidx, oidx + 1, 1)

        X = _subset(self.X, (oidx, vidx)) if not self._adata.isbacked else None

        var = self._var.iloc[vidx]
        new = Raw(self._adata, X=X, var=var)
        if self.varm is not None:
            # Since there is no view of raws
            new.varm = self.varm._view(_RawViewHack(self, vidx), (vidx,)).copy()
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
            X=self.X.copy(),
            var=self.var.copy(),
            varm=None if self._varm is None else self._varm.copy(),
        )

    def to_adata(self) -> AnnData:
        """Create full AnnData object."""
        from anndata import AnnData

        return AnnData(
            X=self.X.copy(),
            var=self.var.copy(),
            varm=None if self._varm is None else self._varm.copy(),
            obs=self._adata.obs.copy(),
            obsm=self._adata.obsm.copy(),
            obsp=self._adata.obsp.copy(),
            uns=self._adata.uns.copy(),
        )

    def _normalize_indices(
        self, packed_index: Index
    ) -> tuple[Index1DNorm | int | np.integer, Index1DNorm | int | np.integer]:
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

    def var_vector(self, k: str) -> np.ndarray:
        # TODO decorator to copy AnnData.var_vector docstring
        return get_vector(self, k, "var", "obs")

    def obs_vector(self, k: str) -> np.ndarray:
        # TODO decorator to copy AnnData.obs_vector docstring
        idx = self._normalize_indices((slice(None), k))
        a = self.X[idx]
        if issparse(a):
            a = a.toarray()
        return np.ravel(a)


# This exists to accommodate AlignedMappings,
# until we implement a proper RawView or get rid of Raw in favor of modes.
class _RawViewHack:
    def __init__(self, raw: Raw, vidx: slice | np.ndarray):
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
        return self.parent_raw.var_names[self.vidx]


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
