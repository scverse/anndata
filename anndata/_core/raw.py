from typing import Union, Mapping, Sequence, Tuple

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse import issparse

from . import anndata
from .index import _normalize_index, _subset, unpack_index
from .aligned_mapping import AxisArrays, AxisArraysView
from .sparse_dataset import SparseDataset


# TODO: Implement views for Raw
class Raw:
    def __init__(
        self,
        adata: "anndata.AnnData",
        X: Union[np.ndarray, sparse.spmatrix, None] = None,
        var: Union[pd.DataFrame, Mapping[str, Sequence], None] = None,
        varm: Union[AxisArrays, Mapping[str, np.ndarray], None] = None,
    ):
        from .anndata import _gen_dataframe

        self._adata = adata
        self._n_obs = adata.n_obs
        # construct manually
        if adata.isbacked == (X is None):
            self._X = X
            self._var = _gen_dataframe(var, self.X.shape[1], ["var_names"])
            self._varm = AxisArrays(self, 1, varm)
        elif X is None:  # construct from adata
            self._X = adata.X.copy()
            self._var = adata.var.copy()
            self._varm = AxisArrays(self, 1, adata.varm.copy())
        elif adata.isbacked:
            raise ValueError("Cannot specify X if adata is backed")

    @property
    def X(self) -> Union[SparseDataset, np.ndarray, sparse.spmatrix]:
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
            raise AttributeError(
                f"Could not find dataset for raw X in file: "
                f"{self._adata.file.filename}."
            )
        if isinstance(X, h5py.Group):
            X = SparseDataset(X)
        # Check if we need to subset
        if self._adata.is_view:
            # TODO: As noted above, implement views of raw
            #       so we can know if we need to subset by var
            return X[self._adata._oidx, slice(None)]
        else:
            return X

    @property
    def shape(self):
        return self.n_obs, self.n_vars

    @property
    def var(self):
        return self._var

    @property
    def n_vars(self):
        return self._var.shape[0]

    @property
    def n_obs(self):
        return self._n_obs

    @property
    def varm(self):
        return self._varm

    @property
    def var_names(self):
        return self.var.index

    @property
    def obs_names(self):
        return self._adata.obs_names

    def __getitem__(self, index):
        oidx, vidx = self._normalize_indices(index)

        # To preserve two dimensional shape
        if isinstance(vidx, (int, np.integer)):
            vidx = slice(vidx, vidx + 1, 1)
        if isinstance(oidx, (int, np.integer)):
            oidx = slice(oidx, oidx + 1, 1)

        if not self._adata.isbacked:
            X = _subset(self.X, (oidx, vidx))
        else:
            X = None

        var = self._var.iloc[vidx]
        new = Raw(self._adata, X=X, var=var)
        if self._varm is not None:
            # Since there is no view of raws
            new._varm = self._varm._view(_RawViewHack(self, vidx), (vidx,)).copy()
        return new

    def __str__(self):
        descr = f"Raw AnnData with n_obs × n_vars = {self.n_obs} × {self.n_vars}"
        for attr in ["var", "varm"]:
            keys = getattr(self, attr).keys()
            if len(keys) > 0:
                descr += f"\n    {attr}: {str(list(keys))[1:-1]}"
        return descr

    def copy(self):
        return Raw(
            self._adata,
            X=self._X.copy(),
            var=self._var.copy(),
            varm=None if self._varm is None else self._varm.copy(),
        )

    def to_adata(self):
        """Create full AnnData object."""
        return anndata.AnnData(
            X=self._X.copy(),
            var=self._var.copy(),
            varm=None if self._varm is None else self._varm.copy(),
            obs=self._adata.obs.copy(),
            obsm=self._adata.obsm.copy(),
            uns=self._adata.uns.copy(),
        )

    def _normalize_indices(self, packed_index):
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
        if k in self.var:
            return self.var[k].values
        else:
            idx = self._normalize_indices((k, slice(None)))
            a = self.X[idx]
        if issparse(a):
            a = a.toarray()
        return np.ravel(a)

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
    def __init__(self, raw: Raw, vidx: Union[slice, np.ndarray]):
        self.parent_raw = raw
        self.vidx = vidx

    @property
    def shape(self) -> Tuple[int, int]:
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
