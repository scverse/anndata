from collections import OrderedDict, abc as cabc
from copy import copy
from functools import cached_property
from pathlib import Path
from typing import (
    Any,
    Iterable,
    Mapping,
    MutableMapping,
    Optional,
    Union,
    List,
    Sequence,
    Tuple,
)
from anndata._core.access import ElementRef
from anndata._core.aligned_mapping import AlignedMapping, Layers, PairwiseArrays
from anndata._core.anndata import StorageType, _check_2d_shape, _gen_dataframe
from anndata._core.anndata_base import AbstractAnnData
from anndata._core.file_backing import AnnDataFileManager
from anndata._core.index import Index, _normalize_indices, _subset
from anndata._core.raw import Raw
from anndata._core.sparse_dataset import sparse_dataset
from anndata._core.views import _resolve_idxs, as_view
from anndata._io.specs.registry import read_elem
from anndata.compat import _move_adj_mtx, _read_attr
from anndata.utils import convert_to_dict

import zarr
import pandas as pd
import numpy as np
from xarray.core.indexing import ExplicitlyIndexedNDArrayMixin

from ..._core import AnnData, AxisArrays
from .. import read_dispatched


class LazyCategoricalArray(ExplicitlyIndexedNDArrayMixin):
    __slots__ = ("codes", "attrs", "_categories", "_categories_cache")

    def __init__(self, group, *args, **kwargs):
        """Class for lazily reading categorical data from formatted zarr group

        Args:
            group (zarr.Group): group containing "codes" and "categories" key as well as "ordered" attr
        """
        self.codes = group["codes"]
        self._categories = group["categories"]
        self._categories_cache = None
        self.attrs = dict(group.attrs)

    @property
    def categories(self):  # __slots__ and cached_property are incompatible
        if self._categories_cache is None:
            self._categories_cache = self._categories[...]
        return self._categories_cache

    @property
    def dtype(self) -> pd.CategoricalDtype:
        return pd.CategoricalDtype(self.categories, self.ordered)

    @property
    def shape(self) -> Tuple[int, ...]:
        return self.codes.shape

    @property
    def ordered(self):
        return bool(self.attrs["ordered"])

    def __getitem__(self, selection) -> pd.Categorical:
        codes = self.codes.oindex[selection]
        if codes.shape == ():  # handle 0d case
            codes = np.array([codes])
        return pd.Categorical.from_codes(
            codes=codes,
            categories=self.categories,
            ordered=self.ordered,
        )

    def __repr__(self) -> str:
        return f"LazyCategoricalArray(codes=..., categories={self.categories}, ordered={self.ordered})"

    def __eq__(self, __o) -> np.ndarray:
        return self[()] == __o

    def __ne__(self, __o) -> np.ndarray:
        return ~(self == __o)


@_subset.register(LazyCategoricalArray)
def _subset_lazy_cat(a: LazyCategoricalArray, subset_idx: Index):
    return a[subset_idx]


@as_view.register(pd.Categorical)
def _subset_lazy_cat(a: pd.Categorical, view_args):
    return a


class AxisArraysRemote(AxisArrays):
    def __getattr__(self, __name: str):
        # If we a method has been accessed that is not here, try the pandas implementation
        if hasattr(pd.DataFrame, __name):
            return self.to_df().__getattribute__(__name)
        return object.__getattribute__(self, __name)

    def to_df(self) -> pd.DataFrame:
        """Convert to pandas dataframe."""
        df = pd.DataFrame(index=self.dim_names)
        for key in self.keys():
            if "index" not in key:
                df[key] = self[key][()]
        return df

    @property
    def iloc(self):
        class IlocDispatch:
            def __getitem__(self_iloc, idx):
                return self._view(self.parent, (idx,))

        return IlocDispatch()


class AnnDataRemote(AbstractAnnData):
    def __init__(
        self,
        X=None,
        obs=None,
        var=None,
        uns=None,
        obsm=None,
        varm=None,
        layers=None,
        raw=None,
        dtype=None,
        shape=None,
        file=None,
        *,
        obsp=None,
        varp=None,
        oidx=None,
        vidx=None,
    ):
        self._is_view = False
        if oidx is not None and vidx is not None:  # and or or?
            self._is_view = True  # hack needed for clean use of views below
        adata_ref = self
        # init from AnnData
        if issubclass(type(X), AbstractAnnData):
            if any((obs, var, uns, obsm, varm, obsp, varp)):
                raise ValueError(
                    "If `X` is a dict no further arguments must be provided."
                )
            if X.is_view:
                prev_oidx, prev_vidx = X._oidx, X._vidx
                self._oidx, self._vidx = _resolve_idxs(
                    (prev_oidx, prev_vidx), (oidx, vidx), X._X
                )
            else:
                self._oidx = oidx
                self._vidx = vidx
            if self._is_view:
                adata_ref = X  # seems to work
            if file is None:
                file = X.file
            X, obs, var, uns, obsm, varm, obsp, varp, layers, raw = (
                X._X,
                X.obs,
                X.var,
                X.uns,
                X.obsm,
                X.varm,
                X.obsp,
                X.varp,
                X.layers,
                X.raw,
            )
        else:
            self._oidx = oidx
            self._vidx = vidx

        if X is not None:
            for s_type in StorageType:
                if isinstance(X, s_type.value):
                    break
            else:
                class_names = ", ".join(c.__name__ for c in StorageType.classes())
                raise ValueError(
                    f"`X` needs to be of one of {class_names}, not {type(X)}."
                )
            if shape is not None:
                raise ValueError("`shape` needs to be `None` if `X` is not `None`.")
            _check_2d_shape(X)
            # if type doesn’t match, a copy is made, otherwise, use a view
            if dtype is not None:
                X = X.astype(dtype)
            # data matrix and shape
            self._X = X
        else:
            self._X = None

        # annotations - need names already for AxisArrays to work.
        self.file = file
        self.obs_names = obs[self.file["obs"].attrs["_index"]]
        if oidx is not None:
            self.obs_names = self.obs_names[oidx]
        self.var_names = var[self.file["var"].attrs["_index"]]
        if vidx is not None:
            self.var_names = self.var_names[vidx]
        self.obs = AxisArraysRemote(adata_ref, 0, vals=convert_to_dict(obs))
        self.var = AxisArraysRemote(adata_ref, 1, vals=convert_to_dict(var))

        self.obsm = AxisArraysRemote(adata_ref, 0, vals=convert_to_dict(obsm))
        self.varm = AxisArraysRemote(adata_ref, 1, vals=convert_to_dict(varm))

        self.obsp = PairwiseArrays(adata_ref, 0, vals=convert_to_dict(obsp))
        self.varp = PairwiseArrays(adata_ref, 1, vals=convert_to_dict(varp))

        self.layers = Layers(layers)
        if self.is_view:
            self.obs = self.obs._view(self, (oidx,))
            self.var = self.var._view(self, (vidx,))
            self.obsm = self.obsm._view(self, (oidx,))
            self.varm = self.varm._view(self, (vidx,))
            self.obsp = self.obsp._view(self, oidx)
            self.varp = self.varp._view(self, vidx)
            self.layers = self.layers._view(self, (oidx, vidx))
        self.uns = uns or OrderedDict()

        if not raw:
            self.raw = None
        elif isinstance(raw, cabc.Mapping):
            self.raw = Raw(self, **raw)
        else:  # is a Raw from another AnnData
            self.raw = Raw(self, raw._X, raw.var, raw.varm)

        self._run_checks()

    def _run_checks(self):
        assert len(self.obs_names) == self.shape[0]
        assert len(self.var_names) == self.shape[1]
        for layer in self.layers:
            assert len(self.obs_names) == layer.shape[0]
            assert len(self.var_names) == layer.shape[1]

    def __getitem__(self, index: Index) -> "AnnData":
        """Returns a sliced view of the object."""
        oidx, vidx = self._normalize_indices(index)
        return AnnDataRemote(self, oidx=oidx, vidx=vidx)

    def _normalize_indices(self, index: Optional[Index]) -> Tuple[slice, slice]:
        return _normalize_indices(index, self.obs_names, self.var_names)

    @property
    def obs_names(self) -> pd.Index:
        return self._obs_names

    @property
    def var_names(self) -> pd.Index:
        return self._var_names

    @obs_names.setter
    def obs_names(self, names: Sequence[str]):
        self._obs_names = names

    @var_names.setter
    def var_names(self, names: Sequence[str]):
        self._var_names = names

    @property
    def X(self):
        if hasattr(self, "_X"):
            if self.is_view:
                return _subset(self._X, (self._oidx, self._vidx))
            return self._X
        return None

    @X.setter
    def X(self, X):
        self._X = X

    @property
    def obs(self):
        if hasattr(self, "_obs"):
            return self._obs
        return None

    @obs.setter
    def obs(self, obs):
        self._obs = obs

    @property
    def obsm(self):
        if hasattr(self, "_obsm"):
            return self._obsm
        return None

    @obsm.setter
    def obsm(self, obsm):
        self._obsm = obsm

    @property
    def obsp(self):
        if hasattr(self, "_obsp"):
            return self._obsp
        return None

    @obsp.setter
    def obsp(self, obsp):
        self._obsp = obsp

    @property
    def var(self):
        if hasattr(self, "_var"):
            return self._var
        return None

    @var.setter
    def var(self, var):
        self._var = var

    @property
    def uns(self):
        if hasattr(self, "_uns"):
            return self._uns
        return None

    @uns.setter
    def uns(self, uns):
        self._uns = uns

    @property
    def varm(self):
        if hasattr(self, "_varm"):
            return self._varm
        return None

    @varm.setter
    def varm(self, varm):
        self._varm = varm

    @property
    def varp(self):
        if hasattr(self, "_varp"):
            return self._varp
        return None

    @varp.setter
    def varp(self, varp):
        self._varp = varp

    @property
    def raw(self):
        if hasattr(self, "_raw"):
            return self._raw
        return None

    @raw.setter
    def raw(self, raw):
        self._raw = raw

    @property
    def raw(self):
        if hasattr(self, "_raw"):
            return self._raw
        return None

    @raw.setter
    def raw(self, raw):
        self._raw = raw

    @property
    def is_view(self) -> bool:
        """`True` if object is view of another AnnData object, `False` otherwise."""
        return self._is_view

    @property
    def n_vars(self) -> int:
        """Number of variables/features."""
        return len(self.var_names)

    @property
    def n_obs(self) -> int:
        """Number of observations."""
        return len(self.obs_names)

    def __repr__(self):
        descr = f"AnnData object with n_obs × n_vars = {self.n_obs} × {self.n_vars}"
        for attr in [
            "obs",
            "var",
            "uns",
            "obsm",
            "varm",
            "layers",
            "obsp",
            "varp",
        ]:
            keys = getattr(self, attr).keys()
            if len(keys) > 0:
                descr += f"\n    {attr}: {str(list(keys))[1:-1]}"
        return descr


def read_remote(store: Union[str, Path, MutableMapping, zarr.Group]) -> AnnData:
    if isinstance(store, Path):
        store = str(store)

    is_consolidated = True
    try:
        f = zarr.open_consolidated(store, mode="r")
    except KeyError:
        is_consolidated = False
        f = zarr.open(store, mode="r")

    def callback(func, elem_name: str, elem, iospec):
        if iospec.encoding_type == "anndata" or elem_name.endswith("/"):
            cols = ["obs", "var", "obsm", "varm", "obsp", "varp", "layers", "X", "raw"]
            iter_object = (
                elem.items()
                if is_consolidated
                else [(k, elem[k]) for k in cols if k in elem]
            )
            return AnnDataRemote(
                **{k: read_dispatched(v, callback) for k, v in iter_object}, file=elem
            )
        elif elem_name.startswith("/raw"):
            return None
        elif elem_name in {"/obs", "/var"}:
            # override to only return AxisArray that will be accessed specially via our special AnnData object
            iter_object = (
                elem.items()
                if is_consolidated
                else [(k, elem[k]) for k in elem.attrs["column-order"]]
                + [(elem.attrs["_index"], elem[elem.attrs["_index"]])]
            )
            return {k: read_dispatched(v, callback) for k, v in iter_object}
        elif iospec.encoding_type == "categorical":
            return LazyCategoricalArray(elem)
        elif iospec.encoding_type in {"array", "string-array"}:
            return elem
        elif iospec.encoding_type in {"csr_matrix", "csc_matrix"}:
            return sparse_dataset(elem).to_backed()
        return func(elem)

    adata = read_dispatched(f, callback=callback)

    return adata
