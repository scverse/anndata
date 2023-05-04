from collections import OrderedDict, abc as cabc
from copy import copy, deepcopy
from enum import Enum
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
from anndata._core.aligned_mapping import (
    AlignedMapping,
    Layers,
    PairwiseArrays,
    AxisArraysView,
)
from anndata._core.anndata import StorageType, _check_2d_shape, _gen_dataframe
from anndata._core.anndata_base import AbstractAnnData
from anndata._core.file_backing import AnnDataFileManager
from anndata._core.index import Index, _normalize_indices, _subset
from anndata._core.raw import Raw
from anndata._core.sparse_dataset import sparse_dataset
from anndata._core.views import _resolve_idx, as_view, _resolve_idxs
from anndata._io.specs.registry import read_elem
from anndata.compat import _move_adj_mtx, _read_attr
from anndata.utils import convert_to_dict

import zarr
import pandas as pd
import numpy as np
from xarray.core.indexing import ExplicitlyIndexedNDArrayMixin
import dask.array as da

from ..._core import AnnData, AxisArrays
from .. import read_dispatched


class MaskedArrayMixIn(ExplicitlyIndexedNDArrayMixin):
    def _resolve_idx(self, new_idx):
        return (
            new_idx
            if self.subset_idx is None
            else _resolve_idx(self.subset_idx, new_idx, self.shape[0])
        )

    @property
    def subset_idx(self):
        return self._subset_idx

    @subset_idx.setter
    def subset_idx(self, new_idx):
        self._subset_idx = self._resolve_idx(new_idx)

    @property
    def shape(self) -> Tuple[int, ...]:
        if self.subset_idx is None:
            return self.values.shape
        if isinstance(self.subset_idx, slice):
            if self.subset_idx == slice(None, None, None):
                return self.values.shape
            return (self.subset_idx.stop - self.subset_idx.start,)
        else:
            return (len(self.subset_idx),)

    def __eq__(self, __o) -> np.ndarray:
        return self[()] == __o

    def __ne__(self, __o) -> np.ndarray:
        return ~(self == __o)


class LazyCategoricalArray(MaskedArrayMixIn):
    __slots__ = ("values", "attrs", "_categories", "_categories_cache", "_subset_idx")

    def __init__(self, group, *args, **kwargs):
        """Class for lazily reading categorical data from formatted zarr group

        Args:
            group (zarr.Group): group containing "codes" and "categories" key as well as "ordered" attr
        """
        self.values = group["codes"]
        self._categories = group["categories"]
        self._categories_cache = None
        self._subset_idx = None
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
    def ordered(self):
        return bool(self.attrs["ordered"])

    def __getitem__(self, selection) -> pd.Categorical:
        idx = self._resolve_idx(selection)
        codes = self.values.oindex[idx]
        if codes.shape == ():  # handle 0d case
            codes = np.array([codes])
        return pd.Categorical.from_codes(
            codes=codes,
            categories=self.categories,
            ordered=self.ordered,
        ).remove_unused_categories()

    def __repr__(self) -> str:
        return f"LazyCategoricalArray(codes=..., categories={self.categories}, ordered={self.ordered})"


class LazyMaskedArray(MaskedArrayMixIn):
    __slots__ = ("mask", "values", "_subset_idx", "_dtype_str")

    def __init__(self, group, dtype_str, *args, **kwargs):
        """Class for lazily reading categorical data from formatted zarr group

        Args:
            group (zarr.Group): group containing "codes" and "categories" key as well as "ordered" attr
            dtype_str (Nullable): group containing "codes" and "categories" key as well as "ordered" attr
        """
        self.values = group["values"]
        self.mask = group["mask"] if "mask" in group else None
        self._subset_idx = None
        self._dtype_str = dtype_str

    @property
    def dtype(self) -> pd.CategoricalDtype:
        if self.mask is not None:
            if self._dtype_str == "nullable-integer":
                return pd.arrays.IntegerArray
            elif self._dtype_str == "nullable-boolean":
                return pd.arrays.BooleanArray
        return pd.array

    def __getitem__(self, selection) -> pd.Categorical:
        idx = self._resolve_idx(selection)
        values = self.values[idx]
        if self.mask is not None:
            mask = self.mask[idx]
            if self._dtype_str == "nullable-integer":
                return pd.arrays.IntegerArray(values, mask=mask)
            elif self._dtype_str == "nullable-boolean":
                return pd.arrays.BooleanArray(values, mask=mask)
        return pd.array(values)

    def __repr__(self) -> str:
        if self._dtype_str == "nullable-integer":
            return "LazyNullableIntegerArray"
        elif self._dtype_str == "nullable-boolean":
            return "LazyNullableBooleanArray"


@_subset.register(MaskedArrayMixIn)
def _subset_masked(a: MaskedArrayMixIn, subset_idx: Index):
    a_copy = deepcopy(a)
    a_copy.subset_idx = subset_idx[0]  # this is a tuple?
    return a_copy


@as_view.register(MaskedArrayMixIn)
def _view_masked(a: MaskedArrayMixIn, view_args):
    return a


@as_view.register(pd.Categorical)
def _view_pd_categorical(a: pd.Categorical, view_args):
    return a


@as_view.register(pd.api.extensions.ExtensionArray)
def _view_pd_array(a: pd.api.extensions.ExtensionArray, view_args):
    return a


@as_view.register(pd.arrays.IntegerArray)
def _view_pd_integer_array(a: pd.arrays.IntegerArray, view_args):
    return a


@as_view.register(pd.arrays.BooleanArray)
def _view_pd_boolean_array(a: pd.arrays.BooleanArray, view_args):
    return a


class AxisArraysRemote(AxisArrays):
    def __getattr__(self, __name: str):
        # If we a method has been accessed that is not here, try the pandas implementation
        if hasattr(pd.DataFrame, __name):
            return self.to_df().__getattribute__(__name)
        return object.__getattribute__(self, __name)

    @property
    def iloc(self):
        class IlocDispatch:
            def __getitem__(self_iloc, idx):
                return self._view(self.parent, (idx,))

        return IlocDispatch()

    @property
    def dim_names(self) -> pd.Index:
        return (self.parent.obs_names, self.parent.var_names)[self._axis].compute()


def to_df_1d_axis_arrays(axis_arrays):
    """Convert to pandas dataframe."""
    df = pd.DataFrame(index=axis_arrays.dim_names)
    for key in axis_arrays.keys():
        if "index" not in key:
            df[key] = axis_arrays[key][()]
    return df


class AxisArrays1dRemote(AxisArraysRemote):
    def to_df(self) -> pd.DataFrame:
        return to_df_1d_axis_arrays(self)


class AxisArraysRemoteView(AxisArraysView):
    def to_df(self) -> pd.DataFrame:
        return to_df_1d_axis_arrays(self)


AxisArrays1dRemote._view_class = AxisArraysRemoteView


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
        self.obs = AxisArrays1dRemote(adata_ref, 0, vals=convert_to_dict(obs))
        self.var = AxisArrays1dRemote(adata_ref, 1, vals=convert_to_dict(var))

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
        elif "nullable" in iospec.encoding_type:
            return LazyMaskedArray(elem, iospec.encoding_type)
        elif iospec.encoding_type in {"array", "string-array"}:
            return da.from_zarr(elem)
        elif iospec.encoding_type in {"csr_matrix", "csc_matrix"}:
            return sparse_dataset(elem).to_backed()
        elif iospec.encoding_type in {"awkward-array"}:
            return read_dispatched(elem, None)
        elif iospec.encoding_type in {"dataframe"}:
            return read_dispatched(elem, None)
        return func(elem)

    adata = read_dispatched(f, callback=callback)

    return adata
