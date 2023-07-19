from collections import OrderedDict, abc as cabc
from pathlib import Path
from typing import (
    MutableMapping,
    Optional,
    Union,
    Sequence,
    Tuple,
)

import h5py
from anndata._core.aligned_mapping import (
    Layers,
    PairwiseArrays,
    AxisArrays
)
from anndata._core.anndata import StorageType, _check_2d_shape
from anndata._core.anndata_base import AbstractAnnData
from anndata._core.index import Index, _normalize_indices, _subset
from anndata._core.raw import Raw
from anndata._core.sparse_dataset import BaseCompressedSparseDataset, sparse_dataset
from anndata._core.views import _resolve_idxs
from anndata._io.h5ad import read_dataset
from anndata.compat import DaskArray
from anndata.utils import asarray, convert_to_dict

import zarr
import xarray as xr
import pandas as pd
import dask.array as da
from scipy import sparse
from ..._core import AnnData
from .. import read_dispatched
from .lazy_arrays import LazyCategoricalArray, LazyMaskedArray

class AnnDataBacked(AbstractAnnData):
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
        self.obs_names = pd.Index(obs['obs_names'].data.compute() if isinstance(obs['obs_names'].data, DaskArray) else obs['obs_names'].data)
        if oidx is not None:
            self.obs_names = self.obs_names[oidx]
        self.var_names = pd.Index(var['var_names'].data.compute() if isinstance(var['var_names'].data, DaskArray) else var['var_names'].data)
        if vidx is not None:
            self.var_names = self.var_names[vidx]
        
        self.obs = xr.Dataset(obs)
        self.var = xr.Dataset(var)

        self.obsm = AxisArrays(adata_ref, 0, vals=convert_to_dict(obsm))
        self.varm = AxisArrays(adata_ref, 1, vals=convert_to_dict(varm))
        self.obsp = PairwiseArrays(adata_ref, 0, vals=convert_to_dict(obsp))
        self.varp = PairwiseArrays(adata_ref, 1, vals=convert_to_dict(varp))
        self.layers = Layers(adata_ref, layers)
        if self.is_view:
            self.obs = self.obs.isel(obs_names=oidx)
            self.var = self.var.isel(var_names=vidx)
            self.obsm = self.obsm._view(self, oidx)
            self.varm = self.varm._view(self, vidx)
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

    def _sanitize(self):
        pass

    def _run_checks(self):
        assert len(self.obs_names) == self.shape[0]
        assert len(self.var_names) == self.shape[1]

    def __getitem__(self, index: Index) -> "AnnData":
        """Returns a sliced view of the object."""
        oidx, vidx = self._normalize_indices(index)
        return AnnDataBacked(self, oidx=oidx, vidx=vidx)

    def _normalize_indices(self, index: Optional[Index]) -> Tuple[slice, slice]:
        return _normalize_indices(
            index,
            pd.Index(self.obs_names.compute() if isinstance(self.obs_names, DaskArray) else self.obs_names), # could be either dask or in-memory from xarray
            pd.Index(self.var_names.compute() if isinstance(self.var_names, DaskArray) else self.var_names),
        )

    def to_memory(self, exclude=[]):
        # handling for AxisArrays
        def backed_dict_to_memory(d, prefix):
            res = {}
            for k, v in d.items():
                full_key = prefix + "/" + k
                if any([full_key == exclude_key for exclude_key in exclude]):
                    continue
                if isinstance(v, DaskArray):
                    res[k] = v.compute()
                elif isinstance(v, BaseCompressedSparseDataset):
                    res[k] = v.to_memory()
                else:
                    res[k] = v
            return res

        # nullable and categoricals need special handling because xarray will convert them to numpy arrays first with dtype object
        def get_nullable_and_categorical_cols(ds):
            cols = []
            for c in ds:
                dtype = ds[c].dtype
                if isinstance(dtype, pd.CategoricalDtype) or dtype == pd.arrays.BooleanArray or dtype == pd.arrays.IntegerArray:
                    cols += [c]
            return cols     
        def to_df(ds):
            nullable_and_categorical_df_cols = get_nullable_and_categorical_cols(ds)
            df = ds.drop_vars(list(set(exclude + nullable_and_categorical_df_cols))).to_dataframe()
            for c in nullable_and_categorical_df_cols:
                df[c] = ds[c].data[()]
            df.index.name = None # matches old AnnData object
            df = df[list(ds.keys())]
            return df
        
        obs = to_df(self.obs)
        var = to_df(self.var)
        obsm = backed_dict_to_memory(convert_to_dict(self.obsm), "obsm")
        varm = backed_dict_to_memory(convert_to_dict(self.varm), "varm")
        varp = backed_dict_to_memory(convert_to_dict(self.varp), "varp")
        obsp = backed_dict_to_memory(convert_to_dict(self.obsp), "obsp")
        layers = backed_dict_to_memory(dict(self.layers), "layers")
        X = None
        if "X" not in exclude:
            if isinstance(self.X, BaseCompressedSparseDataset):
                X = self.X.to_memory()
            else:
                X = self.X.compute()
        return AnnData(
            X=X,
            obs=obs,
            var=var,
            obsm=obsm,
            varm=varm,
            obsp=obsp,
            varp=varp,
            layers=layers,
            uns=self.uns,
        )

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
        descr = (
            f"AnnDataBacked object with n_obs × n_vars = {self.n_obs} × {self.n_vars}"
        )
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


def read_backed(store: Union[str, Path, MutableMapping, zarr.Group, h5py.Dataset]) -> AnnData:
    """Lazily read in on-disk/in-cloud AnnData stores.  A new, but familiar, AnnData object will be returned.
    No array data should need to be read into memory, with exception of non-obs/var dataframes and Awkward Arrays.

    Args:
        store (Union[str, Path, MutableMapping, zarr.Group, h5py.Dataset]): A store-like object to be read in.  If `zarr`, it is best
        for it to be consolidated.

    Returns:
        AnnData: A lazily read-in AnnData object.
    """
    is_h5 = False
    if isinstance(store, Path) or isinstance(store, str):
        store = str(store)
        if store.endswith('h5ad'):
            is_h5 = True

    has_keys = True # true if consolidated or h5ad
    if not is_h5:
        try:
            f = zarr.open_consolidated(store, mode="r")
        except KeyError:
            has_keys = False
            f = zarr.open(store, mode="r")
    else:
         f = h5py.File(store, mode="r")

    def callback(func, elem_name: str, elem, iospec):
        if iospec.encoding_type == "anndata" or elem_name.endswith("/"):
            cols = ["obs", "var", "obsm", "varm", "obsp", "varp", "layers", "X", "raw"]
            iter_object = (
                elem.items()
                if has_keys
                else [(k, elem[k]) for k in cols if k in elem]
            )
            return AnnDataBacked(
                **{k: read_dispatched(v, callback) for k, v in iter_object}, file=elem
            )
        elif elem_name.startswith("/raw"):
            return None
        elif elem_name in {"/obs", "/var"}:
            iter_object = [(k, elem[k]) for k in elem.attrs["column-order"]] + [
                (elem.attrs["_index"], elem[elem.attrs["_index"]])
            ]
            d = {k: read_dispatched(v, callback) for k, v in iter_object}
            d_with_xr = {}
            index_label = f'{elem_name.replace("/", "")}_names'
            for k in d:
                v = d[k]
                if type(v) == DaskArray and k != elem.attrs["_index"]:
                    d_with_xr[k] = xr.DataArray(v, coords=[d[elem.attrs["_index"]]], dims=[index_label], name=k)
                elif (type(v) == LazyCategoricalArray or type(v) == LazyMaskedArray) and k != elem.attrs["_index"]:
                    d_with_xr[k] = xr.DataArray(xr.core.indexing.LazilyIndexedArray(v), coords=[d[elem.attrs["_index"]]], dims=[index_label], name=k)
                elif k == elem.attrs["_index"]:
                    d_with_xr[index_label] = xr.DataArray(v, coords=[v], dims=[index_label], name=index_label)
                else:
                    d_with_xr[k] = v
            return d_with_xr
        elif iospec.encoding_type == "categorical":
            return LazyCategoricalArray(elem["codes"], elem["categories"], elem.attrs)
        elif "nullable" in iospec.encoding_type:
            return LazyMaskedArray(
                elem["values"],
                elem["mask"] if "mask" in elem else None,
                iospec.encoding_type,
            )
        elif iospec.encoding_type in {"array", "string-array"}:
            if is_h5:
                if iospec.encoding_type == "string-array":
                    elem = read_dataset(elem)
                if not hasattr(elem, "chunks") or elem.chunks is None:
                    return da.from_array(elem, chunks=(1000,) * len(elem.shape))
                return da.from_array(elem)
            return da.from_zarr(elem)
        elif iospec.encoding_type in {"csr_matrix", "csc_matrix"}:
            return sparse_dataset(elem)
        elif iospec.encoding_type in {"awkward-array"}:
            return read_dispatched(elem, None)
        elif iospec.encoding_type in {"dataframe"}:
            return read_dispatched(elem, None)
        return func(elem)

    adata = read_dispatched(f, callback=callback)

    return adata
