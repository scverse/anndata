from collections import OrderedDict, abc as cabc
from copy import copy
from functools import cached_property
from pathlib import Path
from typing import Any, MutableMapping, Union, List, Sequence
from anndata._core.aligned_mapping import Layers, PairwiseArrays
from anndata._core.anndata import StorageType, _check_2d_shape, _gen_dataframe
from anndata._core.file_backing import AnnDataFileManager
from anndata._core.index import Index
from anndata._core.raw import Raw
from anndata._core.sparse_dataset import SparseDataset
from anndata._core.views import _resolve_idxs
from anndata._io.specs.registry import read_elem
from anndata.compat import _move_adj_mtx, _read_attr
from anndata.utils import convert_to_dict

import zarr
import pandas as pd
import numpy as np

from ..._core import AnnData, AxisArrays
from .. import read_dispatched


# TODO: Do we really need to subclass the Array class here?  Ryan Abernathy seems to say "no"
# but I don't really want to mess with the methods.  The downside is that (for some reason), it's
# reading the `zarray` of the `codes` path which should not have to happen, but I can't figure out a way around it.
class CategoricalZarrArray(zarr.core.Array):
    def __init__(self, group, *args, **kwargs):
        codes_path = group.path + "/codes"
        super().__init__(group.store.store, codes_path, *args, **kwargs)
        self._categories = group["categories"]
        self._group_attrs = group.attrs

    @cached_property
    def categories(self):
        return self._categories[()]

    @cached_property
    def ordered(self):
        return bool(_read_attr(self._group_attrs, "ordered"))

    def __array__(self, *args):  # may need to override this, copied for now
        a = self[...]
        if args:
            a = a.astype(args[0])
        return a

    def __getitem__(self, selection):
        result = super().__getitem__(selection)
        return pd.Categorical.from_codes(
            codes=result,
            categories=self.categories,
            ordered=self.ordered,
        )


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
            z = self[key]
            if isinstance(z, zarr.Group) and "codes" in z:  # catrgoricql
                value = pd.Categorical.from_codes(
                    codes=read_elem(z["codes"]),
                    categories=read_elem(z["categories"]),
                    ordered=bool(_read_attr(z.attrs, "ordered")),
                )
            else:
                value = z[()]
            df[key] = value
        return df

    @property
    def iloc(self):
        class IlocDispatch:
            def __getitem__(self_iloc, idx):
                return self._view(self.parent, (idx,))

        return IlocDispatch()


class AnnDataRemote(AnnData):
    def _assign_X(self, X, shape, dtype):
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
            self._n_obs, self._n_vars = self._X.shape
        else:
            self._X = None

    def _initialize_indices(self, shape, obs, var):
        # annotations - need names already for AxisArrays to work.
        self.obs_names = pd.Index(
            (obs["index"] if "index" in obs else obs["_index"])[()]
        )
        self.var_names = pd.Index(
            (var["index"] if "index" in var else var["_index"])[()]
        )
        if self._X is not None:
            self._n_obs, self._n_vars = self._X.shape
        else:
            self._n_obs = len([] if obs is None else self.obs_names)
            self._n_vars = len([] if var is None else self.var_names)
            # check consistency with shape
            if shape is not None:
                if self._n_obs == 0:
                    self._n_obs = shape[0]
                else:
                    if self._n_obs != shape[0]:
                        raise ValueError("`shape` is inconsistent with `obs`")
                if self._n_vars == 0:
                    self._n_vars = shape[1]
                else:
                    if self._n_vars != shape[1]:
                        raise ValueError("`shape` is inconsistent with `var`")

    # annotations
    def _assign_obs(self, obs):
        self._obs = AxisArraysRemote(self, 0, vals=convert_to_dict(obs))

    def _assign_var(self, var):
        self._var = AxisArraysRemote(self, 1, vals=convert_to_dict(var))

    def _assign_obsm(self, obsm):
        self._obsm = AxisArraysRemote(self, 0, vals=convert_to_dict(obsm))

    def _assign_varm(self, varm):
        self._varm = AxisArraysRemote(self, 1, vals=convert_to_dict(varm))

    def _run_checks(self):
        pass  # for now

    def _init_as_view(self, adata_ref: "AnnData", oidx: Index, vidx: Index):
        if adata_ref.isbacked and adata_ref.is_view:
            raise ValueError(
                "Currently, you cannot index repeatedly into a backed AnnData, "
                "that is, you cannot make a view of a view."
            )
        self._is_view = True
        if isinstance(oidx, (int, np.integer)):
            if not (-adata_ref.n_obs <= oidx < adata_ref.n_obs):
                raise IndexError(f"Observation index `{oidx}` is out of range.")
            oidx += adata_ref.n_obs * (oidx < 0)
            oidx = slice(oidx, oidx + 1, 1)
        if isinstance(vidx, (int, np.integer)):
            if not (-adata_ref.n_vars <= vidx < adata_ref.n_vars):
                raise IndexError(f"Variable index `{vidx}` is out of range.")
            vidx += adata_ref.n_vars * (vidx < 0)
            vidx = slice(vidx, vidx + 1, 1)
        if adata_ref.is_view:
            prev_oidx, prev_vidx = adata_ref._oidx, adata_ref._vidx
            adata_ref = adata_ref._adata_ref
            oidx, vidx = _resolve_idxs((prev_oidx, prev_vidx), (oidx, vidx), adata_ref)
        # self._adata_ref is never a view
        self._adata_ref = adata_ref
        self._oidx = oidx
        self._vidx = vidx
        # the file is the same as of the reference object
        self.file = adata_ref.file
        # views on attributes of adata_ref
        obs_sub = adata_ref.obs.iloc[oidx]
        var_sub = adata_ref.var.iloc[vidx]
        self._obsm = adata_ref.obsm._view(self, (oidx,))
        self._varm = adata_ref.varm._view(self, (vidx,))
        self._layers = adata_ref.layers._view(self, (oidx, vidx))
        self._obsp = adata_ref.obsp._view(self, oidx)
        self._varp = adata_ref.varp._view(self, vidx)
        # fix categories
        uns = copy(adata_ref._uns)
        self._remove_unused_categories(adata_ref.obs, obs_sub, uns)
        self._remove_unused_categories(adata_ref.var, var_sub, uns)
        # set attributes
        self._obs = adata_ref.obs._view(self, (oidx,))
        self._var = adata_ref.var._view(self, (vidx,))
        self._uns = uns
        self._n_obs = len(
            self.obs["index"] if "index" in self.obs else self.obs["_index"]
        )
        self._n_vars = len(
            self.var["index"] if "index" in self.var else self.var["_index"]
        )

        # set data
        if self.isbacked:
            self._X = None

        # set raw, easy, as it’s immutable anyways...
        if adata_ref._raw is not None:
            # slicing along variables axis is ignored
            self._raw = adata_ref.raw[oidx]
            self._raw._adata = self
        else:
            self._raw = None

    # TODO: this is not quite complete in the original but also here, what do we do about this?
    def __delitem__(self, index: Index):
        obs, var = self._normalize_indices(index)
        # TODO: does this really work?
        if not self.isbacked:
            del self._X[obs, var]
        else:
            X = self.file["X"]
            del X[obs, var]
            self._set_backed("X", X)

    def __getitem__(self, index: Index) -> "AnnData":
        """Returns a sliced view of the object."""
        oidx, vidx = self._normalize_indices(index)
        return AnnDataRemote(self, oidx=oidx, vidx=vidx, asview=True)

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


def read_remote(store: Union[str, Path, MutableMapping, zarr.Group]) -> AnnData:
    if isinstance(store, Path):
        store = str(store)

    f = zarr.open_consolidated(store, mode="r")

    def callback(func, elem_name: str, elem, iospec):
        if iospec.encoding_type == "anndata" or elem_name.endswith("/"):
            return AnnDataRemote(
                **{k: read_dispatched(v, callback) for k, v in elem.items()}
            )
        elif elem_name.startswith("/raw"):
            return None
        elif elem_name in {"/obs", "/var"}:
            # override to only return AxisArray that will be accessed specially via our special AnnData object
            return {k: read_dispatched(v, callback) for k, v in elem.items()}
        elif iospec.encoding_type == "categorical":
            return CategoricalZarrArray(elem)
        elif iospec.encoding_type in {"array", "string_array"}:
            return elem
        elif iospec.encoding_type in {"csr_matrix", "csc_matrix"}:
            return SparseDataset(elem).to_backed()
        return func(elem)

    adata = read_dispatched(f, callback=callback)

    return adata
