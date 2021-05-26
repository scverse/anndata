from collections.abc import Mapping
from functools import reduce
import numpy as np
import pandas as pd

from ..._core.anndata import AnnData
from ..._core.index import _normalize_indices, Index
from ..._core.views import _resolve_idx
from ..._core.merge import concat_arrays, inner_concat_aligned_mapping
from ...logging import anndata_logger as logger

ATTRS = ["obs", "obsm", "layers"]


def _merge(arrs):
    rxers = [lambda x, fill_value, axis: x] * len(arrs)
    return concat_arrays(arrs, rxers)


def _select_convert(key, convert, arr=None):
    key_convert = None

    if callable(convert):
        key_convert = convert
    elif isinstance(convert, dict) and key in convert:
        key_convert = convert[key]

    if arr is not None:
        return key_convert(arr) if key_convert is not None else arr
    else:
        return key_convert


class _ConcatViewMixin:
    def _resolve_idx(self, oidx, vidx):
        reverse = None
        adatas_oidx = []

        old_oidx = getattr(self, "oidx", None)
        if old_oidx is not None:
            oidx = _resolve_idx(old_oidx, oidx, self.limits[-1])

        if isinstance(oidx, slice):
            start, stop, step = oidx.indices(self.limits[-1])
            u_oidx = np.arange(start, stop, step)
        else:
            if self.has_backed and oidx.size > 1 and not np.all(np.diff(oidx) > 0):
                u_oidx, reverse = np.unique(oidx, return_inverse=True)
            else:
                u_oidx = oidx

        for lower, upper in zip([0] + self.limits, self.limits):
            mask = (u_oidx >= lower) & (u_oidx < upper)
            adatas_oidx.append(u_oidx[mask] - lower if mask.any() else None)

        old_vidx = getattr(self, "vidx", None)
        if old_vidx is not None:
            vidx = _resolve_idx(old_vidx, vidx, self.adatas[0].n_vars)

        return adatas_oidx, oidx, vidx, reverse


class MapObsView:
    def __init__(
        self, attr, adatas, keys, adatas_oidx, reverse, vidx=None, convert=None
    ):
        self.adatas = adatas
        self._keys = keys
        self.adatas_oidx = adatas_oidx
        self.reverse = reverse
        self.vidx = vidx
        self.attr = attr
        self.convert = convert

    def __getitem__(self, key, use_convert=True):
        if self._keys is not None and key not in self._keys:
            raise KeyError(f"No {key} in {self.attr} view")

        arrs = []
        for i, oidx in enumerate(self.adatas_oidx):
            if oidx is not None:
                arr = getattr(self.adatas[i], self.attr)[key]
                idx = (oidx, self.vidx) if self.vidx is not None else oidx
                arrs.append(arr[idx])

        _arr = _merge(arrs) if len(arrs) > 1 else arrs[0]
        _arr = _arr[self.reverse] if self.reverse is not None else _arr

        if self.convert is not None and use_convert:
            _arr = _select_convert(key, self.convert, _arr)

        return _arr

    def keys(self):
        if self._keys is not None:
            return self._keys
        else:
            return list(getattr(self.adatas[0], self.attr).keys())

    def to_dict(self):
        dct = {}
        for key in self.keys():
            dct[key] = self.__getitem__(key, False)
        return dct

    def __repr__(self):
        descr = f"View of {self.attr} with keys: {str(self.keys())[1:-1]}"
        return descr


class AnnDataSetView(_ConcatViewMixin):
    def __init__(self, reference, resolved_idx):
        self.reference = reference

        self.adatas = self.reference.adatas
        self.limits = self.reference.limits

        self.adatas_oidx, self.oidx, self.vidx, self.reverse = resolved_idx

        self._view_attrs_keys = self.reference._view_attrs_keys
        self._attrs = self.reference._attrs

        self._layers_view, self._obsm_view, self._obs_view = None, None, None
        self._X = None

        self._convert = None
        self._convert_X = None
        self.convert = reference.convert

    def _lazy_init_attr(self, attr, set_vidx=False):
        if getattr(self, f"_{attr}_view") is not None:
            return
        keys = None
        reverse = None
        if attr in self._view_attrs_keys:
            keys = self._view_attrs_keys[attr]
            if len(keys) == 0:
                return
            adatas = self.adatas
            adatas_oidx = self.adatas_oidx
            reverse = self.reverse
        else:
            adatas = [self.reference]
            adatas_oidx = [self.oidx]
        vidx = self.vidx if set_vidx else None

        attr_convert = None
        if self.convert is not None:
            attr_convert = _select_convert(attr, self.convert)

        setattr(
            self,
            f"_{attr}_view",
            MapObsView(attr, adatas, keys, adatas_oidx, reverse, vidx, attr_convert),
        )

    def _gather_X(self):
        if self._X is not None:
            return self._X

        Xs = []
        for i, oidx in enumerate(self.adatas_oidx):
            if oidx is not None:
                adata = self.adatas[i]
                if adata.is_view and adata.isbacked:
                    adata_ref = adata._adata_ref
                    vidx = _resolve_idx(adata._vidx, self.vidx, adata_ref.n_vars)
                    X = adata_ref.X
                else:
                    vidx = self.vidx
                    X = adata.X

                try:
                    Xs.append(X[oidx, vidx])
                except TypeError:
                    # only one indexing array is allowed for advanced selection
                    # very inefficient, needs to be fixed
                    Xs.append(X[oidx][:, vidx])

        _X = _merge(Xs) if len(Xs) > 1 else Xs[0]
        _X = _X[self.reverse] if self.reverse is not None else _X

        self._X = _X

        return _X

    @property
    def X(self):
        _X = self._gather_X()

        return self._convert_X(_X) if self._convert_X is not None else _X

    @property
    def layers(self):
        self._lazy_init_attr("layers", set_vidx=True)
        return self._layers_view

    @property
    def obsm(self):
        self._lazy_init_attr("obsm")
        return self._obsm_view

    @property
    def obs(self):
        self._lazy_init_attr("obs")
        return self._obs_view

    @property
    def obs_names(self):
        return self.reference.obs_names[self.oidx]

    @property
    def var_names(self):
        return self.reference.var_names[self.vidx]

    @property
    def shape(self):
        return len(self.obs_names), len(self.var_names)

    @property
    def convert(self):
        return self._convert

    @convert.setter
    def convert(self, value):
        self._convert = value
        self._convert_X = _select_convert("X", self._convert)
        for attr in ATTRS:
            setattr(self, f"_{attr}_view", None)

    def __len__(self):
        return len(self.obs_names)

    def __getitem__(self, index: Index):
        oidx, vidx = _normalize_indices(index, self.obs_names, self.var_names)
        resolved_idx = self._resolve_idx(oidx, vidx)

        return AnnDataSetView(self.reference, resolved_idx)

    @property
    def has_backed(self):
        for i, adata in enumerate(self.adatas):
            if adata.isbacked and self.adatas_oidx[i] is not None:
                return True
        return False

    def __repr__(self):
        n_obs, n_vars = self.shape
        descr = f"AnnDataSetView object with n_obs × n_vars = {n_obs} × {n_vars}"
        all_attrs_keys = self._view_attrs_keys.copy()
        for attr in self._attrs:
            all_attrs_keys[attr] = list(getattr(self.reference, attr).keys())
        for attr, keys in all_attrs_keys.items():
            if len(keys) > 0:
                descr += f"\n    {attr}: {str(keys)[1:-1]}"
        return descr

    def to_adata(self, ignore_X=False, ignore_layers=False):
        if ignore_layers or self.layers is None:
            layers = None
        else:
            layers = self.layers.to_dict()
        obsm = None if self.obsm is None else self.obsm.to_dict()
        obs = None if self.obs is None else pd.DataFrame(self.obs.to_dict())

        if ignore_X:
            X = None
            shape = self.shape
        else:
            X = self._gather_X()
            shape = None

        adata = AnnData(X, obs=obs, obsm=obsm, layers=layers, shape=shape)
        adata.obs_names = self.obs_names
        adata.var_names = self.var_names
        return adata


class AnnDataSet(_ConcatViewMixin):
    def __init__(
        self,
        adatas,
        join_obs="inner",
        join_obsm=None,
        join_vars=None,
        label=None,
        keys=None,
        index_unique=None,
        convert=None,
    ):
        if isinstance(adatas, Mapping):
            if keys is not None:
                raise TypeError(
                    "Cannot specify categories in both mapping keys and using `keys`. "
                    "Only specify this once."
                )
            keys, adatas = list(adatas.keys()), list(adatas.values())
        else:
            adatas = list(adatas)

        # check if the variables are the same in all adatas
        vars_names_list = [adata.var_names for adata in adatas]
        vars_eq = all([adatas[0].var_names.equals(vrs) for vrs in vars_names_list[1:]])
        if vars_eq:
            self.var_names = adatas[0].var_names
        elif join_vars == "inner":
            var_names = reduce(pd.Index.intersection, vars_names_list)
            adatas = [adata[:, var_names] for adata in adatas]
            self.var_names = var_names
        else:
            raise ValueError(
                "Adatas have different variables. "
                "Please specify join_vars='inner' for intersection."
            )

        if keys is None:
            keys = np.arange(len(adatas)).astype(str)

        concat_indices = pd.concat(
            [pd.Series(a.obs_names) for a in adatas], ignore_index=True
        )
        label_col = pd.Categorical.from_codes(
            np.repeat(np.arange(len(adatas)), [a.shape[0] for a in adatas]),
            categories=keys,
        )
        if index_unique is not None:
            concat_indices = concat_indices.str.cat(
                label_col.map(str), sep=index_unique
            )
        self.obs_names = pd.Index(concat_indices)

        if not self.obs_names.is_unique:
            logger.info("Observation names are not unique.")

        view_attrs = ATTRS.copy()

        self._attrs = []
        # process obs joins
        if join_obs is not None:
            view_attrs.remove("obs")
            self._attrs.append("obs")
            concat_annot = pd.concat(
                [a.obs for a in adatas], join=join_obs, ignore_index=True
            )
            concat_annot.index = self.obs_names
            self.obs = concat_annot
        else:
            self.obs = pd.DataFrame(index=self.obs_names)
        if label is not None:
            self.obs[label] = label_col

        # process obsm inner join
        self.obsm = None
        if join_obsm == "inner":
            view_attrs.remove("obsm")
            self._attrs.append("obsm")
            self.obsm = inner_concat_aligned_mapping(
                [a.obsm for a in adatas], index=self.obs_names
            )

        # process inner join of views
        self._view_attrs_keys = {}
        for attr in view_attrs:
            self._view_attrs_keys[attr] = list(getattr(adatas[0], attr).keys())

        for a in adatas[1:]:
            if not adatas[0].var_names.equals(a.var_names):
                raise ValueError("Variables in the adatas are different.")
            for attr, keys in self._view_attrs_keys.items():
                ai_attr = getattr(a, attr)
                a0_attr = getattr(adatas[0], attr)
                new_keys = []
                for key in keys:
                    if key in ai_attr.keys():
                        a0_ashape = a0_attr[key].shape
                        ai_ashape = ai_attr[key].shape
                        if len(a0_ashape) < 2 or a0_ashape[1] == ai_ashape[1]:
                            new_keys.append(key)
                self._view_attrs_keys[attr] = new_keys

        self.var_names = adatas[0].var_names

        self.adatas = adatas

        self.limits = [adatas[0].n_obs]
        for i in range(len(adatas) - 1):
            self.limits.append(self.limits[i] + adatas[i + 1].n_obs)

        self.convert = convert

    def __getitem__(self, index: Index):
        oidx, vidx = _normalize_indices(index, self.obs_names, self.var_names)
        resolved_idx = self._resolve_idx(oidx, vidx)

        return AnnDataSetView(self, resolved_idx)

    @property
    def shape(self):
        return self.limits[-1], self.adatas[0].n_vars

    def __len__(self):
        return self.limits[-1]

    def to_adata(self):
        if "obs" in self._view_attrs_keys or "obsm" in self._view_attrs_keys:
            concat_view = self[self.obs_names]

        if "obsm" in self._view_attrs_keys:
            obsm = concat_view.obsm.to_dict() if concat_view.obsm is not None else None
        else:
            obsm = self.obsm.copy()

        obs = self.obs.copy()
        if "obs" in self._view_attrs_keys and concat_view.obs is not None:
            for key, value in concat_view.obs.to_dict().items():
                obs[key] = value

        adata = AnnData(X=None, obs=obs, obsm=obsm, shape=self.shape)
        adata.obs_names = self.obs_names
        adata.var_names = self.var_names
        return adata

    @property
    def has_backed(self):
        return any([adata.isbacked for adata in self.adatas])

    @property
    def attrs_keys(self):
        _attrs_keys = {}
        for attr in self._attrs:
            keys = list(getattr(self, attr).keys())
            _attrs_keys[attr] = keys
        _attrs_keys.update(self._view_attrs_keys)
        return _attrs_keys

    def __repr__(self):
        n_obs, n_vars = self.shape
        descr = f"AnnDataSet object with n_obs × n_vars = {n_obs} × {n_vars}"
        descr += f"\n  constructed from {len(self.adatas)} AnnData objects"
        for attr, keys in self._view_attrs_keys.items():
            if len(keys) > 0:
                descr += f"\n    view of {attr}: {str(keys)[1:-1]}"
        for attr in self._attrs:
            keys = list(getattr(self, attr).keys())
            if len(keys) > 0:
                descr += f"\n    {attr}: {str(keys)[1:-1]}"
        if "obs" in self._view_attrs_keys:
            keys = list(self.obs.keys())
            if len(keys) > 0:
                descr += f"\n    own obs: {str(keys)[1:-1]}"

        return descr
