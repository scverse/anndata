from collections.abc import Mapping
import numpy as np
import pandas as pd

from .index import _normalize_indices, Index
from .views import _resolve_idx
from .merge import concat_arrays
from ..logging import anndata_logger as logger


def _merge(arrs):
    rxers = [lambda x, fill_value, axis: x] * len(arrs)
    return concat_arrays(arrs, rxers)


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
            if oidx.size > 1 and not np.all(np.diff(oidx) > 0):
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
    def __init__(self, attr, adatas, keys, adatas_oidx, reverse, vidx=None):
        self.adatas = adatas
        self.keys = keys
        self.adatas_oidx = adatas_oidx
        self.reverse = reverse
        self.vidx = vidx
        self.attr = attr

    def __getitem__(self, key):
        if key not in self.keys:
            raise KeyError(f'No {key} in {self.attr} view')

        arrs = []
        for i, oidx in enumerate(self.adatas_oidx):
            if oidx is not None:
                arr = getattr(self.adatas[i], self.attr)[key]
                idx = (oidx, self.vidx) if self.vidx is not None else oidx
                arrs.append(arr[idx])

        _arr = _merge(arrs)

        return _arr[self.reverse] if self.reverse is not None else _arr

    def __repr__(self):
        descr = f"View of {self.attr} with keys: {str(self.keys)[1:-1]}"
        return descr


class AnnDataConcatView(_ConcatViewMixin):
    def __init__(self, adatas, limits, resolved_idx, obs_names, var_names, attrs_keys):
        self.adatas = adatas
        self.limits = limits
        self.obs_names = obs_names
        self.var_names = var_names

        self.adatas_oidx, self.oidx, self.vidx, self.reverse = resolved_idx

        for attr, keys in attrs_keys.items():
            set_vidx = self.vidx if attr == "layers" else None
            setattr(
                self,
                attr + "_view",
                MapObsView(
                    attr, adatas, keys, self.adatas_oidx, self.reverse, set_vidx
                )
            )

        self._attrs_keys = attrs_keys

    @property
    def X(self):
        Xs = []
        for i, oidx in enumerate(self.adatas_oidx):
            if oidx is not None:
                Xs.append(self.adatas[i].X[oidx, self.vidx])

        _X = _merge(Xs)

        return _X[self.reverse] if self.reverse is not None else _X

    @property
    def layers(self):
        return self.layers_view

    @property
    def obsm(self):
        return self.obsm_view

    @property
    def obs(self):
        return self.obs_view

    @property
    def shape(self):
        return len(self.oidx), len(self.var_names)

    def __getitem__(self, index: Index):
        oidx, vidx = _normalize_indices(index, self.obs_names, self.var_names)
        resolved_idx = self._resolve_idx(oidx, vidx)

        return AnnDataConcatView(
            self.adatas,
            self.limits,
            resolved_idx,
            self.obs_names[oidx],
            self.var_names[vidx],
            self._attrs_keys,
        )

    def __repr__(self):
        n_obs, n_vars = self.shape
        descr = f"AnnDataConcatView object with n_obs × n_vars = {n_obs} × {n_vars}"
        for attr, keys in self._attrs_keys.items():
            if len(keys) > 0:
                descr += f"\n    {attr}: {str(keys)[1:-1]}"
        return descr


class AnnDataConcatObs(_ConcatViewMixin):
    def __init__(self, adatas, keys=None, index_unique=None):
        if isinstance(adatas, Mapping):
            if keys is not None:
                raise TypeError(
                    "Cannot specify categories in both mapping keys and using `keys`. "
                    "Only specify this once."
                )
            keys, adatas = list(adatas.keys()), list(adatas.values())
        else:
            adatas = list(adatas)

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

        self._attrs_keys = {}
        for attr in ["obs", "obsm", "layers"]:
            self._attrs_keys[attr] = list(getattr(adatas[0], attr).keys())

        for a in adatas[1:]:
            if not adatas[0].var_names.equals(a.var_names):
                raise ValueError("Variables in the adatas are different.")
            for attr, keys in self._attrs_keys.items():
                ai_attr = getattr(a, attr)
                a0_attr = getattr(adatas[0], attr)
                new_keys = []
                for key in keys:
                    if key in ai_attr.keys():
                        a0_ashape = a0_attr[key].shape
                        ai_ashape = ai_attr[key].shape
                        if len(a0_ashape) < 2 or a0_ashape[1] == ai_ashape[1]:
                            new_keys.append(key)
                self._attrs_keys[attr] = new_keys

        self.var_names = adatas[0].var_names

        self.adatas = adatas

        self.limits = [adatas[0].n_obs]
        for i in range(len(adatas) - 1):
            self.limits.append(self.limits[i] + adatas[i + 1].n_obs)

        self.obs = pd.DataFrame(index=self.obs_names)

    def __getitem__(self, index: Index):
        oidx, vidx = _normalize_indices(index, self.obs_names, self.var_names)
        resolved_idx = self._resolve_idx(oidx, vidx)

        return AnnDataConcatView(
            self.adatas,
            self.limits,
            resolved_idx,
            self.obs_names[oidx],
            self.var_names[vidx],
            self._attrs_keys,
        )

    @property
    def shape(self):
        return self.limits[-1], self.adatas[0].n_vars

    def __repr__(self):
        n_obs, n_vars = self.shape
        descr = f"AnnDataConcatObs object with n_obs × n_vars = {n_obs} × {n_vars}"
        for attr, keys in self._attrs_keys.items():
            if len(keys) > 0:
                descr += f"\n    view of {attr}: {str(keys)[1:-1]}"

        obs_keys = list(self.obs.keys())
        if len(obs_keys) > 0:
            descr += f"\n    obs: {str(obs_keys)[1:-1]}"

        return descr
