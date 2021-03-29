from collections.abc import Mapping
import numpy as np
import pandas as pd

from ..._core.index import _normalize_indices, Index
from ..._core.views import _resolve_idx
from ..._core.merge import concat_arrays, inner_concat_aligned_mapping
from ...logging import anndata_logger as logger

ATTRS = ["obs", "obsm", "layers"]


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
        self._keys = keys
        self.adatas_oidx = adatas_oidx
        self.reverse = reverse
        self.vidx = vidx
        self.attr = attr

    def __getitem__(self, key):
        if self._keys is not None and key not in self._keys:
            raise KeyError(f"No {key} in {self.attr} view")

        arrs = []
        for i, oidx in enumerate(self.adatas_oidx):
            if oidx is not None:
                arr = getattr(self.adatas[i], self.attr)[key]
                idx = (oidx, self.vidx) if self.vidx is not None else oidx
                arrs.append(arr[idx])

        _arr = _merge(arrs) if len(arrs) > 1 else arrs[0]

        return _arr[self.reverse] if self.reverse is not None else _arr

    def __repr__(self):
        if self._keys is not None:
            keys = self._keys
        else:
            keys = list(getattr(self.adatas[0], self.attr).keys())
        descr = f"View of {self.attr} with keys: {str(keys)[1:-1]}"
        return descr


class AnnDataConcatView(_ConcatViewMixin):
    def __init__(self, reference, resolved_idx, obs_names, var_names):
        self.reference = reference

        self.adatas = self.reference.adatas
        self.limits = self.reference.limits

        self.obs_names = obs_names
        self.var_names = var_names

        self.adatas_oidx, self.oidx, self.vidx, self.reverse = resolved_idx

        self._view_attrs_keys = self.reference._view_attrs_keys
        self._attrs = self.reference._attrs

        for attr, keys in self._view_attrs_keys.items():
            set_vidx = self.vidx if attr == "layers" else None
            setattr(
                self,
                attr + "_view",
                MapObsView(
                    attr, self.adatas, keys, self.adatas_oidx, self.reverse, set_vidx
                ),
            )
        # non-view attrs
        for attr in self._attrs:
            setattr(
                self,
                attr + "_view",
                MapObsView(attr, [self.reference], None, [self.oidx], None),
            )

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
            self.reference,
            resolved_idx,
            self.obs_names[oidx],
            self.var_names[vidx],
        )

    def __repr__(self):
        n_obs, n_vars = self.shape
        descr = f"AnnDataConcatView object with n_obs × n_vars = {n_obs} × {n_vars}"
        all_attrs_keys = self._view_attrs_keys.copy()
        for attr in self._attrs:
            all_attrs_keys[attr] = list(getattr(self.reference, attr).keys())
        for attr, keys in all_attrs_keys.items():
            if len(keys) > 0:
                descr += f"\n    {attr}: {str(keys)[1:-1]}"
        return descr


class AnnDataConcatObs(_ConcatViewMixin):
    def __init__(
        self,
        adatas,
        join_obs="inner",
        join_obsm=None,
        label=None,
        keys=None,
        index_unique=None,
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

    def __getitem__(self, index: Index):
        oidx, vidx = _normalize_indices(index, self.obs_names, self.var_names)
        resolved_idx = self._resolve_idx(oidx, vidx)

        return AnnDataConcatView(
            self, resolved_idx, self.obs_names[oidx], self.var_names[vidx]
        )

    @property
    def shape(self):
        return self.limits[-1], self.adatas[0].n_vars

    def __repr__(self):
        n_obs, n_vars = self.shape
        descr = f"AnnDataConcatObs object with n_obs × n_vars = {n_obs} × {n_vars}"
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
