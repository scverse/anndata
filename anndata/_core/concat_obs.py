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


class _IndxViewMixin:
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
    def __init__(self, attr, adatas, adatas_oidx, reverse, vidx=None):
        self.adatas = adatas
        self.adatas_oidx = adatas_oidx
        self.reverse = reverse
        self.vidx = vidx
        self.attr = attr

    def __getitem__(self, key):
        arrs = []
        for i, oidx in enumerate(self.adatas_oidx):
            if oidx is not None:
                arr = getattr(self.adatas[i], self.attr)[key]
                idx = (oidx, self.vidx) if self.vidx is not None else oidx
                arrs.append(arr[idx])

        _arr = _merge(arrs)

        return _arr[self.reverse] if self.reverse is not None else _arr


class AnnDataConcatView(_IndxViewMixin):
    def __init__(self, adatas, limits, resolved_idx, obs_names, var_names):
        self.adatas = adatas
        self.limits = limits
        self.obs_names = obs_names
        self.var_names = var_names

        self.adatas_oidx, self.oidx, self.vidx, self.reverse = resolved_idx

        self.layers_view = MapObsView(
            "layers", adatas, self.adatas_oidx, self.reverse, self.vidx
        )
        self.obsm_view = MapObsView("obsm", adatas, self.adatas_oidx, self.reverse)
        self.obs_view = MapObsView("obs", adatas, self.adatas_oidx, self.reverse)

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
        return len(self.oidx), self.adatas[0].n_vars

    def __getitem__(self, index: Index):
        oidx, vidx = _normalize_indices(index, self.obs_names, self.var_names)
        resolved_idx = self._resolve_idx(oidx, vidx)

        return AnnDataConcatView(
            self.adatas,
            self.limits,
            resolved_idx,
            self.obs_names[oidx],
            self.var_names[vidx],
        )


class AnnDataConcatObs(_IndxViewMixin):
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

        for a in adatas[1:]:
            if not adatas[0].var_names.equals(a.var_names):
                raise ValueError("Variables in the adatas are different.")
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
        )

    @property
    def shape(self):
        return self.limits[-1], self.adatas[0].n_vars
