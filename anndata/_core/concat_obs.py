from collections.abc import Mapping
import numpy as np
import pandas as pd

from .index import _normalize_indices, Index
from .merge import concat_arrays
from ..logging import anndata_logger as logger


def merge(arrs):
    rxers = [lambda x, fill_value, axis: x] * len(arrs)
    return concat_arrays(arrs, rxers)


class MapObsView:
    def __init__(self, attr, adatas, adatas_oidx, reverse, vidx):
        self.adatas = adatas
        self.adatas_oidx = adatas_oidx
        self.reverse = reverse
        self.vidx = vidx
        self.attr = attr

    def set_index(adatas_oidx, reverse, vidx):
        self.adatas_oidx = adatas_oidx
        self.reverse = reverse
        self.vidx = vidx

    def __getitem__(self, key):
        print(self.vidx)
        arrs = []
        for i, oidx in enumerate(self.adatas_oidx):
            if oidx is not None:
                arr = getattr(self.adatas[i], self.attr)[key]
                idx = (oidx, self.vidx) if len(arr.shape) > 1 else oidx
                arrs.append(arr[idx])

        _arr = merge(arrs)

        return _arr[self.reverse] if self.reverse is not None else _arr


class AnnDataConcatView:
    def __init__(self, adatas, adatas_oidx, reverse, vidx):
        self.adatas = adatas
        self.adatas_oidx = adatas_oidx
        self.reverse = reverse
        self.vidx = vidx

        self.obsm_view = MapObsView("obsm", adatas, adatas_oidx, reverse, vidx)
        self.obs_view = MapObsView("obs", adatas, adatas_oidx, reverse, vidx)
        self.layers_view = MapObsView("layers", adatas, adatas_oidx, reverse, vidx)

    @property
    def X(self):
        Xs = []
        for i, oidx in enumerate(self.adatas_oidx):
            if oidx is not None:
                Xs.append(self.adatas[i].X[oidx, self.vidx])

        _X = merge(Xs)

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

    def set_index(self, adatas_oidx, reverse, vidx):
        self.adatas_oidx = adatas_oidx
        self.reverse = reverse
        self.vidx = vidx

        self.obs_view.set_index(adatas_oidx, reverse, vidx)
        self.obsm_view.set_index(adatas_oidx, reverse, vidx)
        self.layers_view.set_index(adatas_oidx, reverse, vidx)


class AnnDataConcatObs:
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

        self.concat_view = None

    def __getitem__(self, index: Index):
        reverse = None
        adatas_oidx = []

        oidx, vidx = _normalize_indices(index, self.obs_names, self.var_names)
        if isinstance(oidx, slice):
            start, stop, step = oidx.indices(self.limits[-1])
            u_oidx = np.arange(start, stop, step)
        else:
            if not np.all(np.diff(oidx) > 0):
                u_oidx, reverse = np.unique(oidx, return_inverse=True)
            else:
                u_oidx = oidx
        for lower, upper in zip([0] + self.limits, self.limits):
            mask = (u_oidx >= lower) & (u_oidx < upper)
            adatas_oidx.append(u_oidx[mask] - lower if mask.any() else None)

        if self.concat_view is None:
            self.concat_view = AnnDataConcatView(
                self.adatas, adatas_oidx, reverse, vidx
            )
        else:
            self.concat_view.set_index(adatas_oidx, reverse, vidx)

        return self.concat_view
