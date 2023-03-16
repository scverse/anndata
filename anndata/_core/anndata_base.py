from abc import abstractmethod
from typing import Tuple
import pandas as pd

from .file_backing import AnnDataFileManager
from ..utils import DeprecationMixinMeta


class AbstractAnnData(metaclass=DeprecationMixinMeta):
    def _init_as_actual(
        self,
        X=None,
        obs=None,
        var=None,
        uns=None,
        obsm=None,
        varm=None,
        varp=None,
        obsp=None,
        raw=None,
        layers=None,
        dtype=None,
        shape=None,
        filename=None,
        filemode=None,
    ):
        # view attributes
        self._is_view = False
        self._adata_ref = None
        self._oidx = None
        self._vidx = None
        x_indices = []
        if filename is not None:
            self.file = AnnDataFileManager(self, filename, filemode)
        else:
            self.file = AnnDataFileManager(self, None)
            (
                X,
                obs,
                var,
                uns,
                obsm,
                varm,
                obsp,
                varp,
                layers,
                raw,
                x_indices,
            ) = self._reformat_axes_args_from_X(
                X, obs, var, uns, obsm, varm, obsp, varp, layers, raw
            )
        self._assign_X(X, shape, dtype)

        self._initialize_indices(shape, obs, var)
        assert self.n_obs == self.shape[0]
        assert self.n_vars == self.shape[1]
        if self.X is not None:
            assert self.n_obs == self.X.shape[0]
            assert self.n_vars == self.X.shape[1]

        self._assign_obs(obs)
        self._assign_var(var)
        # now we can verify if indices match!
        for attr_name, x_name, idx in x_indices:
            attr = getattr(self, attr_name)
            if isinstance(attr.index, pd.RangeIndex):
                attr.index = idx
            elif not idx.equals(attr.index):
                raise ValueError(f"Index of {attr_name} must match {x_name} of X.")

        self._assign_uns(uns)
        self._assign_obsm(obsm)
        self._assign_varm(varm)
        self._assign_obsp(obsp)
        self._assign_varp(varp)
        self._assign_layers(layers)
        self._run_checks()
        self._cleanup_raw_and_uns(raw, uns)

    @abstractmethod
    def _init_as_view(self, *args, **kwargs):
        pass

    @abstractmethod
    def _assign_X(self, X, shape, dtype):
        pass

    def _reformat_axes_args_from_X(self, *args):
        return args

    @abstractmethod
    def _initialize_indices(self, *args):
        pass

    @abstractmethod
    def _assign_obs(self, obs):
        pass

    @abstractmethod
    def _assign_var(self, var):
        pass

    @abstractmethod
    def _assign_layers(self, layers):
        pass

    @abstractmethod
    def _assign_uns(self, uns):
        pass

    @abstractmethod
    def _assign_obsm(self, obsm):
        pass

    @abstractmethod
    def _assign_varm(self, varm):
        pass

    @abstractmethod
    def _assign_obsp(self, obsp):
        pass

    @abstractmethod
    def _assign_varp(self, varp):
        pass

    @abstractmethod
    def _assign_layers(self, layers):
        pass

    @abstractmethod
    def _run_checks(self, *args):
        pass

    @abstractmethod
    def _cleanup(self, *args):
        pass

    @property
    @abstractmethod
    def X(self):
        pass

    @property
    @abstractmethod
    def obs(self):
        pass

    @property
    @abstractmethod
    def obsm(self):
        pass

    @property
    @abstractmethod
    def obsp(self):
        pass

    @property
    @abstractmethod
    def var(self):
        pass

    @property
    @abstractmethod
    def uns(self):
        pass

    @property
    @abstractmethod
    def varm(self):
        pass

    @property
    @abstractmethod
    def varp(self):
        pass

    @property
    @abstractmethod
    def n_obs(self) -> int:
        pass

    @property
    @abstractmethod
    def obs_names(self):
        pass

    @property
    @abstractmethod
    def n_vars(self) -> int:
        pass

    @property
    @abstractmethod
    def var_names(self):
        pass

    @property
    @abstractmethod
    def is_view(self) -> bool:
        pass

    @property
    def shape(self) -> Tuple[int, int]:
        """Shape of data matrix (:attr:`n_obs`, :attr:`n_vars`)."""
        return self.n_obs, self.n_vars
