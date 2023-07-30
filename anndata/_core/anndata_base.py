from abc import abstractmethod
from typing import Tuple

from ..utils import DeprecationMixinMeta


class AbstractAnnData(metaclass=DeprecationMixinMeta):
    @property
    @abstractmethod
    def X(self):
        pass

    @X.setter
    @abstractmethod
    def X(self, X):
        pass

    @property
    @abstractmethod
    def obs(self):
        pass

    @obs.setter
    @abstractmethod
    def obs(self, obs):
        pass

    @property
    @abstractmethod
    def obsm(self):
        pass

    @obsm.setter
    @abstractmethod
    def obsm(self, obsm):
        pass

    @property
    @abstractmethod
    def obsp(self):
        pass

    @obsp.setter
    @abstractmethod
    def obsp(self, obsp):
        pass

    @property
    @abstractmethod
    def var(self):
        pass

    @var.setter
    @abstractmethod
    def var(self, var):
        pass

    @property
    @abstractmethod
    def uns(self):
        pass

    @uns.setter
    @abstractmethod
    def uns(self, uns):
        pass

    @property
    @abstractmethod
    def varm(self):
        pass

    @varm.setter
    @abstractmethod
    def varm(self, varm):
        pass

    @property
    @abstractmethod
    def varp(self):
        pass

    @varp.setter
    @abstractmethod
    def varp(self, varp):
        pass

    @property
    @abstractmethod
    def raw(self):
        pass

    @raw.setter
    @abstractmethod
    def raw(self, raw):
        pass

    @property
    def n_obs(self) -> int:
        return len(self.obs_names)

    @property
    @abstractmethod
    def obs_names(self):
        pass

    @property
    def n_obs(self) -> int:
        return len(self.var_names)

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
