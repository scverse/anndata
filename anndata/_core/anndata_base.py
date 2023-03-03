from abc import abstractmethod
from typing import Tuple
from ..utils import DeprecationMixinMeta


class AnnDataBase(metaclass=DeprecationMixinMeta):
    @abstractmethod
    def _init_as_actual(self, *args, **kwargs):
        pass

    @abstractmethod
    def _init_as_view(self, *args, **kwargs):
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
