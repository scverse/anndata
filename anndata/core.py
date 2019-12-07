from warnings import warn

warn("Please only import from anndata, not anndata.core", DeprecationWarning)

from ._core import *
