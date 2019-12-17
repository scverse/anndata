from warnings import warn

warn("Please only import from anndata, not anndata.readwrite", DeprecationWarning)

from ._io import *
