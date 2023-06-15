"""Annotated multivariate observation data."""

try:  # See https://github.com/maresb/hatch-vcs-footgun-example
    from setuptools_scm import get_version

    __version__ = get_version(root="..", relative_to=__file__)
except (ImportError, LookupError):
    try:
        from ._version import __version__
    except ModuleNotFoundError:
        raise RuntimeError(
            "anndata is not correctly installed. Please install it, e.g. with pip."
        )

from ._core.anndata import AnnData
from ._core.merge import concat
from ._core.raw import Raw
from ._io import (
    read_h5ad,
    read_loom,
    read_hdf,
    read_excel,
    read_umi_tools,
    read_csv,
    read_text,
    read_mtx,
    read_zarr,
)
from ._warnings import (
    OldFormatWarning,
    WriteWarning,
    ImplicitModificationWarning,
    ExperimentalFeatureWarning,
)

# backwards compat / shortcut for default format
from ._io import read_h5ad as read
