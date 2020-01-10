from ._core.anndata import AnnData
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

# backwards compat / shortcut for default format
from ._io import read_h5ad as read

__doc__ = """\
API
===

The central class:

.. autosummary::
   :toctree: .

   AnnData


Reading
-------

Reading anndata’s native file format `.h5ad`.

.. autosummary::
   :toctree: .

   read_h5ad

Reading other file formats.

.. autosummary::
   :toctree: .

   read_csv
   read_excel
   read_hdf
   read_loom
   read_mtx
   read_text
   read_umi_tools
   read_zarr


Writing
-------

Writing to anndata’s native file format `.h5ad`.

.. autosummary::
   :toctree: .

   AnnData.write

Writing to other formats.

.. autosummary::
   :toctree: .

   AnnData.write_csvs
   AnnData.write_loom
   AnnData.write_zarr

"""

__author__ = ", ".join(
    ["Philipp Angerer*", "Alex Wolf*", "Isaac Virshup", "Sergei Rybakov"]
)
__email__ = ", ".join(
    [
        "philipp.angerer@helmholtz-muenchen.de",
        "f.alex.wolf@gmx.de",
        # We don’t need all, the main authors are sufficient.
    ]
)

try:
    from setuptools_scm import get_version

    __version__ = get_version(root="..", relative_to=__file__)
    del get_version
except (LookupError, ImportError):
    from .compat import pkg_version

    __version__ = pkg_version(__name__)
    del pkg_version
