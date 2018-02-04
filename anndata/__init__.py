from .base import AnnData, _MAIN_NARRATIVE
from .readwrite import (
    read_h5ad, read_loom, read_hdf,
    read_excel, read_umi_tools,
    read_csv, read_text, read_mtx,
)
from .readwrite import read_h5ad as read  # backwards compat / shortcut for default format

__doc__ = """\
API
===

The central class:

.. autosummary::
  :toctree: .

   AnnData


Reading
-------

Reading anndata's native file format `.h5ad`.

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


Writing
-------

Writing to anndata's native file format `.h5ad`.

.. autosummary::
  :toctree: .

   AnnData.write

Writing to other formats.

.. autosummary::
  :toctree: .

   AnnData.write_csvs
   AnnData.write_loom

h5py
----

Independent of :class:`~anndata.AnnData`, the submodule :class:`anndata.h5py`
provides a thin wrapper of `h5py <http://www.h5py.org/>`_ that is able to handle
sparse matrices in addition to the standard functionality.

.. autosummary::
  :toctree: .

   h5py

""".format(main_narrative=_MAIN_NARRATIVE)

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
