from .base import AnnData
from .base import read_h5ad as read

__doc__ = """\
anndata
=======

The `anndata` package provides functions for :class:`~anndata.AnnData` objects and their native file format `.h5ad`: {main_narrative}

.. autosummary::
  :toctree: .

   AnnData


Reading
-------

Reading anndata's native file format `.h5ad`. Doing this, by default, only
relevant contents of the file are loaded into memory.

.. autosummary::
  :toctree: .

   read

Reading other file formats.

.. autosummary::
  :toctree: .

   read_csv
   read_excel
   read_hdf
   read_loom
   read_mtx
   read_text


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
""".format(main_narrative=AnnData._main_narrative)

from .readwrite import read_loom, \
    read_csv, read_excel, read_text, read_hdf, read_mtx

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
