"""\
Annotated Data
==============

The central class:

.. autosummary::
   :toctree: .

   AnnData

Combining
---------

Combining AnnData objects. See also the section on concatenation.

.. autosummary::
   :toctree: .

   concat

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


Errors and warnings
-------------------

.. autosummary::
   :toctree: .

   ImplicitModificationWarning

"""

from ._metadata import __version__, __author__, __email__, within_flit

if not within_flit():
    del within_flit
    from ._core.anndata import AnnData, ImplicitModificationWarning
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

    # backwards compat / shortcut for default format
    from ._io import read_h5ad as read
