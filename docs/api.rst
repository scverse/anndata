API
===

.. module:: anndata

The central class:

.. autosummary::
   :toctree: generated/

   AnnData

Combining
---------

Combining AnnData objects. See also the section on concatenation.

.. autosummary::
   :toctree: generated/

   concat

Reading
-------

Reading anndata’s native file format `.h5ad`.

.. autosummary::
   :toctree: generated/

   read_h5ad

Reading other file formats.

.. autosummary::
   :toctree: generated/

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
   :toctree: generated/

   AnnData.write

Writing to other formats.

.. autosummary::
   :toctree: generated/

   AnnData.write_csvs
   AnnData.write_loom
   AnnData.write_zarr


Errors and warnings
-------------------

.. autosummary::
   :toctree: generated/

   ImplicitModificationWarning
