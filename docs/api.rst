API
===

.. module:: anndata

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


Experimental development API
----------------------------

.. autosummary::
   :toctree: .

   dev.multi_files.AnnDataSet
   dev.multi_files.AnnDataLoader


Errors and warnings
-------------------

.. autosummary::
   :toctree: .

   ImplicitModificationWarning
