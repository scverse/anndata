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


Experimental API
----------------

.. warning::

   API's in the experimenal module are currently in development and subject to change at any time.

Two classes for working with batched access to collections of many `AnnData` objects or `h5ad` files. In paritcular, for pytorch-based models.

.. autosummary::
   :toctree: .

   experimental.AnnCollection
   experimental.AnnLoader

Errors and warnings
-------------------

.. autosummary::
   :toctree: .

   ImplicitModificationWarning
