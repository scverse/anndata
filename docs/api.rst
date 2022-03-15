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

.. _experimental_api:

Experimental API
----------------

.. warning::

   API's in the experimenal module are currently in development and subject to change at any time.

Two classes for working with batched access to collections of many `AnnData` objects or `h5ad` files. In paritcular, for pytorch-based models.

.. autosummary::
   :toctree: generated/

   experimental.AnnCollection
   experimental.AnnLoader

Low level methods for reading and writing elements of an `AnnData`` object to a store:


.. autosummary::
   :toctree: generated/

   experimental.read_elem
   experimental.write_elem


Errors and warnings
-------------------

.. autosummary::
   :toctree: generated/

   ImplicitModificationWarning
