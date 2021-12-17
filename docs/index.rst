.. include:: ../README.rst
   :end-line: 22

.. image:: _static/img/anndata_schema.svg
   :align: right
   :width: 40%

anndata is a Python package for handling annotated data matrices in memory and on disk, positioned between pandas and xarray. anndata offers a broad range of computationally efficient features including, among others, sparse data support, lazy operations, and a PyTorch interface.

* Discuss development on `GitHub <https://github.com/theislab/anndata>`_.
* Install via `pip install anndata` or `conda install anndata -c conda-forge`.
* Consider reading the `anndata paper <https://github.com/ivirshup/anndata-paper>`__.
* See `Scanpy's documentation <https://scanpy.readthedocs.io/>`__ for usage
  related to single cell data. anndata was initially built for Scanpy.

News
----

.. include:: news.rst

Latest additions
----------------

.. include:: release-latest.rst

.. toctree::
   :maxdepth: 1
   :hidden:

   tutorials
   api
   concatenation
   fileformat-prose
   benchmarks
   release-notes
   references
