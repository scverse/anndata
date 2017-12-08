from .base import AnnData

__doc__ = """\
anndata
=======

The `anndata` package provides functions for the :class:`~anndata.AnnData` data
container, which extends pandas `DataFrames <https://pandas.pydata.org/pandas-docs/stable/dsintro.html#dataframe>`_: {main_narrative}

The AnnData class
-----------------

.. autosummary::
  :toctree: .

   AnnData

Reading
-------

.. autosummary::
  :toctree: .

   read_anndata
   read_csv
   read_excel
   read_hdf
   read_loom
   read_mtx
   read_text

Writing
-------

.. autosummary::
  :toctree: .

   AnnData.write_anndata
   AnnData.write_csvs
   AnnData.write_loom

References
----------

.. [Huber15]
   Huber *et al.* (2015),
   *Orchestrating high-throughput genomic analysis with Bioconductor*,
   `Nature Methods <https://doi.org/10.1038/nmeth.3252>`__.

.. [Wolf17]
   Wolf *et al.* (2017),
   *Scanpy for analysis of large-scale single-cell gene expression data*,
   `bioRxiv <https://doi.org/10.1101/174029>`__.

""".format(main_narrative=AnnData._main_narrative)

from .readwrite import read_anndata, read_loom, \
    read_csv, read_excel, read_text, read_hdf, read_mtx

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
