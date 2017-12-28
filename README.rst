|Docs| |PyPI| |Build Status| |Coverage|

.. |Docs| image:: https://readthedocs.org/projects/scanpy/badge/?version=latest
   :target: https://scanpy.readthedocs.io
.. |PyPI| image:: https://badge.fury.io/py/anndata.svg
   :target: https://pypi.python.org/pypi/anndata
.. |Build Status| image:: https://travis-ci.org/theislab/anndata.svg?branch=master
   :target: https://travis-ci.org/theislab/anndata
.. |Coverage| image:: https://codecov.io/gh/theislab/anndata/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/theislab/anndata

AnnData
=======

At the very basic level: an ``AnnData`` object ``adata`` stores a data matrix (``adata.X``),
dataframe-like annotation of observations (``adata.obs``) and variables (``adata.var``)
and unstructured dict-like annotation (``adata.uns``).

Read the `documentation <http://anndata.readthedocs.org>`_.

Install from `PyPI <https://pypi.python.org/pypi/anndata/>`__ via ``pip install anndata``.

We will slowly support all of AnnData's generic functions (plotting,
preprocessing, reading, writing from disk, hdf5 backing on disk) via the anndata
package.

AnnData has been developed together with Scanpy. We are grateful if you consider citing
our `preprint, soon in Genome Biology <https://doi.org/10.1101/174029>`_.
