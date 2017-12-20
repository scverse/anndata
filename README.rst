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

An ``AnnData`` object ``adata`` stores a data matrix (``adata.data``),
dataframe-like sample (``adata.smp``) and variable (``adata.var``) annotation
and unstructured dict-like annotation (``adata.uns``).

Read the `documentation <http://scanpy.readthedocs.io/en/latest/api/scanpy.api.AnnData.html>`_.

You can install ``anndata`` independently of Scanpy from `PyPI <https://pypi.python.org/pypi/anndata/>`__: ``pip install anndata``.

We will slowly support all of AnnData's generic functions (plotting,
preprocessing, reading, writing from disk, hdf5 backing on disk) via the anndata
package.

AnnData has been developed within Scanpy. We are grateful if you consider citing
our `preprint on Scanpy and AnnData <https://doi.org/10.1101/174029>`_.
