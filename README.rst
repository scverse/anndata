|Stars| |PyPI| |PyPIDownloadsTotal| |PyPIDownloadsMonth| |Conda| |Docs| |Build Status| |Coverage|

.. |Stars| image:: https://img.shields.io/github/stars/theislab/anndata?logo=GitHub&color=yellow
   :target: https://github.com/theislab/anndata/stargazers
.. |PyPI| image:: https://img.shields.io/pypi/v/anndata.svg
   :target: https://pypi.org/project/anndata
.. |PyPIDownloadsTotal| image:: https://pepy.tech/badge/anndata
   :target: https://pepy.tech/project/anndata
.. |PyPIDownloadsMonth| image:: https://img.shields.io/pypi/dm/scanpy?logo=PyPI&color=blue
   :target: https://pypi.org/project/anndata
.. |Conda| image:: https://img.shields.io/conda/vn/conda-forge/anndata.svg
   :target: https://anaconda.org/conda-forge/anndata
.. |Docs| image:: https://readthedocs.com/projects/icb-anndata/badge/?version=latest
   :target: https://anndata.readthedocs.io
.. |Build Status| image:: https://dev.azure.com/theislab/anndata/_apis/build/status/theislab.anndata?branchName=master
   :target: https://dev.azure.com/theislab/anndata/_build
.. |Coverage| image:: https://api.codacy.com/project/badge/Coverage/b92ae35b691141ceb5f2ee74beaf39d3
   :target: https://www.codacy.com/manual/theislab/anndata

anndata - Annotated Data
========================

.. image:: _static/img/anndata_schema.svg
   :align: right
   :width: 350px

anndata is a Python package for handling annotated data matrices in memory and on disk. It is positioned between pandas and xarray by providing structure that organizes data matrix annotations. anndata offers a broad range of computationally efficient features including, among others, sparse data support, lazy operations, and a PyTorch interface.

* Read the `documentation <https://anndata.readthedocs.io>`_.
* Install via ``pip install anndata`` or ``conda install anndata -c conda-forge``.

.. would be nice to have http://falexwolf.de/img/scanpy/anndata.svg also on GitHub, but it’s much too wide there;
.. GitHub doesn’t plan to resolve scaling images: https://github.com/github/markup/issues/295
