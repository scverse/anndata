
**Future plans:** We will slowly support all of AnnData's generic data analysis functions -
currently in Scanpy (plotting, preprocessing) - via the anndata package. All
code will remain backwards compatible but :class:`anndata` will grow in its own
right.

**December 29, 2017**: version 0.4.2

1. fixed text file support (``.csv``, …)
2. read csv files without ``\n`` in the last column name

**December 28, 2017**: version 0.4.1

1. read `UMI tools <https://github.com/CGATOxford/UMI-tools>`_ files: :func:`~anndata.read_umi_tools`
