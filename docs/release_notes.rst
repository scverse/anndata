.. role:: small
.. role:: smaller
.. role:: noteversion


On Master :small:`July 29, 2019`
--------------------------------

.. warning::

    Breaking changes on master include: 

    - Elements of `AnnData` objects don't have their dimensionality reduced when the main object is subset. This is to maintain consistency when subsetting. See discussion `here <https://github.com/theislab/anndata/issues/145>`__.
    - While the `anndata.core` module currently exists, it may be renamed to `anndata._core` `#174 <https://github.com/theislab/anndata/issues/174>`__.

    Currently broken features

    - `sc.pp.normalize_per_cell` doesn't work on dask arrays. It just doesn't modify the matrix.


- Views have been overhauled  `PR #164 <https://github.com/theislab/anndata/pull/164>`__.

   - Indexing into a view no longer keeps a reference to intermediate view, `issue #62 <https://github.com/theislab/anndata/issues/62>`__.
   - Views are now lazy. Elements of view of AnnData are not indexed until they're accessed.
   - Indexing with scalars no longer reduces dimensionality of contained arrays, `issue #145 <https://github.com/theislab/anndata/issues/145>`__.
   - All elements of AnnData should now follow the same rules about how they're subset, `issue #145 <https://github.com/theislab/anndata/issues/145>`__.
   - Can now index by observations and variables at the same time.

Post v0.6 :small:`June 6, 2019`
---------------------------------

- convenience accesors :func:`~anndata.AnnData.obs_vector`, :func:`~anndata.AnnData.var_vector` for 1d arrays, see `here <https://github.com/theislab/anndata/pull/144>`__ :noteversion:`0.6.21` :smaller:`thanks to I Virshup`
- compatibility with Scipy >=1.3 by removing `IndexMixin` dependency, see `here <https://github.com/theislab/anndata/commit/6fb083477bc0b1f3eeccc62e10e4b477ae532346>`__ :noteversion:`0.6.20` :smaller:`thanks to P Angerer`
- bug fix for second-indexing into views, see `here <https://github.com/theislab/anndata/issues/126>`__ :noteversion:`0.6.19` :smaller:`thanks to P Angerer`
- bug fix for reading excel files :noteversion:`0.6.19` :smaller:`thanks to A Wolf`
- :attr:`~anndata.AnnData.layers` inspired by `.loom <http://loompy.org>`__ files allows their information lossless reading via :func:`~anndata.read_loom` :smaller:`thanks to S Rybakov`
- initialization from pandas DataFrames :smaller:`thanks to A Wolf`
- iteration over chunks :func:`~anndata.AnnData.chunked_X` and :func:`~anndata.AnnData.chunk_X` :smaller:`thanks to S Rybakov`
- support for reading zarr files: :func:`~anndata.read_zarr` :smaller:`thanks to T White`
- changed default compression to `None` in :func:`~anndata.AnnData.write_h5ad` to speed up read and write, disk space use is usually less critical :noteversion:`0.6.16`
- maintain dtype upon copy :noteversion:`0.6.13` :smaller:`thanks to A Wolf`


Version 0.6 :small:`May 1, 2018`
--------------------------------

- compatibility with Seurat converter
- tremendous speedup for :func:`~anndata.AnnData.concatenate`
- bug fix for deep copy of unstructured annotation after slicing
- bug fix for reading HDF5 stored single-category annotations
- 'outer join' concatenation: adds zeros for concatenation of sparse data and nans for dense data
- better memory efficiency in loom exports


Version 0.5 :small:`February 9, 2018`
-------------------------------------

- inform about duplicates in :class:`~anndata.AnnData.var_names` and resolve them using :func:`~anndata.AnnData.var_names_make_unique`
- automatically remove unused categories after slicing
- read/write ``.loom`` files using loompy 2
- fixed read/write for a few text file formats
- read `UMI tools <https://github.com/CGATOxford/UMI-tools>`__ files: :func:`~anndata.read_umi_tools`


Version 0.4 :small:`December 23, 2017`
--------------------------------------

- read/write `.loom <http://loompy.org>`__ files
- scalability beyond dataset sizes that fit into memory: see this `blog post <http://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/>`__
- :class:`~anndata.AnnData` has a :class:`~anndata.AnnData.raw` attribute, which simplifies storing the data matrix when you consider it *raw*: see the `clustering tutorial <https://github.com/theislab/scanpy_usage/tree/master/170505_seurat>`__
