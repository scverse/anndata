See all releases `here <https://github.com/theislab/anndata/releases>`_. The following lists selected improvements.


**December 16, 2018**: on GitHub and 0.6.16

1. :func:`~anndata.AnnData.layers` inspired by `.loom <http://loompy.org>`__ files allows their information lossless reading via :func:`~anndata.read_loom`
2. initialatization from pandas DataFrames
3. iteration over chunks :func:`~anndata.AnnData.chunked_X` and :func:`~anndata.AnnData.chunk_X`
4. support for reading zarr files: :func:`~anndata.read_zarr`
5. changed default compression to `None` in :func:`~anndata.AnnData.write_h5ad` to speed up read and write, disk space use is usually less critical
      

**May 1, 2018**: version 0.6

1. compatibility with Seurat converter
2. tremendous speedup for :func:`~anndata.AnnData.concatenate`
3. bug fix for deep copy of unstructured annotation after slicing
4. bug fix for reading HDF5 stored single-category annotations
5. 'outer join' concatenation: adds zeros for concatenation of sparse data and nans for dense data
6. better memory efficiency in loom exports


**February 9, 2018**: version 0.5

1. inform about duplicates in :class:`~anndata.AnnData.var_names` and resolve them using :func:`~anndata.AnnData.var_names_make_unique`
2. automatically remove unused categories after slicing
3. read/write `.loom` files using loompy 2
4. fixed read/write for a few text file formats
5. read `UMI tools <https://github.com/CGATOxford/UMI-tools>`__ files: :func:`~anndata.read_umi_tools`


**December 23, 2017**: version 0.4

1. read/write `.loom <http://loompy.org>`__ files
2. scalability beyond dataset sizes that fit into memory: see this
   `blog post
   <http://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/>`__
3. :class:`~anndata.AnnData` has a :class:`~anndata.AnnData.raw` attribute
   that simplifies storing the data matrix when you consider it "raw": see the
   `clustering tutorial
   <https://github.com/theislab/scanpy_usage/tree/master/170505_seurat>`__
