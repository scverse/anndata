.. role:: small
.. role:: smaller
.. role:: noteversion


.. note::
   Upcoming changes:

   - :attr:`~anndata.AnnData.layers` and :attr:`~anndata.AnnData.X` will be unified.
   - :attr:`~anndata.AnnData.filename` and :attr:`~anndata.AnnData.isbacked` will be unified under a new name.
   - The types of :attr:`~anndata.AnnData.raw`, :attr:`~anndata.AnnData.layers`, :attr:`~anndata.AnnData.obsm`,
     :attr:`~anndata.AnnData.varm`, :attr:`~anndata.AnnData.obsp` and :attr:`~anndata.AnnData.varp` will be exported.
   - Square matrices in :attr:`~anndata.AnnData.uns` will no longer be sliced (use `.{obs,var}p` instead).


0.7 :small:`January 22, 2020`
-----------------------------

.. warning::
   Breaking changes introduced between 0.6.22.post1 and 0.7:

   - Elements of :class:`~anndata.AnnData`\ s don’t have their dimensionality reduced when the main object is subset.
     This is to maintain consistency when subsetting. See discussion in :issue:`145`.
   - Internal modules like `anndata.core` are private and their contents are not stable: See :issue:`174`.
   - The old deprecated attributes `.smp*`. `.add` and `.data` have been removed.

   Currently broken features

   - `sc.pp.normalize_per_cell` doesn’t work on dask arrays. It just doesn’t modify the matrix.


View overhaul – PR :pr:`164`
  - Indexing into a view no longer keeps a reference to intermediate view, see :issue:`62`.
  - Views are now lazy. Elements of view of AnnData are not indexed until they’re accessed.
  - Indexing with scalars no longer reduces dimensionality of contained arrays, see :issue:`145`.
  - All elements of AnnData should now follow the same rules about how they’re subset, see :issue:`145`.
  - Can now index by observations and variables at the same time.


IO overhaul – PR :pr:`167`
  - Reading and writing has been overhauled for simplification and speed.

    - Time and memory usage can be half of previous in typical use cases

  - Zarr backend now supports sparse arrays, and generally is closer to having the same features as HDF5.
  - Backed mode should see significant speed and memory improvements for access along compressed dimensions and IO. PR :pr:`241`.
  - :class:`~pandas.Categorical`\ s can now be ordered (PR :pr:`230`) and written to disk with a large number of categories (PR :pr:`217`).


Mapping attributes overhaul :smaller:`(obsm, varm, layers, …)`
  - New attributes :attr:`~anndata.AnnData.obsp` and :attr:`~anndata.AnnData.varp` have been added for two dimensional arrays where each axis corresponds to a single axis of the AnnData object. PR :pr:`207`.

    - These are intended to store values like cell-by-cell graphs, which are currently stored in :attr:`~anndata.AnnData.uns`.

  - Sparse arrays are now allowed as values in all mapping attributes.
  - DataFrames are now allowed as values in :attr:`~anndata.AnnData.obsm` and :attr:`~anndata.AnnData.varm`.
  - All mapping attributes now share an implementation and will have the same behaviour. PR :pr:`164`.


Miscellaneous improvements
  - Mapping attributes now have ipython tab completion (e.g. `adata.obsm["\\t` can provide suggestions) PR :pr:`183`.
  - :class:`~anndata.AnnData` attributes are now delete-able (e.g. `del adata.raw`) PR :pr:`242`.
  - Many many bug fixes


Versions 0.6.*
--------------

- better support for aligned mappings (obsm, varm, layers)
  :noteversion:`0.6.22` :pr:`155` :smaller:`thanks to I Virshup`
- convenience accesors :func:`~anndata.AnnData.obs_vector`, :func:`~anndata.AnnData.var_vector` for 1d arrays.
  :noteversion:`0.6.21` :pr:`144` :smaller:`thanks to I Virshup`
- compatibility with Scipy >=1.3 by removing `IndexMixin` dependency.
  :noteversion:`0.6.20` :pr:`151` :smaller:`thanks to P Angerer`
- bug fix for second-indexing into views.
  :noteversion:`0.6.19` :commit:`0ab553f368a93c52923f8cc700a066440824e8d8` :smaller:`thanks to P Angerer`
- bug fix for reading excel files.
  :noteversion:`0.6.19` :commit:`90bea2c1721d5dbfad20975b14809c63cc126ae8` :smaller:`thanks to A Wolf`
- changed default compression to `None` in :func:`~anndata.AnnData.write_h5ad` to speed up read and write, disk space use is usually less critical.
  :noteversion:`0.6.16` :commit:`21d8033dc560794b8eb8b58a693e30f4d154554e` :smaller:`thanks to A Wolf`
- maintain dtype upon copy.
  :noteversion:`0.6.13` :commit:`534bea4b04a542d33743050a63c8b7dbff8b4d9a` :smaller:`thanks to A Wolf`
- :attr:`~anndata.AnnData.layers` inspired by `.loom`_ files allows their information lossless reading via :func:`~anndata.read_loom`.
  :noteversion:`0.6.7`–:noteversion:`0.6.9` :pr:`46` & :pr:`48` :smaller:`thanks to S Rybakov`
- support for reading zarr files: :func:`~anndata.read_zarr`
  :noteversion:`0.6.7` :pr:`38` :smaller:`thanks to T White`
- initialization from pandas DataFrames
  :noteversion:`0.6.` :commit:`648bcc8a33f645de1e483bd6f9f5a3cb34ff43a3` :smaller:`thanks to A Wolf`
- iteration over chunks :func:`~anndata.AnnData.chunked_X` and :func:`~anndata.AnnData.chunk_X`
  :noteversion:`0.6.1` :pr:`20` :smaller:`thanks to S Rybakov`

Version 0.6 :small:`May 1, 2018`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- compatibility with Seurat converter
- tremendous speedup for :func:`~anndata.AnnData.concatenate`
- bug fix for deep copy of unstructured annotation after slicing
- bug fix for reading HDF5 stored single-category annotations
- “outer join” concatenation: adds zeros for concatenation of sparse data and nans for dense data
- better memory efficiency in loom exports


Version 0.5 :small:`February 9, 2018`
-------------------------------------

- inform about duplicates in :class:`~anndata.AnnData.var_names` and resolve them using :func:`~anndata.AnnData.var_names_make_unique`
- automatically remove unused categories after slicing
- read/write `.loom`_ files using loompy 2
- fixed read/write for a few text file formats
- read `UMI tools`_ files: :func:`~anndata.read_umi_tools`

.. _UMI tools: https://github.com/CGATOxford/UMI-tools


Version 0.4 :small:`December 23, 2017`
--------------------------------------

- read/write `.loom`_ files
- scalability beyond dataset sizes that fit into memory: see this `blog post`_
- :class:`~anndata.AnnData` has a :class:`~anndata.AnnData.raw` attribute, which simplifies storing the data matrix when you consider it *raw*: see the `clustering tutorial`_

.. _.loom: http://loompy.org
.. _blog post: http://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/
.. _clustering tutorial: https://github.com/theislab/scanpy_usage/tree/master/170505_seurat
