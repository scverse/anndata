On-disk format
--------------

.. note::
   These docs are written for anndata 0.7.
   Files written before this version may differ in some conventions,
   but will still be read by newer versions of the library.

AnnData objects are saved on disk to hierarchichal array stores like HDF5_
(via :doc:`H5py <h5py:index>`) and :doc:`zarr:index`.
This allows us to have very similar structures in disk and on memory.

In general, `AnnData` objects can hold three kinds of values:
(1) dense arrays, (2) sparse arrays, and (3) DataFrames.
As an example we’ll look into a typical `.h5ad` object that’s been through an analysis.
This structure should be largely equivalent to Zarr structure, though there are a few minor differences.

.. _HDF5: https://en.wikipedia.org/wiki/Hierarchical_Data_Format
.. I’ve started using h5py since I couldn’t figure out a nice way to print attributes from bash.

>>> import h5py
>>> f = h5py.File("02_processed.h5ad", "r")
>>> list(f.keys())
['X', 'layers', 'obs', 'obsm', 'uns', 'var', 'varm']

.. .. code:: bash

..    $ h5ls 02_processed.h5ad
..    X                        Group
..    layers                   Group
..    obs                      Group
..    obsm                     Group
..    uns                      Group
..    var                      Group
..    varm                     Group

Dense arrays
~~~~~~~~~~~~

Dense arrays have the most simple representation on disk,
as they have native equivalents in H5py :doc:`h5py:high/dataset` and Zarr :ref:`Arrays <zarr:tutorial_create>`.
We can see an example of this with dimensionality reductions stored in the `obsm` group:

>>> f["obsm"].visititems(print)
X_pca <HDF5 dataset "X_pca": shape (38410, 50), type "<f4">
X_umap <HDF5 dataset "X_umap": shape (38410, 2), type "<f4">

.. .. code:: bash

..    $ h5ls 02_processed.h5ad/obsm
..    X_pca                    Dataset {38410, 50}
..    X_umap                   Dataset {38410, 2}

Sparse arrays
~~~~~~~~~~~~~

Sparse arrays don’t have a native representations in HDF5 or Zarr,
so we've defined our own based on their in-memory structure.
Currently two sparse data formats are supported by `AnnData` objects, CSC and CSR
(corresponding to :class:`scipy.sparse.csc_matrix` and :class:`scipy.sparse.csr_matrix` respectivley).
These formats represent a two-dimensional sparse array with
three one-dimensional arrays, `indptr`, `indices`, and `data`.

.. note::
   A full description of these formats is out of scope for this document,
   but are `easy to find`_.

.. _easy to find: https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)

We represent a sparse array as a `Group` on-disk,
where the kind and shape of the sparse array is defined in the `Group`'s attributes:

>>> dict(f["X"].attrs)
{'encoding-type': 'csr_matrix',
 'encoding-version': '0.1.0',
 'shape': array([38410, 27899])}

Inside the group are the three constituent arrays:

>>> f["X"].visititems(print)
data <HDF5 dataset "data": shape (41459314,), type "<f4">
indices <HDF5 dataset "indices": shape (41459314,), type "<i4">
indptr <HDF5 dataset "indptr": shape (38411,), type "<i4">

.. .. code:: bash

..    $ h5ls 02_processed.h5ad/X
..    data                     Dataset {41459314/Inf}
..    indices                  Dataset {41459314/Inf}
..    indptr                   Dataset {38411/Inf}

DataFrames
~~~~~~~~~~

DataFrames are saved as a columnar format in a group, so each column of a DataFrame gets its own dataset.
To maintain efficiency with categorical values, only the numeric codes are stored for each row,
while categories values are saved in a reserved subgroup `__categories`.

Dataframes can be identified from other groups by their attributes:

>>> dict(f["obs"].attrs)
{'_index': 'Cell',
 'column-order': array(['sample', 'cell_type', 'n_genes_by_counts',
       'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts',
       'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes',
       'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes',
       'total_counts_mito', 'log1p_total_counts_mito', 'pct_counts_mito',
       'label_by_score'], dtype=object),
 'encoding-type': 'dataframe',
 'encoding-version': '0.1.0'}

These attributes identify the column used as an index,
the order of the original columns, and some type information.

>>> f["obs"].visititems(print)
Cell <HDF5 dataset "Cell": shape (38410,), type "|O">
__categories <HDF5 group "/obs/__categories" (3 members)>
__categories/cell_type <HDF5 dataset "cell_type": shape (22,), type "|O">
__categories/label_by_score <HDF5 dataset "label_by_score": shape (16,), type "|O">
__categories/sample <HDF5 dataset "sample": shape (41,), type "|O">
cell_type <HDF5 dataset "cell_type": shape (38410,), type "|i1">
label_by_score <HDF5 dataset "label_by_score": shape (38410,), type "|i1">
log1p_n_genes_by_counts <HDF5 dataset "log1p_n_genes_by_counts": shape (38410,), type "<f8">
[...]

Categorical Series can be identified by the presence of the attribute `"categories"`,
which contains a pointer to the categories' values:

>>> dict(f["obs/cell_type"].attrs)
{'categories': <HDF5 object reference>}

.. note::
   In `zarr`, as there are no reference objects, the `categories` attribute
   is an absolute path to the category values.

Other values
~~~~~~~~~~~~

Mappings
^^^^^^^^

Mappings are simply stored as `Group` s on disk.
These are distinct from DataFrames and sparse arrays since they don’t have any special attributes.
A `Group` is created for any `Mapping` in the AnnData object,
including the standard `obsm`, `varm`, `layers`, and `uns`.
Notably, this definition is used recursively within `uns`:

>>> f["uns"].visititems(print)
[...]
pca <HDF5 group "/uns/pca" (2 members)>
pca/variance <HDF5 dataset "variance": shape (50,), type "<f4">
pca/variance_ratio <HDF5 dataset "variance_ratio": shape (50,), type "<f4">
[...]

Scalars
^^^^^^^

Zero dimensional arrays are used for scalar values (i.e. single values like strings, numbers or booleans).
These should only occur inside of `uns`, and are common inside of saved parameters:

>>> f["uns/neighbors/params"].visititems(print)
method <HDF5 dataset "method": shape (), type "|O">
metric <HDF5 dataset "metric": shape (), type "|O">
n_neighbors <HDF5 dataset "n_neighbors": shape (), type "<i8">
>>> f["uns/neighbors/params/metric"][()]
'euclidean'
