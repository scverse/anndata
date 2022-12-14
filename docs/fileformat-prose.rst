On-disk format
--------------

.. note::
   These docs are written for anndata 0.8.
   Files written before this version may differ in some conventions,
   but will still be read by newer versions of the library.

AnnData objects are saved on disk to hierarchical array stores like HDF5_
(via :doc:`H5py <h5py:index>`) and :doc:`zarr:index`.
This allows us to have very similar structures in disk and on memory.

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

In general, `AnnData` objects are comprised of a various types of elements.
Each element is encoded as either an Array (or Dataset in hdf5 terminology) or a collection of elements (e.g. Group) in the store.
We record the type of an element using the `encoding-type` and `encoding-version` keys in it's attributes.
For example, we can this file represents an `AnnData` object from this metadata:

>>> dict(f.attrs)
{'encoding-type': 'anndata', 'encoding-version': '0.1.0'}

Using this information, we're able to dispatch onto readers for the different element types that you'd find in an anndata.

Dense arrays
~~~~~~~~~~~~

Dense numeric arrays have the most simple representation on disk,
as they have native equivalents in H5py :doc:`h5py:high/dataset` and Zarr :ref:`Arrays <zarr:tutorial_create>`.
We can see an example of this with dimensionality reductions stored in the `obsm` group:

>>> f["obsm"].visititems(print)
X_pca <HDF5 dataset "X_pca": shape (38410, 50), type "<f4">
X_umap <HDF5 dataset "X_umap": shape (38410, 2), type "<f4">

>>> dict(f["obsm"]["X_pca"].attrs)
{'encoding-type': 'array', 'encoding-version': '0.2.0'}

.. .. code:: bash

..    $ h5ls 02_processed.h5ad/obsm
..    X_pca                    Dataset {38410, 50}
..    X_umap                   Dataset {38410, 2}

Sparse arrays
~~~~~~~~~~~~~

Sparse arrays don’t have a native representations in HDF5 or Zarr,
so we've defined our own based on their in-memory structure.
Currently two sparse data formats are supported by `AnnData` objects, CSC and CSR
(corresponding to :class:`scipy.sparse.csc_matrix` and :class:`scipy.sparse.csr_matrix` respectively).
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

AwkwardArrays
~~~~~~~~~~~~~

Ragged arrays are supported in ``anndata`` through the `Awkward
Array <https://awkward-array.org/>`__ library. For storage on disk, we
break down the awkward array into it’s constituent arrays using
```ak.to_buffers`` <https://awkward-array.readthedocs.io/en/latest/_auto/ak.to_buffers.html>`__,
then writing these arrays using ``anndata``\ ’s methods.

The container of arrays is stored in a group called ``"container"``

.. code:: python

   >>> import zarr
   >>> z = zarr.open("airr.zarr", "r")
   >>> awkward_group = z["obsm/airr"]
   >>> awkward_group.tree()

::

   airr
    └── container
        ├── node0-offsets (17,) int64
        ├── node2-offsets (40,) int64
        ├── node3-data (117,) uint8
        ├── node4-offsets (40,) int64
        └── node5-data (117,) uint8

The length of the array is saved to it’s own ``"length"`` attribute,
while metadata for the array structure is serialized and saved to the
“form” attribute.

.. code:: python

   >>> dict(awkward_group.attrs)

.. code:: json

   {'encoding-type': 'awkward-array',
    'encoding-version': '0.1.0',
    'form': '{"class": "ListOffsetArray", "offsets": "i64", "content": {"class": '
            '"RecordArray", "contents": {"locus": {"class": "ListOffsetArray", '
            '"offsets": "i64", "content": {"class": "NumpyArray", "primitive": '
            '"uint8", "inner_shape": [], "has_identifier": false, "parameters": '
            '{"__array__": "char"}, "form_key": "node3"}, "has_identifier": '
            'false, "parameters": {"__array__": "string"}, "form_key": "node2"}, '
            '"junction_aa": {"class": "ListOffsetArray", "offsets": "i64", '
            '"content": {"class": "NumpyArray", "primitive": "uint8", '
            '"inner_shape": [], "has_identifier": false, "parameters": '
            '{"__array__": "char"}, "form_key": "node5"}, "has_identifier": '
            'false, "parameters": {"__array__": "string"}, "form_key": "node4"}}, '
            '"has_identifier": false, "parameters": {}, "form_key": "node1"}, '
            '"has_identifier": false, "parameters": {}, "form_key": "node0"}',
    'length': 16}

These can be read back as awkward arrays using the
```ak.from_buffers`` <https://awkward-array.readthedocs.io/en/latest/_auto/ak.from_buffers.html>`__
function:

.. code:: python

   >>> import awkward._v2 as ak
   >>> from anndata.experimental import read_elem
   >>> ak.from_buffers(
   ...     awkward_group.attrs["form"],
   ...     awkward_group.attrs["length"],
   ...     {k: read_elem(v) for k, v in awkward_group["container"].items()}
   ... )

::

   <Array [[], [...], ..., [{locus: 'TRD', ...}]] type='16 * var * {locus: str...'>

DataFrames
~~~~~~~~~~

DataFrames are saved as a columnar format in a group, so each column of a DataFrame is saved as a separate array.
We save a little more information in the attributes here.

>>> dict(f["obs"].attrs)
{'_index': 'Cell',
 'column-order': array(['sample', 'cell_type', 'n_genes_by_counts',
        'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts',
        'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes',
        'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes',
        'total_counts_mito', 'log1p_total_counts_mito', 'pct_counts_mito',
        'label_by_score'], dtype=object),
 'encoding-type': 'dataframe',
 'encoding-version': '0.2.0'}

These attributes identify the index of the dataframe, as well as the original order of the columns.
Each column in this dataframe is encoded as it's own array.

>>> dict(f["obs"]["total_counts"].attrs)
{'encoding-type': 'array', 'encoding-version': '0.2.0'}

>>> dict(f["obs"]["cell_type"].attrs)
{'encoding-type': 'categorical', 'encoding-version': '0.2.0', 'ordered': False}

Mappings
~~~~~~~~

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
~~~~~~~

Zero dimensional arrays are used for scalar values (i.e. single values like strings, numbers or booleans).
These should only occur inside of `uns`, and are common inside of saved parameters:

>>> f["uns/neighbors/params"].visititems(print)
method <HDF5 dataset "method": shape (), type "|O">
metric <HDF5 dataset "metric": shape (), type "|O">
n_neighbors <HDF5 dataset "n_neighbors": shape (), type "<i8">
>>> f["uns/neighbors/params/metric"][()]
'euclidean'
>>> dict(f["uns/neighbors/params/metric"].attrs)
{'encoding-type': 'string', 'encoding-version': '0.2.0'}

Categorical arrays
~~~~~~~~~~~~~~~~~~

>>> categorical = f["obs"]["cell_type"]
>>> dict(categorical.attrs)
{'encoding-type': 'categorical', 'encoding-version': '0.2.0', 'ordered': False}

Discrete labels can be efficiently represented with categorical arrays (similar to `factors` in `R`).
These arrays encode the labels as small width integers (`codes`), which map to the original label set (`categories`).
We store these two arrays separately.

>>> categorical.visititems(print)
categories <HDF5 dataset "categories": shape (22,), type "|O">
codes <HDF5 dataset "codes": shape (38410,), type "|i1">

String arrays
~~~~~~~~~~~~~

Arrays of strings are handled differently than numeric arrays since numpy doesn't really have a good way of representing arrays of unicode strings.
`anndata` assumes strings are text like data, so are variable length.

>>> dict(categorical["categories"].attrs)
{'encoding-type': 'string-array', 'encoding-version': '0.2.0'}

Nullable integers and booleans
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We support IO with Pandas nullable integer and boolean arrays.
We represent these on disk similar to `numpy` masked arrays, `julia` nullable arrays, or `arrow` validity bitmaps (see :issue:`504` for more discussion).
That is, we store a indicator array (or mask) of null values alongside the array of all values.

>>> h5_file = h5py.File("anndata_format.h5", "a")
>>> int_array = pd.array([1, None, 3, 4])
>>> int_array
<IntegerArray>
[1, <NA>, 3, 4]
Length: 4, dtype: Int64
>>> write_elem(h5_file, "nullable_integer", int_array)

>>> h5_file["nullable_integer"].visititems(print)
mask <HDF5 dataset "mask": shape (4,), type "|b1">
values <HDF5 dataset "values": shape (4,), type "<i8">

>>> dict(h5_file["nullable_integer"].attrs)
{'encoding-type': 'nullable-integer', 'encoding-version': '0.1.0'}
