# On-disk format

```{note}
These docs are written for anndata 0.8+.
Files written before this version may differ in some conventions,
but will still be read by newer versions of the library.
```

AnnData objects are saved on disk to hierarchical array stores like [HDF5]
(via {doc}`H5py <h5py:index>`) and {doc}`zarr:index`.
This allows us to have very similar structures in disk and on memory.

As an example we’ll look into a typical `.h5ad`/ `.zarr` object that’s been through an analysis.
The structures are largely equivalent, though there are a few minor differences when it comes to type encoding.

## Elements

 <!-- I’ve started using h5py since I couldn’t figure out a nice way to print attributes from bash. -->


`````{tab-set}

````{tab-item} HDF5
:sync: hdf5

```python
>>> import h5py
>>> store = h5py.File("for-ondisk-docs/cart-164k-processed.h5ad", mode="r")
>>> list(store.keys())
['X', 'layers', 'obs', 'obsm', 'obsp', 'uns', 'var', 'varm', 'varp']
```

````

````{tab-item} Zarr
:sync: zarr

```python
>>> import zarr
>>> store = zarr.open("for-ondisk-docs/cart-164k-processed.zarr", mode="r")
>>> list(store.keys())
['X', 'layers', 'obs', 'obsm', 'obsp', 'uns', 'var', 'varm', 'varp']
```

````

`````

<!-- ```bash
$ h5ls 02_processed.h5ad
X                        Group
layers                   Group
obs                      Group
obsm                     Group
uns                      Group
var                      Group
varm                     Group
``` -->

In general, `AnnData` objects are comprised of various types of elements.
Each element is encoded as either an Array (or Dataset in hdf5 terminology) or a collection of elements (e.g. Group) in the store.
We record the type of an element using the `encoding-type` and `encoding-version` keys in its attributes.
For example, we can see that this file represents an `AnnData` object from its metadata:

```python
>>> dict(store.attrs)
{'encoding-type': 'anndata', 'encoding-version': '0.1.0'}
```

Using this information, we're able to dispatch onto readers for the different element types that you'd find in an anndata.

### Element Specification

* An element can be any object within the storage hierarchy (typically an array or group) with associated metadata
* An element MUST have a string-valued field `"encoding-type"` in its metadata
* An element MUST have a string-valued field `"encoding-version"` in its metadata that can be evaluated to a version

### AnnData specification (v0.1.0)

* An `AnnData` object MUST be a group.
* The group's metadata MUST include entries: `"encoding-type": "anndata"`, `"encoding-version": "0.1.0"`.
* An `AnnData` group MUST contain entries `"obs"` and `"var"`, which MUST be dataframes (though this may only have an index with no columns).
* The group MAY contain an entry `X`, which MUST be either a dense or sparse array and whose shape MUST be (`n_obs`, `n_var`)
* The group MAY contain a mapping `layers`. Entries in `layers` MUST be dense or sparse arrays which have shapes (`n_obs`, `n_var`)
* The group MAY contain a mapping `obsm`. Entries in `obsm` MUST be sparse arrays, dense arrays, or dataframes. These entries MUST have a first dimension of size `n_obs`
* The group MAY contain a mapping `varm`. Entries in `varm` MUST be sparse arrays, dense arrays, or dataframes. These entries MUST have a first dimension of size `n_var`
* The group MAY contain a mapping `obsp`. Entries in `obsp` MUST be sparse or dense arrays. The entries first two dimensions MUST be of size `n_obs`
* The group MAY contain a mapping `varp`. Entries in `varp` MUST be sparse or dense arrays. The entries first two dimensions MUST be of size `n_var`
* The group MAY contain a mapping `uns`. Entries in `uns` MUST be an anndata encoded type.

## Dense arrays

Dense numeric arrays have the most simple representation on disk,
as they have native equivalents in H5py {doc}`h5py:high/dataset` and Zarr {doc}`Arrays <zarr:user-guide/arrays>`.
We can see an example of this with dimensionality reductions stored in the `obsm` group:

`````{tab-set}

````{tab-item} HDF5
:sync: hdf5

```python
>>> store["obsm/X_pca"]
<HDF5 dataset "X_pca": shape (164114, 50), type "<f4">
```

````

````{tab-item} Zarr
:sync: zarr

```python
>>> store["obsm/X_pca"]
<zarr.core.Array '/obsm/X_pca' (164114, 50) float32 read-only>
```

````

`````

```python
>>> dict(store["obsm"]["X_pca"].attrs)
{'encoding-type': 'array', 'encoding-version': '0.2.0'}
```

<!-- ```bash
$ h5ls 02_processed.h5ad/obsm
X_pca                    Dataset {38410, 50}
X_umap                   Dataset {38410, 2}
``` -->

### Dense arrays specification (v0.2.0)

* Dense arrays MUST be stored in an Array object
* Dense arrays MUST have the entries `'encoding-type': 'array'` and `'encoding-version': '0.2.0'` in their metadata

## Sparse arrays

Sparse arrays don’t have a native representations in HDF5 or Zarr,
so we've defined our own based on their in-memory structure.
Currently two sparse data formats are supported by `AnnData` objects, CSC and CSR
(corresponding to {class}`scipy.sparse.csc_matrix` and {class}`scipy.sparse.csr_matrix` respectively).
These formats represent a two-dimensional sparse array with
three one-dimensional arrays, `indptr`, `indices`, and `data`.

```{note}
A full description of these formats is out of scope for this document,
but are [easy to find].
```

We represent a sparse array as a `Group` on-disk,
where the kind and shape of the sparse array is defined in the `Group`'s attributes:

```python
>>> dict(store["X"].attrs)
{'encoding-type': 'csr_matrix',
 'encoding-version': '0.1.0',
 'shape': [164114, 40145]}
```

The group contains three arrays:

`````{tab-set}

````{tab-item} HDF5
:sync: hdf5

```python
>>> store["X"].visititems(print)
data <HDF5 dataset "data": shape (495079432,), type "<f4">
indices <HDF5 dataset "indices": shape (495079432,), type "<i4">
indptr <HDF5 dataset "indptr": shape (164115,), type "<i4">
```

````

````{tab-item} Zarr
:sync: zarr

```python
>>> store["X"].visititems(print)
data <zarr.core.Array '/X/data' (495079432,) float32 read-only>
indices <zarr.core.Array '/X/indices' (495079432,) int32 read-only>
indptr <zarr.core.Array '/X/indptr' (164115,) int32 read-only>
```

````

`````

### Sparse array specification (v0.1.0)

* Each sparse array MUST be its own group
* The group MUST contain arrays `indices`, `indptr`, and `data`
* The group's metadata MUST contain:
    * `"encoding-type"`, which is set to `"csr_matrix"` or `"csc_matrix"` for compressed sparse row and compressed sparse column, respectively.
    * `"encoding-version"`, which is set to `"0.1.0"`
    * `"shape"` which is an integer array of length 2 whose values are the sizes of the array's dimensions

## DataFrames

DataFrames are saved as a columnar format in a group, so each column of a DataFrame is saved as a separate array.
We save a little more information in the attributes here.

```python
>>> dict(store["var"].attrs)
{'_index': 'ensembl_id',
 'column-order': ['highly_variable',
  'means',
  'variances',
  'variances_norm',
  'feature_is_filtered',
  'feature_name',
  'feature_reference',
  'feature_biotype',
  'mito'],
 'encoding-type': 'dataframe',
 'encoding-version': '0.2.0'}
```

These attributes identify the index of the dataframe, as well as the original order of the columns.
Each column in this dataframe is encoded as its own array.

`````{tab-set}

````{tab-item} HDF5
:sync: hdf5

```python
>>> store["var"].visititems(print)
ensembl_id <HDF5 dataset "ensembl_id": shape (40145,), type "|O">
feature_biotype <HDF5 group "/var/feature_biotype" (2 members)>
feature_biotype/categories <HDF5 dataset "categories": shape (1,), type "|O">
feature_biotype/codes <HDF5 dataset "codes": shape (40145,), type "|i1">
feature_is_filtered <HDF5 dataset "feature_is_filtered": shape (40145,), type "|b1">
...
```

````

````{tab-item} Zarr
:sync: zarr

```python
>>> store["var"].visititems(print)
ensembl_id <zarr.core.Array '/var/ensembl_id' (40145,) object read-only>
feature_biotype <zarr.hierarchy.Group '/var/feature_biotype' read-only>
feature_biotype/categories <zarr.core.Array '/var/feature_biotype/categories' (1,) object read-only>
feature_biotype/codes <zarr.core.Array '/var/feature_biotype/codes' (40145,) int8 read-only>
feature_is_filtered <zarr.core.Array '/var/feature_is_filtered' (40145,) bool read-only>
...
```

````

`````

```python
>>> dict(store["var"]["feature_name"].attrs)
{'encoding-type': 'categorical', 'encoding-version': '0.2.0', 'ordered': False}

>>> dict(store["var"]["feature_is_filtered"].attrs)
{'encoding-type': 'array', 'encoding-version': '0.2.0'}
```

### Dataframe Specification (v0.2.0)

* A dataframe MUST be stored as a group
* The group's metadata:
    * MUST contain the field `"_index"`, whose value is the key of the array to be used as an index/ row labels
    * MUST contain encoding metadata `"encoding-type": "dataframe"`, `"encoding-version": "0.2.0"`
    * MUST contain `"column-order"` an array of strings denoting the order of column entries
* The group MUST contain an array for the index
* Each entry in the group MUST correspond to an array with equivalent first dimensions
* Each entry SHOULD share chunk sizes (in the HDF5 or zarr container)

## Mappings

Mappings are simply stored as `Group`s on disk.
These are distinct from DataFrames and sparse arrays since they don’t have any special attributes.
A `Group` is created for any `Mapping` in the AnnData object,
including the standard `obsm`, `varm`, `layers`, and `uns`.
Notably, this definition is used recursively within `uns`:

`````{tab-set}

````{tab-item} HDF5
:sync: hdf5

```python
>>> store["uns"].visititems(print)
[...]
pca <HDF5 group "/uns/pca" (3 members)>
pca/variance <HDF5 dataset "variance": shape (50,), type "<f8">
pca/variance_ratio <HDF5 dataset "variance_ratio": shape (50,), type "<f8">
[...]
```

````

````{tab-item} Zarr
:sync: zarr

```python
>>> store["uns"].visititems(print)
[...]
pca <zarr.hierarchy.Group '/uns/pca' read-only>
pca/variance <zarr.core.Array '/uns/pca/variance' (50,) float64 read-only>
pca/variance_ratio <zarr.core.Array '/uns/pca/variance_ratio' (50,) float64 read-only>
[...]
```

````

`````



### Mapping specifications (v0.1.0)

* Each mapping MUST be its own group
* The group's metadata MUST contain the encoding metadata `"encoding-type": "dict"`, `"encoding-version": "0.1.0"`

## Scalars

Zero dimensional arrays are used for scalar values (i.e. single values like strings, numbers or booleans).
These should only occur inside of `uns`, and are commonly saved parameters:

`````{tab-set}

````{tab-item} HDF5
:sync: hdf5

```python
>>> store["uns/neighbors/params"].visititems(print)
method <HDF5 dataset "method": shape (), type "|O">
metric <HDF5 dataset "metric": shape (), type "|O">
n_neighbors <HDF5 dataset "n_neighbors": shape (), type "<i8">
random_state <HDF5 dataset "random_state": shape (), type "<i8">
```

````

````{tab-item} Zarr
:sync: zarr

```python
>>> store["uns/neighbors/params"].visititems(print)
method <zarr.core.Array '/uns/neighbors/params/method' () <U4 read-only>
metric <zarr.core.Array '/uns/neighbors/params/metric' () <U9 read-only>
n_neighbors <zarr.core.Array '/uns/neighbors/params/n_neighbors' () int64 read-only>
random_state <zarr.core.Array '/uns/neighbors/params/random_state' () int64 read-only>
```

````

`````

```python
>>> store["uns/neighbors/params/metric"][()]
'euclidean'
>>> dict(store["uns/neighbors/params/metric"].attrs)
{'encoding-type': 'string', 'encoding-version': '0.2.0'}
```

### Scalar specification (v0.2.0)

* Scalars MUST be written as a 0 dimensional array
* Numeric scalars
    * MUST have `"encoding-type": "numeric-scalar"`, `"encoding-version": "0.2.0"` in their metadata
    * MUST be a single numeric value, including boolean, unsigned integer, signed integer,  floating point, or complex floating point
* String scalars
    * MUST have `"encoding-type": "string"`, `"encoding-version": "0.2.0"` in their metadata
    * In zarr, scalar strings MUST be stored as a fixed length unicode dtype
    * In HDF5, scalar strings MUST be stored as a variable length utf-8 encoded string dtype

## Categorical arrays

```python
>>> categorical = store["obs"]["development_stage"]
>>> dict(categorical.attrs)
{'encoding-type': 'categorical', 'encoding-version': '0.2.0', 'ordered': False}
```

Discrete values can be efficiently represented with categorical arrays (similar to `factors` in `R`).
These arrays encode the values as small width integers (`codes`), which map to the original label set (`categories`).
Each entry in the `codes` array is the zero-based index of the encoded value in the `categories` array.
To represent a missing value, a code of `-1` is used.
We store these two arrays separately.

`````{tab-set}

````{tab-item} HDF5
:sync: hdf5

```python
>>> categorical.visititems(print)
categories <HDF5 dataset "categories": shape (7,), type "|O">
codes <HDF5 dataset "codes": shape (164114,), type "|i1">
```

````

````{tab-item} Zarr
:sync: zarr

```python
>>> categorical.visititems(print)
categories <zarr.core.Array '/obs/development_stage/categories' (7,) object read-only>
codes <zarr.core.Array '/obs/development_stage/codes' (164114,) int8 read-only>
```

````

`````

### Categorical array specification (v0.2.0)

* Categorical arrays MUST be stored as a group
* The group's metadata MUST contain the encoding metadata `"encoding-type": "categorical"`, `"encoding-version": "0.2.0"`
* The group's metadata MUST contain the boolean valued field `"ordered"`, which indicates whether the categories are ordered
* The group MUST contain an integer valued array named `"codes"` whose maximum value is the number of categories - 1
    * The `"codes"` array MAY contain signed integer values. If so, the code `-1` denotes a missing value
* The group MUST contain an array called `"categories"`

## String arrays

Arrays of strings are handled differently than numeric arrays since numpy doesn't really have a good way of representing arrays of unicode strings.
`anndata` assumes strings are text-like data, so it uses a variable length encoding.

`````{tab-set}

````{tab-item} HDF5
:sync: hdf5

```python
>>> store["var"][store["var"].attrs["_index"]]
<HDF5 dataset "ensembl_id": shape (40145,), type "|O">
```

````

````{tab-item} Zarr
:sync: zarr

```python
>>> store["var"][store["var"].attrs["_index"]]
<zarr.core.Array '/var/ensembl_id' (40145,) object read-only>
```

````

`````

```python
>>> dict(categorical["categories"].attrs)
{'encoding-type': 'string-array', 'encoding-version': '0.2.0'}
```

### String array specifications (v0.2.0)

* String arrays MUST be stored in arrays
* The arrays's metadata MUST contain the encoding metadata `"encoding-type": "string-array"`, `"encoding-version": "0.2.0"`
* In `zarr`, string arrays MUST be stored using `numcodecs`' `VLenUTF8` codec
* In `HDF5`, string arrays MUST be stored using the variable length string data type, with a utf-8 encoding

## Nullable integers and booleans

We support IO with Pandas nullable integer and boolean arrays.
We represent these on disk similar to `numpy` masked arrays, `julia` nullable arrays, or `arrow` validity bitmaps (see {issue}`504` for more discussion).
That is, we store an indicator array (or mask) of null values alongside the array of all values.

`````{tab-set}

````{tab-item} HDF5
:sync: hdf5

```python
>>> from anndata import write_elem
>>> null_store = h5py.File("tmp.h5", mode="w")
>>> int_array = pd.array([1, None, 3, 4])
>>> int_array
<IntegerArray>
[1, <NA>, 3, 4]
Length: 4, dtype: Int64

>>> write_elem(null_store, "nullable_integer", int_array)

>>> null_store.visititems(print)
nullable_integer <HDF5 group "/nullable_integer" (2 members)>
nullable_integer/mask <HDF5 dataset "mask": shape (4,), type "|b1">
nullable_integer/values <HDF5 dataset "values": shape (4,), type "<i8">
```

````

````{tab-item} Zarr
:sync: zarr

```python
>>> from anndata import write_elem
>>> null_store = zarr.open()
>>> int_array = pd.array([1, None, 3, 4])
>>> int_array
<IntegerArray>
[1, <NA>, 3, 4]
Length: 4, dtype: Int64

>>> write_elem(null_store, "nullable_integer", int_array)

>>> null_store.visititems(print)
nullable_integer <zarr.hierarchy.Group '/nullable_integer'>
nullable_integer/mask <zarr.core.Array '/nullable_integer/mask' (4,) bool>
nullable_integer/values <zarr.core.Array '/nullable_integer/values' (4,) int64>
```

````

`````

```python
>>> dict(null_store["nullable_integer"].attrs)
{'encoding-type': 'nullable-integer', 'encoding-version': '0.1.0'}
```

### Nullable integer specifications (v0.1.0)

* Nullable integers MUST be stored as a group
* The group's attributes MUST have contain the encoding metadata `"encoding-type": "nullable-integer"`, `"encoding-version": "0.1.0"`
* The group MUST contain an integer valued array under the key `"values"`
* The group MUST contain an boolean valued array under the key `"mask"`

### Nullable boolean specifications (v0.1.0)

* Nullable booleans MUST be stored as a group
* The group's attributes MUST have contain the encoding metadata `"encoding-type": "nullable-boolean"`, `"encoding-version": "0.1.0"`
* The group MUST contain an boolean valued array under the key `"values"`
* The group MUST contain an boolean valued array under the key `"mask"`
* The `"values"` and `"mask"` arrays MUST be the same shape

## AwkwardArrays

```{warning}
**Experimental**

Support for ragged arrays via awkward array is considered experimental under the 0.9.0 release series.
Please direct feedback on it's implementation to [https://github.com/scverse/anndata](https://github.com/scverse/anndata).
```

Ragged arrays are supported in `anndata` through the [Awkward
Array](https://awkward-array.org/) library. For storage on disk, we
break down the awkward array into it’s constituent arrays using
[`ak.to_buffers`](https://awkward-array.readthedocs.io/en/latest/_auto/ak.to_buffers.html)
then writing these arrays using `anndata`’s methods.

`````{tab-set}

````{tab-item} HDF5
:sync: hdf5

```python
>>> store["varm/transcript"].visititems(print)
node1-mask <HDF5 dataset "node1-mask": shape (5019,), type "|u1">
node10-data <HDF5 dataset "node10-data": shape (250541,), type "<i8">
node11-mask <HDF5 dataset "node11-mask": shape (5019,), type "|u1">
node12-offsets <HDF5 dataset "node12-offsets": shape (40146,), type "<i8">
node13-mask <HDF5 dataset "node13-mask": shape (250541,), type "|i1">
node14-data <HDF5 dataset "node14-data": shape (250541,), type "<i8">
node16-offsets <HDF5 dataset "node16-offsets": shape (40146,), type "<i8">
node17-data <HDF5 dataset "node17-data": shape (602175,), type "|u1">
node2-offsets <HDF5 dataset "node2-offsets": shape (40146,), type "<i8">
node3-data <HDF5 dataset "node3-data": shape (600915,), type "|u1">
node4-mask <HDF5 dataset "node4-mask": shape (5019,), type "|u1">
node5-offsets <HDF5 dataset "node5-offsets": shape (40146,), type "<i8">
node6-data <HDF5 dataset "node6-data": shape (59335,), type "|u1">
node7-mask <HDF5 dataset "node7-mask": shape (5019,), type "|u1">
node8-offsets <HDF5 dataset "node8-offsets": shape (40146,), type "<i8">
node9-mask <HDF5 dataset "node9-mask": shape (250541,), type "|i1">
```

````

````{tab-item} Zarr
:sync: zarr

```python
>>> store["varm/transcript"].visititems(print)
node1-mask <zarr.core.Array '/varm/transcript/node1-mask' (5019,) uint8 read-only>
node10-data <zarr.core.Array '/varm/transcript/node10-data' (250541,) int64 read-only>
node11-mask <zarr.core.Array '/varm/transcript/node11-mask' (5019,) uint8 read-only>
node12-offsets <zarr.core.Array '/varm/transcript/node12-offsets' (40146,) int64 read-only>
node13-mask <zarr.core.Array '/varm/transcript/node13-mask' (250541,) int8 read-only>
node14-data <zarr.core.Array '/varm/transcript/node14-data' (250541,) int64 read-only>
node16-offsets <zarr.core.Array '/varm/transcript/node16-offsets' (40146,) int64 read-only>
node17-data <zarr.core.Array '/varm/transcript/node17-data' (602175,) uint8 read-only>
node2-offsets <zarr.core.Array '/varm/transcript/node2-offsets' (40146,) int64 read-only>
node3-data <zarr.core.Array '/varm/transcript/node3-data' (600915,) uint8 read-only>
node4-mask <zarr.core.Array '/varm/transcript/node4-mask' (5019,) uint8 read-only>
node5-offsets <zarr.core.Array '/varm/transcript/node5-offsets' (40146,) int64 read-only>
node6-data <zarr.core.Array '/varm/transcript/node6-data' (59335,) uint8 read-only>
node7-mask <zarr.core.Array '/varm/transcript/node7-mask' (5019,) uint8 read-only>
node8-offsets <zarr.core.Array '/varm/transcript/node8-offsets' (40146,) int64 read-only>
node9-mask <zarr.core.Array '/varm/transcript/node9-mask' (250541,) int8 read-only>
```

````

`````



The length of the array is saved to it’s own `"length"` attribute,
while metadata for the array structure is serialized and saved to the
`“form”` attribute.

```python
>>> dict(store["varm/transcript"].attrs)
{'encoding-type': 'awkward-array',
 'encoding-version': '0.1.0',
 'form': '{"class": "RecordArray", "fields": ["tx_id", "seq_name", '
         '"exon_seq_start", "exon_seq_end", "ensembl_id"], "contents": '
         '[{"class": "BitMaskedArray", "mask": "u8", "valid_when": true, '
         '"lsb_order": true, "content": {"class": "ListOffsetArray", '
         '"offsets": "i64", "content": {"class": "NumpyArray", "primitive": '
         '"uint8", "inner_shape": [], "parameters": {"__array__": "char"}, '
         '"form_key": "node3"}, "parameters": {"__array__": "string"}, '
         '"form_key": "node2"}, "parameters": {}, "form_key": "node1"}, '
        ...
 'length': 40145}
```

These can be read back as awkward arrays using the
[`ak.from_buffers`](https://awkward-array.readthedocs.io/en/latest/_auto/ak.from_buffers.html)
function:

```python
>>> import awkward as ak
>>> from anndata.io import read_elem
>>> awkward_group = store["varm/transcript"]
>>> ak.from_buffers(
...     awkward_group.attrs["form"],
...     awkward_group.attrs["length"],
...     {k: read_elem(v) for k, v in awkward_group.items()}
... )
>>> transcript_models[:5]
[{tx_id: 'ENST00000450305', seq_name: '1', exon_seq_start: [...], ...},
 {tx_id: 'ENST00000488147', seq_name: '1', exon_seq_start: [...], ...},
 {tx_id: 'ENST00000473358', seq_name: '1', exon_seq_start: [...], ...},
 {tx_id: 'ENST00000477740', seq_name: '1', exon_seq_start: [...], ...},
 {tx_id: 'ENST00000495576', seq_name: '1', exon_seq_start: [...], ...}]
-----------------------------------------------------------------------
type: 5 * {
    tx_id: ?string,
    seq_name: ?string,
    exon_seq_start: option[var * ?int64],
    exon_seq_end: option[var * ?int64],
    ensembl_id: ?string
}
>>> transcript_models[0]
{tx_id: 'ENST00000450305',
 seq_name: '1',
 exon_seq_start: [12010, 12179, 12613, 12975, 13221, 13453],
 exon_seq_end: [12057, 12227, 12697, 13052, 13374, 13670],
 ensembl_id: 'ENSG00000223972'}
------------------------------------------------------------
type: {
    tx_id: ?string,
    seq_name: ?string,
    exon_seq_start: option[var * ?int64],
    exon_seq_end: option[var * ?int64],
    ensembl_id: ?string
}
```


[easy to find]: https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)
[hdf5]: https://en.wikipedia.org/wiki/Hierarchical_Data_Format
