# On-disk format

```{note}
These docs are written for anndata 0.8.
Files written before this version may differ in some conventions,
but will still be read by newer versions of the library.
```

AnnData objects are saved on disk to hierarchical array stores like [HDF5]
(via {doc}`H5py <h5py:index>`) and {doc}`zarr:index`.
This allows us to have very similar structures in disk and on memory.

As an example we’ll look into a typical `.h5ad` object that’s been through an analysis.
This structure should be largely equivalent to Zarr structure, though there are a few minor differences.

## Elements


 <!-- I’ve started using h5py since I couldn’t figure out a nice way to print attributes from bash. -->

```python
>>> import h5py
>>> f = h5py.File("02_processed.h5ad", "r")
>>> list(f.keys())
['X', 'layers', 'obs', 'obsm', 'uns', 'var', 'varm']
```

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

In general, `AnnData` objects are comprised of a various types of elements.
Each element is encoded as either an Array (or Dataset in hdf5 terminology) or a collection of elements (e.g. Group) in the store.
We record the type of an element using the `encoding-type` and `encoding-version` keys in it's attributes.
For example, we can this file represents an `AnnData` object from this metadata:

```python
>>> dict(f.attrs)
{'encoding-type': 'anndata', 'encoding-version': '0.1.0'}
```

Using this information, we're able to dispatch onto readers for the different element types that you'd find in an anndata.

### Element Specification

* An element can be any object within the storage hierarchy (typically an array or group) with associated metadata
* An element MUST have a string valued field `"encoding-type"` in it's metadata
* An element MUST have a string valued field `"encoding-version"` in it's metadata that can be evaluated to a version

### AnnData specification (v0.1.0)

* An `AnnData` object MUST be a group. 
* The group's metadata MUST include entries: `"encoding-type": "anndata"`, `"encoding-version": "0.1.0"`.
* An `AnnData` group MUST contain entries `"obs"` and `"var"`, which MUST be dataframes.
* The group MAY contain an entry `X`, which MUST be either a dense or sparse array and whose shape MUST be (`n_obs`, `n_var`)
* The group MAY contain a mapping `layers`. Entries in `layers` MUST be dense or sparse arrays which have shapes (`n_obs`, `n_var`)
* The group MAY contain a mapping `obsm`. Entries in `obsm` MUST be sparse arrays, dense arrays, or dataframes. These entries MUST have a first dimension of size `n_obs`
* The group MAY contain a mapping `varm`. Entries in `varm` MUST be sparse arrays, dense arrays, or dataframes. These entries MUST have a first dimension of size `n_var`
* The group MAY contain a mapping `obsp`. Entries in `obsp` MUST be sparse or dense arrays. The entries first two dimensions MUST be of size `n_obs`
* The group MAY contain a mapping `varp`. Entries in `varp` MUST be sparse or dense arrays. The entries first two dimensions MUST be of size `n_var`
* The group MAY contain a mapping `uns`. Entries in `uns` MUST be an anndata encoded type.

## Dense arrays

Dense numeric arrays have the most simple representation on disk,
as they have native equivalents in H5py {doc}`h5py:high/dataset` and Zarr {ref}`Arrays <zarr:tutorial_create>`.
We can see an example of this with dimensionality reductions stored in the `obsm` group:

```python
>>> f["obsm"].visititems(print)
X_pca <HDF5 dataset "X_pca": shape (38410, 50), type "<f4">
X_umap <HDF5 dataset "X_umap": shape (38410, 2), type "<f4">

>>> dict(f["obsm"]["X_pca"].attrs)
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
>>> dict(f["X"].attrs)
{'encoding-type': 'csr_matrix',
 'encoding-version': '0.1.0',
 'shape': array([38410, 27899])}
```

Inside the group are the three constituent arrays:

```python
>>> f["X"].visititems(print)
data <HDF5 dataset "data": shape (41459314,), type "<f4">
indices <HDF5 dataset "indices": shape (41459314,), type "<i4">
indptr <HDF5 dataset "indptr": shape (38411,), type "<i4">
```

<!-- ```bash
$ h5ls 02_processed.h5ad/X
data                     Dataset {41459314/Inf}
indices                  Dataset {41459314/Inf}
indptr                   Dataset {38411/Inf}
``` -->

### Sparse array specification (v0.1.0)

* Each sparse array MUST be it's own group
* The group MUST contain arrays `indices`, `indptr`, and `data`
* The group's metadata MUST contain:
    * `"encoding-type"`, which is set to `"csr_matrix"` or `"csc_matrix"` for compressed sparse row and compressed sparse column, respectively.
    * `"encoding-version"`, which is set to `"0.1.0"`
    * `"shape"` which is an integer array of length 2 whose values are the sizes of the array's dimensions

## DataFrames

DataFrames are saved as a columnar format in a group, so each column of a DataFrame is saved as a separate array.
We save a little more information in the attributes here.

```python
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
```

These attributes identify the index of the dataframe, as well as the original order of the columns.
Each column in this dataframe is encoded as it's own array.

```python
>>> dict(f["obs"]["total_counts"].attrs)
{'encoding-type': 'array', 'encoding-version': '0.2.0'}

>>> dict(f["obs"]["cell_type"].attrs)
{'encoding-type': 'categorical', 'encoding-version': '0.2.0', 'ordered': False}
```

### Dataframe Specification (v0.2.0)

* A dataframe MUST be stored as a group
* The group's metadata
    * MUST contain the field `"_index"`, which is the name of the column to be used as an index
    * MUST contain encoding metadata `"encoding-type: "dataframe"`, `"encoding-version": "0.2.0"`
    * MUST contain `"column-order"` a string valued array denoting the order of column entries
* The group MUST contain a column for the index
* Each entry in the group MUST correspond to an array with equivalent first dimensions
* Each entry SHOULD share chunk sizes

## Mappings

Mappings are simply stored as `Group` s on disk.
These are distinct from DataFrames and sparse arrays since they don’t have any special attributes.
A `Group` is created for any `Mapping` in the AnnData object,
including the standard `obsm`, `varm`, `layers`, and `uns`.
Notably, this definition is used recursively within `uns`:

```python
>>> f["uns"].visititems(print)
[...]
pca <HDF5 group "/uns/pca" (2 members)>
pca/variance <HDF5 dataset "variance": shape (50,), type "<f4">
pca/variance_ratio <HDF5 dataset "variance_ratio": shape (50,), type "<f4">
[...]
```

### Mapping specifications (v0.1.0)

* Each mapping MUST be it's own group
* The groups metadata MUST contain the encoding metadata `"encoding-type": "dict"`, `"encoding-version": "0.1.0"`

## Scalars

Zero dimensional arrays are used for scalar values (i.e. single values like strings, numbers or booleans).
These should only occur inside of `uns`, and are commonly saved parameters:

```python
>>> f["uns/neighbors/params"].visititems(print)
method <HDF5 dataset "method": shape (), type "|O">
metric <HDF5 dataset "metric": shape (), type "|O">
n_neighbors <HDF5 dataset "n_neighbors": shape (), type "<i8">
>>> f["uns/neighbors/params/metric"][()]
'euclidean'
>>> dict(f["uns/neighbors/params/metric"].attrs)
{'encoding-type': 'string', 'encoding-version': '0.2.0'}
```

### Scalar specification (v0.2.0)

* Scalars MUST be written as a 0 dimensional array
* Numeric scalars
    * MUST have `"encoding-type": "numeric-scalar"`, `"encoding-version": "0.2.0"` in their metadata
    * MUST be a single numeric value, including boolean, unsigned integer, signed integer,  floating point, or complex floating point.
* String scalars
    * MUST have `"encoding-type": "string"`, `"encoding-version": "0.2.0"` in their metadata
    * In zarr, scalar strings MUST be stored as a fixed length unicode dtype
    * In HDF5, scalar strings MUST be stored as a variable length utf-8 encoded string dtype

## Categorical arrays

```python
>>> categorical = f["obs"]["cell_type"]
>>> dict(categorical.attrs)
{'encoding-type': 'categorical', 'encoding-version': '0.2.0', 'ordered': False}
```

Discrete labels can be efficiently represented with categorical arrays (similar to `factors` in `R`).
These arrays encode the labels as small width integers (`codes`), which map to the original label set (`categories`).
We store these two arrays separately.

```python
>>> categorical.visititems(print)
categories <HDF5 dataset "categories": shape (22,), type "|O">
codes <HDF5 dataset "codes": shape (38410,), type "|i1">
```

### Categorical array specification (v0.2.0)

* Categorical arrays MUST be stored as a group
* The group's metadata MUST contain the encoding metadata `"encoding-type": "categorical"`, `"encoding-version": "0.2.0"`
* The group's metadata MUST contain the boolean valued field `"ordered"`, which indicates whether the categories are ordered
* The group MUST contain an integer valued array named `"codes"`
    * The `"codes"` array MAY contain signed integer values. If so, the code `-1` denotes a missing value.
* The group MUST contain an array called `"categories"`

## String arrays

Arrays of strings are handled differently than numeric arrays since numpy doesn't really have a good way of representing arrays of unicode strings.
`anndata` assumes strings are text like data, so are variable length.

```python
>>> dict(categorical["categories"].attrs)
{'encoding-type': 'string-array', 'encoding-version': '0.2.0'}
```

### String array specifications (v0.2.0)

* String arrays MUST be stored in arrays
* The arrays's metadata MUST contain the encoding metadata `"encoding-type": "string-array"`, `"encoding-version": "0.2.0"`
* In `zarr`, string arrays MUST be stored using `numcodecs`' `VLenUTF8` codec.
* In `HDF5`, string arrays MUST be stored using the variable length string data type, with a utf-8 encoding.

## Nullable integers and booleans

We support IO with Pandas nullable integer and boolean arrays.
We represent these on disk similar to `numpy` masked arrays, `julia` nullable arrays, or `arrow` validity bitmaps (see {issue}`504` for more discussion).
That is, we store a indicator array (or mask) of null values alongside the array of all values.

```python
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
* The `"values"` and `"mask"` array MUST be the same shape


[easy to find]: https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)
[hdf5]: https://en.wikipedia.org/wiki/Hierarchical_Data_Format