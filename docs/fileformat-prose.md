## Disk representation

*Note: these docs are written for anndata 0.7. Files written before this version may differ in some conventions, but will still be read by newer versions of the library.*

AnnData objects are saved on disk to hierarchichal array stores like `HDF5` and `zarr`. This allows us to have very similar structures in disk and on memory.

AnnData objects can hold three kinds of objects in it's dimensioned mappings (i.e. `X`, `obsm`, `layers` etc.). These are (1) dense arrays, (2) sparse arrays, and (3) data frames. As an example we'll look into a typical `.h5ad` object that's been through an analysis.

**NOTE:** I've started using h5py since I couldn't figure out a nice way to print attributes from bash

```python
>>> import h5py
>>> f = h5py.File("02_processed.h5ad", "r")
>>> list(f.keys())
['X', 'layers', 'obs', 'obsm', 'uns', 'var', 'varm'] 
```

```bash
$ h5ls 02_processed.h5ad
X                        Group
layers                   Group
obs                      Group
obsm                     Group
uns                      Group
var                      Group
varm                     Group
```

### Dense arrays

Dense arrays have the most simple representation on disk, as they have clear equivalents in hdf5 and zarr. Any dense array will be written to the file as a `dataset`. We can see an example of this for principle components stored in the `obsm` group:

```python
>>> f["obsm"].visititems(lambda k, v: print(f"{k}: {v}")) 
X_pca: <HDF5 dataset "X_pca": shape (38410, 50), type "<f4">
X_umap: <HDF5 dataset "X_umap": shape (38410, 2), type "<f4">
```



```bash
$ h5ls 02_processed.h5ad/obsm
X_pca                    Dataset {38410, 50}
X_umap                   Dataset {38410, 2}
```

### Sparse arrays

Sparse arrays don't have a native  representations in hdf5 or zarr, so we use a representation as close to the in memory datasets as we can. Currently two sparse data formats are supported by `AnnData` objects. These are CSC and CSR formats, where a two dimensional sparse array is represented by three one dimensional arrays, `indptr`, `indices`, and `data`. A full description of this format is out of scope for this document, but are widley available ([wikipedia description](https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)))

On disk we represent a sparse array as a group. These the kind and shape of sparse array can be identified by their attributes:

```python
>>> dict(f["X"].attrs)
{'encoding-type': 'csr_matrix',
 'encoding-version': '0.1.0',
 'shape': array([38410, 27899])}
```

Inside the group are the three constituent arrays:

```python
>>> f["X"].visititems(lambda k, v: print(f"{k}: {v}")) 
data: <HDF5 dataset "data": shape (41459314,), type "<f4">
indices: <HDF5 dataset "indices": shape (41459314,), type "<i4">
indptr: <HDF5 dataset "indptr": shape (38411,), type "<i4">
```

```bash
$ h5ls 02_processed.h5ad/X
data                     Dataset {41459314/Inf}
indices                  Dataset {41459314/Inf}
indptr                   Dataset {38411/Inf}
```

### DataFrames

Data frames are saved as a columnar format in a group, so each column of a dataframe gets it's own dataset. To maintain the efficiency of categorical values they are stored by as their numeric codes with their values saved in a reserved subgroup `__categories`.

Dataframes can be identified from other groups by their attributes:

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
 'encoding-version': '0.1.0'}
```

These attributes identify the column used as an index, the order of the original columns, and some type information.

```python
>>> f["obs"].visititems(lambda k, v: print(f"{k}: {v}"))
Cell: <HDF5 dataset "Cell": shape (38410,), type "|O">
__categories: <HDF5 group "/obs/__categories" (3 members)>
__categories/cell_type: <HDF5 dataset "cell_type": shape (22,), type "|O">
__categories/label_by_score: <HDF5 dataset "label_by_score": shape (16,), type "|O">
__categories/sample: <HDF5 dataset "sample": shape (41,), type "|O">
cell_type: <HDF5 dataset "cell_type": shape (38410,), type "|i1">
label_by_score: <HDF5 dataset "label_by_score": shape (38410,), type "|i1">
log1p_n_genes_by_counts: <HDF5 dataset "log1p_n_genes_by_counts": shape (38410,), type "<f8">
log1p_total_counts: <HDF5 dataset "log1p_total_counts": shape (38410,), type "<f4">
log1p_total_counts_mito: <HDF5 dataset "log1p_total_counts_mito": shape (38410,), type "<f4">
n_genes_by_counts: <HDF5 dataset "n_genes_by_counts": shape (38410,), type "<i4">
pct_counts_in_top_100_genes: <HDF5 dataset "pct_counts_in_top_100_genes": shape (38410,), type "<f8">
pct_counts_in_top_200_genes: <HDF5 dataset "pct_counts_in_top_200_genes": shape (38410,), type "<f8">
pct_counts_in_top_500_genes: <HDF5 dataset "pct_counts_in_top_500_genes": shape (38410,), type "<f8">
pct_counts_in_top_50_genes: <HDF5 dataset "pct_counts_in_top_50_genes": shape (38410,), type "<f8">
pct_counts_mito: <HDF5 dataset "pct_counts_mito": shape (38410,), type "<f4">
sample: <HDF5 dataset "sample": shape (38410,), type "|i1">
total_counts: <HDF5 dataset "total_counts": shape (38410,), type "<f4">
total_counts_mito: <HDF5 dataset "total_counts_mito": shape (38410,), type "<f4">
```

Categorical series can be identified by the presence of the attribute `"categories"`, which contains a pointer to their categorical values:

*Note:* as `zarr` does not have reference objects, in zarr files the `categories` attribute is an absolute path to the category values.

```python
>>> dict(f["obs/cell_type"].attrs)
{'categories': <HDF5 object reference>}
```

## Other values:

### Mappings

Mappings are stored as native groups in an `h5ad` file. These can be identified as being seperate from dataframes and sparse arrays since they don't have any special attributes. These are used for any `Mapping` in the AnnData object, including the default `obsm`, `varm`, `layers`, and `uns`. This definition is used recursivley within `uns`:

```python
>>> f["uns"].visititems(print) 
...
pca <HDF5 group "/uns/pca" (2 members)>
pca/variance <HDF5 dataset "variance": shape (50,), type "<f4">
pca/variance_ratio <HDF5 dataset "variance_ratio": shape (50,), type "<f4">
...
```

### Scalars

Zero dimensional arrays are used for scalar values (i.e. single values like strings, numbers or booleans). These should only occur inside of `uns`, and are common inside of saved parameters:

```python
>>> f["uns/neighbors/params"].visititems(print)
method <HDF5 dataset "method": shape (), type "|O">
metric <HDF5 dataset "metric": shape (), type "|O">
n_neighbors <HDF5 dataset "n_neighbors": shape (), type "<i8">

>>> f["uns/neighbors/params/metric"][()]   
'euclidean'
```
