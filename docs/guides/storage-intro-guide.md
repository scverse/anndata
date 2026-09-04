# How to work with with stored AnnData objects

By using array storage libraries for it's storage format, `AnnData` allows efficient partial access both within python, and beyond. This article will briefly introduce a few APIs for working with stored `AnnData` objects.

## IO with individual elements

An AnnData object is made up of a set of constituent elements. These elements and their on-disk formats are described in {doc}`file format </fileformat-prose>`.

Each of these elements can be individually read and written with {func}`anndata.experimental.read_elem` and {func}`anndata.experimental.write_elem`. These functions work on all storage backends that `AnnData` supports and all the types.

For example, given an stored anndata object:

```python
import s3fs
import zarr
import anndata as ad
from anndata.experimental import read_elem, write_elem

store = zarr.open_consolidated("...")

obs = read_elem(store["obs"])
obs
```

You can use this to load minimal amounts of information into your anndata object at once. E.g.:

```python
adata = ad.AnnData(
    obs=read_elem(store["obs"]),
    var=read_elem(store["var"]),
    obsm={"X_umap": read_elem(store["obsm/X_umap"])},
    uns=read_elem(store["uns"]),
)
adata
```

```python
import scanpy as sc

sc.pl.umap(adata, color="cell_type")
```

## Partial access to stored elements

Some elements can be partially accessed. All of these element types can be referenced in an AnnData object.

### `zarr.Array` and `h5py.Dataset`

The simplest access can be done with {class}`zarr.core.Array` and {class}`h5py.Dataset`.

```python
embedding = store["obsm/X_embedding"]
embedding
# zarr.Array
```

This storage array can be indexed into, read, or written to.

```python
```

### `sparse_dataset`

A slightly more limited set of accessor classes is defined in the `anndata` library. These classes allow backed partial access to sparse arrays in CSC or CSR format.

```{note}
It's very important to access sparse matrices on their major dimension. That means selecting rows from a CSR matrix and columns from a CSC matrix.

While this is pretty important in memory, it's hugely important when these objects are in slower storage. If you know that you'll need both column and row based access into your matrices, we suggest saving a CSC and CSR copy of the matrix.
```

```python
from anndata.experimental import sparse_dataset
```

```python
%%time
X = sparse_dataset(store["X"])
X
```

```python
%%time
X[:1000]
```

### Dask

For more complex interactions with out of core data, you can use a framework like [`dask`](https://www.dask.org).  Dask works best with dense arrays, but has limited support for sparse arrays. This is an active area of upstream development work by scverse devs.

Unlike the `sparse_dataset`s and `HDF5` and `zarr` classes, you can perform more complex operations with `dask` arrays. In addition, you can index into them more than once.

#### Dask with dense arrays

```python
SparseDataset.write_dense("", chunks=())
```

#### Dask with sparse arrays

We're not quite comfortable with the stability of this to recommend it in the docs. For the brave of heart, you can find an example of working with dask with sparse chunks in this [`gist`](https://gist.github.com/ivirshup/3fbe634b648304978ea77469b5d88961).

## Combining stored objects


## Further reading:

* [`dask.array` documentation](https://docs.dask.org/en/stable/array.html)
* Documentation for [`h5py`](https://docs.h5py.org/en/stable/) and [`zarr`](https://zarr.readthedocs.io/en/stable/index.html)
* {func}`anndata.experimental.concat_on_disk`
* {func}`anndata.experimental.read_elem`, {func}`anndata.experimental.write_elem`
