# API

```{eval-rst}
.. module:: anndata
```

The central class:

```{eval-rst}
.. autosummary::
   :toctree: generated/

   AnnData
```

## Combining

Combining AnnData objects. See also the section on concatenation.

```{eval-rst}
.. autosummary::
   :toctree: generated/

   concat
```

## Reading

Reading anndata’s native file format `.h5ad`.

```{eval-rst}
.. autosummary::
   :toctree: generated/

   read_h5ad
```

Reading other file formats.

```{eval-rst}
.. autosummary::
   :toctree: generated/

   read_csv
   read_excel
   read_hdf
   read_loom
   read_mtx
   read_text
   read_umi_tools
   read_zarr

```

## Writing

Writing to anndata’s native file format `.h5ad`.

```{eval-rst}
.. autosummary::
   :toctree: generated/

   AnnData.write
```

Writing to other formats.

```{eval-rst}
.. autosummary::
   :toctree: generated/

   AnnData.write_csvs
   AnnData.write_loom
   AnnData.write_zarr
```

(experimental-api)=

## Experimental API

```{warning}
API's in the experimental module are currently in development and subject to change at any time.
```

Two classes for working with batched access to collections of many `AnnData` objects or `h5ad` files. In paritcular, for pytorch-based models.

```{eval-rst}
.. autosummary::
   :toctree: generated/

   experimental.AnnCollection
   experimental.AnnLoader
```

Interface for accessing on-disk sparse data:

```{eval-rst}
.. autosummary::
   :toctree: generated/

   experimental.sparse_dataset
   experimental.CSRDataset
   experimental.CSCDataset
```

Out of core concatenation

```{eval-rst}
.. autosummary::
   :toctree: generated/

   experimental.concat_on_disk
```

Low level methods for reading and writing elements of an `` AnnData` `` object to a store:

```{eval-rst}
.. autosummary::
   :toctree: generated/

   experimental.read_elem
   experimental.write_elem
```

Utilities for customizing the IO process:

```{eval-rst}
.. autosummary::
   :toctree: generated/

   experimental.read_dispatched
   experimental.write_dispatched
   experimental.IOSpec
   experimental.InMemoryReadElem
   experimental.InMemoryArrayOrScalarType
   experimental.Reader
   experimental.Writer
   experimental.read_callback
   experimental.write_callback
   experimental.StorageType
```

## Errors and warnings

```{eval-rst}
.. autosummary::
   :toctree: generated/

   ImplicitModificationWarning
```

## Settings

```{eval-rst}
.. autosummary::
   :toctree: generated/

   settings
   settings.override
```
