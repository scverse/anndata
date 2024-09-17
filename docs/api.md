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

Combining {class}`AnnData` objects.
See also the section on concatenation.

```{eval-rst}
.. autosummary::
   :toctree: generated/

   concat
```

## Reading

Reading anndata’s native formats `.h5ad` and `zarr`.

```{eval-rst}
.. autosummary::
   :toctree: generated/

   read_h5ad
   read_zarr
```

Reading individual portions ({attr}`~AnnData.obs`, {attr}`~AnnData.varm` etc.) of the {class}`AnnData` object.

```{eval-rst}
.. autosummary::
   :toctree: generated/

   read_elem
   sparse_dataset
```

Reading file formats that cannot represent all aspects of {class}`AnnData` objects.

```{tip}
You might have more success by assembling the {class}`AnnData` object yourself from the individual parts.
```

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
```

## Writing

Writing a complete {class}`AnnData` object to disk in anndata’s native formats `.h5ad` and `zarr`.

```{eval-rst}
.. autosummary::
   :toctree: generated/

   AnnData.write
   AnnData.write_zarr
```

Writing individual portions ({attr}`~AnnData.obs`, {attr}`~AnnData.varm` etc.) of the {class}`AnnData` object.

```{eval-rst}
.. autosummary::
   :toctree: generated/

   write_elem
```

Writing formats that cannot represent all aspects of {class}`AnnData` objects.

```{eval-rst}
.. autosummary::
   :toctree: generated/

   AnnData.write_csvs
   AnnData.write_loom
```

(experimental-api)=

## Experimental API

```{warning}
APIs in the experimental module are currently in development and subject to change at any time.
```

Two classes for working with batched access to collections of many {class}`AnnData` objects or `.h5ad` files.
In particular, for pytorch-based models.

```{eval-rst}
.. autosummary::
   :toctree: generated/

   experimental.AnnCollection
   experimental.AnnLoader
```

Out of core concatenation

```{eval-rst}
.. autosummary::
   :toctree: generated/

   experimental.concat_on_disk
```

Low level methods for reading and writing elements of an {class}`AnnData` object to a store:

```{eval-rst}
.. autosummary::
   :toctree: generated/

   experimental.read_elem_as_dask
```

Utilities for customizing the IO process:

```{eval-rst}
.. autosummary::
   :toctree: generated/

   experimental.read_dispatched
   experimental.write_dispatched
```

Types used by the former:

```{eval-rst}
.. autosummary::
   :toctree: generated/

   experimental.IOSpec
   experimental.Read
   experimental.Write
   experimental.ReadCallback
   experimental.WriteCallback
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

## Custom Types/Classes for Readable/Writeable Elements

```{eval-rst}
.. autosummary::
   :toctree: generated/

   abc.CSRDataset
   abc.CSCDataset
   typing.Index
   typing.AxisStorable
   typing.RWAble
```
