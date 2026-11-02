# Installation

anndata needs Python 3.13 or newer.

## PyPI

```console
pip install anndata
```

## conda-forge

```console
conda install -c conda-forge anndata
```

## Development version

```console
pip install git+https://github.com/scverse/anndata.git
```

To work on anndata itself, see {doc}`contributing`.

## Optional dependencies

`anndata[dask]` lets an `AnnData` hold {doc}`Dask arrays <tutorials/notebooks/anndata_dask_array>`.

`anndata[lazy]` adds Dask, xarray and what reading remote zarr stores needs, for {func}`~anndata.experimental.read_lazy` and {func}`~anndata.experimental.read_elem_lazy`.

`anndata[cu11]`, `anndata[cu12]` or `anndata[cu13]` installs the CuPy wheel for your CUDA version, so an `AnnData` can hold GPU arrays.
`anndata[gpu]` installs CuPy from source instead.
