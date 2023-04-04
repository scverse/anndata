# Interoperability

The on-disk representation of anndata files can be read from other
languages. Here we list interfaces for working with AnnData from your
language of choice:

## R

- [zellkonverter](https://bioconductor.org/packages/release/bioc/html/zellkonverter.html) zellkonverter provides basilisk based tooling for loading from `h5ad` files to `SingleCellExperiment`
- [anndata](https://anndata.dynverse.org) provides an R implementation of `AnnData` as well as IO for the HDF5 format.
- [MuData](https://bioconductor.org/packages/release/bioc/html/MuData.html) provides IO for `AnnData` and `MuData` stored in HDF5 to Bioconductor's `SingleCellExperiment` and `MultiAssayExperiment` objects.
- [MuDataSeurat](https://pmbio.github.io/MuDataSeurat/) provides IO from `AnnData` and `MuData` stored in HDF5 to `Seurat` objects.

## Julia

- [Muon.jl](https://docs.juliahub.com/Muon/QfqCh/0.1.1/objects/) provides Julia implementations of `AnnData` and `MuData` objects, as well as IO for the HDF5 format
- [scVI.jl](https://maren-ha.github.io/scVI.jl/index.html) provides a Julia implementation of `AnnData` as well as IO for the HDF5 format.

## Javascript

- [Vitessce](https://github.com/vitessce/vitessce) contains loaders from `AnnData`s stored as Zarr, and uses this to provide interactive visualization

## Rust

- [anndata-rs](https://github.com/kaizhang/anndata-rs) provides a Rust implementation of `AnnData` as well as advanced IO support for the HDF5 storage format.
