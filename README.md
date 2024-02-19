[![Build Status](https://dev.azure.com/scverse/anndata/_apis/build/status/scverse.anndata?branchName=main)](https://dev.azure.com/scverse/anndata/_build)
[![Conda](https://img.shields.io/conda/vn/conda-forge/anndata.svg)](https://anaconda.org/conda-forge/anndata)
[![Coverage](https://codecov.io/gh/scverse/anndata/branch/main/graph/badge.svg?token=IN1mJN1Wi8)](https://codecov.io/gh/scverse/anndata)
[![Docs](https://readthedocs.com/projects/icb-anndata/badge/?version=latest)](https://anndata.readthedocs.io)
[![PyPI](https://img.shields.io/pypi/v/anndata.svg)](https://pypi.org/project/anndata)
[![Downloads](https://static.pepy.tech/badge/anndata/month)](https://pepy.tech/project/anndata)
[![Downloads](https://static.pepy.tech/badge/anndata)](https://pepy.tech/project/anndata)
[![Stars](https://img.shields.io/github/stars/scverse/anndata?logo=GitHub&color=yellow)](https://github.com/scverse/anndata/stargazers)
[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](http://numfocus.org)

<img
  src="https://raw.githubusercontent.com/scverse/anndata/main/docs/_static/img/anndata_schema.svg"
  class="dark-light" align="right" width="350" alt="image"
/>

# anndata - Annotated data

anndata is a Python package for handling annotated data matrices in memory and on disk, positioned between pandas and xarray. anndata offers a broad range of computationally efficient features including, among others, sparse data support, lazy operations, and a PyTorch interface.

- Discuss development on [GitHub](https://github.com/scverse/anndata).
- Read the [documentation](https://anndata.readthedocs.io).
- Ask questions on the [scverse Discourse](https://discourse.scverse.org).
- Install via `pip install anndata` or `conda install anndata -c conda-forge`.
- See [Scanpy's documentation](https://scanpy.readthedocs.io/) for usage related to single cell data. anndata was initially built for Scanpy.

[//]: # (numfocus-fiscal-sponsor-attribution)

anndata is part of the scverse project ([website](https://scverse.org), [governance](https://scverse.org/about/roles)) and is fiscally sponsored by [NumFOCUS](https://numfocus.org/).
Please consider making a tax-deductible [donation](https://numfocus.org/donate-to-scverse) to help the project pay for developer time, professional services, travel, workshops, and a variety of other needs.


<a href="https://numfocus.org/project/scverse">
  <img
    src="https://raw.githubusercontent.com/numfocus/templates/master/images/numfocus-logo.png"
    width="200"
  >
</a>

## Citation

If you use `anndata` in your work, please cite the `anndata` pre-print as follows:

> **anndata: Annotated data**
>
> Isaac Virshup, Sergei Rybakov, Fabian J. Theis, Philipp Angerer, F. Alexander Wolf
>
> _bioRxiv_ 2021 Dec 19. doi: [10.1101/2021.12.16.473007](https://doi.org/10.1101/2021.12.16.473007).

You can cite the scverse publication as follows:

> **The scverse project provides a computational ecosystem for single-cell omics data analysis**
>
> Isaac Virshup, Danila Bredikhin, Lukas Heumos, Giovanni Palla, Gregor Sturm, Adam Gayoso, Ilia Kats, Mikaela Koutrouli, Scverse Community, Bonnie Berger, Dana Pe’er, Aviv Regev, Sarah A. Teichmann, Francesca Finotello, F. Alexander Wolf, Nir Yosef, Oliver Stegle & Fabian J. Theis
>
> _Nat Biotechnol._ 2023 Apr 10. doi: [10.1038/s41587-023-01733-8](https://doi.org/10.1038/s41587-023-01733-8).
