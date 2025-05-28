[![Tests](https://github.com/scverse/anndata/actions/workflows/test-cpu.yml/badge.svg)](https://github.com/scverse/anndata/actions)
[![Conda](https://img.shields.io/conda/vn/conda-forge/anndata.svg)](https://anaconda.org/conda-forge/anndata)
[![Coverage](https://codecov.io/gh/scverse/anndata/branch/main/graph/badge.svg?token=IN1mJN1Wi8)](https://codecov.io/gh/scverse/anndata)
[![Docs](https://readthedocs.com/projects/icb-anndata/badge/?version=latest)](https://anndata.readthedocs.io)
[![PyPI](https://img.shields.io/pypi/v/anndata.svg)](https://pypi.org/project/anndata)
[![Downloads](https://static.pepy.tech/badge/anndata/month)](https://pepy.tech/project/anndata)
[![Downloads](https://static.pepy.tech/badge/anndata)](https://pepy.tech/project/anndata)
[![Stars](https://img.shields.io/github/stars/scverse/anndata?style=flat&logo=github&color=yellow)](https://github.com/scverse/anndata/stargazers)
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

anndata is part of the scverse® project ([website](https://scverse.org), [governance](https://scverse.org/about/roles)) and is fiscally sponsored by [NumFOCUS](https://numfocus.org/).
If you like scverse® and want to support our mission, please consider making a tax-deductible [donation](https://numfocus.org/donate-to-scverse) to help the project pay for developer time, professional services, travel, workshops, and a variety of other needs.

<div align="center">
<a href="https://numfocus.org/project/scverse">
  <img
    src="https://raw.githubusercontent.com/numfocus/templates/master/images/numfocus-logo.png"
    width="200"
  >
</a>
</div>

## Public API

Our public API is documented in the [API section][] of these docs.
We cannot guarantee the stability of our internal APIs, whether it's the location of a function, its arguments, or something else.
In other words, we do not officially support (or encourage users to do) something like `from anndata._core import AnnData` as `_core` is both not documented and contains a [leading underscore][].
However, we are aware that [many users do use these internal APIs][] and thus encourage them to [open an issue][] or migrate to the public API.
That is, if something is missing from our public API as documented, for example a feature you wish to be exported publicly, please open an issue.

[api section]: https://anndata.readthedocs.io/en/stable/api.html
[leading underscore]: https://peps.python.org/pep-0008/#public-and-internal-interfaces
[many users do use these internal APIs]: https://github.com/search?q=%22anndata._io%22&type=code
[open an issue]: https://github.com/scverse/anndata/issues/new/choose


## Citation

If you use `anndata` in your work, please cite the `anndata` publication as follows:

> **anndata: Annotated data**
>
> Isaac Virshup, Sergei Rybakov, Fabian J. Theis, Philipp Angerer, F. Alexander Wolf
>
> _JOSS_ 2024 Sep 16. doi: [10.21105/joss.04371](https://doi.org/10.21105/joss.04371).

You can cite the scverse publication as follows:

> **The scverse project provides a computational ecosystem for single-cell omics data analysis**
>
> Isaac Virshup, Danila Bredikhin, Lukas Heumos, Giovanni Palla, Gregor Sturm, Adam Gayoso, Ilia Kats, Mikaela Koutrouli, Scverse Community, Bonnie Berger, Dana Pe’er, Aviv Regev, Sarah A. Teichmann, Francesca Finotello, F. Alexander Wolf, Nir Yosef, Oliver Stegle & Fabian J. Theis
>
> _Nat Biotechnol._ 2023 Apr 10. doi: [10.1038/s41587-023-01733-8](https://doi.org/10.1038/s41587-023-01733-8).
