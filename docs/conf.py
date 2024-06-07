from __future__ import annotations

import sys
from datetime import datetime
from importlib import metadata
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sphinx.application import Sphinx

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE / "extensions")]


# -- General configuration ------------------------------------------------

# General information
project = "anndata"
author = f"{project} developers"
copyright = f"{datetime.now():%Y}, {author}"
release = version = metadata.version("anndata")

# default settings
templates_path = ["_templates"]
html_static_path = ["_static"]
source_suffix = [".rst", ".md"]
master_doc = "index"
default_role = "literal"
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
    "tutorials/notebooks/*.rst",
]
pygments_style = "sphinx"

extensions = [
    "myst_parser",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx_autodoc_typehints",  # needs to be after napoleon
    "sphinx_issues",
    "sphinx_design",
    "sphinx_search.extension",
    "sphinxext.opengraph",
    "scanpydoc",  # needs to be before linkcode
    "sphinx.ext.linkcode",
    "nbsphinx",
    "IPython.sphinxext.ipython_console_highlighting",
]
myst_enable_extensions = [
    "html_image",  # So README.md can be used on github and sphinx docs
]

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = "bysource"
issues_github_path = "scverse/anndata"
# autodoc_default_flags = ['members']
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
typehints_defaults = "braces"
todo_include_todos = False
nitpicky = True  # Report broken links
nitpick_ignore = [
    ("py:class", "scipy.sparse.base.spmatrix"),
    ("py:meth", "pandas.DataFrame.iloc"),
    ("py:meth", "pandas.DataFrame.loc"),
    ("py:class", "anndata._core.views.ArrayView"),
    ("py:class", "anndata._core.raw.Raw"),
    *[
        ("py:class", f"anndata._core.aligned_mapping.{cls}{kind}")
        for cls in "Layers AxisArrays PairwiseArrays".split()
        for kind in ["", "View"]
    ],
    # TODO: sphinx’ builtin autodoc.typehints extension isn’t handled by `qualname_overrides` yet
    # https://github.com/theislab/scanpydoc/issues/140
    ("py:class", "h5py._hl.group.Group"),
]
suppress_warnings = [
    "ref.citation",
    "myst.header",  # https://github.com/executablebooks/MyST-Parser/issues/262
]


def setup(app: Sphinx):
    # Don’t allow broken links. DO NOT CHANGE THIS LINE, fix problems instead.
    app.warningiserror = True


intersphinx_mapping = dict(
    h5py=("https://docs.h5py.org/en/latest/", None),
    hdf5plugin=("https://hdf5plugin.readthedocs.io/en/latest/", None),
    loompy=("https://linnarssonlab.org/loompy/", None),
    numpy=("https://numpy.org/doc/stable/", None),
    pandas=("https://pandas.pydata.org/pandas-docs/stable/", None),
    python=("https://docs.python.org/3", None),
    scipy=("https://docs.scipy.org/doc/scipy/", None),
    sklearn=("https://scikit-learn.org/stable/", None),
    zarr=("https://zarr.readthedocs.io/en/stable/", None),
    xarray=("https://xarray.pydata.org/en/stable/", None),
)
qualname_overrides = {
    "h5py._hl.group.Group": "h5py.Group",
    "h5py._hl.files.File": "h5py.File",
    "h5py._hl.dataset.Dataset": "h5py.Dataset",
    "anndata._core.anndata.AnnData": "anndata.AnnData",
}

# -- Social cards ---------------------------------------------------------

ogp_site_url = "https://anndata.readthedocs.io/"
ogp_image = "https://anndata.readthedocs.io/en/latest/_static/img/anndata_schema.svg"

# -- Options for HTML output ----------------------------------------------


# The theme is sphinx-book-theme, with patches for readthedocs-sphinx-search
html_theme = "scanpydoc"
html_theme_options = dict(
    use_repository_button=True,
    repository_url="https://github.com/scverse/anndata",
    repository_branch="main",
    navigation_with_keys=False,  # https://github.com/pydata/pydata-sphinx-theme/issues/1492
)
html_logo = "_static/img/anndata_schema.svg"
issues_github_path = "scverse/anndata"
html_show_sphinx = False


# -- Options for other output formats ------------------------------------------


htmlhelp_basename = f"{project}doc"
doc_title = f"{project} Documentation"
latex_documents = [(master_doc, f"{project}.tex", doc_title, author, "manual")]
man_pages = [(master_doc, project, doc_title, [author], 1)]
texinfo_documents = [
    (
        master_doc,
        project,
        doc_title,
        author,
        project,
        "One line description of project.",
        "Miscellaneous",
    )
]
