from __future__ import annotations

import sys
from datetime import datetime
from functools import partial
from importlib import metadata
from pathlib import Path, PurePosixPath
from typing import TYPE_CHECKING

from docutils import nodes

if TYPE_CHECKING:
    from sphinx.application import Sphinx

HERE = Path(__file__).parent
_extension_dir = HERE / "extensions"
sys.path[:0] = [str(_extension_dir)]


# -- General configuration ------------------------------------------------

# General information
project = "anndata"
author = f"{project} developers"
copyright = f"{datetime.now():%Y}, {author}"
release = version = metadata.version("anndata")

# default settings
templates_path = ["_templates"]
html_static_path = ["_static"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
master_doc = "index"
default_role = "literal"
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
    "tutorials/notebooks/*.rst",
    # exclude all 0.x.y.md files, but not index.md
    "release-notes/[!i]*.md",
    "news.md",  # is `include`d into index.md
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
    "patch_sphinx_toolbox_autoprotocol",  # internal extension
    "sphinx_toolbox.more_autodoc.autoprotocol",
    *(  # other internal extensions
        p.stem
        for p in _extension_dir.glob("*.py")
        if p.stem != "patch_sphinx_toolbox_autoprotocol"
    ),
]
myst_enable_extensions = [
    "html_image",  # So README.md can be used on github and sphinx docs
]
myst_heading_anchors = 3

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = "bysource"
issues_github_path = "scverse/anndata"
rtd_links_prefix = PurePosixPath("src")
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
    ("py:class", "h5py._hl.dataset.Dataset"),
    # for experimental callback exports
    ("py:class", "anndata.compat.ZappyArray"),
    ("py:class", "anndata.compat.DaskArray"),
    ("py:class", "anndata.compat.CupyArray"),
    ("py:class", "anndata.compat.CupySparseMatrix"),
    ("py:class", "numpy.ma.core.MaskedArray"),
    ("py:class", "dask.array.core.Array"),
    ("py:class", "awkward.highlevel.Array"),
    ("py:class", "anndata._core.sparse_dataset.BaseCompressedSparseDataset"),
    ("py:obj", "numpy._typing._array_like._ScalarType_co"),
    # https://github.com/sphinx-doc/sphinx/issues/10974
    ("py:class", "numpy.int64"),
]


def setup(app: Sphinx):
    app.add_generic_role("small", partial(nodes.inline, classes=["small"]))
    app.add_generic_role("smaller", partial(nodes.inline, classes=["smaller"]))


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
    dask=("https://docs.dask.org/en/stable/", None),
)
qualname_overrides = {
    "h5py._hl.group.Group": "h5py.Group",
    "h5py._hl.files.File": "h5py.File",
    "h5py._hl.dataset.Dataset": "h5py.Dataset",
    "anndata._core.anndata.AnnData": "anndata.AnnData",
    "anndata._types.ReadCallback": "anndata.experimental.ReadCallback",
    "anndata._types.WriteCallback": "anndata.experimental.WriteCallback",
    "anndata._types.Read": "anndata.experimental.Read",
    "anndata._types.Write": "anndata.experimental.Write",
}
autodoc_type_aliases = dict(
    NDArray=":data:`~numpy.typing.NDArray`",
    AxisStorable=":data:`~anndata.typing.AxisStorable`",
    **{
        f"{v}variantRWAble": ":data:`~anndata.typing.RWAble`"
        for v in ["In", "Co", "Contra"]
    },
)

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
