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
copyright = f"{datetime.now():%Y}, scverse"
release = version = metadata.version("anndata")

# default settings
templates_path = ["_templates"]
html_static_path = ["_static"]
source_suffix = {".rst": "restructuredtext", ".md": "myst-nb"}
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
    "myst_nb",
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
    "sphinxext.opengraph",
    "scanpydoc",  # needs to be before linkcode
    "sphinx.ext.linkcode",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinx_toolbox.more_autodoc.autoprotocol",
    *(p.stem for p in _extension_dir.glob("*.py")),
]
myst_enable_extensions = [
    "html_image",  # So README.md can be used on github and sphinx docs
    "colon_fence",
    "dollarmath",
]
myst_heading_anchors = 3
nb_execution_mode = "off"

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_mock_imports = ["torch"]
# autodoc_default_flags = ['members']
issues_github_path = "scverse/anndata"
rtd_links_prefix = PurePosixPath("src")
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
typehints_defaults = "braces"
todo_include_todos = False
nitpicky = True  # Report broken links
nitpick_ignore = [  # APIs without an intersphinx entry
    # This API isn’t actually documented
    ("py:class", "anndata._core.raw.Raw"),
    # TODO: remove zappy support; the zappy repo is archived
    ("py:class", "anndata.compat.ZappyArray"),
]


def setup(app: Sphinx):
    app.add_generic_role("small", partial(nodes.inline, classes=["small"]))
    app.add_generic_role("smaller", partial(nodes.inline, classes=["smaller"]))


intersphinx_mapping = dict(
    awkward=("https://awkward-array.org/doc/stable", None),
    cupy=("https://docs.cupy.dev/en/stable", None),
    dask=("https://docs.dask.org/en/stable", None),
    h5py=("https://docs.h5py.org/en/latest", None),
    hdf5plugin=("https://hdf5plugin.readthedocs.io/en/latest", None),
    loompy=("https://linnarssonlab.org/loompy", None),
    numpy=("https://numpy.org/doc/stable", None),
    pandas=("https://pandas.pydata.org/pandas-docs/stable", None),
    python=("https://docs.python.org/3", None),
    scipy=("https://docs.scipy.org/doc/scipy", None),
    sklearn=("https://scikit-learn.org/stable", None),
    # TODO: move back to stable once `ObjectStore` is released
    zarr=("https://zarr.readthedocs.io/en/latest/", None),
    xarray=("https://docs.xarray.dev/en/stable", None),
    obstore=("https://developmentseed.org/obstore/latest/", None),
    kvikio=("https://docs.rapids.ai/api/kvikio/stable/", None),
    zarrs=("https://zarrs-python.readthedocs.io/en/stable/", None),
)

qualname_overrides = {
    "h5py._hl.group.Group": "h5py.Group",
    "h5py._hl.files.File": "h5py.File",
    "h5py._hl.dataset.Dataset": "h5py.Dataset",
    "anndata._core.anndata.AnnData": "anndata.AnnData",
    **{
        f"anndata._core.aligned_mapping.{cls}{kind}": "collections.abc.Mapping"
        for cls in ["Layers", "AxisArrays", "PairwiseArrays"]
        for kind in ["", "View"]
    },
    "anndata._types.ReadCallback": "anndata.experimental.ReadCallback",
    "anndata._types.WriteCallback": "anndata.experimental.WriteCallback",
    "anndata._types.Read": "anndata.experimental.Read",
    "anndata._types.Write": "anndata.experimental.Write",
    "zarr.core.array.Array": "zarr.Array",
    "zarr.core.group.Group": "zarr.Group",
    # Buffer is not yet exported, so the buffer class registry is the closest thing
    "zarr.core.buffer.core.Buffer": "zarr.registry.Registry",
    "zarr.storage._common.StorePath": "zarr.storage.StorePath",
    "anndata.compat.DaskArray": "dask.array.Array",
    "anndata.compat.CupyArray": "cupy.ndarray",
    "anndata.compat.CupySparseMatrix": "cupyx.scipy.sparse.spmatrix",
    "anndata.compat.XDataArray": "xarray.DataArray",
    "awkward.highlevel.Array": "ak.Array",
    "numpy.int64": ("py:attr", "numpy.int64"),
    "pandas.DataFrame.iloc": ("py:attr", "pandas.DataFrame.iloc"),
    "pandas.DataFrame.loc": ("py:attr", "pandas.DataFrame.loc"),
    # should be fixed soon: https://github.com/tox-dev/sphinx-autodoc-typehints/pull/516
    "types.EllipsisType": ("py:data", "types.EllipsisType"),
    "pathlib._local.Path": "pathlib.Path",
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
