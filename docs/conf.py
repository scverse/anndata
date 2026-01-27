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
always_use_bars_union = True  # use `|`, not `Union` in types even when on Python ≤3.14
todo_include_todos = False


def setup(app: Sphinx):
    app.add_generic_role("small", partial(nodes.inline, classes=["small"]))
    app.add_generic_role("smaller", partial(nodes.inline, classes=["smaller"]))

    # TODO: move to scanpydoc
    if TYPE_CHECKING:
        from docutils.nodes import TextElement, reference
        from sphinx.addnodes import pending_xref
        from sphinx.environment import BuildEnvironment

    def res(
        app: Sphinx, env: BuildEnvironment, node: pending_xref, contnode: TextElement
    ) -> reference | None:
        return env.domains["py"].resolve_xref(
            env,
            node["refdoc"],
            app.builder,
            node["reftype"],
            node["reftarget"],
            node,
            contnode,
        )

    app.connect("missing-reference", res, priority=502)


intersphinx_mapping = dict(
    awkward=("https://awkward-array.org/doc/stable", None),
    cupy=("https://docs.cupy.dev/en/stable", None),
    dask=("https://docs.dask.org/en/stable", None),
    fsspec=("https://filesystem-spec.readthedocs.io/en/stable/", None),
    h5py=("https://docs.h5py.org/en/latest", None),
    hdf5plugin=("https://hdf5plugin.readthedocs.io/en/latest", None),
    kvikio=("https://docs.rapids.ai/api/kvikio/stable/", None),
    loompy=("https://linnarssonlab.org/loompy", None),
    numpy=("https://numpy.org/doc/stable", None),
    obstore=("https://developmentseed.org/obstore/latest/", None),
    pandas=("https://pandas.pydata.org/pandas-docs/stable", None),
    python=("https://docs.python.org/3", None),
    scipy=("https://docs.scipy.org/doc/scipy", None),
    sklearn=("https://scikit-learn.org/stable", None),
    xarray=("https://docs.xarray.dev/en/stable", None),
    zarr=("https://zarr.readthedocs.io/en/stable/", None),
    zarrs=("https://zarrs-python.readthedocs.io/en/stable/", None),
)

# Fix mis-documented types. Use `anndata.utils.set_module` for ours instead.
qualname_overrides = {
    #### stdlib
    "types.EllipsisType": ("py:data", "Ellipsis"),
    #### anndata
    **{
        f"anndata._core.aligned_mapping.{cls}{kind}": "collections.abc.Mapping"
        for cls in ["Layers", "AxisArrays", "PairwiseArrays"]
        for kind in ["", "View"]
    },
    # Can’t use `set_module` for `type`s. When moving out of .experimental, define in actual location.
    "anndata.compat.Index": "anndata.typing.Index",
    "anndata._types.StorageType": "anndata.experimental.StorageType",
    #### h5py
    "h5py._hl.group.Group": "h5py.Group",
    "h5py._hl.files.File": "h5py.File",
    "h5py._hl.dataset.Dataset": "h5py.Dataset",
    #### arrays
    "awkward.highlevel.Array": "ak.Array",
    "numpy.int64": ("py:attr", "numpy.int64"),
    "numpy.dtypes.StringDType": ("py:attr", "numpy.dtypes.StringDType"),
    "pandas.DataFrame.iloc": ("py:attr", "pandas.DataFrame.iloc"),
    "pandas.DataFrame.loc": ("py:attr", "pandas.DataFrame.loc"),
}
# Sphinx consults this {alias → name} mapping when rendering types
# sphinx-autodoc-typehints uses when importing types to resolve them
autodoc_type_aliases = dict()
# if nothing else helps, modify `nitpick_ignore`
nitpicky = True  # Report broken links, this stays on
nitpick_ignore = [  # APIs without an intersphinx entry
    # These APIs aren’t actually documented
    ("py:class", "anndata._core.raw.Raw"),
    ("py:class", "pandas.api.typing.NAType"),
    # TODO: remove zappy support; the zappy repo is archived
    ("py:class", "anndata.compat.ZappyArray"),
]

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
