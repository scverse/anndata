import sys
import logging
from pathlib import Path
from datetime import datetime

from sphinx.application import Sphinx

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent), str(HERE / "extensions")]
import anndata  # noqa


logger = logging.getLogger(__name__)

for generated in HERE.glob("anndata.*.rst"):
    generated.unlink()


# -- General configuration ------------------------------------------------


needs_sphinx = "1.7"  # autosummary bugfix

# General information
project = "anndata"
author = anndata.__author__
copyright = f"{datetime.now():%Y}, {author}."
version = anndata.__version__.replace(".dirty", "")
release = version

# default settings
templates_path = ["_templates"]
source_suffix = ".rst"
master_doc = "index"
default_role = "literal"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
pygments_style = "sphinx"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx_autodoc_typehints",  # needs to be after napoleon
    "scanpydoc",
    *[p.stem for p in (HERE / "extensions").glob("*.py")],
]

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = "bysource"
# autodoc_default_flags = ['members']
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
todo_include_todos = False
nitpicky = True  # Report broken links
suppress_warnings = ["ref.citation"]


def work_around_issue_6785():
    """See https://github.com/sphinx-doc/sphinx/issues/6785"""
    from docutils.parsers.rst import directives
    from sphinx.ext import autodoc
    from sphinx.domains.python import PyAttribute

    # check if the code changes on the sphinx side and we can remove this
    assert autodoc.PropertyDocumenter.directivetype == "method"
    autodoc.PropertyDocumenter.directivetype = "attribute"

    def get_signature_prefix(self, sig: str) -> str:
        # TODO: abstract attributes
        return "property " if "property" in self.options else ""

    PyAttribute.option_spec["property"] = directives.flag
    PyAttribute.get_signature_prefix = get_signature_prefix


def setup(app: Sphinx):
    work_around_issue_6785()
    # Don’t allow broken links. DO NOT CHANGE THIS LINE, fix problems instead.
    app.warningiserror = True


intersphinx_mapping = dict(
    h5py=("http://docs.h5py.org/en/latest/", None),
    loompy=("https://linnarssonlab.org/loompy/", None),
    numpy=("https://docs.scipy.org/doc/numpy/", None),
    pandas=("https://pandas.pydata.org/pandas-docs/stable/", None),
    python=("https://docs.python.org/3", None),
    scipy=("https://docs.scipy.org/doc/scipy/reference/", None),
    sklearn=("https://scikit-learn.org/stable/", None),
    zarr=("https://zarr.readthedocs.io/en/stable/", None),
    xarray=("http://xarray.pydata.org/en/stable/", None),
)
qualname_overrides = {
    "anndata._core.anndata.AnnData": "anndata.AnnData",
    # Temporarily
    "anndata._core.raw.Raw": "anndata.AnnData",
    "anndata._core.views.ArrayView": "numpy.ndarray",
    **{
        f"anndata._core.aligned_mapping.{cls}{kind}": "typing.Mapping"
        for cls in "Layers AxisArrays PairwiseArrays".split()
        for kind in ["", "View"]
    },
}

# -- Options for HTML output ----------------------------------------------


html_theme = "scanpydoc"
html_theme_options = dict(navigation_depth=4)
html_context = dict(
    display_github=True,  # Integrate GitHub
    github_user="theislab",  # Username
    github_repo="anndata",  # Repo name
    github_version="master",  # Version
    conf_py_path="/docs/",  # Path in the checkout to the docs root
)
issues_github_path = "{github_user}/{github_repo}".format_map(html_context)
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
