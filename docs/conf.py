from __future__ import annotations

import sys
from datetime import datetime
from functools import partial
from importlib import metadata
from pathlib import Path
from typing import TYPE_CHECKING

import scverse_doc.config
from docutils import nodes
from scverse_doc.registry import intersphinx

if TYPE_CHECKING:
    from sphinx.application import Sphinx

HERE = Path(__file__).parent
_extension_dir = HERE / "extensions"
acc_schema = HERE.parent / "src/anndata/acc/acc-schema-v1.json"
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
html_extra_path = [str(acc_schema)]
exclude_patterns = [
    *scverse_doc.config.DEFAULTS.exclude_patterns,
    "tutorials/notebooks/*.rst",
    # exclude all 0.x.y.md files, but not index.md
    "release-notes/[!i]*.md",
    "news.md",  # is `include`d into index.md
]
pygments_style = "sphinx"

extensions = [
    "scverse_doc",  # sets up napoleon, so needs to be before sphinx_autodoc_typehints
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx_autodoc_typehints",
    "sphinx_issues",
    "scverse_misc.sphinx_ext",
    "sphinx.ext.linkcode",
    "IPython.sphinxext.ipython_console_highlighting",
    *(p.stem for p in _extension_dir.glob("*.py")),
]
myst_heading_anchors = 3

# Generate the API documentation when building
autodoc_member_order = "bysource"
autodoc_mock_imports = ["torch"]
# autodoc_default_options = {}
issues_github_path = "scverse/anndata"
source_repository = "https://github.com/scverse/anndata"
napoleon_custom_sections = [("Params", "Parameters")]
always_use_bars_union = True  # use `|`, not `Union` in types even when on Python ≤3.14


def setup(app: Sphinx) -> None:
    app.add_generic_role("small", partial(nodes.inline, classes=["small"]))  # type: ignore[arg-type]
    app.add_generic_role("smaller", partial(nodes.inline, classes=["smaller"]))  # type: ignore[arg-type]

    # https://github.com/mkdocstrings/python/issues/339
    from sphinx.domains import ObjType
    from sphinx.domains.python import PythonDomain

    PythonDomain.object_types.setdefault(
        "typealias", ObjType("type alias", "type", "obj")
    )

    # TODO: move to scverse-doc
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


intersphinx_mapping = intersphinx("mudata", core=False) | dict(
    array_api=("https://array-api.readthedocs.io/en/stable", None),
    awkward=("https://awkward-array.org/doc/stable", None),
    cupy=("https://docs.cupy.dev/en/stable", None),
    dask=("https://docs.dask.org/en/stable", None),
    fsspec=("https://filesystem-spec.readthedocs.io/en/stable/", None),
    h5py=("https://docs.h5py.org/en/latest", None),
    hdf5plugin=("https://hdf5plugin.readthedocs.io/en/latest", None),
    kvikio=("https://docs.rapids.ai/api/kvikio/stable/", None),
    loompy=("https://linnarssonlab.org/loompy", None),
    obstore=("https://developmentseed.org/obstore/latest/", None),
    scverse_misc=("https://scverse-misc.readthedocs.io/stable", None),
    sklearn=("https://scikit-learn.org/stable", None),
    xarray=("https://docs.xarray.dev/en/stable", None),
    zarr=("https://zarr.readthedocs.io/en/stable/", None),
    zarrs=("https://zarrs-python.readthedocs.io/en/stable/", None),
    annbatch=("https://annbatch.readthedocs.io/en/stable/", None),
)

# Fix mis-documented types. Use `anndata.utils.set_module` for ours instead.
xref_overrides = {
    #### stdlib
    "types.EllipsisType": ("py:data", "Ellipsis"),
    #### anndata
    **{
        f"anndata._core.aligned_mapping.{cls}{kind}": "collections.abc.Mapping"
        for cls in ["Layers", "AxisArrays", "PairwiseArrays"]
        for kind in ["", "View"]
    },
    # Can’t use `set_module` for `type`s. When moving out of .experimental, define in actual location.
    "anndata._types.StorageType": "anndata.experimental.StorageType",
    #### zarr
    "zarr.core.group.StoreLike": "zarr.storage.StoreLike",
    #### h5py
    "h5py._hl.group.Group": "h5py.Group",
    "h5py._hl.files.File": "h5py.File",
    "h5py._hl.dataset.AsTypeView": ("py:meth", "h5py.Dataset.astype"),
    "h5py._hl.dataset.Dataset": "h5py.Dataset",
    #### arrays
    "awkward.highlevel.Array": "ak.Array",
    "numpy.int64": ("py:attr", "numpy.int64"),
    "numpy.dtypes.StringDType": ("py:attr", "numpy.dtypes.StringDType"),
    "numpy._typing._array_like.NDArray": ("py:data", "numpy.typing.NDArray"),
    "pandas.DataFrame.iloc": ("py:attr", "pandas.DataFrame.iloc"),
    "pandas.DataFrame.loc": ("py:attr", "pandas.DataFrame.loc"),
}
# Sphinx consults this {alias → name} mapping when rendering types
# sphinx-autodoc-typehints uses when importing types to resolve them
autodoc_type_aliases = dict()
# if nothing else helps, modify `nitpick_ignore`
nitpick_ignore = [  # APIs without an intersphinx entry
    # These APIs aren’t actually documented
    ("py:class", "anndata._core.raw.Raw"),
    ("py:class", "pandas.api.typing.NAType"),
    # TODO: remove zappy support; the zappy repo is archived
    ("py:class", "anndata.compat.ZappyArray"),
    ("py:class", "anndata.utils.Default"),
    # this happens when a `type` or `class` is generic
    ("py:class", "anndata.acc.GenericAlias"),
    ("py:obj", "typing.R"),
    #### type variables
    *(("py:class", v) for v in "DIMR"),
]

# -- Social cards ---------------------------------------------------------

ogp_site_url = "https://anndata.scverse.org/"
ogp_image = f"{ogp_site_url}page/_static/img/anndata_schema.svg"

# -- Options for HTML output ----------------------------------------------

html_logo = "_static/img/anndata_schema.svg"
html_show_sphinx = False


# -- Options for other output formats ------------------------------------------


htmlhelp_basename = f"{project}doc"
doc_title = f"{project} Documentation"
latex_documents = [("index", f"{project}.tex", doc_title, author, "manual")]
man_pages = [("index", project, doc_title, [author], 1)]
texinfo_documents = [
    (
        "index",
        project,
        doc_title,
        author,
        project,
        "One line description of project.",
        "Miscellaneous",
    )
]
