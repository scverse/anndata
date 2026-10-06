"""Fix attribute autodoc in presence of `set_module`.

Autodoc uses `cls.__module__` to decide which source file to scan for attribute docstrings.
Since `set_module` points that at the public re-export module instead of the module the class body actually lives in,
autodoc scans the wrong file and finds nothing.
Work around it by merging the attribute docs found in each class’s real defining module
into the analyzer results for its public module.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from sphinx.errors import PycodeError
from sphinx.pycode import ModuleAnalyzer
from sphinx.util.typing import ExtensionMetadata

if TYPE_CHECKING:
    from sphinx.application import Sphinx


# Classes decorated with `anndata.utils.set_module` for a cleaner public API
# (e.g. so docs/reprs read `anndata.AnnData` instead of
# `anndata._core.anndata.AnnData`) that also have class-body attribute
# docstrings (plain string literals following e.g.
# `obsm: ... = AlignedMappingProperty(...)`).
_ALIASED_CLASSES = ["anndata.AnnData"]


def merge_attr_docs(app: Sphinx) -> None:
    from importlib import import_module

    for dotted_name in _ALIASED_CLASSES:
        modname, _, clsname = dotted_name.rpartition(".")
        cls = getattr(import_module(modname), clsname)
        try:
            public = ModuleAnalyzer.for_module(cls.__module__)
            private = ModuleAnalyzer.for_module(cls.__source_module__)
            public.analyze()
            private.analyze()
        except PycodeError:
            continue
        public.attr_docs.update(private.attr_docs)


def setup(app: Sphinx) -> ExtensionMetadata:
    app.connect("builder-inited", merge_attr_docs)
    return ExtensionMetadata(parallel_read_safe=True)
