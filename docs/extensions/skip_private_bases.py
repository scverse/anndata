from __future__ import annotations

from typing import TYPE_CHECKING, Generic, get_origin

from sphinx.util.typing import ExtensionMetadata

if TYPE_CHECKING:
    from sphinx.application import Sphinx


def skip_private_bases(
    app: Sphinx, name: str, obj: type, _unused, bases: list[type]
) -> None:
    bases[:] = [
        b
        for b in bases
        if b is not object
        if get_origin(b) is not Generic
        if not b.__name__.startswith("_")
    ]


def setup(app: Sphinx) -> ExtensionMetadata:
    app.connect("autodoc-process-bases", skip_private_bases)
    return ExtensionMetadata(parallel_read_safe=True)
