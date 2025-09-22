"""Sphinx extension to not skip abstract methods."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Literal

    from sphinx.application import Sphinx
    from sphinx.ext.autodoc import Options


def autodoc_skip_member(  # noqa: PLR0917
    app: Sphinx,
    what: Literal["module", "class", "exception", "function", "method", "attribute"],
    name: str,
    obj: object,
    skip: bool,  # noqa: FBT001
    options: Options,
):
    if what == "method" and getattr(obj, "__isabstractmethod__", False):
        return False
    return None


def setup(app: Sphinx):
    app.connect("autodoc-skip-member", autodoc_skip_member)
