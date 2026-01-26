"""Sphinx extension to not skip abstract/protocol methods."""

from __future__ import annotations

import abc
import sys
from traceback import walk_stack
from typing import TYPE_CHECKING

if sys.version_info >= (3, 13):
    from typing import get_protocol_members, is_protocol
else:
    from typing_extensions import get_protocol_members, is_protocol

if TYPE_CHECKING:
    from sphinx.application import Sphinx
    from sphinx.ext.autodoc import Options
    from sphinx.ext.autodoc._property_types import _AutodocObjType


def no_skip_abc_members(  # noqa: PLR0917
    app: Sphinx,
    what: _AutodocObjType,
    name: str,
    obj: object,
    skip: bool,  # noqa: FBT001
    options: Options,
) -> bool | None:
    # Find parent class
    for frame, _ in walk_stack(None):
        if frame.f_code.co_name == "_get_members" and frame.f_code.co_filename.endswith(
            "/generate.py"
        ):
            parent = frame.f_locals["obj"]
            if not isinstance(parent, type):
                return None
            break
    else:
        return None

    # Don’t skip abstract methods or properties
    if issubclass(type(parent), abc.ABCMeta) and getattr(
        obj, "__isabstractmethod__", False
    ):
        return False

    # If we’re documenting a protocol attribute, include it
    if is_protocol(parent) and name in get_protocol_members(parent):
        return False

    return None


def setup(app: Sphinx) -> None:
    app.connect("autodoc-skip-member", no_skip_abc_members)
