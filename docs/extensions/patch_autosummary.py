"""Extension to
- skip inherited methods and properties in autosummary.
- include all methods and properties in protocols
"""

from __future__ import annotations

import sys
from traceback import walk_stack
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Literal

    from sphinx.application import Sphinx
    from sphinx.ext.autodoc import Options

if sys.version_info >= (3, 13):
    from typing import get_protocol_members, is_protocol
else:
    from typing_extensions import get_protocol_members, is_protocol


def skip_member(  # noqa: PLR0917
    app: Sphinx,
    what: Literal[
        "module", "class", "exception", "function", "method", "attribute", "property"
    ],
    name: str,
    obj: object,
    skip: bool,  # noqa: FBT001
    options: Options | dict[str, object],
) -> bool | None:
    """Include protocol attributes and skip inherited members."""
    # Skip `getdoc` property
    if what == "method" and name == "getdoc":
        return True

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

    # If we’re documenting a protocol attribute, include it
    if is_protocol(parent) and name in get_protocol_members(parent):
        return False

    # Skip if it’s not a member of the parent class
    typ = parent
    while typ is not type:
        if name in typ.__dict__:
            break
        typ = type(typ)
    else:
        return True

    return None  # Leave decision up to autosummary


def setup(app: Sphinx) -> None:
    """App setup hook."""
    app.connect("autodoc-skip-member", skip_member)
