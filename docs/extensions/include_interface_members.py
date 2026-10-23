"""Sphinx extension to list ``_``-prefixed members that are part of a class’s interface.

Autosummary’s ``methods`` and ``attributes`` template variables, which our class templates
iterate over, only contain public names (and ``__init__``, which the templates drop).
That hides members that start with ``_`` but are still part of what a class offers:

- abstract members of an :class:`abc.ABC`, e.g. ``CSRDataset.__getitem__``,
- members of a :class:`typing.Protocol`, e.g. ``SupportsArrayApi.__array_namespace__``,
- special methods with a docstring of their own, e.g. ``AnnData.__getitem__``.

Inherited special methods (e.g. ``object.__repr__``) and undocumented ones stay hidden.
"""

from __future__ import annotations

import abc
from traceback import walk_stack
from typing import TYPE_CHECKING, get_protocol_members, is_protocol

from sphinx.util.typing import ExtensionMetadata

if TYPE_CHECKING:
    from sphinx.application import Sphinx
    from sphinx.ext.autodoc import Options
    from sphinx.ext.autodoc._property_types import _AutodocObjType


def include_interface_members(  # noqa: PLR0917
    app: Sphinx,
    what: _AutodocObjType,
    name: str,
    obj: object,
    skip: bool,  # noqa: FBT001
    options: Options,
) -> bool | None:
    if what in {"module", "class", "exception", "decorator"}:
        return None  # speed up by not getting parent of non-class members

    # The event doesn’t pass the class, so get it from autosummary’s member collection.
    # This also means we don’t interfere when the event is emitted from elsewhere.
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

    is_abstract = issubclass(type(parent), abc.ABCMeta) and getattr(
        obj, "__isabstractmethod__", False
    )
    is_protocol_member = is_protocol(parent) and name in get_protocol_members(parent)
    # `callable` excludes e.g. `__module__`, whose `__doc__` is `str`’s docstring
    is_documented_special = (
        name.startswith("__")
        and name.endswith("__")
        and name in vars(parent)
        and callable(obj)
        and obj.__doc__ is not None
    )
    # autosummary treats members as public when a handler returns `False`
    return False if is_abstract or is_protocol_member or is_documented_special else None


def setup(app: Sphinx) -> ExtensionMetadata:
    app.connect("autodoc-skip-member", include_interface_members)
    return ExtensionMetadata(parallel_read_safe=True)
