from __future__ import annotations

from typing import TYPE_CHECKING

from sphinx.ext.autodoc import ObjectMember
from sphinx_toolbox.more_autodoc.autoprotocol import ProtocolDocumenter

if TYPE_CHECKING:
    from typing import Self

    from sphinx.application import Sphinx


def patch_sphinx_toolbox_autoprotocol():
    """Compat hack: https://github.com/sphinx-toolbox/sphinx-toolbox/issues/168"""

    class ObjectMemberCompat(ObjectMember):
        @classmethod
        def from_other(cls, other: ObjectMember) -> Self:
            return cls(
                other.__name__,
                other.object,
                docstring=other.docstring,
                class_=other.class_,
                skipped=other.skipped,
            )

        def __iter__(self):
            return iter([self.__name__, self.object])

    filter_orig = ProtocolDocumenter.filter_members

    def filter_members(
        self, members: list[ObjectMember], want_all: bool
    ) -> list[tuple[str, object, bool]]:
        member_tuples = [ObjectMemberCompat.from_other(m) for m in members]
        return filter_orig(self, member_tuples, want_all)

    ProtocolDocumenter.filter_members = filter_members


def setup(_app: Sphinx) -> None:
    patch_sphinx_toolbox_autoprotocol()
