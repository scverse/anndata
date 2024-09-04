from __future__ import annotations

import itertools
import re
from pathlib import Path
from typing import TYPE_CHECKING

from docutils import nodes
from packaging.version import Version
from sphinx.util.docutils import SphinxDirective

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from typing import ClassVar

    from myst_parser.mdit_to_docutils.base import DocutilsRenderer
    from sphinx.application import Sphinx


FULL_VERSION_RE = re.compile(r"^(\d+)\.(\d+)\.(\d+)(?:(\.dev.*)|(rc.*))?$")


class ReleaseNotes(SphinxDirective):
    required_arguments: ClassVar = 1

    def run(self) -> Sequence[nodes.Node]:
        dir_ = Path(self.arguments[0])
        # resolve relative dir
        if not dir_.is_absolute():
            src_file = Path(self.get_source_info()[0])
            if not src_file.is_file():
                msg = f"Cannot find relative path to: {src_file}"
                raise self.error(msg)
            dir_ = src_file.parent / self.arguments[0]
        if not dir_.is_dir():
            msg = f"Not a directory: {dir_}"
            raise self.error(msg)

        versions = sorted(
            (
                (Version(f.stem), f)
                for f in dir_.iterdir()
                if FULL_VERSION_RE.match(f.stem)
            ),
            reverse=True,  # descending
        )
        version_groups = itertools.groupby(
            versions, key=lambda vf: (vf[0].major, vf[0].minor)
        )
        for (major, minor), versions in version_groups:
            self.render_version_group(major, minor, versions)
        return []

    def render_version_group(
        self, major: int, minor: int, versions: Iterable[tuple[Version, Path]]
    ) -> None:
        target = nodes.target(
            ids=[f"v{major}-{minor}"],
            names=[f"v{major}.{minor}"],
        )
        section = nodes.section(
            "",
            nodes.title("", f"Version {major}.{minor}"),
            ids=[],
            names=[f"version {major}.{minor}"],
        )
        self.state.document.note_implicit_target(section)
        self.state.document.note_explicit_target(target)
        # append target and section to parent
        self.renderer.current_node.append(target)
        self.renderer.update_section_level_state(section, 2)
        # append children to section
        with self.renderer.current_node_context(section):
            for _, p in versions:
                self.render_include(p)

    def render_include(self, path: Path) -> None:
        # hacky solution because of https://github.com/executablebooks/MyST-Parser/issues/967
        from docutils.parsers.rst.directives.misc import Include
        from myst_parser.mocking import MockIncludeDirective

        srcfile, lineno = self.get_source_info()
        parent_dir = Path(srcfile).parent

        d = MockIncludeDirective(
            renderer=self.renderer,
            name=type(self).__name__,
            klass=Include,  # type: ignore  # wrong type hint
            arguments=[str(path.relative_to(parent_dir))],
            options={},
            body=[],
            lineno=lineno,
        )
        d.run()

        # TODO: replace the above with this once the above mentioned bug is fixed
        # from sphinx.util.parsing import nested_parse_to_nodes
        # return nested_parse_to_nodes(
        #    self.state,
        #    path.read_text(),
        #    source=str(path),
        #    offset=self.content_offset,
        # )

    @property
    def renderer(self) -> DocutilsRenderer:
        return self.state._renderer


def setup(app: Sphinx) -> None:
    app.add_directive("release-notes", ReleaseNotes)
