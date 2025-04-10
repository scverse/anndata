"""See <https://github.com/sloria/sphinx-issues/issues/174>"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from sphinx_issues import IssueRole

if TYPE_CHECKING:
    from sphinx.application import Sphinx


def setup(app: Sphinx) -> None:
    IssueRole.EXTERNAL_REPO_REGEX = re.compile(
        rf"^(.+)/(.+)([{IssueRole.ELEMENT_SEPARATORS}])(\w+)$"
    )
