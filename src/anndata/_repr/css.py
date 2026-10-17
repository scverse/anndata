"""CSS styles for AnnData HTML representation."""

from __future__ import annotations

import re
from functools import cache
from importlib.resources import files

_CSS_COMMENT = re.compile(r"/\*.*?\*/", re.DOTALL)


def _minify_css(css: str) -> str:
    """Drop comments, indentation and blank lines.

    The stylesheet is inlined into every repr output, so this keeps notebooks
    smaller. Line breaks are kept, which keeps the result readable and safe.
    """
    lines = (line.strip() for line in _CSS_COMMENT.sub("", css).splitlines())
    return "\n".join(line for line in lines if line)


@cache
def get_css() -> str:
    """Get the complete CSS for the HTML representation.

    Dark/light theming is handled entirely in CSS via ``light-dark()``
    and ``color-scheme`` — no Python-side substitution needed.
    """
    css = files("anndata._repr.static").joinpath("repr.css").read_text(encoding="utf-8")
    return f"<style>\n{_minify_css(css)}\n</style>"
