"""
Minimal JavaScript-based markdown parser for README rendering.

This module previously provided an inline markdown parser. The implementation
has been moved to static/markdown-parser.js for better maintainability.

Note
----
This is a temporary solution to avoid adding markdown parsing dependencies.
In the future, this could be replaced by a proper markdown library such as:
- marked.js (client-side): https://marked.js.org/
- markdown-it (client-side): https://github.com/markdown-it/markdown-it
- Python's markdown package with pre-rendering (server-side)

The current implementation supports:
- Headers (h1-h4)
- Bold and italic text
- Inline code and code blocks
- Ordered and unordered lists
- Links (with target="_blank" for security)
- Auto-linking raw URLs (https://...)
- Blockquotes
- Tables (basic, no alignment)
- Paragraphs

Limitations:
- No nested lists
- No table column alignment
- No images
- No horizontal rules
- Limited edge case handling

To replace this with an external library:
1. Update static/markdown-parser.js with the library's code
2. Ensure the parser function is named `parseMarkdown` and accepts a string
3. Update the CSS in static/repr.css if the library produces different HTML structure
"""

from __future__ import annotations

from functools import lru_cache
from importlib.resources import files


@lru_cache(maxsize=1)
def get_markdown_parser_js() -> str:
    """
    Get the JavaScript markdown parser function.

    Returns
    -------
    str
        JavaScript code defining the `parseMarkdown(text)` function.
        The function converts markdown text to HTML.

    Note
    ----
    This minimal implementation avoids external dependencies but has limitations.
    See module docstring for details on replacing with a full-featured library.
    """
    return (
        files("anndata._repr.static")
        .joinpath("markdown-parser.js")
        .read_text(encoding="utf-8")
    )
