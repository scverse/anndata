"""
JavaScript for AnnData HTML representation interactivity.

Provides:
- Section folding/unfolding
- Search/filter functionality across all levels
- Copy to clipboard
- Nested content expansion
- README modal with plain text display

The JavaScript is loaded from static/repr.js and wrapped in an IIFE
that scopes it to a specific container element.
"""

from __future__ import annotations

import json
from functools import cache
from importlib.resources import files

from .utils import get_anndata_version


def _minify_js(js: str) -> str:
    """Drop indentation, blank lines and full-line ``//`` comments.

    Line breaks are kept, so automatic semicolon insertion (the script uses no
    semicolons) and the meaning of the code are unaffected.
    """
    lines = (line.strip() for line in js.splitlines())
    return "\n".join(line for line in lines if line and not line.startswith("//"))


@cache
def _load_js_content() -> str:
    """Load main JS content from static file (cached)."""
    js = files("anndata._repr.static").joinpath("repr.js").read_text(encoding="utf-8")
    return _minify_js(js)


def get_javascript(container_id: str) -> str:
    """
    Get the JavaScript code for a specific container.

    Each rendered repr ships the full source so that any cell is
    self-sufficient (surviving deletion, reorder, or notebook reopen),
    but only the first to execute installs the ``init`` function, keyed by
    anndata version in ``window.anndataRepr``. Subsequent cells reuse it for
    their own container, while outputs from a different anndata version (e.g.
    after an upgrade and kernel restart without reloading the page) install
    and use their own.

    Parameters
    ----------
    container_id
        Unique ID for the container element

    Returns
    -------
    JavaScript code wrapped in script tags
    """
    js_content = _load_js_content()
    # json.dumps produces valid JS string literals; "</" cannot appear in them
    version = json.dumps(get_anndata_version()).replace("</", "<\\/")
    container = json.dumps(container_id).replace("</", "<\\/")
    return f"""<script>
(function() {{
const container = document.getElementById({container});
if (!container) return;
const registry = (window.anndataRepr ??= {{}});
const version = {version};
registry[version] ??= function(container) {{
{js_content}
}};
registry[version](container);
}})();
</script>"""
