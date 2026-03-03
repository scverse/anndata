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

from functools import lru_cache
from importlib.resources import files


@lru_cache(maxsize=1)
def _load_js_content() -> str:
    """Load main JS content from static file (cached)."""
    return files("anndata._repr.static").joinpath("repr.js").read_text(encoding="utf-8")


def get_javascript(container_id: str) -> str:
    """
    Get the JavaScript code for a specific container.

    Parameters
    ----------
    container_id
        Unique ID for the container element

    Returns
    -------
    JavaScript code wrapped in script tags
    """
    js_content = _load_js_content()
    return f"""<script>
(function() {{
    // Scoped to avoid conflicts with multiple repr instances
    const container = document.getElementById('{container_id}');
    if (!container) return;

    {js_content}
}})();
</script>"""
