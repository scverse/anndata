"""CSS styles for AnnData HTML representation."""

from __future__ import annotations

from functools import lru_cache
from importlib.resources import files


@lru_cache(maxsize=1)
def get_css() -> str:
    """Get the complete CSS for the HTML representation."""
    css = files("anndata._repr.static").joinpath("repr.css").read_text(encoding="utf-8")
    return f"<style>\n{css}\n</style>"
