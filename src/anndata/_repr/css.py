"""
CSS styles for AnnData HTML representation.

Provides inline CSS with:
- Light and dark mode support
- Jupyter notebook theme detection
- Responsive design
- Accessible color contrasts

The CSS is loaded from static/repr.css and processed to substitute
dark mode variables into the appropriate places in the template.
"""

from __future__ import annotations

from functools import lru_cache
from importlib.resources import files

# Dark mode CSS variables (shared by media query and theme detection selectors)
# These values are substituted into the CSS template at the /* __DARK_MODE_VARS__ */ placeholder.
# Using a CSS comment as placeholder keeps the CSS valid for linters/formatters.
# The same values are used in:
# 1. @media (prefers-color-scheme: dark) - OS-level dark mode
# 2. Theme detection selectors - Jupyter, VS Code, Sphinx (Furo/sphinx-book-theme), manual class
_DARK_MODE_VARS = """--anndata-bg-primary: #1e1e1e;
    --anndata-bg-secondary: #252526;
    --anndata-bg-tertiary: #2d2d2d;
    --anndata-highlight: #264f78;
    --anndata-text-primary: #e0e0e0;
    --anndata-text-secondary: #a0a0a0;
    --anndata-text-muted: #707070;
    --anndata-border-color: #404040;
    --anndata-border-light: #333333;
    --anndata-accent-color: #58a6ff;
    --anndata-accent-hover: #79b8ff;
    --anndata-warning-color: #d29922;
    --anndata-warning-bg: #3d3200;
    --anndata-error-color: #f85149;
    --anndata-error-bg: #3d1a1a;
    --anndata-success-color: #3fb950;
    --anndata-info-color: #58a6ff;
    --anndata-link-color: #58a6ff;
    --anndata-code-bg: #2d2d2d;
    --anndata-shadow: 0 1px 3px rgba(0, 0, 0, 0.3);
    /* Dtype colors - Dark mode */
    --anndata-dtype-category: #d2a8ff;
    --anndata-dtype-int: #79c0ff;
    --anndata-dtype-float: #79c0ff;
    --anndata-dtype-bool: #ff7b72;
    --anndata-dtype-string: #a5d6ff;
    --anndata-dtype-sparse: #7ee787;
    --anndata-dtype-array: #79c0ff;
    --anndata-dtype-dataframe: #d2a8ff;
    --anndata-dtype-anndata: #ff7b72;
    --anndata-dtype-dask: #ffc168;
    --anndata-dtype-gpu: #a0db63;
    --anndata-dtype-awkward: #ff9d76;
    --anndata-dtype-array-api: #e6c400;
"""


@lru_cache(maxsize=1)
def _load_css_template() -> str:
    """Load CSS template from static file (cached)."""
    return (
        files("anndata._repr.static").joinpath("repr.css").read_text(encoding="utf-8")
    )


def get_css() -> str:
    """Get the complete CSS for the HTML representation."""
    # Build CSS by substituting dark mode variables into the template
    css = _load_css_template().replace("/* __DARK_MODE_VARS__ */", _DARK_MODE_VARS)
    return f"<style>\n{css}\n</style>"
