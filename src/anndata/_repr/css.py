"""CSS styles for AnnData HTML representation."""

from __future__ import annotations

from functools import lru_cache
from importlib.resources import files

# Dark mode CSS variables — single source of truth.
# Substituted into /* @dark-vars */ placeholders in repr.css (Tier 1 and Tier 3).
_DARK_VARS = """\
        --anndata-bg-primary: #1e1e1e;
        --anndata-bg-secondary: #252526;
        --anndata-bg-tertiary: #2d2d2d;
        --anndata-highlight: #264f78;
        --anndata-text-primary: #e0e0e0;
        --anndata-text-secondary: #a0a0a0;
        --anndata-text-muted: #707070;
        --anndata-border-color: #404040;
        --anndata-border-light: #333333;
        --anndata-accent-color: #58a6ff;
        --anndata-warning-color: #d29922;
        --anndata-warning-bg: #3d3200;
        --anndata-error-color: #f85149;
        --anndata-error-bg: #3d1a1a;
        --anndata-success-color: #3fb950;
        --anndata-info-color: #58a6ff;
        --anndata-link-color: #58a6ff;
        --anndata-code-bg: #2d2d2d;
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
        --anndata-dtype-tpu: #67e8f9;
        --anndata-dtype-awkward: #ff9d76;
        --anndata-dtype-array-api: #e6c400;"""

# Light mode CSS variables — single source of truth.
# Substituted into /* @light-vars */ placeholder in repr.css (Tier 2).
_LIGHT_VARS = """\
        --anndata-bg-primary: #ffffff;
        --anndata-bg-secondary: #f8f9fa;
        --anndata-bg-tertiary: #e9ecef;
        --anndata-highlight: #e7f1ff;
        --anndata-text-primary: #212529;
        --anndata-text-secondary: #6c757d;
        --anndata-text-muted: #adb5bd;
        --anndata-border-color: #dee2e6;
        --anndata-border-light: #e9ecef;
        --anndata-accent-color: #0d6efd;
        --anndata-warning-color: #ffc107;
        --anndata-warning-bg: #fff3cd;
        --anndata-error-color: #dc3545;
        --anndata-error-bg: #f8d7da;
        --anndata-success-color: #198754;
        --anndata-info-color: #0dcaf0;
        --anndata-link-color: #0d6efd;
        --anndata-code-bg: #f8f9fa;
        /* Dtype colors - Light mode */
        --anndata-dtype-category: #8250df;
        --anndata-dtype-int: #0550ae;
        --anndata-dtype-float: #0550ae;
        --anndata-dtype-bool: #cf222e;
        --anndata-dtype-string: #0a3069;
        --anndata-dtype-object: #6e7781;
        --anndata-dtype-sparse: #1a7f37;
        --anndata-dtype-array: #0550ae;
        --anndata-dtype-dataframe: #8250df;
        --anndata-dtype-anndata: #cf222e;
        --anndata-dtype-unknown: #6e7781;
        --anndata-dtype-extension: #8250df;
        --anndata-dtype-dask: #fb8500;
        --anndata-dtype-gpu: #76b900;
        --anndata-dtype-tpu: #0891b2;
        --anndata-dtype-awkward: #e85d04;
        --anndata-dtype-array-api: #9a6700;"""


@lru_cache(maxsize=1)
def get_css() -> str:
    """Get the complete CSS for the HTML representation.

    Replaces ``/* @dark-vars */`` and ``/* @light-vars */`` placeholders
    in ``repr.css`` with the actual variable blocks, so dark/light variables
    are each defined once (in this module) and reused across tiers.
    """
    css = files("anndata._repr.static").joinpath("repr.css").read_text(encoding="utf-8")
    css = css.replace("/* @dark-vars */", _DARK_VARS)
    css = css.replace("/* @light-vars */", _LIGHT_VARS)
    return f"<style>\n{css}\n</style>"
