"""
Shared fixtures for repr tests.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from anndata import AnnData

# Import HTML validation utilities from separate module
from .html_validator import (
    HTMLValidator,
    StrictHTMLParser,
    validate_html5_strict,
)

# Re-export validation utilities for use by other test modules
__all__ = ["HTMLValidator", "StrictHTMLParser", "validate_html5_strict"]


def pytest_configure(config):
    """Suppress ImplicitModificationWarning for all repr tests.

    This warning is expected when AnnData transforms indices internally
    during singledispatch in functools.
    """
    import warnings

    from anndata._warnings import ImplicitModificationWarning

    warnings.filterwarnings("ignore", category=ImplicitModificationWarning)


# =============================================================================
# Optional Dependencies
# =============================================================================

try:
    import dask.array as da  # noqa: F401

    HAS_DASK = True
except ImportError:
    HAS_DASK = False

try:
    import cupy as cp  # noqa: F401

    HAS_CUPY = True
except ImportError:
    HAS_CUPY = False

try:
    import awkward as ak  # noqa: F401

    HAS_AWKWARD = True
except ImportError:
    HAS_AWKWARD = False

try:
    import xarray  # noqa: F401

    HAS_XARRAY = True
except ImportError:
    HAS_XARRAY = False


# =============================================================================
# HTML5 Validation Fixture
# =============================================================================


@pytest.fixture
def validate_html5():
    """
    Fixture for strict HTML5 validation using Nu Html Checker.

    Usage:
        def test_valid_html5(validate_html5):
            html = adata._repr_html_()
            errors = validate_html5(html)
            assert not errors, f"HTML5 validation errors: {errors}"

    Skips validation if vnu is not installed.
    Install via: brew install vnu OR pip install html5-validator
    """

    def _validate(html: str) -> list[str]:
        return validate_html5_strict(html)

    return _validate


# =============================================================================
# Optional JavaScript Validation
# =============================================================================

# Check for esprima (pure Python JS parser) availability
try:
    import esprima  # noqa: F401

    HAS_ESPRIMA = True
except ImportError:
    HAS_ESPRIMA = False


def validate_javascript_syntax(html: str) -> list[str]:
    """
    Validate JavaScript syntax in HTML script tags.

    Returns list of syntax errors.
    Requires: pip install esprima

    This is a lightweight alternative to ESLint that doesn't require Node.js.
    """
    if not HAS_ESPRIMA:
        return []

    import esprima

    errors = []
    script_pattern = r"<script[^>]*>(.*?)</script>"
    scripts = re.findall(script_pattern, html, re.DOTALL | re.I)

    for i, script in enumerate(scripts):
        if not script.strip():
            continue
        try:
            esprima.parseScript(script, tolerant=True)
        except esprima.Error as e:
            errors.append(f"Script {i + 1}: {e}")

    return errors


@pytest.fixture
def validate_js():
    """
    Fixture for JavaScript syntax validation.

    Usage:
        def test_valid_js(validate_js):
            html = adata._repr_html_()
            errors = validate_js(html)
            assert not errors, f"JavaScript errors: {errors}"

    Skips validation if esprima is not installed.
    Install via: pip install esprima
    """

    def _validate(html: str) -> list[str]:
        return validate_javascript_syntax(html)

    return _validate


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def adata():
    """Basic AnnData for testing."""
    return AnnData(
        np.random.randn(100, 50).astype(np.float32),
        obs=pd.DataFrame(
            {"batch": ["A", "B"] * 50}, index=[f"cell_{i}" for i in range(100)]
        ),
        var=pd.DataFrame(
            {"gene_name": [f"gene_{i}" for i in range(50)]},
            index=[f"gene_{i}" for i in range(50)],
        ),
    )


@pytest.fixture
def adata_full():
    """AnnData with all attributes populated."""
    n_obs, n_vars = 100, 50
    adata = AnnData(
        sp.random(n_obs, n_vars, density=0.1, format="csr", dtype=np.float32),
        obs=pd.DataFrame({
            "batch": pd.Categorical(["A", "B"] * (n_obs // 2)),
            "n_counts": np.random.randint(1000, 10000, n_obs),
            "cell_type": pd.Categorical(
                ["T", "B", "NK"] * (n_obs // 3) + ["T"] * (n_obs % 3)
            ),
        }),
        var=pd.DataFrame({
            "gene_name": [f"gene_{i}" for i in range(n_vars)],
            "highly_variable": np.random.choice([True, False], n_vars),
        }),
    )
    adata.uns["neighbors"] = {"params": {"n_neighbors": 15}}
    adata.uns["batch_colors"] = ["#FF0000", "#00FF00"]
    adata.obsm["X_pca"] = np.random.randn(n_obs, 50).astype(np.float32)
    adata.obsm["X_umap"] = np.random.randn(n_obs, 2).astype(np.float32)
    adata.varm["PCs"] = np.random.randn(n_vars, 50).astype(np.float32)
    adata.layers["raw"] = sp.random(n_obs, n_vars, density=0.1, format="csr")
    adata.obsp["distances"] = sp.random(n_obs, n_obs, density=0.01, format="csr")
    adata.varp["gene_corr"] = sp.random(n_vars, n_vars, density=0.1, format="csr")
    return adata


@pytest.fixture
def adata_with_colors():
    """AnnData with color annotations."""
    adata = AnnData(np.zeros((10, 5)))
    adata.obs["cluster"] = pd.Categorical(["A", "B", "C"] * 3 + ["A"])
    adata.uns["cluster_colors"] = ["#FF0000", "#00FF00", "#0000FF"]
    return adata


@pytest.fixture
def adata_with_nested():
    """AnnData with nested AnnData in uns."""
    inner = AnnData(np.zeros((5, 3)))
    outer = AnnData(np.zeros((10, 5)))
    outer.uns["nested_adata"] = inner
    return outer


@pytest.fixture
def adata_with_special_chars():
    """AnnData with special characters in names."""
    adata = AnnData(np.zeros((10, 5)))
    adata.obs["col<script>"] = list(range(10))  # XSS attempt
    adata.uns["key&value"] = "test"
    adata.uns["quotes\"test'"] = "value"
    return adata


@pytest.fixture
def adata_empty():
    """Empty AnnData."""
    return AnnData()


@pytest.fixture
def adata_sparse_csr():
    """AnnData with sparse CSR X."""
    return AnnData(sp.random(1000, 500, density=0.1, format="csr"))


@pytest.fixture
def adata_sparse_csc():
    """AnnData with sparse CSC X."""
    return AnnData(sp.random(1000, 500, density=0.1, format="csc"))


@pytest.fixture
def validate_html():
    """Factory fixture for creating HTMLValidator instances.

    Usage:
        def test_something(validate_html):
            html = adata._repr_html_()
            v = validate_html(html)
            v.assert_section_exists("obs")
            v.assert_text_visible("batch")
    """

    def _create_validator(html: str) -> HTMLValidator:
        return HTMLValidator(html)

    return _create_validator
