from __future__ import annotations

import re
from string import ascii_letters

import numpy as np
import pandas as pd
import pytest

import anndata as ad

ADATA_ATTRS = ("obs", "var", "varm", "obsm", "layers", "obsp", "varp", "uns")


@pytest.fixture
def adata():
    return ad.AnnData(
        np.zeros((20, 10)),
        obs=pd.DataFrame(
            dict(obs_key=list(ascii_letters[:20])),
            index=[f"cell{i}" for i in range(20)],
        ),
        var=pd.DataFrame(
            dict(var_key=np.arange(10)), index=[f"gene{i}" for i in range(10)]
        ),
        varm=dict(varm_key=np.zeros((10, 20))),
        obsm=dict(obsm_key=np.zeros((20, 20))),
        layers=dict(layers_key=np.zeros((20, 10))),
        obsp=dict(obsp_key=np.zeros((20, 20))),
        varp=dict(varp_key=np.zeros((10, 10))),
        uns=dict(uns_key=dict(zip("abc", range(3), strict=True))),
    )


@pytest.fixture(params=ADATA_ATTRS)
def adata_attr(request):
    return request.param


def test_anndata_repr(adata):
    assert f"{adata.n_obs} × {adata.n_vars}" in repr(adata)

    for idxr in [
        (slice(10, 20), 9),
        (12, 9),
        (["cell1", "cell2"], slice(10, 15)),
    ]:
        v = adata[idxr]
        v_repr = repr(v)
        assert f"{v.n_obs} × {v.n_vars}" in v_repr
        assert "View of" in v_repr
        for attr in ADATA_ATTRS:
            assert re.search(
                rf"^\s+{attr}:[^$]*{attr}_key.*$", v_repr, flags=re.MULTILINE
            )


def test_removal(adata, adata_attr):
    attr = adata_attr
    assert re.search(rf"^\s+{attr}:.*$", repr(adata), flags=re.MULTILINE)
    delattr(adata, attr)
    assert re.search(rf"^\s+{attr}:.*$", repr(adata), flags=re.MULTILINE) is None


# ---------------------------------------------------------------------------
# HTML repr tests
# ---------------------------------------------------------------------------


def test_repr_html_basic(adata):
    """HTML repr returns a string with expected sections."""
    html = adata._repr_html_()
    assert html is not None
    assert "anndata-repr" in html
    for section in ("obs", "var", "obsm", "varm", "obsp", "varp", "layers", "uns"):
        assert f'data-section="{section}"' in html


def test_repr_html_section_formatter_get_entries(adata):
    """Custom SectionFormatter with get_entries produces a foldable section."""
    from anndata._repr import (
        FormattedEntry,
        FormattedOutput,
        SectionFormatter,
        register_formatter,
    )

    @register_formatter
    class TestEntriesSection(SectionFormatter):
        section_names = ("_test_entries",)

        @property
        def section_name(self):
            return "_test_entries"

        @property
        def after_section(self):
            return "obs"

        def get_entries(self, obj, context):
            return [
                FormattedEntry(
                    key="test_key",
                    output=FormattedOutput(type_name="str", preview="hello"),
                )
            ]

    try:
        html = adata._repr_html_()
        assert "test_key" in html
        assert "hello" in html
    finally:
        from anndata._repr.registry import formatter_registry

        formatter_registry._section_formatters.pop("_test_entries", None)


def test_repr_html_section_formatter_render_html(adata):
    """Custom SectionFormatter with render_html produces raw HTML (no <details>)."""
    from anndata._repr import (
        SectionFormatter,
        register_formatter,
    )

    @register_formatter
    class TestRawHtmlSection(SectionFormatter):
        section_names = ("_test_raw_html",)

        @property
        def section_name(self):
            return "_test_raw_html"

        @property
        def after_section(self):
            return "obs"

        def get_entries(self, obj, context):
            return []

        def render_html(self, obj, context):
            return '<div class="test-raw">custom_raw_content</div>'

    try:
        html = adata._repr_html_()
        assert "custom_raw_content" in html
        assert "test-raw" in html
    finally:
        from anndata._repr.registry import formatter_registry

        formatter_registry._section_formatters.pop("_test_raw_html", None)


def test_repr_html_section_formatter_render_html_escaping(adata):
    """render_html output is inserted as-is — formatter must escape values."""
    from anndata._repr import (
        SectionFormatter,
        register_formatter,
    )
    from anndata._repr.utils import escape_html

    @register_formatter
    class TestEscapedSection(SectionFormatter):
        section_names = ("_test_escaped",)

        @property
        def section_name(self):
            return "_test_escaped"

        @property
        def after_section(self):
            return "obs"

        def get_entries(self, obj, context):
            return []

        def render_html(self, obj, context):
            # Properly escaped — this is the formatter's responsibility
            val = '<script>alert("xss")</script>'
            return f"<div>{escape_html(val)}</div>"

    try:
        html = adata._repr_html_()
        # The raw script tag must NOT appear
        assert '<script>alert("xss")</script>' not in html
        # The escaped version should appear
        assert "&lt;script&gt;" in html
    finally:
        from anndata._repr.registry import formatter_registry

        formatter_registry._section_formatters.pop("_test_escaped", None)


def test_repr_html_section_formatter_render_html_crash_fallback(adata):
    """Crashing render_html falls back to get_entries."""
    import warnings

    from anndata._repr import (
        FormattedEntry,
        FormattedOutput,
        SectionFormatter,
        register_formatter,
    )

    @register_formatter
    class TestCrashFallbackSection(SectionFormatter):
        section_names = ("_test_crash_fb",)

        @property
        def section_name(self):
            return "_test_crash_fb"

        @property
        def after_section(self):
            return "obs"

        def get_entries(self, obj, context):
            return [
                FormattedEntry(
                    key="fallback_key",
                    output=FormattedOutput(type_name="str", preview="fallback_value"),
                )
            ]

        def render_html(self, obj, context):
            raise RuntimeError("intentional crash")

    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            html = adata._repr_html_()
            # Should still produce valid HTML
            assert html is not None
            assert "anndata-repr" in html
            # Should warn about the crash
            crash_warnings = [x for x in w if "render_html failed" in str(x.message)]
            assert len(crash_warnings) == 1
            # Should fall back to get_entries
            assert "fallback_key" in html
            assert "fallback_value" in html
    finally:
        from anndata._repr.registry import formatter_registry

        formatter_registry._section_formatters.pop("_test_crash_fb", None)
