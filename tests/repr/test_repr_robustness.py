"""
Adversarial robustness tests for the HTML repr module.

These tests verify that the repr system handles malformed, broken, and
adversarial objects gracefully without crashing. The design philosophy is
"report what is there and what could not be done" - errors should be visible
in the output, not hidden or causing crashes.

Test categories:
- XSS injection attempts (must be escaped)
- Unicode edge cases (emoji, CJK, RTL, Zalgo)
- Huge data (large strings, many categories, deep nesting)
- Broken objects (properties that raise, missing attributes)
- Type confusion and lying hasattr
- Circular references
- Thread safety

All tests use the HTMLValidator to ensure proper HTML output and error reporting.
"""

# ruff: noqa: EM101, RUF003
# EM101: Exception string literals are used intentionally in test fixtures to create
# identifiable error messages that can be verified in test assertions.
# RUF003: Unicode lookalike characters in comments are intentional - we're testing
# that the repr handles confusable characters correctly (e.g., Cyrillic 'а' vs Latin 'a').

from __future__ import annotations

import re
import threading
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from typing import Any
import pandas as pd
import pytest
import scipy.sparse as sp

import anndata as ad
from anndata import AnnData
from anndata._repr.core import render_x_entry
from anndata._repr.html import generate_repr_html
from anndata._repr.lazy import is_lazy_adata
from anndata._repr.registry import FormatterContext, formatter_registry
from anndata._repr.utils import (
    _get_categories_from_column,
    get_backing_info,
    is_backed,
    is_serializable,
    is_view,
)

# =============================================================================
# Evil object fixtures - the most absurd data and behavior
# =============================================================================


class PropertyBomb:
    """Object where every property access raises an exception."""

    @property
    def X(self) -> None:
        raise RuntimeError("X exploded")

    @property
    def obs(self) -> None:
        raise MemoryError("obs exploded")

    @property
    def is_view(self) -> None:
        raise MemoryError("is_view exploded")

    @property
    def isbacked(self) -> None:
        raise RuntimeError("isbacked exploded")

    @property
    def shape(self) -> None:
        raise TypeError("shape exploded")

    @property
    def dtype(self) -> None:
        raise TypeError("dtype exploded")

    def __len__(self) -> int:
        raise MemoryError("len exploded")


class LyingHasattr:
    """Object where hasattr returns True but getattr fails."""

    def __getattribute__(self, name: str) -> Any:
        if name in ("X", "obs", "var", "uns"):
            msg = f"Gotcha! {name} doesn't really exist"
            raise AttributeError(msg)
        return object.__getattribute__(self, name)


class BrokenRepr:
    """Object where __repr__ and __str__ crash."""

    def __repr__(self) -> str:
        msg = "repr is broken"
        raise ValueError(msg)

    def __str__(self) -> str:
        msg = "str is broken"
        raise TypeError(msg)


class RecursiveDict(dict):
    """Dict that contains itself."""

    def __init__(self) -> None:
        super().__init__()
        self["self"] = self
        self["deeper"] = {"even_deeper": self}


class BrokenCategories:
    """Object with broken categorical accessor."""

    @property
    def cat(self) -> Any:
        class FakeCat:
            @property
            def categories(self) -> None:
                raise RuntimeError("categories exploded")

        return FakeCat()


class ZalgoText:
    """Generator for Zalgo (heavily combined) text."""

    @staticmethod
    def generate(base: str = "EVIL") -> str:
        """Generate Zalgo text with many combining characters."""
        combiners = [
            "\u0300",  # grave
            "\u0301",  # acute
            "\u0302",  # circumflex
            "\u0303",  # tilde
            "\u0304",  # macron
            "\u0305",  # overline
            "\u0306",  # breve
            "\u0307",  # dot above
            "\u0308",  # diaeresis
            "\u0309",  # hook above
            "\u030a",  # ring above
            "\u030b",  # double acute
            "\u030c",  # caron
            "\u030d",  # vertical line above
            "\u030e",  # double vertical line above
            "\u030f",  # double grave
        ]
        result = ""
        for char in base:
            result += char
            # Add random combiners
            for combiner in combiners[:8]:
                result += combiner
        return result


# =============================================================================
# Tests for XSS prevention
# =============================================================================


class TestXSSPrevention:
    """Test that XSS injection attempts are properly escaped."""

    def test_script_tag_in_column_name(self, validate_html):
        """Script tags in column names must be escaped."""
        adata = AnnData(np.zeros((3, 3)))
        adata.obs['<script>alert("XSS")</script>'] = [1, 2, 3]

        html = adata._repr_html_()
        v = validate_html(html)

        # Must not contain raw script tag
        v.assert_no_raw_xss()
        v.assert_html_well_formed()
        # Escaped version should be visible
        v.assert_text_visible("&lt;script&gt;")

        # Explicit negative check: raw payload must not appear unescaped
        # (outside of legitimate script blocks)
        content_after_style = html.split("</style>")[-1]
        assert (
            'alert("XSS")'
            not in content_after_style
            .replace("&quot;", '"')
            .replace("&lt;", "<")
            .replace("&gt;", ">")
            .split("<script")[0]
        )

    def test_img_onerror_in_uns_key(self, validate_html):
        """Image onerror handlers in uns keys must be escaped."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["<img onerror=alert(1)>"] = "xss_test"

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_no_raw_xss()
        v.assert_html_well_formed()
        v.assert_section_exists("uns")

        # Explicit negative check: no raw img tag with onerror handler
        # Check for actual img tag with onerror attribute (not escaped text)
        assert not re.search(r"<img[^>]+onerror\s*=", html, re.I), (
            "Raw img onerror found - XSS!"
        )

    def test_onclick_in_var_column(self, validate_html):
        """Onclick handlers in var column names must be escaped."""
        adata = AnnData(np.zeros((3, 5)))
        adata.var['<div onclick="evil()">click</div>'] = range(5)

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_no_raw_xss()
        v.assert_html_well_formed()

        # Explicit negative check: no unescaped onclick handler
        # Check that onclick="evil()" doesn't appear as a real attribute
        # This pattern matches onclick as an actual HTML attribute
        assert not re.search(r'<[^>]+\sonclick\s*=\s*["\']?evil', html, re.I), (
            "Raw onclick handler found - XSS!"
        )

    def test_xss_in_category_values(self, validate_html):
        """XSS attempts in categorical values must be escaped in preview."""
        adata = AnnData(np.zeros((6, 3)))
        # Category VALUES (not column names) with XSS payloads
        adata.obs["evil_cats"] = pd.Categorical([
            '<script>alert("cat_xss")</script>',
            "<img onerror=alert(1)>",
            '<svg onload="evil()">',
            "normal_value",
            "javascript:alert(1)",
            "<div onclick=bad()>click</div>",
        ])

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_no_raw_xss()
        v.assert_html_well_formed()
        v.assert_section_contains_entry("obs", "evil_cats")

        # Escaped versions should appear as text (in category preview)
        assert "&lt;script&gt;" in html, "Script tag should be escaped"

        # Raw payloads must NOT appear unescaped
        assert '<script>alert("cat_xss")</script>' not in html
        assert not re.search(r"<img[^>]+onerror\s*=", html, re.I)
        assert not re.search(r"<svg[^>]+onload\s*=", html, re.I)

    def test_xss_in_dataframe_column_names(self, validate_html):
        """XSS in DataFrame column names (obsm preview) must be escaped."""
        adata = AnnData(np.zeros((5, 3)))
        # DataFrame with XSS payloads in column names (shown in obsm preview)
        evil_df = pd.DataFrame(
            {
                '<script>alert("df_col")</script>': np.random.rand(5),
                "<img onerror=alert(1)>": np.random.rand(5),
                "normal_col": np.random.rand(5),
            },
            index=adata.obs_names,
        )
        adata.obsm["X_evil"] = evil_df

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_no_raw_xss()
        v.assert_html_well_formed()
        v.assert_section_contains_entry("obsm", "X_evil")

        # Escaped versions should appear in column preview
        assert "&lt;script&gt;" in html, "Script tag in column name should be escaped"

        # Raw payloads must NOT appear
        assert '<script>alert("df_col")</script>' not in html
        assert not re.search(r"<img[^>]+onerror\s*=", html, re.I)

    def test_javascript_url_in_readme(self, validate_html):
        """JavaScript URLs in README must be escaped."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["README"] = 'Click <a href="javascript:alert(1)">here</a>'

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_no_raw_xss()
        v.assert_element_exists(".anndata-readme__icon")

        # Explicit negative check: javascript: URL must not appear as real href
        assert not re.search(r'href\s*=\s*["\']?\s*javascript:', html, re.I), (
            "Raw javascript: URL found - XSS!"
        )

    def test_markdown_link_attribute_injection(self, validate_html):
        """Markdown links with quotes in URL must be escaped to prevent XSS.

        This tests that the markdown parser properly escapes quotes in URLs.
        Attack: [text](https://evil.com" onclick="alert(1)) would break out
        of the href attribute if quotes aren't escaped.
        """
        adata = AnnData(np.zeros((3, 3)))
        # Markdown link with quote injection attempt
        adata.uns["README"] = '[Click me](https://evil.com" onclick="alert(1)" x=")'

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_element_exists(".anndata-readme__icon")

        # The README content is stored in data-readme attribute (HTML-escaped)
        # Check that quote injection is escaped in the data attribute
        readme_match = re.search(r'data-readme="([^"]*)"', html)
        assert readme_match, "README icon should have data-readme attribute"

        readme_content = readme_match.group(1)
        # Quotes should be escaped as &quot;
        assert "onclick" in readme_content, "Attack payload should be present"
        assert "&quot;" in readme_content, "Quotes should be HTML-escaped"
        # Raw unescaped quote followed by onclick should NOT appear
        assert '" onclick=' not in readme_content, (
            "Unescaped quote before onclick - XSS risk!"
        )

    def test_css_breakout_attempt(self, validate_html):
        """CSS breakout attempts must be escaped."""
        adata = AnnData(np.zeros((3, 3)))
        adata.var["</style><script>bad()</script>"] = range(3)

        html = adata._repr_html_()
        v = validate_html(html)

        # Should only have 1 closing style tag (the legitimate one)
        style_count = html.lower().count("</style>")
        assert style_count == 1, f"Found {style_count} </style> tags (expected 1)"
        v.assert_html_well_formed()
        v.assert_no_raw_xss()

        # Explicit check: the injected script must be escaped, not executable
        assert "<script>bad()</script>" not in html, "Raw script tag found - XSS!"
        # The escaped version should appear as visible text
        assert "&lt;script&gt;" in html or "&#60;script&#62;" in html

    def test_html_div_breakout_attempt(self, validate_html):
        """HTML div breakout attempts must be escaped."""
        adata = AnnData(np.zeros((3, 3)))
        adata.obs["</div></div></div><script>pwned</script>"] = [1, 2, 3]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_no_raw_xss()
        v.assert_html_well_formed()

        # Explicit negative check: injected script must be escaped
        assert "<script>pwned</script>" not in html, "Raw script tag found - XSS!"

    def test_xss_via_exception_class_name(self, validate_html):
        """XSS via malicious exception __name__ in error messages must be escaped."""
        # This tests a real vulnerability: if an object's property raises an exception
        # with a malicious __name__, it could be rendered unescaped in error messages

        class XSSException(Exception):
            pass

        XSSException.__name__ = "<img src=x onerror=alert(1)>"

        class MaliciousObject:
            @property
            def shape(self):
                raise XSSException()

        adata = AnnData(np.zeros((5, 3)))
        adata.uns["attack"] = MaliciousObject()

        with pytest.warns(UserWarning, match="Formatter.*:"):
            html = adata._repr_html_()
        v = validate_html(html)

        v.assert_no_raw_xss()
        v.assert_html_well_formed()

        # The XSS payload must be escaped
        assert "<img src=x onerror=alert(1)>" not in html, (
            "Unescaped XSS via exception name!"
        )
        # Escaped version should appear
        assert "&lt;img" in html, "Exception name should appear escaped"

    def test_xss_via_categorical_error(self, validate_html):
        """XSS via exception in categorical preview generation must be escaped."""
        # This tests the CategoricalFormatter error path specifically
        # (different from FallbackFormatter which handles uns objects)

        class XSSCatException(Exception):
            pass

        XSSCatException.__name__ = "<script>alert('cat')</script>"

        # Create a fake categorical-like object that raises during preview
        # Must have CategoricalDtype to trigger CategoricalFormatter
        class FakeCategorical:
            dtype = pd.CategoricalDtype(categories=["a", "b"])
            shape = (5,)

            @property
            def categories(self):
                raise XSSCatException("boom")

            def __len__(self):
                return 5

        adata = AnnData(np.zeros((5, 3)))
        # Put it in uns since it accepts arbitrary objects
        adata.uns["fake_cat"] = FakeCategorical()

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_no_raw_xss()
        v.assert_html_well_formed()

        # Raw script tag must not appear
        assert "<script>alert('cat')</script>" not in html, (
            "Unescaped XSS via categorical error!"
        )

    def test_xss_via_type_name(self, validate_html):
        """XSS via malicious type __name__ must be escaped."""

        class MaliciousType:
            pass

        MaliciousType.__name__ = "<script>alert('type')</script>"

        adata = AnnData(np.zeros((5, 3)))
        adata.uns["evil_type"] = MaliciousType()

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_no_raw_xss()
        v.assert_html_well_formed()

        # Raw script tag must not appear
        assert "<script>alert('type')</script>" not in html, (
            "Unescaped XSS via type name!"
        )


# =============================================================================
# Tests for Unicode edge cases
# =============================================================================


class TestUnicodeEdgeCases:
    """Test handling of Unicode edge cases."""

    def test_emoji_in_column_names(self, validate_html):
        """Emoji in column names should work."""
        adata = AnnData(np.zeros((3, 3)))
        adata.obs["emoji_\U0001f4a9_poop"] = [1, 2, 3]
        adata.obs["\U0001f600\U0001f601\U0001f602"] = [4, 5, 6]  # Multiple emoji

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("obs")
        # Emoji column should be visible (emoji preserved or the word "emoji")
        assert "emoji" in html.lower() or "\U0001f4a9" in html

    def test_cjk_characters(self, validate_html):
        """CJK characters in data should work."""
        adata = AnnData(np.zeros((4, 3)))
        adata.obs["chinese_中文"] = ["猫", "狗", "鸟", "魚"]
        adata.var["日本語"] = ["遺伝子", "発現", "解析"]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("obs")
        v.assert_section_exists("var")

    def test_rtl_override_character(self, validate_html):
        """RTL override characters should be handled."""
        adata = AnnData(np.zeros((3, 3)))
        # U+202E is right-to-left override
        adata.obs["rtl_\u202eEVIL\u202c_text"] = [1, 2, 3]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("obs")

    def test_zalgo_text(self, validate_html):
        """Zalgo (heavily combined) text should not crash."""
        adata = AnnData(np.zeros((3, 3)))
        zalgo = ZalgoText.generate("EVIL")
        adata.obs[f"zalgo_{zalgo}"] = [1, 2, 3]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("obs")

    def test_null_byte_in_string(self, validate_html):
        """Null bytes in strings should be handled."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["null_byte"] = "before\x00after"

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")

    def test_mixed_unicode_categories(self, validate_html):
        """Mixed unicode in categorical should work."""
        adata = AnnData(np.zeros((6, 3)))
        adata.obs["mixed"] = pd.Categorical([
            "English",
            "日本語",
            "العربية",
            "עברית",
            "emoji\U0001f600",
            "Ελληνικά",
        ])

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("obs", "mixed")

    def test_zero_width_characters_in_names(self, validate_html):
        """Zero-width chars create identical-looking but different columns."""
        adata = AnnData(np.zeros((3, 3)))
        # These columns look identical but are different keys
        adata.obs["gene"] = [1, 2, 3]
        adata.obs["gene\u200b"] = [4, 5, 6]  # Zero-width space
        adata.obs["gene\u200d"] = [7, 8, 9]  # Zero-width joiner

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("obs")
        # All three columns should be shown as distinct entries
        # They look the same but are different keys
        assert html.count("gene") >= 3, "All 3 'gene' columns should be shown"


# =============================================================================
# Tests for huge/large data handling
# =============================================================================


class TestHugeDataHandling:
    """Test handling of extremely large data."""

    def test_huge_categorical_truncated(self, validate_html):
        """Categoricals with many categories should be truncated."""
        adata = AnnData(np.zeros((50, 3)))
        # Use 500 categories (still triggers truncation, but faster than 10000)
        cats = [f"category_{i}" for i in range(500)]
        adata.obs["huge_cat"] = pd.Categorical(
            np.random.choice(cats[:50], size=50), categories=cats
        )

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("obs", "huge_cat")
        v.assert_truncation_indicator()
        # Should not show all 500 categories
        category_count = html.count("category_")
        assert category_count < 200, f"Too many categories shown: {category_count}"

    def test_giant_string_in_uns_truncated(self, validate_html):
        """Giant strings (100KB) in uns should be truncated."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["giant"] = "x" * 100_000  # 100KB string (faster than 1MB)

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")
        # HTML should not be bloated (base template is ~50KB CSS/JS)
        # 100KB string should be truncated to add only ~10-20KB content
        assert len(html) < 80_000, f"HTML too large: {len(html)} chars"
        # The string should be truncated (not all 100K x's present)
        assert html.count("x") < 10_000, "Giant string should be truncated"

    def test_deeply_nested_uns(self, validate_html):
        """Deeply nested structures (20 levels) should be handled."""
        adata = AnnData(np.zeros((3, 3)))
        nested: dict = {}
        current = nested
        # 20 levels is enough to test depth limiting (reduced from 100)
        for i in range(20):
            current["level"] = {"depth": i}
            current = current["level"]
        adata.uns["deep"] = nested

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")

    def test_many_uns_keys_truncated(self, validate_html):
        """Many uns keys (300) should be truncated."""
        adata = AnnData(np.zeros((3, 3)))
        # Use batch update instead of individual assignments (faster)
        adata.uns.update({f"key_{i}": i for i in range(300)})

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")
        v.assert_truncation_indicator()
        # Should not show all 300 keys
        key_count = html.count("key_")
        assert key_count < 1500, f"Way too many keys shown: {key_count}"

    def test_wide_array_in_obsm(self, validate_html):
        """Wide array (500 columns) in obsm should be handled."""
        adata = AnnData(np.zeros((10, 5)))
        # 500 columns is enough to test (reduced from 1000)
        adata.obsm["wide"] = np.random.rand(10, 500)

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("obsm")
        v.assert_section_contains_entry("obsm", "wide")


# =============================================================================
# Tests for broken/adversarial objects
# =============================================================================


class TestBrokenObjects:
    """Test handling of objects with broken attributes."""

    @pytest.fixture
    def context(self) -> FormatterContext:
        return FormatterContext()

    def test_is_view_with_exploding_property(self) -> None:
        """is_view should return False when property raises."""
        result = is_view(PropertyBomb())
        assert result is False

    def test_is_backed_with_exploding_property(self) -> None:
        """is_backed should return False when property raises."""
        result = is_backed(PropertyBomb())
        assert result is False

    def test_get_backing_info_with_exploding_property(self) -> None:
        """get_backing_info should return default when property raises."""
        result = get_backing_info(PropertyBomb())
        assert result == {"backed": False}

    def test_is_lazy_adata_with_exploding_obs(self) -> None:
        """is_lazy_adata should return False when .obs raises."""
        result = is_lazy_adata(PropertyBomb())
        assert result is False

    def test_get_categories_with_exploding_cat_accessor(self) -> None:
        """_get_categories_from_column should return [] when .cat raises."""
        with pytest.warns(UserWarning, match="Failed to extract categories"):
            result = _get_categories_from_column(BrokenCategories())
        assert result == []

    def test_is_serializable_with_circular_reference(self) -> None:
        """is_serializable should detect circular references."""
        result = is_serializable(RecursiveDict())
        assert isinstance(result, tuple)
        assert result[0] is False

    def test_render_x_entry_missing_x(self, context, validate_html) -> None:
        """render_x_entry should show error for missing X attribute."""

        class NoX:
            pass

        result = render_x_entry(NoX(), context)
        v = validate_html(result)

        v.assert_error_shown("AttributeError")

    def test_render_x_entry_x_raises(self, context, validate_html) -> None:
        """render_x_entry should show error when X property raises."""
        result = render_x_entry(PropertyBomb(), context)
        v = validate_html(result)

        v.assert_error_shown("RuntimeError")

    def test_fallback_formatter_len_raises(self, context) -> None:
        """FallbackFormatter should handle __len__ raising."""

        class LenRaises:
            def __len__(self) -> int:
                raise MemoryError("len exploded")

        output = formatter_registry.format_value(LenRaises(), context)
        # Errors are now in output.error, not output.warnings
        assert output.error is not None
        assert "len()" in output.error

    def test_fallback_formatter_shape_raises(self, context) -> None:
        """FallbackFormatter should handle .shape raising."""

        class ShapeRaises:
            @property
            def shape(self) -> None:
                raise TypeError("shape exploded")

        with pytest.warns(UserWarning, match="shape exploded"):
            output = formatter_registry.format_value(ShapeRaises(), context)
        # Errors are now in output.error, not output.warnings
        assert output.error is not None
        assert ".shape" in output.error

    def test_fallback_formatter_dtype_raises(self, context) -> None:
        """FallbackFormatter should handle .dtype raising."""

        class DtypeRaises:
            shape = (3, 3)

            @property
            def dtype(self) -> None:
                raise TypeError("dtype exploded")

        with pytest.warns(UserWarning, match="dtype exploded"):
            output = formatter_registry.format_value(DtypeRaises(), context)
        # Errors are now in output.error, not output.warnings
        assert output.error is not None
        assert ".dtype" in output.error

    def test_fallback_formatter_broken_repr(self, context) -> None:
        """FallbackFormatter should handle broken __repr__."""
        output = formatter_registry.format_value(BrokenRepr(), context)
        assert output.type_name == "BrokenRepr"

    def test_object_with_failing_repr_in_uns(self, validate_html) -> None:
        """Objects with failing __repr__ in uns should show type name."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["broken"] = BrokenRepr()

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("uns", "broken")
        v.assert_text_visible("BrokenRepr")

    def test_object_with_failing_sizeof_in_uns(self, validate_html) -> None:
        """Objects with failing __sizeof__ should still render."""

        class FailingSizeof:
            def __sizeof__(self):
                msg = "Sizeof failed"
                raise RuntimeError(msg)

        adata = AnnData(np.zeros((3, 3)))
        adata.uns["failing_size"] = FailingSizeof()

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")


# =============================================================================
# Tests for concurrent access (thread safety)
# =============================================================================


class TestThreadSafety:
    """Test that repr generation doesn't crash under concurrent access.

    Note: This tests crash resistance, not full thread safety with synchronization.
    Concurrent modification of AnnData while generating repr may raise exceptions,
    but should never cause memory corruption or segfaults.
    """

    def test_concurrent_repr_same_object(self) -> None:
        """Multiple threads generating reprs of the SAME object should work."""
        adata = ad.AnnData(
            X=np.random.rand(50, 20),
            obs=pd.DataFrame({"cat": pd.Categorical(["a", "b"] * 25)}),
            var=pd.DataFrame({"gene": [f"g{i}" for i in range(20)]}),
        )
        adata.obsm["X_pca"] = np.random.rand(50, 10)
        adata.uns["params"] = {"k": 10}

        errors: list[Exception] = []
        results: list[str] = []

        def generate_reprs() -> None:
            try:
                for _ in range(5):
                    html = generate_repr_html(adata)
                    results.append(html)
            except Exception as e:  # noqa: BLE001
                errors.append(e)

        threads = [threading.Thread(target=generate_reprs) for _ in range(4)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=30)

        assert len(errors) == 0, f"Errors in concurrent repr: {errors}"
        # All results should be valid HTML with expected content
        assert len(results) == 20  # 4 threads × 5 iterations
        for html in results:
            assert "anndata-repr" in html
            assert "50 × 20" in html or "50</span> ×" in html  # shape in header

    def test_concurrent_repr_different_objects(self) -> None:
        """Multiple threads generating reprs of different objects should work."""
        errors: list[Exception] = []

        def generate_reprs() -> None:
            try:
                adata = ad.AnnData(X=np.random.rand(10, 10))
                for _ in range(5):
                    generate_repr_html(adata)
            except Exception as e:  # noqa: BLE001
                errors.append(e)

        threads = [threading.Thread(target=generate_reprs) for _ in range(4)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=30)

        assert len(errors) == 0, f"Errors in concurrent repr: {errors}"

    def test_concurrent_repr_with_modifications(self) -> None:
        """Concurrent repr generation while modifying should not crash."""
        adata = ad.AnnData(X=np.random.rand(10, 10))
        errors: list[Exception] = []
        successful_reprs: list[str] = []
        stop_flag = threading.Event()

        def modify_adata() -> None:
            """Modify AnnData while repr is being generated."""
            i = 0
            # Limit iterations to avoid infinite loop if stop_flag never set
            while not stop_flag.is_set() and i < 100:
                try:
                    adata.obs[f"col_{i % 5}"] = np.random.rand(10)
                    i += 1
                except Exception:  # noqa: BLE001
                    pass  # Expected during concurrent access

        def generate_reprs() -> None:
            for _ in range(2):  # Reduced from 5 to 2
                try:
                    html = generate_repr_html(adata)
                    if html:
                        successful_reprs.append(html)
                except Exception as e:  # noqa: BLE001
                    errors.append(e)

        modifier = threading.Thread(target=modify_adata)
        generators = [threading.Thread(target=generate_reprs) for _ in range(2)]

        modifier.start()
        for g in generators:
            g.start()
        for g in generators:
            g.join(timeout=5)  # Reduced from 30 to 5
        stop_flag.set()
        modifier.join(timeout=2)

        # Critical errors that indicate memory corruption or crashes - must never happen
        critical_errors = [
            e
            for e in errors
            if isinstance(e, (SystemError, SegmentationError, MemoryError))
        ]
        assert not critical_errors, (
            f"Critical errors during concurrent repr: {critical_errors}"
        )

        # At least some reprs should succeed even with concurrent modification
        assert len(successful_reprs) > 0, (
            f"No successful reprs generated. Errors: {errors}"
        )

        # Successful reprs should be valid HTML (not corrupted)
        for html in successful_reprs:
            # HTML can start with <style> or <div> depending on structure
            assert html.startswith(("<style", "<div")), (
                f"Corrupted HTML start: {html[:100]}"
            )
            assert "</div>" in html, "HTML missing closing div"
            assert "anndata-repr" in html, "Missing anndata-repr class"


class SegmentationError(Exception):
    """Placeholder for segfault (should never happen)."""


# =============================================================================
# Ultimate evil AnnData - combines ALL adversarial scenarios
# =============================================================================


class TestUltimateEvilAnnData:
    """
    The most absurd, evil AnnData combining ALL attack vectors.

    This AnnData has:
    - XSS injection attempts in column names
    - Unicode bombs (emojis, Chinese, RTL override, Zalgo)
    - Null bytes in strings
    - Giant strings (25KB in uns)
    - Huge categoricals (1000 categories)
    - Deeply nested structures (6 levels)
    - HTML/CSS breakout attempts
    - Many uns keys (100 items)

    The repr MUST:
    1. NOT crash
    2. Escape all XSS payloads
    3. Show errors visibly for things that couldn't be shown
    4. Produce valid, well-formed HTML
    5. Truncate large data appropriately
    """

    @pytest.fixture
    def evil_adata(self) -> AnnData:
        """Create the ultimate evil AnnData."""
        adata = AnnData(X=np.random.rand(50, 30).astype(np.float32))

        # XSS in obs column names
        adata.obs["normal_column"] = np.random.choice(["A", "B", "C"], size=50)
        adata.obs['<script>alert("XSS")</script>'] = np.random.randint(0, 10, size=50)
        adata.obs['<img onerror="alert(1)">'] = np.random.rand(50)
        adata.obs['onclick="evil()"'] = np.random.rand(50)

        # Unicode edge cases
        adata.obs["emoji_\U0001f4a9_poop"] = np.random.rand(50)
        adata.obs["中文_chinese"] = pd.Categorical(
            np.random.choice(["猫", "狗", "鸟", "魚"], size=50)
        )
        adata.obs["rtl_\u202eEVIL\u202c_text"] = np.random.rand(50)
        adata.obs[f"zalgo_{ZalgoText.generate('X')}"] = np.random.rand(50)

        # HTML/CSS breakout in var
        adata.var["gene_name"] = [f"gene_{i}" for i in range(30)]
        adata.var["</style><script>bad()</script>"] = np.random.rand(30)
        adata.var["</div></div></div>breakout"] = np.random.rand(30)

        # uns with various evil content
        adata.uns["normal"] = {"key": "value", "number": 42}
        adata.uns["nested_deep"] = {
            "l1": {"l2": {"l3": {"l4": {"l5": {"l6": "deep"}}}}}
        }
        adata.uns["giant_string"] = "x" * 10000  # Reduced from 25000
        adata.uns["<script>uns_xss</script>"] = "xss_attempt"
        adata.uns["null_byte"] = "before\x00after"

        # Many uns keys - use batch update (faster than loop)
        adata.uns.update({f"spam_key_{i}": {"value": i} for i in range(50)})

        # Categorical with many categories (reduced from 1000 to 300)
        big_cats = [f"category_{i}" for i in range(300)]
        adata.obs["huge_categorical"] = pd.Categorical(
            np.random.choice(big_cats[:50], size=50), categories=big_cats
        )

        # obsm/varm/layers
        adata.obsm["X_pca"] = np.random.rand(50, 10)
        adata.obsm["X_umap"] = np.random.rand(50, 2)
        adata.layers["raw"] = sp.random(50, 30, density=0.1, format="csr")
        adata.obsp["distances"] = sp.random(50, 50, density=0.05, format="csr")

        return adata

    def test_evil_adata_renders_without_crash(self, evil_adata, validate_html) -> None:
        """The evil AnnData should render without crashing."""
        html = evil_adata._repr_html_()
        v = validate_html(html)

        v.assert_element_exists(".anndata-repr")
        v.assert_shape_displayed(50, 30)

    def test_evil_adata_html_well_formed(self, evil_adata, validate_html) -> None:
        """The evil AnnData should produce well-formed HTML."""
        html = evil_adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()

    def test_evil_adata_xss_escaped(self, evil_adata, validate_html) -> None:
        """All XSS payloads should be escaped."""
        html = evil_adata._repr_html_()
        v = validate_html(html)

        v.assert_no_raw_xss()

        # Verify script tags don't appear raw
        assert '<script>alert("XSS")</script>' not in html
        assert "<script>alert" not in html
        assert "<img onerror" not in html

    def test_evil_adata_all_sections_present(self, evil_adata, validate_html) -> None:
        """All standard sections should be present."""
        html = evil_adata._repr_html_()
        v = validate_html(html)

        v.assert_section_exists("obs")
        v.assert_section_exists("var")
        v.assert_section_exists("uns")
        v.assert_section_exists("obsm")
        v.assert_section_exists("layers")
        v.assert_section_exists("obsp")

    def test_evil_adata_truncation_applied(self, evil_adata, validate_html) -> None:
        """Large data should be truncated."""
        html = evil_adata._repr_html_()
        v = validate_html(html)

        # Huge categorical should be truncated
        v.assert_truncation_indicator()

        # 300 categories should NOT all be shown (truncated to ~100)
        # (Note: categories may appear in multiple places like tooltips)
        category_count = html.count("category_")
        assert category_count < 300, f"Too many categories: {category_count}"

        # 100 spam keys - may appear multiple times in HTML
        spam_count = html.count("spam_key_")
        assert spam_count < 1000, f"Way too many spam keys: {spam_count}"

        # Giant string (25KB) should not make HTML huge
        assert len(html) < 500000, f"HTML too large: {len(html)}"

    def test_evil_adata_unicode_survives(self, evil_adata, validate_html) -> None:
        """Unicode should be preserved or escaped but not crash."""
        html = evil_adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # At least check we have obs section with entries
        v.assert_section_exists("obs")

    def test_evil_adata_css_breakout_prevented(self, evil_adata, validate_html) -> None:
        """CSS breakout attempts should be prevented."""
        html = evil_adata._repr_html_()

        # Should only have 1 legitimate </style> tag
        style_count = html.lower().count("</style>")
        assert style_count == 1, f"Found {style_count} </style> tags"


# =============================================================================
# Tests for arbitrary object types
# =============================================================================


class TestArbitraryObjects:
    """Test that functions handle completely arbitrary objects."""

    @pytest.fixture
    def context(self) -> FormatterContext:
        return FormatterContext()

    @pytest.mark.parametrize(
        "obj",
        [
            None,
            42,
            3.14,
            True,
            "string",
            b"bytes",
            [],
            [1, 2, 3],
            {},
            {"a": 1},
            (),
            (1, 2),
            set(),
            {1, 2},
            lambda x: x,
            type("Empty", (), {})(),
        ],
    )
    def test_is_view_arbitrary(self, obj: Any) -> None:
        """is_view should handle arbitrary objects without crash."""
        result = is_view(obj)
        assert result is False

    @pytest.mark.parametrize(
        "obj",
        [None, 42, "string", [], {}, type("Empty", (), {})()],
    )
    def test_is_backed_arbitrary(self, obj: Any) -> None:
        """is_backed should handle arbitrary objects without crash."""
        result = is_backed(obj)
        assert result is False

    @pytest.mark.parametrize(
        "obj",
        [None, 42, "string", [], {}, type("Empty", (), {})()],
    )
    def test_is_lazy_adata_arbitrary(self, obj: Any) -> None:
        """is_lazy_adata should handle arbitrary objects without crash."""
        result = is_lazy_adata(obj)
        assert result is False

    @pytest.mark.parametrize(
        "obj",
        [None, 42, "string", [], {}, type("Empty", (), {})()],
    )
    def test_get_categories_arbitrary(self, obj: Any) -> None:
        """_get_categories_from_column should return [] for arbitrary objects."""
        result = _get_categories_from_column(obj)
        assert result == []

    @pytest.mark.parametrize(
        "obj",
        [
            None,
            42,
            "string",
            [],
            {},
            np.array([1, 2, 3]),
            pd.DataFrame({"a": [1, 2, 3]}),
        ],
    )
    def test_format_value_arbitrary(self, obj: Any, context: FormatterContext) -> None:
        """format_value should handle arbitrary objects without crash."""
        output = formatter_registry.format_value(obj, context)
        assert output.type_name is not None


# =============================================================================
# Tests combining errors with real AnnData
# =============================================================================


class TestRealAnnDataWithErrors:
    """Test real AnnData objects with various problematic data."""

    def test_circular_reference_in_uns(self, validate_html) -> None:
        """Circular references in uns should not crash."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["circular"] = {"ref": adata.uns}

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")

    def test_anndata_self_reference_in_uns(self, validate_html) -> None:
        """AnnData that contains itself in uns should not crash."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["self"] = adata  # AnnData containing itself

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")
        # Should show nested AnnData info
        assert "AnnData" in html

    def test_none_values_in_obs(self, validate_html) -> None:
        """None values in obs columns should be handled."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["with_none"] = [None, "a", None, "b", None]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("obs", "with_none")

    def test_nan_and_inf_in_obs(self, validate_html) -> None:
        """NaN and inf values in obs columns should be handled."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["with_nan"] = [np.nan, 1.0, np.nan, 2.0, np.nan]
        adata.obs["with_inf"] = [np.inf, -np.inf, 0, 1, 2]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("obs", "with_nan")
        v.assert_section_contains_entry("obs", "with_inf")

    def test_empty_categorical(self, validate_html) -> None:
        """Empty categorical (all None values) should be handled."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["empty_cat"] = pd.Categorical([None] * 5, categories=["a", "b", "c"])

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("obs", "empty_cat")

    def test_zero_size_array_in_obsm(self, validate_html) -> None:
        """Zero-size arrays in obsm should be handled."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obsm["empty"] = np.zeros((5, 0))

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("obsm", "empty")

    def test_special_chars_in_keys(self, validate_html) -> None:
        """Special characters in keys should be escaped."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["key<with>special&chars"] = "value"
        adata.uns["quotes\"and'apostrophes"] = "value"

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")

    def test_mixed_type_list_in_uns(self, validate_html) -> None:
        """Mixed type lists in uns should be handled."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["mixed"] = [1, "string", None, 3.14, True, [1, 2], {"a": 1}]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("uns", "mixed")


# =============================================================================
# Tests for error visibility (crashing objects should show error messages)
# =============================================================================


class TestErrorVisibility:
    """Test that errors from crashing objects are visible in the repr.

    Error messages should be visible in the HTML output, not just in tooltips.
    This is important for users to understand why their data might not display correctly.
    """

    def test_crashing_repr_shows_error_visibly(self, validate_html) -> None:
        """Objects with crashing __repr__ should show error info visibly in the preview."""

        class ExplodingRepr:
            def __repr__(self):
                raise RuntimeError("BOOM! __repr__ exploded")

        adata = AnnData(np.zeros((3, 3)))
        adata.uns["exploding"] = ExplodingRepr()

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("uns", "exploding")
        # Type name should still be shown (as fallback)
        v.assert_text_visible("ExplodingRepr")
        # Error should be VISIBLE in the text (not just tooltip) with red styling
        v.assert_text_visible("RuntimeError")
        v.assert_text_visible("repr()")
        # Should use error text CSS class (red color)
        v.assert_element_exists(".anndata-text--error")

    def test_crashing_len_shows_error_visibly(self, validate_html) -> None:
        """Objects with crashing __len__ should show error info visibly."""

        class ExplodingLen:
            def __len__(self):
                raise MemoryError("BOOM! __len__ exploded")

        adata = AnnData(np.zeros((3, 3)))
        adata.uns["exploding_len"] = ExplodingLen()

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("uns", "exploding_len")
        # Type name should still be shown
        v.assert_text_visible("ExplodingLen")
        # Error should be VISIBLE in the text (not just tooltip)
        v.assert_text_visible("MemoryError")
        v.assert_text_visible("len()")

    def test_crashing_shape_shows_error_visibly(self, validate_html) -> None:
        """Objects with crashing .shape property should show error info visibly."""

        class ExplodingShape:
            @property
            def shape(self):
                raise TypeError("BOOM! shape exploded")

        adata = AnnData(np.zeros((3, 3)))
        adata.uns["exploding_shape"] = ExplodingShape()

        with pytest.warns(UserWarning, match="BOOM! shape exploded"):
            html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("uns", "exploding_shape")
        v.assert_text_visible("ExplodingShape")
        # Error should be VISIBLE in the text (not just tooltip)
        v.assert_text_visible("TypeError")
        v.assert_text_visible(".shape")

    def test_very_long_error_message_not_in_html(self, validate_html) -> None:
        """Very long error messages should NOT appear in HTML (only exception type)."""

        class VeryLongError:
            @property
            def shape(self):
                # Create a very long error message (2KB+)
                raise TypeError(
                    "LONG_ERROR_MSG " * 100 + "This is additional context. " * 50
                )

        adata = AnnData(np.zeros((3, 3)))
        adata.uns["long_error"] = VeryLongError()

        with pytest.warns(UserWarning, match="LONG_ERROR_MSG"):
            html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("uns", "long_error")
        v.assert_text_visible("VeryLongError")
        # Only the exception TYPE should be in HTML, not the full message
        # The error display shows ".shape raised TypeError", not the message content
        assert ".shape raised TypeError" in html, "Error indicator should be visible"
        # The full error message content should NOT be in HTML
        assert "LONG_ERROR_MSG" not in html, (
            "Full error message should not appear in HTML (only type name)"
        )

    def test_section_error_truncation_shows_ellipsis(self, validate_html) -> None:
        """Section rendering errors with long messages should show '...' truncation."""
        from anndata._repr.sections import _render_error_entry

        # Create a very long error message (>200 chars)
        long_error = "X" * 300

        html = _render_error_entry("test_section", long_error)

        # Should be truncated with "..."
        assert "..." in html, "Truncated error should show '...' indicator"
        # Should not contain the full 300 X's
        assert "X" * 250 not in html, "Error message should be truncated"
        # Should contain some of the error
        assert "X" * 50 in html, "Some error content should be visible"


# =============================================================================
# Tests for section truncation with many entries
# =============================================================================


class TestSectionTruncation:
    """Test that sections with many entries are truncated properly."""

    def test_varp_with_many_entries_truncated(self, validate_html) -> None:
        """varp with many entries should show truncation indicator."""
        adata = AnnData(np.zeros((30, 30)))

        # Create one valid entry, then populate internal store directly
        tiny_sparse = sp.csr_matrix(([1.0], ([0], [0])), shape=(30, 30))
        adata.varp["varp_000"] = tiny_sparse
        # Add more entries directly to bypass validation
        for i in range(1, 250):
            adata.varp._data[f"varp_{i:03d}"] = tiny_sparse

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("varp")
        # Should show truncation (default max_items=200)
        v.assert_truncation_indicator()
        # Should show "(250 items)" or similar count
        assert "250" in html or "items" in html.lower()

    def test_uns_with_many_keys_truncated(self, validate_html) -> None:
        """uns with many keys should be truncated."""
        adata = AnnData(np.zeros((3, 3)))
        # Use batch update instead of loop (much faster)
        adata.uns.update({f"key_{i:04d}": i for i in range(300)})

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")
        v.assert_truncation_indicator()
        # Should not show all 300 keys (allowing 4x for key appearing in multiple places:
        # data-key, data-copy, visible text, tooltip)
        key_count = html.count("key_")
        assert key_count < 1200, f"Too many keys in HTML: {key_count}"
        # Also verify that key_0299 (the last one) is NOT shown since it's beyond max_items=200
        assert "key_0299" not in html, "Last key should not be shown (truncation)"


# =============================================================================
# Tests for nested object visibility
# =============================================================================


class TestNestedObjectVisibility:
    """Test that nested/deeply nested objects are visible in the repr."""

    def test_nested_dict_visible(self, validate_html) -> None:
        """Nested dicts in uns should be visible."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["nested"] = {
            "level1": {
                "level2": {
                    "level3": "deep_value",
                }
            }
        }

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("uns", "nested")
        # The nested structure should be visible (type dict shown)
        assert "dict" in html.lower()

    def test_multiple_references_to_same_anndata(self, validate_html) -> None:
        """Multiple references to same AnnData should work."""
        adata = AnnData(np.zeros((3, 3)))
        shared = AnnData(np.zeros((2, 2)))
        # Same AnnData referenced multiple times
        adata.uns["ref1"] = shared
        adata.uns["ref2"] = shared
        adata.uns["nested"] = {"ref3": shared}

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")


# =============================================================================
# Additional crashing object tests
# =============================================================================


class LyingObject:
    """Object where shape/dtype/len/str all raise via explicit properties.

    Note: Using explicit properties because Python's special method lookup
    bypasses __getattr__. This simulates objects that claim to have properties
    but crash when accessed.
    """

    @property
    def shape(self):
        raise AttributeError("I have no shape")

    @property
    def dtype(self):
        raise AttributeError("I have no dtype")

    def __len__(self):
        raise AttributeError("I have no length")

    def __repr__(self):
        return "LyingObject(all properties lie)"

    def __str__(self):
        raise AttributeError("I have no str")


# =============================================================================
# Bad color array tests
# =============================================================================


class TestBadColorArrays:
    """Tests for malformed color arrays in uns."""

    def test_too_many_colors(self, validate_html) -> None:
        """More colors than categories should be handled."""
        adata = AnnData(np.zeros((10, 3)))
        adata.obs["cat"] = pd.Categorical(np.random.choice(["A", "B", "C"], size=10))
        # 6 colors for 3 categories
        adata.uns["cat_colors"] = ["red", "green", "blue", "yellow", "purple", "orange"]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("obs", "cat")

    def test_too_few_colors(self, validate_html) -> None:
        """Fewer colors than categories should be handled."""
        adata = AnnData(np.zeros((10, 3)))
        adata.obs["cat"] = pd.Categorical(
            np.random.choice(["X", "Y", "Z", "W"], size=10)
        )
        # 1 color for 4 categories
        adata.uns["cat_colors"] = ["red"]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("obs", "cat")

    def test_invalid_color_strings(self, validate_html) -> None:
        """Invalid color strings should be handled gracefully."""
        adata = AnnData(np.zeros((10, 3)))
        adata.obs["cat"] = pd.Categorical(np.random.choice(["alpha", "beta"], size=10))
        adata.uns["cat_colors"] = ["not_a_color", "also_invalid"]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("obs", "cat")

    def test_strange_color_formats(self, validate_html) -> None:
        """Various color formats (hex, rgb, rgba) should work."""
        adata = AnnData(np.zeros((10, 3)))
        adata.obs["cat"] = pd.Categorical(
            np.random.choice(["one", "two", "three"], size=10)
        )
        adata.uns["cat_colors"] = [
            "#FF0000",  # Valid hex
            "rgb(0,255,0)",  # Valid RGB
            "rgba(0,0,255,0.5)",  # Valid RGBA
        ]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("obs", "cat")

    def test_empty_colors_array(self, validate_html) -> None:
        """Empty colors array should be handled."""
        adata = AnnData(np.zeros((10, 3)))
        adata.obs["cat"] = pd.Categorical(np.random.choice(["p", "q"], size=10))
        adata.uns["cat_colors"] = []

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("obs", "cat")

    def test_css_injection_in_colors_blocked(self, validate_html) -> None:
        """CSS injection attempts in color values must be blocked."""
        adata = AnnData(np.zeros((10, 3)))
        adata.obs["batch"] = pd.Categorical(["a", "b", "c"] * 3 + ["a"])

        # Attempt CSS injection through "color" values
        adata.uns["batch_colors"] = [
            "#ff0000",  # Valid color (passes is_color_list check)
            "blue; } .anndata-section__table { display:none } .x {",  # CSS injection!
            "red; background-image: url(https://evil.com/steal)",  # Data exfil
        ]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()

        # Extract all style attributes to check for CSS injection
        style_attrs = re.findall(r'style="([^"]*)"', html)
        for style in style_attrs:
            # These patterns in a style attribute would indicate CSS injection
            assert ".anndata-section__table" not in style, (
                f"CSS injection in style: {style}"
            )
            assert "url(https://evil" not in style, f"URL injection in style: {style}"
            assert "background-image" not in style, (
                f"background-image in style: {style}"
            )

        # Valid color should still work
        assert "background:#ff0000" in html, "Valid color should be rendered"

        # Malicious strings should appear escaped in title attributes (which is safe)
        # Title format is "Invalid color: '{color}'"
        assert "Invalid color: &#x27;blue;" in html or "Invalid color: 'blue;" in html

    def test_css_injection_via_semicolon_blocked(self, validate_html) -> None:
        """Colors with semicolons (CSS property separator) must be blocked."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["cat"] = pd.Categorical(["x", "y"] * 2 + ["x"])
        adata.uns["cat_colors"] = [
            "red;padding:100px",  # Semicolon injection
            "blue;font-size:100px",  # Another injection
        ]

        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_html_well_formed()

        # Extract all style attributes to check for CSS injection
        style_attrs = re.findall(r'style="([^"]*)"', html)
        for style in style_attrs:
            # Semicolon-based injection patterns should not appear in styles
            assert "padding:100px" not in style, f"padding injection: {style}"
            assert "font-size:100px" not in style, f"font-size injection: {style}"

    def test_css_url_injection_blocked(self, validate_html) -> None:
        """CSS url() and expression() must be blocked."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["cat"] = pd.Categorical(["x", "y"] * 2 + ["x"])
        adata.uns["cat_colors"] = [
            "url(https://evil.com/track)",  # URL injection
            "expression(alert(1))",  # IE expression injection
        ]

        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_html_well_formed()

        style_attrs = re.findall(r'style="([^"]*)"', html)
        for style in style_attrs:
            assert "url(" not in style, f"url() in style: {style}"
            assert "expression(" not in style, f"expression() in style: {style}"

    def test_css_escape_sequences_blocked(self, validate_html) -> None:
        """CSS escape sequences and special chars must be blocked."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["cat"] = pd.Categorical(["a", "b"] * 2 + ["a"])
        adata.uns["cat_colors"] = [
            "red\\;background:url(evil)",  # Backslash escape
            "blue\nbackground:url(x)",  # Newline injection
        ]

        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_html_well_formed()

        style_attrs = re.findall(r'style="([^"]*)"', html)
        for style in style_attrs:
            assert "url(" not in style, f"url via escape in style: {style}"

    def test_very_long_color_rejected(self, validate_html) -> None:
        """Very long color strings must be rejected (DoS protection)."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["cat"] = pd.Categorical(["a", "b"] * 2 + ["a"])
        # Create a very long "color" string
        long_color = "red" + "x" * 1000
        adata.uns["cat_colors"] = [long_color, "blue"]

        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_html_well_formed()

        # The long string should not appear in any style attribute
        style_attrs = re.findall(r'style="([^"]*)"', html)
        for style in style_attrs:
            assert "xxxx" not in style, f"Long color in style: {style[:100]}"

    def test_hsl_var_blocked(self, validate_html) -> None:
        """CSS hsl() and var() must be blocked (not whitelisted)."""
        adata = AnnData(np.zeros((5, 3)))
        adata.obs["cat"] = pd.Categorical(["a", "b"] * 2 + ["a"])
        adata.uns["cat_colors"] = [
            "hsl(120, 100%, 50%)",  # HSL not supported
            "var(--user-color)",  # CSS variables not supported
        ]

        html = adata._repr_html_()
        v = validate_html(html)
        v.assert_html_well_formed()

        style_attrs = re.findall(r'style="([^"]*)"', html)
        for style in style_attrs:
            assert "hsl(" not in style, f"hsl() in style: {style}"
            assert "var(" not in style, f"var() in style: {style}"


# =============================================================================
# Nested AnnData with errors tests
# =============================================================================


class TestNestedAnnDataWithErrors:
    """Tests for nested AnnData objects containing broken objects."""

    def test_nested_anndata_with_broken_objects(self, validate_html) -> None:
        """Nested AnnData with broken objects should render gracefully."""
        parent = AnnData(np.zeros((5, 5)))
        child = AnnData(np.zeros((3, 3)))

        # Add broken objects to child
        child.uns["broken_repr"] = BrokenRepr()
        child.uns["lying"] = LyingObject()

        parent.uns["nested"] = child

        html = parent._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")
        # Should show the nested AnnData
        assert "AnnData" in html


# =============================================================================
# README modal with adversarial content
# =============================================================================


class TestEvilReadme:
    """Tests for README modal with malicious content."""

    def test_xss_in_readme_escaped(self, validate_html) -> None:
        """XSS attempts in README should be escaped."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["README"] = """# Evil README
<script>alert('XSS!')</script>
<img src=x onerror="alert('img')">
<svg onload="alert('svg')">
"""

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_no_raw_xss()
        v.assert_element_exists(".anndata-readme__icon")

    def test_html_injection_in_readme(self, validate_html) -> None:
        """HTML injection attempts in README should be escaped."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["README"] = """# Break Out
</div></div></div></table>
<style>body { display: none !important; }</style>
"""

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # Should only have legitimate style tag
        assert html.lower().count("</style>") == 1

    def test_large_readme_handled(self, validate_html) -> None:
        """Large README (50KB+) should be handled without bloating HTML."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["README"] = "# Big README\n" + "A" * 50000

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_element_exists(".anndata-readme__icon")
        # HTML should not be bloated excessively
        # (The 50KB content will be in data attribute but shouldn't break HTML)

    def test_unicode_chaos_in_readme(self, validate_html) -> None:
        """Unicode edge cases in README should be handled."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["README"] = """# Unicode Chaos
RTL: \u202eSIHT DAER\u202c
Null: before\x00after
Zalgo: H̸̡̪̯ͨ͊̽̅̾ḛ̫̞̜̹̙̈́͊̓̑̄̏ c̷̶̻̠̜̲̗̠̪o̶̜̹̠̺̗m̴̨̙̝̯͕̥̞̥͉̲e̴͕̫͉̮͇̣̮̼̱̤s̵̨͖̖̱̻̣͙̥̱͓
Emoji: 💀💀💀💀💀
"""

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_element_exists(".anndata-readme__icon")


class TestCSSAttacks:
    """Tests for CSS-based attacks that could break visualization."""

    def test_css_filter_attack(self, validate_html) -> None:
        """CSS filters that could blur/hide page content."""
        adata = AnnData(np.zeros((3, 3)))
        adata.obs["<style>*{filter:blur(10px)!important}</style>"] = [1, 2, 3]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_no_raw_xss()
        # The style tag in column name should be escaped
        assert "&lt;style&gt;" in html, "Style tag should be escaped to &lt;style&gt;"


class TestEncodingAttacks:
    """Tests for encoding-related attacks."""

    def test_utf7_injection(self, validate_html) -> None:
        """UTF-7 encoding that could bypass filters in some contexts."""
        adata = AnnData(np.zeros((3, 3)))
        # UTF-7 encoded <script> - only dangerous if charset is wrong
        adata.uns["utf7"] = "+ADw-script+AD4-alert(1)+ADw-/script+AD4-"

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # Should appear as literal text, not execute
        v.assert_no_raw_xss()
        # UTF-7 encoded content should appear as-is
        assert "+ADw-" in html, "UTF-7 content should appear as literal text"

    def test_byte_order_marks(self, validate_html) -> None:
        """BOM characters in unexpected places."""
        adata = AnnData(np.zeros((3, 3)))
        # UTF-8 BOM: EF BB BF
        adata.obs["\ufeffbom_prefix"] = [1, 2, 3]
        adata.obs["bom_suffix\ufeff"] = [4, 5, 6]
        adata.obs["bom\ufeffmiddle"] = [7, 8, 9]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # All three columns should be shown
        assert html.count("bom") >= 3, "All BOM columns should be visible"
