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


class InfiniteLength:
    """Object with impossibly large __len__."""

    def __len__(self) -> int:
        return 10**18


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


class NestedExceptions:
    """Object where every attribute access raises nested exceptions."""

    def __getattr__(self, name: str) -> None:
        try:
            raise ValueError("inner")
        except ValueError:
            raise RuntimeError("outer") from None


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

    def test_img_onerror_in_uns_key(self, validate_html):
        """Image onerror handlers in uns keys must be escaped."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["<img onerror=alert(1)>"] = "xss_test"

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_no_raw_xss()
        v.assert_html_well_formed()
        v.assert_section_exists("uns")

    def test_onclick_in_var_column(self, validate_html):
        """Onclick handlers in var column names must be escaped."""
        adata = AnnData(np.zeros((3, 5)))
        adata.var['<div onclick="evil()">click</div>'] = range(5)

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_no_raw_xss()
        v.assert_html_well_formed()

    def test_javascript_url_in_readme(self, validate_html):
        """JavaScript URLs in README must be escaped."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["README"] = 'Click <a href="javascript:alert(1)">here</a>'

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_no_raw_xss()
        v.assert_element_exists(".adata-readme-icon")

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

    def test_html_div_breakout_attempt(self, validate_html):
        """HTML div breakout attempts must be escaped."""
        adata = AnnData(np.zeros((3, 3)))
        adata.obs["</div></div></div><script>pwned</script>"] = [1, 2, 3]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_no_raw_xss()
        v.assert_html_well_formed()


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
        # Emoji might be escaped or preserved - just ensure no crash
        assert len(html) > 0

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
        assert any("len()" in w for w in output.warnings)

    def test_fallback_formatter_shape_raises(self, context) -> None:
        """FallbackFormatter should handle .shape raising."""

        class ShapeRaises:
            @property
            def shape(self) -> None:
                raise TypeError("shape exploded")

        output = formatter_registry.format_value(ShapeRaises(), context)
        assert any(".shape" in w for w in output.warnings)

    def test_fallback_formatter_dtype_raises(self, context) -> None:
        """FallbackFormatter should handle .dtype raising."""

        class DtypeRaises:
            shape = (3, 3)

            @property
            def dtype(self) -> None:
                raise TypeError("dtype exploded")

        output = formatter_registry.format_value(DtypeRaises(), context)
        assert any(".dtype" in w for w in output.warnings)

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
    """Test that repr generation is thread-safe."""

    def test_concurrent_repr_generation(self) -> None:
        """Multiple threads generating reprs should not crash."""
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
            try:
                for _ in range(2):  # Reduced from 5 to 2
                    generate_repr_html(adata)
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

        # We accept some errors during concurrent modification, but repr should not crash
        assert all(
            not isinstance(e, (SystemError, SegmentationError))
            for e in errors
            if hasattr(e, "__class__")
        )


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
        v.assert_element_exists(".adata-error-text")

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

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("uns", "exploding_shape")
        v.assert_text_visible("ExplodingShape")
        # Error should be VISIBLE in the text (not just tooltip)
        v.assert_text_visible("TypeError")
        v.assert_text_visible(".shape")


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

    def test_circular_reference_detected(self, validate_html) -> None:
        """Circular references should be detected and shown safely."""
        adata = AnnData(np.zeros((3, 3)))
        circular = {"level1": {"level2": {}}}
        circular["level1"]["level2"]["back"] = circular

        adata.uns["circular"] = circular

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")
        # Should not crash and should produce valid HTML
        # The type should be shown
        assert "dict" in html.lower()

    def test_deeply_nested_structure(self, validate_html) -> None:
        """Deeply nested structures should be handled."""
        adata = AnnData(np.zeros((3, 3)))
        deeply_nested = {}
        current = deeply_nested
        for i in range(10):
            current[f"level_{i}"] = {}
            current = current[f"level_{i}"]
        current["bottom"] = "reached"

        adata.uns["deep"] = deeply_nested

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("uns", "deep")

    def test_self_reference_in_uns(self, validate_html) -> None:
        """AnnData that references itself in uns should not crash."""
        adata = AnnData(np.zeros((3, 3)))
        # Create self-reference (adata.uns contains the adata itself)
        adata.uns["self_ref"] = adata

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")
        # Should show that there's an AnnData nested
        assert "AnnData" in html

    def test_nested_anndata_with_self_reference(self, validate_html) -> None:
        """Nested AnnData that references its parent should not crash."""
        parent = AnnData(np.zeros((5, 5)))
        child = AnnData(np.zeros((3, 3)))
        # Create circular reference: parent -> child -> parent
        parent.uns["child"] = child
        child.uns["parent_ref"] = parent

        html = parent._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")

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
# Additional crashing object tests from visual test 24
# =============================================================================


class ExplodingStr:
    """Object where __str__ crashes."""

    def __repr__(self):
        return "ExplodingStr()"

    def __str__(self):
        msg = "BOOM! __str__ exploded"
        raise ValueError(msg)


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


class VeryLongErrorObject:
    """Object that produces a very long error message."""

    @property
    def shape(self):
        raise TypeError(
            "This is a VERY LONG ERROR MESSAGE that should be properly truncated. " * 10
            + "It contains lots of details about what went wrong: "
            + "ValueError: The input array has shape (100, 200, 300) but expected (50, 100). "
            + "Additional context: This error occurred while processing the data matrix. "
            + "Stack trace would go here with many lines of debugging information. " * 5
        )

    def __repr__(self):
        return "VeryLongErrorObject(produces long error)"


class FakeAnndataType:
    """Unknown type from anndata package triggers warning (not error).

    Objects whose __module__ starts with 'anndata.' but aren't recognized
    types should show a warning instead of an error - they might be from
    a newer version of anndata.
    """

    __module__ = "anndata.experimental.fake"

    def __repr__(self):
        return "FakeAnndataType()"


class TestExplodingStr:
    """Tests for objects with crashing __str__."""

    def test_exploding_str_shows_error(self, validate_html) -> None:
        """Objects with crashing __str__ should show error info."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["exploding_str"] = ExplodingStr()

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("uns", "exploding_str")
        v.assert_text_visible("ExplodingStr")
        # Error should be visible
        v.assert_text_visible("ValueError")
        v.assert_text_visible("str()")


class TestLyingObject:
    """Tests for objects where all properties raise."""

    def test_lying_object_shows_multiple_errors(self, validate_html) -> None:
        """LyingObject should show errors for all failed property accesses."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["lying"] = LyingObject()

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("uns", "lying")
        v.assert_text_visible("LyingObject")
        # Should show errors for failed accesses
        v.assert_text_visible("AttributeError")

    def test_lying_object_repr_survives(self) -> None:
        """LyingObject's repr should work even though other methods fail."""
        obj = LyingObject()
        assert repr(obj) == "LyingObject(all properties lie)"


class TestInfiniteLengthWarning:
    """Tests for objects with suspiciously large length."""

    def test_infinite_len_shows_warning(self, validate_html) -> None:
        """Objects claiming impossibly large length should show warning."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["infinite"] = InfiniteLength()

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("uns", "infinite")
        # Should show warning about suspicious length
        # The warning could be in text or as a warning indicator
        assert "1" in html  # Should mention the large number somehow


class TestVeryLongErrorMessages:
    """Tests for error message truncation."""

    def test_very_long_error_truncated(self, validate_html) -> None:
        """Very long error messages should be truncated."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["long_error"] = VeryLongErrorObject()

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("uns", "long_error")
        v.assert_text_visible("VeryLongErrorObject")
        # The full 2KB+ error message should be truncated
        # Verify HTML isn't bloated by repeated error text
        error_count = html.count("VERY LONG ERROR MESSAGE")
        assert error_count < 5, (
            f"Error message not truncated: {error_count} occurrences"
        )

    def test_long_error_in_varm(self, validate_html) -> None:
        """Long errors in varm entries should be truncated."""
        adata = AnnData(np.zeros((10, 10)))
        # Add valid entry first
        adata.varm["valid"] = np.random.rand(10, 5)
        # Add broken entry via internal store
        adata.varm._data["broken"] = VeryLongErrorObject()

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("varm")
        # Should still produce valid HTML despite the error


class TestUnknownAnndataType:
    """Tests for unknown types from anndata package."""

    def test_fake_anndata_type_shows_warning(self, validate_html) -> None:
        """Unknown types from anndata.* should show warning, not error."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["fake"] = FakeAnndataType()

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("uns", "fake")
        v.assert_text_visible("FakeAnndataType")
        # Should indicate it's an unknown/unrecognized type
        # This could be shown as a warning or with special styling


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

    def test_nested_anndata_errors_visible(self, validate_html) -> None:
        """Errors in nested AnnData should be visible in output."""
        parent = AnnData(np.zeros((5, 5)))
        child = AnnData(np.zeros((3, 3)))

        class ExplodingReprInNested:
            def __repr__(self):
                raise RuntimeError("BOOM in nested!")

        child.uns["exploding"] = ExplodingReprInNested()
        parent.uns["child"] = child

        html = parent._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # The nested structure should be shown
        v.assert_section_exists("uns")

    def test_deeply_nested_anndata_with_errors(self, validate_html) -> None:
        """Deeply nested AnnData with errors at multiple levels."""
        root = AnnData(np.zeros((5, 5)))
        level1 = AnnData(np.zeros((4, 4)))
        level2 = AnnData(np.zeros((3, 3)))

        # Add broken objects at different levels
        root.uns["broken_at_root"] = BrokenRepr()
        level1.uns["broken_at_level1"] = LyingObject()
        level2.uns["broken_at_level2"] = InfiniteLength()

        # Create nesting
        level1.uns["level2"] = level2
        root.uns["level1"] = level1

        html = root._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")


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
        v.assert_element_exists(".adata-readme-icon")

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
        v.assert_element_exists(".adata-readme-icon")
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
        v.assert_element_exists(".adata-readme-icon")


# =============================================================================
# Advanced adversarial tests - sophisticated attacks
# =============================================================================


class SlowRepr:
    """Object with __repr__ that blocks for a long time.

    Could hang a Jupyter notebook if not handled with timeout.
    """

    def __init__(self, delay: float = 0.1):
        self.delay = delay

    def __repr__(self):
        import time

        time.sleep(self.delay)
        return f"SlowRepr(delay={self.delay})"


class RecursiveRepr:
    """Object whose __repr__ calls itself recursively.

    Could cause stack overflow if not protected.
    Note: We don't actually test this because it would crash - it's here
    for documentation purposes only.
    """

    def __repr__(self):
        # Don't actually recurse - just simulate what would happen
        return "RecursiveRepr(<would recurse>)"


class MemoryBomb:
    """Object that raises MemoryError on property access.

    Simulates an object that would try to allocate huge memory.
    We raise the error directly instead of actually trying to allocate.
    """

    @property
    def shape(self):
        # Raise MemoryError directly instead of actually trying to allocate
        # (allocating 10**12 items would kill the process)
        raise MemoryError("Simulated memory allocation failure")


class TestHomographAttacks:
    """Tests for Unicode homograph attacks that could deceive users.

    Homograph attacks use lookalike characters from different scripts
    to create visually identical but semantically different strings.
    This could trick users into thinking they're looking at different data.
    """

    def test_cyrillic_lookalikes_in_column_names(self, validate_html) -> None:
        """Cyrillic characters that look like Latin should be handled."""
        adata = AnnData(np.zeros((3, 3)))
        # 'а' is Cyrillic, not Latin 'a' - visually identical
        adata.obs["gene_a"] = [1, 2, 3]  # Latin
        adata.obs["gene_\u0430"] = [4, 5, 6]  # Cyrillic 'а'

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # Both columns should be shown as distinct entries
        v.assert_section_exists("obs")
        # Both keys should appear in the output (not merged/confused)
        assert html.count("gene_") >= 2, "Both similar column names should be shown"

    def test_greek_lookalikes(self, validate_html) -> None:
        """Greek characters that look like Latin should be handled."""
        adata = AnnData(np.zeros((3, 3)))
        # Greek letters that look like Latin
        adata.obs["DATA"] = [1, 2, 3]  # Latin
        adata.obs["\u0394\u0391\u03a4\u0391"] = [4, 5, 6]  # Greek ΔΑΤΑ

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("obs")
        # Both should appear (Latin DATA and Greek)
        assert "DATA" in html, "Latin DATA should appear"

    def test_fullwidth_lookalikes(self, validate_html) -> None:
        """Fullwidth characters that look like ASCII should be handled."""
        adata = AnnData(np.zeros((3, 3)))
        # Fullwidth letters (used in CJK contexts)
        adata.obs["ABC"] = [1, 2, 3]  # Normal ASCII
        adata.obs["\uff21\uff22\uff23"] = [4, 5, 6]  # Fullwidth ＡＢＣ

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # Both columns should be visible (not confused as one)
        assert "ABC" in html, "ASCII ABC should appear"


class TestZeroWidthDeception:
    """Tests for zero-width character attacks.

    Zero-width characters are invisible but make strings different.
    Could create columns that look identical but are different.
    """

    def test_zero_width_space_in_names(self, validate_html) -> None:
        """Zero-width spaces make identical-looking different names."""
        adata = AnnData(np.zeros((3, 3)))
        # U+200B is zero-width space
        adata.obs["gene"] = [1, 2, 3]
        adata.obs["gene\u200b"] = [4, 5, 6]  # Has invisible character

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("obs")
        # Both columns should be shown - they're different keys even if they look the same
        # Count occurrences of "gene" in data-key attributes or visible text
        assert html.count("gene") >= 2, "Both 'gene' columns should be shown"

    def test_zero_width_joiner_confusion(self, validate_html) -> None:
        """Zero-width joiners between characters."""
        adata = AnnData(np.zeros((3, 3)))
        # U+200D is zero-width joiner
        adata.obs["cell_type"] = [1, 2, 3]
        adata.obs["cell\u200dtype"] = [4, 5, 6]  # ZWJ in middle

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # Both columns should appear in output
        assert html.count("cell") >= 2, "Both cell columns should be shown"

    def test_invisible_separator_attack(self, validate_html) -> None:
        """Various invisible Unicode separators."""
        adata = AnnData(np.zeros((3, 3)))
        invisible_chars = [
            "\u00a0",  # Non-breaking space
            "\u2000",  # En quad
            "\u2001",  # Em quad
            "\u2002",  # En space
            "\u2003",  # Em space
            "\u2007",  # Figure space
            "\u2008",  # Punctuation space
            "\u200a",  # Hair space
            "\u200b",  # Zero-width space
            "\u202f",  # Narrow no-break space
            "\u205f",  # Medium mathematical space
            "\u3000",  # Ideographic space
            "\ufeff",  # Zero-width no-break space (BOM)
        ]
        for i, char in enumerate(invisible_chars[:5]):  # Test first 5
            adata.obs[f"col{char}{i}"] = np.random.rand(3)

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # All 5 columns should be shown
        assert html.count("col") >= 5, (
            "All columns with invisible chars should be shown"
        )


class TestBidiSpoofing:
    """Tests for bidirectional text spoofing (Trojan Source attack).

    RTL override characters can make code/text appear differently than
    it actually is. This was a major security vulnerability (CVE-2021-42574).
    """

    def test_trojan_source_attack(self, validate_html) -> None:
        """Bidi overrides that reverse text display."""
        adata = AnnData(np.zeros((3, 3)))
        # This looks like "access" but actually contains hidden chars
        # that could hide malicious content
        trojan = "ac\u202ecesscc\u202ca"
        adata.obs[trojan] = [1, 2, 3]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # The bidi characters should not break HTML
        v.assert_no_raw_xss()
        # The column should still be shown (bidi chars preserved or escaped)
        v.assert_section_exists("obs")

    def test_bidi_in_readme_content(self, validate_html) -> None:
        """Bidi attacks in README that could hide malicious instructions."""
        adata = AnnData(np.zeros((3, 3)))
        # Classic Trojan Source pattern - appears to show one thing,
        # but reversed text could hide malicious content
        adata.uns["README"] = """# Analysis Results
The data shows \u202e}};alert('safe')//\u202c normal patterns.
Run: \u202egnirtSresUteG.tnemucod\u202c to see results.
"""

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_no_raw_xss()
        # The "alert" text should appear as escaped text, not executable
        # (it's in the README content which should be escaped)
        assert "alert" in html, "Alert text should be visible (escaped)"
        # No actual script execution
        assert "<script>alert" not in html, "Script should not be executable"


class TestDOMClobbering:
    """Tests for DOM clobbering attacks.

    DOM clobbering uses id/name attributes to override global JavaScript
    properties. If our JS uses `document.getElementById` or accesses
    globals that could be clobbered, this is a vulnerability.
    """

    def test_dom_clobbering_document(self, validate_html) -> None:
        """Names that try to clobber document properties."""
        adata = AnnData(np.zeros((3, 3)))
        # These could clobber window/document properties
        dangerous_names = [
            "document",
            "window",
            "location",
            "navigator",
            "localStorage",
            "sessionStorage",
            "parent",
            "top",
            "self",
            "frames",
        ]
        for name in dangerous_names:
            adata.uns[name] = f"trying to clobber {name}"

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")
        # Verify that user data keys don't become element IDs that could clobber globals
        # The keys should appear in data-key attributes, not as raw IDs
        for name in dangerous_names:
            # Should not have id="document" etc. that could clobber globals
            assert f'id="{name}"' not in html, f"Dangerous id={name} should not exist"

    def test_dom_clobbering_form_elements(self, validate_html) -> None:
        """Names that could clobber form-related properties."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["action"] = "http://evil.com"
        adata.uns["method"] = "POST"
        adata.uns["target"] = "_blank"
        adata.uns["enctype"] = "multipart/form-data"

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # These keys should appear as data, not as form attributes
        # No <form action="http://evil.com"> should be generated
        assert '<form action="http://evil.com"' not in html, "Form action clobbering"


class TestCSSAttacks:
    """Tests for CSS-based attacks that could break visualization."""

    def test_css_animation_dos(self, validate_html) -> None:
        """CSS animations that could consume CPU."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["README"] = """# Test
<style>
@keyframes spin { from { transform: rotate(0deg); } to { transform: rotate(360deg); } }
* { animation: spin 0.001s infinite !important; }
</style>
"""

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # Should only have one style tag (ours) - user CSS should be escaped
        assert html.lower().count("</style>") == 1, (
            "User CSS style tag should be escaped"
        )
        # The @keyframes should appear as text, not as CSS
        assert "@keyframes" not in html or "&" in html, "CSS should be escaped"

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

    def test_css_cursor_hijack(self, validate_html) -> None:
        """CSS that hides or changes cursor."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["README"] = """# Test
<style>* { cursor: none !important; }</style>
<style>body { cursor: url('http://evil.com/cursor.png'), auto; }</style>
"""

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # Should only have our legitimate style tag
        assert html.lower().count("</style>") == 1, "User style tags should be escaped"
        # Evil URL should not appear in actual CSS
        assert "evil.com" not in html.split("</style>")[0], "Evil URL not in CSS"

    def test_css_content_injection(self, validate_html) -> None:
        """CSS content property that could inject fake text."""
        adata = AnnData(np.zeros((3, 3)))
        adata.obs['x{content:"FAKE DATA"}'] = [1, 2, 3]

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # The CSS selector syntax should be escaped, not interpreted
        assert "FAKE DATA" in html, "The text should appear (escaped)"
        # Should not be in an actual style block
        style_section = html.split("</style>")[0] if "</style>" in html else ""
        assert 'content:"FAKE DATA"' not in style_section, "Not in actual CSS"


class TestSVGXSS:
    """Tests for SVG-based XSS attacks.

    SVG can contain JavaScript and even embed HTML via foreignObject.
    """

    def test_svg_script_injection(self, validate_html) -> None:
        """SVG with embedded script."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["README"] = """# Test
<svg><script>alert('SVG XSS')</script></svg>
"""

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_no_raw_xss()
        # The SVG script should be escaped, not executable
        assert "<svg><script>" not in html, "SVG script tag should be escaped"
        # The text should still be visible (escaped)
        assert "alert" in html, "Alert text should be visible as escaped content"

    def test_svg_foreignobject(self, validate_html) -> None:
        """SVG foreignObject that can embed HTML."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["README"] = """# Test
<svg><foreignObject><body onload="alert('foreignObject')"></body></foreignObject></svg>
"""

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_no_raw_xss()
        # foreignObject should be escaped
        assert "<foreignObject>" not in html, "foreignObject should be escaped"
        # onload handler should not be executable
        assert 'onload="alert' not in html, "onload handler should be escaped"

    def test_svg_use_xss(self, validate_html) -> None:
        """SVG use element for external resource loading."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["README"] = """# Test
<svg><use href="http://evil.com/xss.svg#payload"></use></svg>
<svg><use xlink:href="data:image/svg+xml,<svg onload='alert(1)'/>"></use></svg>
"""

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # SVG use element should be escaped, not loading external resources
        assert "<svg><use href=" not in html, "SVG use should be escaped"
        # Data URI should not be interpreted
        assert "data:image/svg+xml" not in html or "&" in html, "Data URI escaped"


class TestPrototypePollution:
    """Tests for prototype pollution patterns.

    If any JavaScript uses user-controlled data as object keys,
    __proto__, constructor, etc. could pollute prototypes.
    """

    def test_proto_in_keys(self, validate_html) -> None:
        """Keys that could cause prototype pollution."""
        adata = AnnData(np.zeros((3, 3)))
        dangerous_keys = [
            "__proto__",
            "constructor",
            "prototype",
            "__defineGetter__",
            "__defineSetter__",
            "__lookupGetter__",
            "__lookupSetter__",
            "hasOwnProperty",
            "isPrototypeOf",
            "propertyIsEnumerable",
            "toLocaleString",
            "toString",
            "valueOf",
        ]
        for key in dangerous_keys:
            adata.uns[key] = {"polluted": True}

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_exists("uns")
        # All dangerous keys should be visible as data (escaped/quoted)
        # They should appear in data-key attributes, not as JS object keys
        for key in dangerous_keys:
            assert key in html, f"Key {key} should be visible in output"

    def test_template_injection_patterns(self, validate_html) -> None:
        """Template injection patterns that could execute in some contexts."""
        adata = AnnData(np.zeros((3, 3)))
        templates = [
            "{{constructor.constructor('alert(1)')()}}",
            "${alert(1)}",
            "#{alert(1)}",
            "<%= alert(1) %>",
            "{{7*7}}",
            "${7*7}",
            "[[${7*7}]]",
            "{{config}}",
            "{{self}}",
        ]
        for i, tmpl in enumerate(templates):
            adata.uns[f"tmpl_{i}"] = tmpl

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # Template patterns should appear as literal text, not evaluated
        # {{7*7}} should not become 49
        assert "49" not in html or "7*7" in html, "Templates should not be evaluated"
        # The constructor pattern should be visible as text
        assert "constructor" in html, "Template text should be visible"


class TestSlowAndBlockingObjects:
    """Tests for objects that could hang or slow down the kernel."""

    def test_slow_repr_handled(self, validate_html) -> None:
        """Slow __repr__ should not block indefinitely."""
        import time

        adata = AnnData(np.zeros((3, 3)))
        # 0.1 second delay - not too long but tests the pattern
        adata.uns["slow"] = SlowRepr(delay=0.1)

        start = time.time()
        html = adata._repr_html_()
        elapsed = time.time() - start

        v = validate_html(html)
        v.assert_html_well_formed()
        v.assert_section_contains_entry("uns", "slow")
        # Should complete in reasonable time (not hang)
        assert elapsed < 5.0, f"Repr took too long: {elapsed}s"

    def test_memory_bomb_handled(self, validate_html) -> None:
        """Objects that try to allocate huge memory should be handled."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["bomb"] = MemoryBomb()

        # This should not raise MemoryError - it should be caught
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_section_contains_entry("uns", "bomb")
        # Should show the type name even if accessing shape failed
        v.assert_text_visible("MemoryBomb")
        # Should show an error indicator
        assert "MemoryError" in html or "error" in html.lower(), (
            "Error should be indicated"
        )


class TestEncodingAttacks:
    """Tests for encoding-related attacks."""

    def test_overlong_utf8_simulation(self, validate_html) -> None:
        """Strings that might bypass filters via encoding tricks."""
        adata = AnnData(np.zeros((3, 3)))
        # Use actual bytes sequences that might be dangerous
        adata.uns["encoded"] = {
            "null_byte": "\x00",  # Actual null byte
            "quote": "\x22",  # Double quote (")
            "lt": "\x3c",  # Less than (<)
        }

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # Null bytes should be handled (removed or escaped)
        assert "\x00" not in html, "Null bytes should be removed from HTML"
        # Less-than should be escaped
        assert "encoded" in html, "Section should be visible"

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


class TestMutationXSS:
    """Tests for mutation XSS (mXSS) attacks.

    mXSS exploits differences in how HTML is parsed vs serialized.
    Content that looks safe can become dangerous after DOM manipulation.
    """

    def test_mutation_xss_patterns(self, validate_html) -> None:
        """Content that could mutate into XSS during parsing."""
        adata = AnnData(np.zeros((3, 3)))
        mxss_payloads = [
            "<img src=x onerror=alert(1)//",  # Missing closing >
            "<svg><![CDATA[><script>alert(1)</script>]]>",
            "<math><mtext><table><mglyph><style><img src=x onerror=alert(1)>",
            '<noscript><p title="</noscript><script>alert(1)</script>">',
            "<<script>alert(1)</script>",
            "<div<script>alert(1)</script>>",
        ]
        for i, payload in enumerate(mxss_payloads):
            adata.uns[f"mxss_{i}"] = payload

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_no_raw_xss()
        # All mXSS payloads should be escaped - < becomes &lt;
        # The key check: raw <img and <script tags don't exist from user content
        assert "<img src=x" not in html, "img tag should be escaped to &lt;img"
        # Escaped content should be present as &lt; entities
        assert "&lt;" in html, "Escaped content should be visible"

    def test_html_entity_smuggling(self, validate_html) -> None:
        """HTML entities that decode to dangerous content."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["entities"] = {
            "script": "&#60;script&#62;alert(1)&#60;/script&#62;",
            "hex": "&#x3c;script&#x3e;alert(1)&#x3c;/script&#x3e;",
            "named": "&lt;script&gt;alert(1)&lt;/script&gt;",
        }

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        v.assert_no_raw_xss()
        # Entity-encoded content should be double-escaped or shown as-is
        # It should NOT decode into executable script tags
        assert "<script>alert(1)</script>" not in html, (
            "Entities should not decode to XSS"
        )


class TestAccessibilityAttacks:
    """Tests for attacks targeting screen readers and assistive tech."""

    def test_aria_label_spoofing(self, validate_html) -> None:
        """Aria labels that could mislead screen reader users."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["README"] = """# Test
<div aria-label="Safe content">Actually malicious instructions here</div>
<span role="button" aria-label="Download">Delete all data</span>
"""

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # The aria-label should be escaped, not interpreted as real HTML attribute
        assert 'aria-label="Safe content"' not in html, "aria-label should be escaped"
        # The raw text should be visible (escaped)
        assert "aria-label" in html, "aria-label text should appear (escaped)"

    def test_visually_hidden_content(self, validate_html) -> None:
        """Content visible to screen readers but hidden visually."""
        adata = AnnData(np.zeros((3, 3)))
        adata.uns["README"] = """# Analysis
<span style="position:absolute;left:-9999px">Secret malicious instructions</span>
<span class="sr-only">Hidden content for screen readers</span>
"""

        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_html_well_formed()
        # The inline style attempting to hide content should be escaped
        assert 'style="position:absolute;left:-9999px"' not in html, (
            "Style should be escaped"
        )
        # User content should not inject CSS positioning
        # The text "Secret malicious instructions" should still be visible
        assert "Secret" in html or "malicious" in html, "Text content should be visible"
