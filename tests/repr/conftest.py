"""
Shared fixtures for repr tests.
"""

from __future__ import annotations

import re
from html.parser import HTMLParser

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from anndata import AnnData


def pytest_configure(config):
    """Suppress ImplicitModificationWarning for all repr tests.

    This warning is expected when AnnData transforms indices internally
    during singledispatch in functools.
    """
    import warnings

    from anndata._warnings import ImplicitModificationWarning

    warnings.filterwarnings("ignore", category=ImplicitModificationWarning)


# =============================================================================
# HTML Testing Helpers
# =============================================================================


class HTMLValidator:
    """
    DOM-based HTML validator for testing repr output.

    Provides structured assertions about HTML content, not just string matching.
    Uses regex-based parsing (no external dependencies like BeautifulSoup).

    Usage:
        validator = HTMLValidator(html)
        validator.assert_element_exists(".adata-badge-view")
        validator.assert_text_in_element(".adata-shape", "100")
        validator.assert_element_has_class("div", "anndata-sec")
    """

    def __init__(self, html: str):
        self.html = html
        self._parsed = False

    def assert_element_exists(
        self, selector: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert an element matching the CSS-like selector exists.

        Supports: tag, .class, #id, tag.class, [attr], [attr=value]
        """
        pattern = self._selector_to_pattern(selector)
        if not re.search(pattern, self.html):
            raise AssertionError(
                msg or f"Element matching '{selector}' not found in HTML"
            )
        return self

    def assert_element_not_exists(
        self, selector: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert no element matching the selector exists."""
        pattern = self._selector_to_pattern(selector)
        if re.search(pattern, self.html):
            raise AssertionError(
                msg or f"Element matching '{selector}' unexpectedly found in HTML"
            )
        return self

    def assert_text_visible(
        self, text: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert text appears in visible content (not in style/script tags)."""
        # Remove style and script content
        visible_html = re.sub(
            r"<(style|script)[^>]*>.*?</\1>", "", self.html, flags=re.DOTALL | re.I
        )
        # Remove HTML tags to get text content
        text_content = re.sub(r"<[^>]+>", " ", visible_html)
        if text not in text_content:
            raise AssertionError(
                msg or f"Text '{text}' not found in visible HTML content"
            )
        return self

    def assert_text_in_element(
        self, selector: str, text: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert text appears within elements matching selector."""
        # Find elements matching selector and check their content
        pattern = self._selector_to_content_pattern(selector)
        matches = re.findall(pattern, self.html, re.DOTALL)
        for match in matches:
            if text in match:
                return self
        raise AssertionError(
            msg or f"Text '{text}' not found in elements matching '{selector}'"
        )

    def assert_attribute_value(
        self,
        selector: str,
        attr: str,
        expected: str,
        *,
        msg: str | None = None,
    ) -> HTMLValidator:
        """Assert element has attribute with expected value."""
        # Find element and check attribute
        pattern = self._selector_to_pattern(selector)
        match = re.search(pattern, self.html)
        if not match:
            msg = f"Element matching '{selector}' not found"
            raise AssertionError(msg)

        element_html = match.group(0)
        attr_pattern = rf'{attr}=["\']([^"\']*)["\']'
        attr_match = re.search(attr_pattern, element_html)
        if not attr_match:
            raise AssertionError(
                msg or f"Attribute '{attr}' not found on element matching '{selector}'"
            )
        if attr_match.group(1) != expected:
            raise AssertionError(
                msg
                or f"Attribute '{attr}' has value '{attr_match.group(1)}', expected '{expected}'"
            )
        return self

    def assert_class_present(
        self, selector: str, class_name: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert element has specified CSS class."""
        pattern = self._selector_to_pattern(selector)
        match = re.search(pattern, self.html)
        if not match:
            msg = f"Element matching '{selector}' not found"
            raise AssertionError(msg)

        element_html = match.group(0)
        class_pattern = r'class=["\']([^"\']*)["\']'
        class_match = re.search(class_pattern, element_html)
        if not class_match:
            raise AssertionError(
                msg or f"No class attribute on element matching '{selector}'"
            )

        classes = class_match.group(1).split()
        if class_name not in classes:
            raise AssertionError(
                msg
                or f"Class '{class_name}' not found on element matching '{selector}'. "
                f"Found: {classes}"
            )
        return self

    def assert_section_exists(
        self, section_name: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert a data section (obs, var, uns, etc.) exists."""
        pattern = rf'data-section=["\']?{re.escape(section_name)}["\']?'
        if not re.search(pattern, self.html):
            raise AssertionError(msg or f"Section '{section_name}' not found in HTML")
        return self

    def assert_section_contains_entry(
        self, section_name: str, entry_name: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert a section contains a specific entry.

        Checks that the entry appears in the section by looking for:
        1. A data-key attribute with the entry name, OR
        2. The entry name as text content

        The entry must appear after the section start marker.

        Note: Handles nested AnnData structures by checking ALL sections
        with the given name, not just the first one.
        """
        # Find ALL sections with this name
        section_pattern = rf'data-section=["\']?{re.escape(section_name)}["\']?'
        section_matches = list(re.finditer(section_pattern, self.html))

        if not section_matches:
            msg = f"Section '{section_name}' not found"
            raise AssertionError(msg)

        # Check each section for the entry (handles nested AnnData case)
        entry_pattern = rf'data-key=["\']?{re.escape(entry_name)}["\']?'

        for section_match in section_matches:
            # Get HTML from this section start
            section_start = section_match.start()
            html_from_section = self.html[section_start:]

            # Find next section (any section) or end of document
            next_section = re.search(
                r'data-section=["\'][^"\']+["\']',
                html_from_section[len(section_match.group(0)) :],
            )
            if next_section:
                section_html = html_from_section[
                    : len(section_match.group(0)) + next_section.start()
                ]
            else:
                section_html = html_from_section

            # Check for entry by data-key attribute or text content
            if re.search(entry_pattern, section_html) or entry_name in section_html:
                return self  # Found!

        raise AssertionError(
            msg or f"Entry '{entry_name}' not found in section '{section_name}'"
        )

    def assert_badge_shown(
        self, badge_type: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert a badge (View, Backed, Lazy) is shown."""
        pattern = rf'class=["\'][^"\']*adata-badge-{re.escape(badge_type.lower())}[^"\']*["\']'
        if not re.search(pattern, self.html, re.I):
            raise AssertionError(msg or f"Badge '{badge_type}' not shown in HTML")
        return self

    def assert_badge_not_shown(
        self, badge_type: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert a badge is NOT shown."""
        pattern = rf'class=["\'][^"\']*adata-badge-{re.escape(badge_type.lower())}[^"\']*["\']'
        if re.search(pattern, self.html, re.I):
            raise AssertionError(
                msg or f"Badge '{badge_type}' unexpectedly shown in HTML"
            )
        return self

    def assert_color_swatch(
        self, color: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert a color swatch with specific color is shown."""
        # Color could be in style attribute or data attribute
        color_lower = color.lower()
        if color_lower not in self.html.lower():
            raise AssertionError(msg or f"Color '{color}' not found in HTML")
        return self

    def assert_warning_indicator(self, *, msg: str | None = None) -> HTMLValidator:
        """Assert a warning indicator (⚠ or warning class) is shown."""
        has_warning = (
            "⚠" in self.html or "warning" in self.html.lower() or "(!)" in self.html
        )
        if not has_warning:
            raise AssertionError(msg or "No warning indicator found in HTML")
        return self

    def assert_shape_displayed(
        self, n_obs: int, n_vars: int, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert shape values are displayed."""
        # Check both values appear (not necessarily together)
        if str(n_obs) not in self.html:
            raise AssertionError(msg or f"n_obs={n_obs} not found in HTML")
        if str(n_vars) not in self.html:
            raise AssertionError(msg or f"n_vars={n_vars} not found in HTML")
        return self

    def assert_dtype_displayed(
        self, dtype: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert dtype is displayed."""
        if dtype.lower() not in self.html.lower():
            raise AssertionError(msg or f"dtype '{dtype}' not found in HTML")
        return self

    # =========================================================================
    # CSS/JavaScript validation methods
    # =========================================================================

    def assert_css_variable_defined(
        self, var_name: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert a CSS variable is defined in style tags or inline styles.

        Usage:
            v.assert_css_variable_defined("--anndata-name-col-width")
        """
        var_pattern = rf"{re.escape(var_name)}\s*:"

        # Check in style tags
        style_pattern = r"<style[^>]*>(.*?)</style>"
        styles = re.findall(style_pattern, self.html, re.DOTALL | re.I)
        for style in styles:
            if re.search(var_pattern, style):
                return self

        # Check in inline styles (style="...")
        inline_pattern = r'style=["\']([^"\']*)["\']'
        inline_styles = re.findall(inline_pattern, self.html)
        for inline in inline_styles:
            if re.search(var_pattern, inline):
                return self

        raise AssertionError(
            msg or f"CSS variable '{var_name}' not defined in styles or inline styles"
        )

    def assert_css_variable_value(
        self, var_name: str, expected_value: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert a CSS variable has expected value.

        Usage:
            v.assert_css_variable_value("--anndata-type-col-width", "150px")
        """
        var_pattern = rf"{re.escape(var_name)}\s*:\s*([^;}}\"']+)"

        # Check in style tags
        style_pattern = r"<style[^>]*>(.*?)</style>"
        styles = re.findall(style_pattern, self.html, re.DOTALL | re.I)
        for style in styles:
            match = re.search(var_pattern, style)
            if match:
                actual = match.group(1).strip()
                if actual == expected_value:
                    return self
                raise AssertionError(
                    msg
                    or f"CSS variable '{var_name}' has value '{actual}', "
                    f"expected '{expected_value}'"
                )

        # Check in inline styles
        inline_pattern = r'style=["\']([^"\']*)["\']'
        inline_styles = re.findall(inline_pattern, self.html)
        for inline in inline_styles:
            match = re.search(var_pattern, inline)
            if match:
                actual = match.group(1).strip()
                if actual == expected_value:
                    return self
                raise AssertionError(
                    msg
                    or f"CSS variable '{var_name}' has value '{actual}', "
                    f"expected '{expected_value}'"
                )

        raise AssertionError(
            msg or f"CSS variable '{var_name}' not found in styles or inline styles"
        )

    def assert_style_contains(
        self, css_fragment: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert style tags contain specific CSS.

        Usage:
            v.assert_style_contains(".collapsed")
            v.assert_style_contains("display: none")
        """
        style_pattern = r"<style[^>]*>(.*?)</style>"
        styles = re.findall(style_pattern, self.html, re.DOTALL | re.I)

        for style in styles:
            if css_fragment in style:
                return self

        raise AssertionError(
            msg or f"CSS fragment '{css_fragment}' not found in styles"
        )

    def assert_script_contains(
        self, js_fragment: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert script tags contain specific JavaScript.

        Usage:
            v.assert_script_contains("toggleCollapse")
            v.assert_script_contains("clipboard")
        """
        script_pattern = r"<script[^>]*>(.*?)</script>"
        scripts = re.findall(script_pattern, self.html, re.DOTALL | re.I)

        for script in scripts:
            if js_fragment in script:
                return self

        raise AssertionError(
            msg or f"JavaScript fragment '{js_fragment}' not found in scripts"
        )

    def assert_collapse_functionality_present(
        self, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert collapse/expand functionality is present in HTML.

        Checks for CSS classes and JavaScript that handle collapsing.

        Usage:
            v.assert_collapse_functionality_present()
        """
        # Check for collapsed CSS class definition and JS toggle functionality
        has_css = ".collapsed" in self.html or "collapsed" in self.html
        has_js = "classList.toggle" in self.html or "classList.add" in self.html

        if not (has_css and has_js):
            raise AssertionError(
                msg or "Collapse functionality (CSS/JS) not found in HTML"
            )
        return self

    def assert_section_initially_collapsed(
        self, section_name: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert a section is marked to start in collapsed state.

        Usage:
            v.assert_section_initially_collapsed("obs")
        """
        # Check data-should-collapse attribute on section
        # Use quotes to ensure exact match (avoid "obs" matching "obsm")
        pattern = (
            rf'data-section="{re.escape(section_name)}"[^>]*data-should-collapse="true"'
        )
        alt_pattern = (
            rf'data-should-collapse="true"[^>]*data-section="{re.escape(section_name)}"'
        )
        if not (re.search(pattern, self.html) or re.search(alt_pattern, self.html)):
            raise AssertionError(
                msg or f"Section '{section_name}' not marked for collapse"
            )
        return self

    def assert_section_not_initially_collapsed(
        self, section_name: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert a section is NOT marked to start in collapsed state.

        Usage:
            v.assert_section_not_initially_collapsed("obs")
        """
        # Check section exists
        self.assert_section_exists(section_name)

        # Check data-should-collapse is false (or not true)
        # Use quotes to ensure exact match (avoid "obs" matching "obsm")
        pattern = (
            rf'data-section="{re.escape(section_name)}"[^>]*data-should-collapse="true"'
        )
        alt_pattern = (
            rf'data-should-collapse="true"[^>]*data-section="{re.escape(section_name)}"'
        )
        if re.search(pattern, self.html) or re.search(alt_pattern, self.html):
            raise AssertionError(
                msg or f"Section '{section_name}' unexpectedly marked for collapse"
            )
        return self

    def assert_has_event_handler(
        self, event: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert elements have specific event handlers.

        Usage:
            v.assert_has_event_handler("click")
            v.assert_has_event_handler("mouseover")
        """
        # Check for onclick, addEventListener('click'), etc.
        inline_pattern = rf"on{event}\s*="
        listener_pattern = rf"addEventListener\s*\(\s*['\"]?{event}['\"]?"

        has_handler = re.search(inline_pattern, self.html, re.I) or re.search(
            listener_pattern, self.html
        )

        if not has_handler:
            raise AssertionError(
                msg or f"Event handler for '{event}' not found in HTML"
            )
        return self

    def assert_has_data_attribute(
        self,
        attr_name: str,
        expected_value: str | None = None,
        *,
        msg: str | None = None,
    ) -> HTMLValidator:
        """Assert data-* attributes exist with optional value check.

        Usage:
            v.assert_has_data_attribute("copy")  # checks data-copy exists
            v.assert_has_data_attribute("section", "obs")  # checks data-section="obs"
        """
        if expected_value is not None:
            pattern = (
                rf'data-{re.escape(attr_name)}=["\']?{re.escape(expected_value)}["\']?'
            )
        else:
            pattern = rf"data-{re.escape(attr_name)}="

        if not re.search(pattern, self.html):
            if expected_value:
                raise AssertionError(
                    msg or f"data-{attr_name}='{expected_value}' not found in HTML"
                )
            raise AssertionError(msg or f"data-{attr_name} attribute not found in HTML")
        return self

    def assert_truncation_indicator(self, *, msg: str | None = None) -> HTMLValidator:
        """Assert truncation is indicated (ellipsis, +N more, etc.)."""
        has_truncation = (
            "..." in self.html
            or "…" in self.html  # Unicode ellipsis
            or re.search(r"\+\s*\d+", self.html)  # +N pattern
            or "more" in self.html.lower()
        )
        if not has_truncation:
            raise AssertionError(msg or "No truncation indicator found in HTML")
        return self

    def assert_accessibility_attribute(
        self, attr: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert accessibility attributes exist (aria-*, role, etc.).

        Usage:
            v.assert_accessibility_attribute("aria-label")
            v.assert_accessibility_attribute("role")
        """
        pattern = rf"{re.escape(attr)}="
        if not re.search(pattern, self.html):
            raise AssertionError(
                msg or f"Accessibility attribute '{attr}' not found in HTML"
            )
        return self

    def count_elements(self, selector: str) -> int:
        """Count elements matching selector."""
        pattern = self._selector_to_pattern(selector)
        return len(re.findall(pattern, self.html))

    def get_text_content(self) -> str:
        """Get visible text content (stripped of tags)."""
        visible_html = re.sub(
            r"<(style|script)[^>]*>.*?</\1>", "", self.html, flags=re.DOTALL | re.I
        )
        return re.sub(r"<[^>]+>", " ", visible_html)

    def _selector_to_pattern(self, selector: str) -> str:  # noqa: PLR0911
        """Convert CSS-like selector to regex pattern."""
        # Handle common selectors
        if selector.startswith("#"):
            # ID selector: #my-id
            id_val = selector[1:]
            return rf'<[^>]+id=["\']?{re.escape(id_val)}["\']?[^>]*>'
        elif selector.startswith("."):
            # Class selector: .my-class
            class_val = selector[1:]
            return (
                rf'<[^>]+class=["\'][^"\']*\b{re.escape(class_val)}\b[^"\']*["\'][^>]*>'
            )
        elif "[" in selector:
            # Attribute selector: [attr] or [attr=value]
            match = re.match(r"\[(\w+)(?:=([^\]]+))?\]", selector)
            if match:
                attr, val = match.groups()
                if val:
                    val = val.strip("\"'")
                    return rf'<[^>]+{attr}=["\']?{re.escape(val)}["\']?[^>]*>'
                return rf"<[^>]+{attr}=[^>]*>"
        elif "." in selector:
            # Tag with class: div.my-class
            tag, class_val = selector.split(".", 1)
            return rf'<{tag}[^>]+class=["\'][^"\']*\b{re.escape(class_val)}\b[^"\']*["\'][^>]*>'
        else:
            # Tag selector: div
            return rf"<{selector}[^>]*>"

        return selector  # Fallback

    def _selector_to_content_pattern(self, selector: str) -> str:
        """Convert selector to pattern that captures element content."""
        if selector.startswith("."):
            class_val = selector[1:]
            return (
                rf'<[^>]+class=["\'][^"\']*\b{re.escape(class_val)}\b[^"\']*["\'][^>]*>'
                rf"(.*?)</[^>]+>"
            )
        elif selector.startswith("#"):
            id_val = selector[1:]
            return rf'<[^>]+id=["\']?{re.escape(id_val)}["\']?[^>]*>(.*?)</[^>]+>'
        else:
            return rf"<{selector}[^>]*>(.*?)</{selector}>"


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
# HTML Validation Helpers
# =============================================================================

VOID_ELEMENTS = frozenset({
    "area",
    "base",
    "br",
    "col",
    "embed",
    "hr",
    "img",
    "input",
    "link",
    "meta",
    "param",
    "source",
    "track",
    "wbr",
})


class StrictHTMLParser(HTMLParser):
    """Validates HTML structure and catches common errors."""

    def __init__(self):
        super().__init__()
        self.tag_stack = []
        self.errors = []
        self.ids_seen = set()

    def handle_starttag(self, tag, attrs):
        if tag not in VOID_ELEMENTS:
            self.tag_stack.append(tag)

        # Check for duplicate IDs
        for name, value in attrs:
            if name == "id":
                if value in self.ids_seen:
                    self.errors.append(f"Duplicate ID: {value}")
                self.ids_seen.add(value)

    def handle_endtag(self, tag):
        if tag not in VOID_ELEMENTS:
            if not self.tag_stack or self.tag_stack[-1] != tag:
                self.errors.append(f"Mismatched tag: </{tag}>")
            else:
                self.tag_stack.pop()


# =============================================================================
# Optional W3C HTML5 Validation
# =============================================================================

# Check for html5validator availability
try:
    from html5_parser import parse  # noqa: F401

    HAS_HTML5_PARSER = True
except ImportError:
    HAS_HTML5_PARSER = False

try:
    import subprocess

    # Check if vnu (Nu Html Checker) is available
    result = subprocess.run(
        ["vnu", "--version"],
        check=False,
        capture_output=True,
        text=True,
        timeout=5,
    )
    HAS_VNU = result.returncode == 0
except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
    HAS_VNU = False


def validate_html5_strict(html: str) -> list[str]:
    """
    Validate HTML5 using Nu Html Checker (vnu) if available.

    Returns list of validation errors/warnings.
    Requires: pip install html5-validator OR system vnu installation.

    This provides:
    - Full W3C HTML5 spec validation
    - ARIA accessibility validation
    - Proper attribute checking
    """
    if not HAS_VNU:
        return []  # Skip if not available

    import json
    import tempfile

    with tempfile.NamedTemporaryFile(mode="w", suffix=".html", delete=False) as f:
        # Wrap in minimal HTML5 document
        full_html = f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="utf-8"><title>Test</title></head>
<body>{html}</body>
</html>"""
        f.write(full_html)
        f.flush()

        try:
            result = subprocess.run(
                ["vnu", "--format", "json", f.name],
                check=False,
                capture_output=True,
                text=True,
                timeout=30,
            )
            if result.stderr:
                data = json.loads(result.stderr)
                return [
                    f"{msg['type']}: {msg['message']}"
                    for msg in data.get("messages", [])
                ]
        except (subprocess.TimeoutExpired, json.JSONDecodeError):
            pass
        finally:
            from pathlib import Path

            Path(f.name).unlink()

    return []


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
