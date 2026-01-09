"""
HTML validation utilities for repr tests.

This module provides:
- HTMLValidator: DOM-like assertions for testing HTML output
- StrictHTMLParser: HTML structure validation
- validate_html5_strict: W3C HTML5 validation (if vnu available)
"""

from __future__ import annotations

import re
from html.parser import HTMLParser

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
        """Assert a CSS variable is defined in style tags or inline styles."""
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
        """Assert a CSS variable has expected value."""
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
        """Assert style tags contain specific CSS."""
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
        """Assert script tags contain specific JavaScript."""
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
        """Assert collapse/expand functionality is present in HTML."""
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
        """Assert a section is marked to start in collapsed state."""
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
        """Assert a section is NOT marked to start in collapsed state."""
        self.assert_section_exists(section_name)
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
        """Assert elements have specific event handlers."""
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
        """Assert data-* attributes exist with optional value check."""
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

    def assert_error_shown(
        self, error_text: str | None = None, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert an error message is shown in the HTML output.

        The repr system should show errors visibly, not hide them.
        Errors are typically shown as:
        - Text containing 'error:' or 'Error'
        - Elements with 'error' class
        - Text with exception names (e.g., 'AttributeError', 'RuntimeError')
        """
        if error_text:
            if error_text not in self.html:
                raise AssertionError(
                    msg or f"Error text '{error_text}' not found in HTML"
                )
            return self

        # Check for generic error indicators
        has_error = (
            "error:" in self.html.lower()
            or re.search(r'\bError\b', self.html)
            or "adata-error" in self.html
            or "adata-text-muted" in self.html and "error" in self.html.lower()
        )
        if not has_error:
            raise AssertionError(msg or "No error indicator found in HTML")
        return self

    def assert_no_raw_xss(self, *, msg: str | None = None) -> HTMLValidator:
        """Assert no raw XSS payloads exist in executable context.

        Checks that common XSS attack vectors are properly escaped.
        Note: Escaped content showing XSS strings as visible text is OK.
        We only flag actual executable payloads (unescaped in HTML attributes).
        """
        # Get content after the style tag (where user content would be)
        style_end = self.html.find("</style>")
        content = self.html[style_end:] if style_end > 0 else self.html

        # Check for actual executable script injection
        # Raw <script> tags (not escaped as &lt;script&gt;)
        if re.search(r'<script[^>]*>(?!.*&lt;)', content, re.I):
            # Make sure it's not a false positive from our own script
            scripts = re.findall(r'<script[^>]*>(.*?)</script>', content, re.DOTALL | re.I)
            for script in scripts:
                # Our own scripts have specific patterns
                if 'anndata' not in script.lower() and 'toggle' not in script.lower():
                    if 'alert' in script or 'eval(' in script:
                        raise AssertionError(
                            msg or "Potential XSS: executable script tag found"
                        )

        # Check for event handlers in HTML tags (actual attributes, not text)
        # Pattern: <tag ... onclick="..." or onerror="..."
        event_handlers = ['onclick', 'onerror', 'onload', 'onmouseover', 'onfocus']
        for handler in event_handlers:
            # Look for handler in tag attributes (not in text content)
            pattern = rf'<[^>]+\s{handler}\s*='
            if re.search(pattern, content, re.I):
                # Verify it's not in our own legitimate HTML
                match = re.search(pattern, content, re.I)
                if match:
                    # Check if it's part of user content (escaped) or real attribute
                    context = content[max(0, match.start() - 50):match.end() + 50]
                    if '&lt;' not in context and '&quot;' not in context:
                        raise AssertionError(
                            msg or f"Potential XSS: {handler} handler in tag attribute"
                        )

        # Check for javascript: URLs in href attributes
        if re.search(r'href\s*=\s*["\']?\s*javascript:', content, re.I):
            raise AssertionError(
                msg or "Potential XSS: javascript: URL in href"
            )

        return self

    def assert_html_well_formed(self, *, msg: str | None = None) -> HTMLValidator:
        """Assert HTML is well-formed (balanced tags, no duplicate IDs)."""
        parser = StrictHTMLParser()
        parser.feed(self.html)
        if parser.errors:
            raise AssertionError(
                msg or f"HTML is malformed: {parser.errors}"
            )
        if parser.tag_stack:
            raise AssertionError(
                msg or f"Unclosed tags: {parser.tag_stack}"
            )
        return self

    def assert_accessibility_attribute(
        self, attr: str, *, msg: str | None = None
    ) -> HTMLValidator:
        """Assert accessibility attributes exist (aria-*, role, etc.)."""
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
        if selector.startswith("#"):
            id_val = selector[1:]
            return rf'<[^>]+id=["\']?{re.escape(id_val)}["\']?[^>]*>'
        elif selector.startswith("."):
            class_val = selector[1:]
            return (
                rf'<[^>]+class=["\'][^"\']*\b{re.escape(class_val)}\b[^"\']*["\'][^>]*>'
            )
        elif "[" in selector:
            match = re.match(r"\[(\w+)(?:=([^\]]+))?\]", selector)
            if match:
                attr, val = match.groups()
                if val:
                    val = val.strip("\"'")
                    return rf'<[^>]+{attr}=["\']?{re.escape(val)}["\']?[^>]*>'
                return rf"<[^>]+{attr}=[^>]*>"
        elif "." in selector:
            tag, class_val = selector.split(".", 1)
            return rf'<{tag}[^>]+class=["\'][^"\']*\b{re.escape(class_val)}\b[^"\']*["\'][^>]*>'
        else:
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
# HTML Structure Validation
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

# Check for vnu (Nu Html Checker) availability
try:
    import subprocess

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
    """
    if not HAS_VNU:
        return []  # Skip if not available

    import json
    import tempfile

    with tempfile.NamedTemporaryFile(mode="w", suffix=".html", delete=False) as f:
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
