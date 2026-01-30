"""
Tests for HTMLValidator and examples of proper HTML testing patterns.

These tests demonstrate the recommended approach for testing HTML repr output:
- Use structured assertions instead of string matching
- Validate element presence and attributes
- Check text appears in correct elements
- Verify sections and entries are properly rendered
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from anndata import AnnData

from .conftest import HTMLValidator


def _get_top_level_selectors(css: str) -> list[str]:
    """Extract only top-level CSS selectors (brace depth 0).

    With native CSS nesting, nested selectors inherit scope from
    their parent and don't need individual scope checking.
    """
    selectors = []
    depth = 0
    current: list[str] = []
    for char in css:
        if char == "{":
            if depth == 0:
                selector = "".join(current).strip()
                if selector:
                    selectors.append(selector)
            depth += 1
            current = []
        elif char == "}":
            depth -= 1
            depth = max(depth, 0)
            current = []
        elif depth == 0:
            current.append(char)
    return selectors


class TestHTMLValidatorBasics:
    """Tests for HTMLValidator functionality."""

    def test_assert_element_exists_by_class(self):
        """Test finding elements by class selector."""
        html = '<div class="my-class">content</div>'
        v = HTMLValidator(html)
        v.assert_element_exists(".my-class")

    def test_assert_element_exists_by_id(self):
        """Test finding elements by ID selector."""
        html = '<div id="my-id">content</div>'
        v = HTMLValidator(html)
        v.assert_element_exists("#my-id")

    def test_assert_element_exists_by_tag(self):
        """Test finding elements by tag selector."""
        html = "<span>content</span>"
        v = HTMLValidator(html)
        v.assert_element_exists("span")

    def test_assert_element_exists_by_attribute(self):
        """Test finding elements by attribute selector."""
        html = '<div data-section="obs">content</div>'
        v = HTMLValidator(html)
        v.assert_element_exists("[data-section=obs]")

    def test_assert_element_not_exists(self):
        """Test asserting element doesn't exist."""
        html = '<div class="other">content</div>'
        v = HTMLValidator(html)
        v.assert_element_not_exists(".missing-class")

    def test_assert_element_exists_fails(self):
        """Test assertion fails when element missing."""
        html = "<div>content</div>"
        v = HTMLValidator(html)
        with pytest.raises(AssertionError, match="not found"):
            v.assert_element_exists(".missing")

    def test_assert_text_visible(self):
        """Test asserting text is visible (not in style/script)."""
        html = "<style>.foo{}</style><div>visible text</div><script>hidden</script>"
        v = HTMLValidator(html)
        v.assert_text_visible("visible text")
        with pytest.raises(AssertionError):
            v.assert_text_visible("hidden")

    def test_assert_text_in_element(self):
        """Test asserting text appears in specific element."""
        html = '<div class="target">expected text</div><div class="other">other</div>'
        v = HTMLValidator(html)
        v.assert_text_in_element(".target", "expected text")

    def test_assert_section_exists(self):
        """Test asserting data section exists."""
        html = '<div data-section="obs">content</div>'
        v = HTMLValidator(html)
        v.assert_section_exists("obs")

    def test_assert_section_contains_entry(self):
        """Test asserting section contains entry."""
        html = '<div data-section="obs">batch cell_type</div>'
        v = HTMLValidator(html)
        v.assert_section_contains_entry("obs", "batch")

    def test_assert_badge_shown(self):
        """Test asserting badge is displayed."""
        html = '<span class="anndata-badge--view">View</span>'
        v = HTMLValidator(html)
        v.assert_badge_shown("view")

    def test_assert_badge_not_shown(self):
        """Test asserting badge is not displayed."""
        html = '<span class="other">content</span>'
        v = HTMLValidator(html)
        v.assert_badge_not_shown("view")

    def test_assert_shape_displayed(self):
        """Test asserting shape values are displayed."""
        html = "<div>100 obs × 50 var</div>"
        v = HTMLValidator(html)
        v.assert_shape_displayed(100, 50)

    def test_assert_dtype_displayed(self):
        """Test asserting dtype is displayed."""
        html = '<span class="dtype">float32</span>'
        v = HTMLValidator(html)
        v.assert_dtype_displayed("float32")

    def test_count_elements(self):
        """Test counting matching elements."""
        html = '<div class="item">1</div><div class="item">2</div><div class="other">3</div>'
        v = HTMLValidator(html)
        assert v.count_elements(".item") == 2

    def test_chaining(self):
        """Test method chaining works."""
        html = '<div class="a" data-section="obs">text</div>'
        v = HTMLValidator(html)
        (
            v
            .assert_element_exists(".a")
            .assert_section_exists("obs")
            .assert_text_visible("text")
        )


class TestHTMLValidatorWithAnnData:
    """Tests demonstrating HTMLValidator with actual AnnData repr."""

    def test_basic_structure(self, validate_html):
        """Test basic AnnData repr structure."""
        adata = AnnData(np.zeros((100, 50)))
        html = adata._repr_html_()
        v = validate_html(html)

        # Validate structure, not just string presence
        v.assert_element_exists(".anndata-repr")
        v.assert_shape_displayed(100, 50)

    def test_sections_present(self, validate_html):
        """Test all populated sections are present."""
        adata = AnnData(
            np.zeros((10, 5)),
            obs=pd.DataFrame({"batch": ["A", "B"] * 5}),
            var=pd.DataFrame({"gene": range(5)}),
        )
        adata.uns["key"] = "value"
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_section_exists("obs")
        v.assert_section_exists("var")
        v.assert_section_exists("uns")

    def test_obs_entries_in_correct_section(self, validate_html):
        """Test obs column names appear in obs section."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["cell_type"] = pd.Categorical(["A", "B"] * 5)
        adata.obs["n_counts"] = range(10)
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_section_contains_entry("obs", "cell_type")
        v.assert_section_contains_entry("obs", "n_counts")

    def test_view_badge_displayed_correctly(self, validate_html):
        """Test View badge is shown for views."""
        adata = AnnData(np.zeros((10, 5)))
        view = adata[0:5, :]
        html = view._repr_html_()
        v = validate_html(html)

        v.assert_badge_shown("view")
        v.assert_shape_displayed(5, 5)  # View shape, not original

    def test_view_badge_not_shown_for_non_view(self, validate_html):
        """Test View badge is NOT shown for non-views."""
        adata = AnnData(np.zeros((10, 5)))
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_badge_not_shown("view")

    def test_sparse_matrix_dtype_displayed(self, validate_html):
        """Test sparse matrix shows dtype."""
        X = sp.random(100, 50, density=0.1, format="csr", dtype=np.float32)
        adata = AnnData(X)
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_dtype_displayed("float32")
        v.assert_text_visible("csr")

    def test_categorical_with_colors(self, validate_html):
        """Test categorical with colors shows color values."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["cluster"] = pd.Categorical(["A", "B"] * 5)
        adata.uns["cluster_colors"] = ["#FF0000", "#00FF00"]
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_section_contains_entry("obs", "cluster")
        v.assert_color_swatch("#FF0000")
        v.assert_color_swatch("#00FF00")

    def test_warning_for_unserializable(self, validate_html):
        """Test warning indicator for unserializable objects."""

        class CustomClass:
            pass

        adata = AnnData(np.zeros((5, 3)))
        adata.uns["custom"] = CustomClass()
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_section_contains_entry("uns", "custom")
        v.assert_warning_indicator()

    def test_backed_badge(self, validate_html, tmp_path):
        """Test backed badge is shown for backed AnnData with sparse X."""
        from scipy import sparse

        import anndata as ad

        # Use sparse matrix to cover BackedSparseDatasetFormatter
        adata = AnnData(sparse.random(100, 50, density=0.1, format="csr"))
        path = tmp_path / "test.h5ad"
        adata.write_h5ad(path)

        backed = ad.read_h5ad(path, backed="r")
        html = backed._repr_html_()
        v = validate_html(html)

        v.assert_badge_shown("backed")
        # Verify sparse matrix is shown with "on disk" indicator
        v.assert_text_visible("on disk")
        backed.file.close()

    def test_raw_section_visible(self, validate_html):
        """Test raw section is visible when raw is set."""
        adata = AnnData(np.zeros((10, 20)))
        adata.raw = adata.copy()
        adata = adata[:, :5]
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_section_exists("raw")
        # Raw should show original var count
        v.assert_text_visible("20")

    def test_layers_section(self, validate_html):
        """Test layers section shows layer names."""
        adata = AnnData(np.zeros((10, 5)))
        adata.layers["counts"] = np.ones((10, 5))
        adata.layers["normalized"] = np.zeros((10, 5))
        html = adata._repr_html_()
        v = validate_html(html)

        v.assert_section_exists("layers")
        v.assert_section_contains_entry("layers", "counts")
        v.assert_section_contains_entry("layers", "normalized")


class TestMigrationExamples:
    """Examples showing how to migrate old-style tests to HTMLValidator.

    OLD STYLE (string matching):
        html = adata._repr_html_()
        assert "batch" in html
        assert "View" in html

    NEW STYLE (structured validation):
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "batch")
        v.assert_badge_shown("view")
    """

    def test_old_vs_new_section_check(self, validate_html):
        """Compare old vs new section checking."""
        adata = AnnData(np.zeros((10, 5)))
        adata.obs["batch"] = ["A", "B"] * 5
        html = adata._repr_html_()

        # OLD: String in HTML (could match CSS class name, comment, etc.)
        assert "batch" in html  # Weak assertion

        # NEW: Entry in correct section
        v = validate_html(html)
        v.assert_section_contains_entry("obs", "batch")  # Strong assertion

    def test_old_vs_new_badge_check(self, validate_html):
        """Compare old vs new badge checking."""
        adata = AnnData(np.zeros((10, 5)))
        view = adata[0:5, :]
        html = view._repr_html_()

        # OLD: String matching (could match "View" in documentation, comments, etc.)
        assert "View" in html  # Weak assertion

        # NEW: Check for actual badge element
        v = validate_html(html)
        v.assert_badge_shown("view")  # Strong assertion - checks CSS class

    def test_old_vs_new_shape_check(self, validate_html):
        """Compare old vs new shape checking."""
        adata = AnnData(np.zeros((123, 456)))
        html = adata._repr_html_()

        # OLD: Numbers anywhere in HTML
        assert "123" in html  # Could match anything
        assert "456" in html

        # NEW: Explicit shape check
        v = validate_html(html)
        v.assert_shape_displayed(123, 456)


class TestJupyterNotebookCompatibility:
    """Tests for Jupyter Notebook/Lab HTML compatibility.

    These tests ensure the HTML repr works correctly when embedded
    in Jupyter output cells, including:
    - CSS scoping (no style leakage)
    - JavaScript isolation (no global pollution)
    - Multiple cell compatibility (unique IDs)
    - Jupyter theming support
    """

    def test_css_scoped_to_anndata_repr(self, adata_full):
        """Test CSS rules are scoped to .anndata-repr container.

        Unscoped CSS could affect other notebook cells or UI elements.
        With native CSS nesting, only top-level selectors need checking
        since nested selectors inherit scope from their parent.
        """
        import re

        html = adata_full._repr_html_()
        style_match = re.search(r"<style[^>]*>(.*?)</style>", html, re.DOTALL)

        if style_match:
            css = style_match.group(1)

            # Remove CSS comments
            css_clean = re.sub(r"/\*.*?\*/", "", css, flags=re.DOTALL)

            # Extract only top-level selectors (brace depth 0).
            # With native CSS nesting, nested selectors inherit scope
            # from their parent and don't need individual checking.
            selectors = _get_top_level_selectors(css_clean)

            for selector in selectors:
                # Skip :root (used for CSS variables)
                if ":root" in selector:
                    continue
                # Skip Jupyter theme selectors (these are intentionally global)
                if "[data-jp-theme" in selector or ".jp-" in selector:
                    continue
                # Skip @media / @keyframes (parsed separately)
                if selector.startswith("@"):
                    continue

                # All other selectors should be scoped to anndata-repr
                assert ".anndata-repr" in selector or "anndata" in selector.lower(), (
                    f"CSS selector '{selector}' is not scoped to .anndata-repr. "
                    "This could affect other Jupyter cells."
                )

    def test_no_global_element_selectors(self, adata_full):
        """Test no unscoped element selectors like 'div', 'span', etc.

        Global element selectors would style ALL divs/spans in the notebook.
        With native CSS nesting, only top-level selectors are checked since
        nested element selectors (e.g., `td` inside `.anndata-repr`) are scoped.
        """
        import re

        html = adata_full._repr_html_()
        style_match = re.search(r"<style[^>]*>(.*?)</style>", html, re.DOTALL)

        if style_match:
            css = style_match.group(1)
            css_clean = re.sub(r"/\*.*?\*/", "", css, flags=re.DOTALL)

            # Only check top-level selectors (brace depth 0)
            selectors = _get_top_level_selectors(css_clean)

            bare_elements = {
                "div",
                "span",
                "table",
                "tr",
                "td",
                "th",
                "ul",
                "li",
                "p",
                "a",
                "button",
            }
            global_elements = [
                s
                for s in selectors
                if s.strip().split(",")[0].strip().split()[0] in bare_elements
            ]
            assert not global_elements, (
                f"Found global element selectors: {global_elements}. "
                "These would affect the entire notebook."
            )

    def test_javascript_uses_iife_or_closure(self, adata_full):
        """Test JavaScript is wrapped to avoid global scope pollution."""
        import re

        html = adata_full._repr_html_()
        script_matches = re.findall(
            r"<script[^>]*>(.*?)</script>", html, re.DOTALL | re.I
        )

        for script in script_matches:
            script = script.strip()
            if not script:
                continue

            # Check for IIFE pattern: (function() { ... })() or (() => { ... })()
            has_iife = bool(
                re.search(r"\(\s*function\s*\([^)]*\)\s*\{", script)
                or re.search(r"\(\s*\([^)]*\)\s*=>\s*\{", script)
            )

            # Check for block scope (const/let at top level within braces)
            has_block_scope = bool(re.search(r"\{\s*(const|let)\s+", script))

            # Check for event handler inline (scoped to element)
            is_event_handler = bool(re.search(r"\.addEventListener\s*\(", script))

            # Should use one of these isolation patterns
            assert has_iife or has_block_scope or is_event_handler, (
                "JavaScript should use IIFE, block scope, or event handlers "
                "to avoid polluting global scope in Jupyter."
            )

    def test_no_global_function_declarations(self, adata_full):
        """Test no unscoped 'function name()' declarations.

        Named function declarations at top level would pollute global scope.
        """
        import re

        html = adata_full._repr_html_()
        script_matches = re.findall(
            r"<script[^>]*>(.*?)</script>", html, re.DOTALL | re.I
        )

        for script in script_matches:
            # Look for "function name(" not inside another function/IIFE
            # This is a simplified check - looks for function at start of line
            global_funcs = re.findall(
                r"^\s*function\s+(\w+)\s*\(", script, re.MULTILINE
            )

            # Filter out functions that are clearly inside IIFEs
            # (script starts with "(function" or "(() =>")
            if script.strip().startswith("("):
                continue  # Inside IIFE, OK

            assert not global_funcs, (
                f"Found global function declarations: {global_funcs}. "
                "Use const name = function() or wrap in IIFE."
            )

    def test_unique_ids_per_render(self, adata):
        """Test each render produces unique element IDs.

        Multiple AnnData cells in same notebook must not have ID collisions.
        """
        import re

        html1 = adata._repr_html_()
        html2 = adata._repr_html_()

        ids1 = set(re.findall(r'id=["\']([^"\']+)["\']', html1))
        ids2 = set(re.findall(r'id=["\']([^"\']+)["\']', html2))

        # If both have IDs, they should be different (use unique prefixes)
        if ids1 and ids2:
            # At least the container IDs should be different
            overlap = ids1 & ids2
            # Allow empty overlap or ensure IDs use unique suffixes
            for id_val in overlap:
                # IDs like "anndata-repr" without unique suffix are problematic
                assert re.search(r"[a-f0-9]{6,}|_\d+$", id_val), (
                    f"ID '{id_val}' appears in both renders without unique suffix. "
                    "This could cause conflicts in Jupyter notebooks with multiple cells."
                )

    def test_jupyter_dark_mode_support(self, adata, validate_html):
        """Test dark mode CSS uses Jupyter-compatible selectors."""
        html = adata._repr_html_()

        # Should support at least one Jupyter dark mode detection method
        has_jp_theme = "[data-jp-theme-light" in html
        has_prefers_color = "prefers-color-scheme: dark" in html
        has_jp_dark_class = ".jp-Theme-Dark" in html

        assert has_jp_theme or has_prefers_color or has_jp_dark_class, (
            "HTML should support Jupyter dark mode via "
            "[data-jp-theme-light], prefers-color-scheme, or .jp-Theme-Dark"
        )

    def test_no_document_level_operations(self, adata_full):
        """Test JavaScript doesn't use document-level operations unsafely.

        Operations like document.querySelector without scoping could
        affect elements in other cells.
        """
        import re

        html = adata_full._repr_html_()
        script_matches = re.findall(
            r"<script[^>]*>(.*?)</script>", html, re.DOTALL | re.I
        )

        for script in script_matches:
            # Check for document.querySelector without container scoping
            # OK: container.querySelector, element.querySelector
            # Risky: document.querySelector(".class") without context

            # This is a heuristic check - look for document.querySelector
            # that doesn't immediately follow a variable assignment from
            # a scoped search
            if "document.querySelector" in script:
                # Should be immediately scoped, e.g.:
                # const container = document.getElementById('unique-id')
                # container.querySelector(...)
                has_scoped_search = bool(
                    re.search(r"getElementById\s*\(['\"][\w-]+['\"]", script)
                )
                assert has_scoped_search, (
                    "document.querySelector should be scoped to a container "
                    "obtained via getElementById with unique ID"
                )

    def test_html_valid_as_fragment(self, adata_full, validate_html5):
        """Test HTML is valid when embedded in a typical Jupyter output div.

        Jupyter wraps output in: <div class="output_subarea">...</div>
        """
        html = adata_full._repr_html_()

        # Wrap in typical Jupyter output structure
        jupyter_wrapper = f"""
        <div class="output_wrapper">
            <div class="output">
                <div class="output_area">
                    <div class="output_subarea output_html rendered_html">
                        {html}
                    </div>
                </div>
            </div>
        </div>
        """

        errors = validate_html5(jupyter_wrapper)
        # Filter expected fragment-related warnings
        critical = [
            e
            for e in errors
            if not e.startswith("info:")
            and "style" not in e.lower()
            and "script" not in e.lower()
            # vnu's CSS parser doesn't support native CSS nesting
            and "css: parse error" not in e.lower()
        ]
        assert not critical, "HTML invalid in Jupyter context:\n" + "\n".join(critical)

    def test_multiple_cells_valid_html(self, validate_html5):
        """Test multiple AnnData reprs together produce valid HTML.

        Simulates having multiple AnnData output cells in one notebook view.
        """
        # Create different AnnData objects
        adata1 = AnnData(np.zeros((10, 5)))
        adata1.obs["batch"] = ["A", "B"] * 5

        adata2 = AnnData(np.zeros((20, 8)))
        adata2.uns["key"] = "value"

        adata3 = AnnData(sp.random(100, 50, density=0.1, format="csr"))

        html1 = adata1._repr_html_()
        html2 = adata2._repr_html_()
        html3 = adata3._repr_html_()

        # Combine as if in multiple notebook cells
        combined = f"""
        <div class="cell output">{html1}</div>
        <div class="cell output">{html2}</div>
        <div class="cell output">{html3}</div>
        """

        errors = validate_html5(combined)
        critical = [
            e
            for e in errors
            if not e.startswith("info:")
            and "style" not in e.lower()
            and "script" not in e.lower()
            # Duplicate styles are OK in fragments
            and "duplicate" not in e.lower()
            # vnu's CSS parser doesn't support native CSS nesting
            and "css: parse error" not in e.lower()
        ]
        assert not critical, "Combined cells invalid:\n" + "\n".join(critical)

    def test_no_id_collisions_multiple_cells(self):
        """Test no ID collisions when multiple cells are rendered."""
        import re

        adata1 = AnnData(np.zeros((10, 5)))
        adata2 = AnnData(np.zeros((20, 8)))
        adata3 = AnnData(np.zeros((15, 6)))

        html1 = adata1._repr_html_()
        html2 = adata2._repr_html_()
        html3 = adata3._repr_html_()

        ids1 = re.findall(r'id=["\']([^"\']+)["\']', html1)
        ids2 = re.findall(r'id=["\']([^"\']+)["\']', html2)
        ids3 = re.findall(r'id=["\']([^"\']+)["\']', html3)

        all_ids = ids1 + ids2 + ids3
        if all_ids:
            # Check for duplicates
            seen = set()
            duplicates = []
            for id_val in all_ids:
                if id_val in seen:
                    duplicates.append(id_val)
                seen.add(id_val)

            assert not duplicates, (
                f"Duplicate IDs across cells: {duplicates}. "
                "Each cell must have unique IDs."
            )

    def test_css_variables_prefixed(self, adata):
        """Test CSS variables use anndata prefix to avoid conflicts."""
        import re

        html = adata._repr_html_()

        # Find all CSS variable definitions (not BEM modifiers like --copied)
        # CSS vars are defined as: --name: value; (note single colon)
        # BEM modifiers are: .block--modifier (in class names)
        css_vars = re.findall(r"(?<![.\w])--([\w-]+)\s*:", html)

        for var in css_vars:
            # Should use anndata prefix OR be a standard Jupyter variable
            is_anndata_prefixed = var.startswith("anndata")
            is_jupyter_var = var.startswith("jp-")

            assert is_anndata_prefixed or is_jupyter_var, (
                f"CSS variable '--{var}' should be prefixed with 'anndata' "
                "to avoid conflicts with other Jupyter extensions."
            )

    def test_inline_styles_for_isolation(self, adata):
        """Test critical styles use inline or scoped approach."""
        html = adata._repr_html_()

        # The HTML should either:
        # 1. Use inline styles for critical layout
        # 2. Use scoped <style> with unique ID prefix
        # 3. Use CSS variables that cascade properly

        # Check that the main container uses scoped class
        assert 'class="anndata-repr' in html or "class='anndata-repr" in html, (
            "Main container should use .anndata-repr class for CSS scoping"
        )

    def test_works_without_jupyter_css_variables(self, adata, validate_html):
        """Test repr works even without Jupyter CSS variables defined.

        The repr should have sensible defaults when --jp-* variables
        are not available (e.g., when viewed outside Jupyter).
        """
        html = adata._repr_html_()
        v = validate_html(html)

        # Should still render all key elements
        v.assert_element_exists(".anndata-repr")
        v.assert_shape_displayed(100, 50)

        # Should have fallback colors defined (not relying solely on --jp-* vars)
        # Check that colors are defined in the CSS, not just referenced
        assert "#" in html or "rgb" in html.lower(), (
            "Should define fallback colors for non-Jupyter environments"
        )
