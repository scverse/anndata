"""
Visual inspection script for AnnData HTML representation.

Run this script to generate an HTML file that can be opened in a browser
to visually inspect the _repr_html_ output.

Usage:
    python tests/visual_inspect_repr_html.py

Then open tests/repr_html_visual_test.html in your browser.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from anndata import AnnData
import anndata as ad
from anndata._repr import (
    register_formatter,
    TypeFormatter,
    FormattedOutput,
    extract_uns_type_hint,
)

# Check optional dependencies
try:
    import dask.array as da
    HAS_DASK = True
except ImportError:
    HAS_DASK = False

try:
    import networkx as nx
    from treedata import TreeData
    from anndata._repr import (
        register_formatter,
        SectionFormatter,
        FormattedEntry,
        FormattedOutput,
        FormatterContext,
    )
    HAS_TREEDATA = True

    def _render_tree_svg(tree: nx.DiGraph, max_leaves: int = 30, width: int = 300, height: int = 150) -> str:
        """Render a tree as an SVG visualization.

        Uses a simple top-down layout similar to pycea's approach but generates SVG.
        """
        # Find root and leaves
        roots = [n for n in tree.nodes() if tree.in_degree(n) == 0]
        if not roots:
            return "<span style='color:#888;font-size:11px;'>Invalid tree (no root)</span>"
        root = roots[0]

        leaves = [n for n in tree.nodes() if tree.out_degree(n) == 0]
        n_leaves = len(leaves)

        # Truncate if too many leaves
        if n_leaves > max_leaves:
            return (
                f"<span style='color:#888;font-size:11px;'>"
                f"Tree with {n_leaves} leaves (too large to preview)</span>"
            )

        # Compute depths using BFS
        depths = {root: 0}
        queue = [root]
        while queue:
            node = queue.pop(0)
            for child in tree.successors(node):
                depths[child] = depths[node] + 1
                queue.append(child)

        max_depth = max(depths.values()) if depths else 0
        if max_depth == 0:
            return "<span style='color:#888;font-size:11px;'>Single node tree</span>"

        # Assign y-coordinates (leaves get sequential positions)
        y_coords = {}
        leaf_idx = 0
        def assign_y(node):
            nonlocal leaf_idx
            children = list(tree.successors(node))
            if not children:  # leaf
                y_coords[node] = leaf_idx
                leaf_idx += 1
            else:
                for child in children:
                    assign_y(child)
                # Internal node: average of children
                y_coords[node] = sum(y_coords[c] for c in children) / len(children)

        assign_y(root)

        # Scale coordinates
        margin = 15
        x_scale = (width - 2 * margin) / max_depth if max_depth > 0 else 1
        y_scale = (height - 2 * margin) / (n_leaves - 1) if n_leaves > 1 else 1

        def get_x(node):
            return margin + depths[node] * x_scale

        def get_y(node):
            return margin + y_coords[node] * y_scale

        # Generate SVG
        svg_parts = [
            f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg" '
            f'style="background:#fafafa;border-radius:4px;border:1px solid #e0e0e0;">'
        ]

        # Draw branches (parent -> child)
        for parent, child in tree.edges():
            px, py = get_x(parent), get_y(parent)
            cx, cy = get_x(child), get_y(child)
            # Draw elbow connector (horizontal then vertical)
            svg_parts.append(
                f'<path d="M{px:.1f},{py:.1f} L{cx:.1f},{py:.1f} L{cx:.1f},{cy:.1f}" '
                f'fill="none" stroke="#666" stroke-width="1.5"/>'
            )

        # Draw nodes
        for node in tree.nodes():
            x, y = get_x(node), get_y(node)
            is_leaf = tree.out_degree(node) == 0
            r = 3 if is_leaf else 4
            fill = "#4a90d9" if is_leaf else "#333"
            svg_parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}"/>')

        svg_parts.append('</svg>')
        return "".join(svg_parts)

    # TreeData documentation URL
    TREEDATA_DOCS = "https://treedata.readthedocs.io/en/latest/"

    # Register TreeData section formatters
    @register_formatter
    class ObstSectionFormatter(SectionFormatter):
        """Section formatter for obst (observation trees)."""

        @property
        def section_name(self) -> str:
            return "obst"

        @property
        def after_section(self) -> str:
            return "obsm"

        @property
        def doc_url(self) -> str:
            return TREEDATA_DOCS

        @property
        def tooltip(self) -> str:
            return "Tree annotation of observations (TreeData)"

        def should_show(self, obj) -> bool:
            return hasattr(obj, "obst") and len(obj.obst) > 0

        def get_entries(self, obj, context: FormatterContext) -> list[FormattedEntry]:
            entries = []
            for key, tree in obj.obst.items():
                n_nodes = tree.number_of_nodes()
                n_leaves = sum(1 for n in tree.nodes() if tree.out_degree(n) == 0)
                # Generate SVG preview
                svg_html = _render_tree_svg(tree)
                output = FormattedOutput(
                    type_name=f"DiGraph ({n_nodes} nodes, {n_leaves} leaves)",
                    css_class="dtype-tree",
                    tooltip=f"Phylogenetic tree with {n_nodes} total nodes",
                    html_content=svg_html,
                    is_expandable=True,
                )
                entries.append(FormattedEntry(key=key, output=output))
            return entries

    @register_formatter
    class VartSectionFormatter(SectionFormatter):
        """Section formatter for vart (variable trees)."""

        @property
        def section_name(self) -> str:
            return "vart"

        @property
        def after_section(self) -> str:
            return "varm"

        @property
        def doc_url(self) -> str:
            return TREEDATA_DOCS

        @property
        def tooltip(self) -> str:
            return "Tree annotation of variables (TreeData)"

        def should_show(self, obj) -> bool:
            return hasattr(obj, "vart") and len(obj.vart) > 0

        def get_entries(self, obj, context: FormatterContext) -> list[FormattedEntry]:
            entries = []
            for key, tree in obj.vart.items():
                n_nodes = tree.number_of_nodes()
                n_leaves = sum(1 for n in tree.nodes() if tree.out_degree(n) == 0)
                # Generate SVG preview
                svg_html = _render_tree_svg(tree)
                output = FormattedOutput(
                    type_name=f"DiGraph ({n_nodes} nodes, {n_leaves} leaves)",
                    css_class="dtype-tree",
                    tooltip=f"Phylogenetic tree with {n_nodes} total nodes",
                    html_content=svg_html,
                    is_expandable=True,
                )
                entries.append(FormattedEntry(key=key, output=output))
            return entries

except ImportError:
    HAS_TREEDATA = False


def create_test_treedata():
    """Create a TreeData object with observation and variable trees."""
    if not HAS_TREEDATA:
        return None

    np.random.seed(42)
    n_obs = 24  # Small enough for SVG preview (< 30 leaves)
    n_vars = 45  # Large enough to trigger "too large to preview" (> 30 leaves)
    obs_names = [f"cell_{i}" for i in range(n_obs)]
    var_names = [f"gene_{i}" for i in range(n_vars)]

    # Create observation tree (phylogenetic-like structure)
    obs_tree = nx.DiGraph()
    obs_tree.add_edges_from([
        ("root", "clade_A"), ("root", "clade_B"),
        ("clade_A", "subA1"), ("clade_A", "subA2"),
        ("clade_B", "subB1"), ("clade_B", "subB2"),
    ])
    for i, name in enumerate(obs_names):
        parent = ["subA1", "subA2", "subB1", "subB2"][i % 4]
        obs_tree.add_edge(parent, name)

    # Create variable tree (gene ontology-like structure, >30 leaves)
    var_tree = nx.DiGraph()
    var_tree.add_edges_from([
        ("all_genes", "pathway_X"), ("all_genes", "pathway_Y"), ("all_genes", "pathway_Z"),
        ("pathway_X", "module_1"), ("pathway_X", "module_2"),
        ("pathway_Y", "module_3"), ("pathway_Y", "module_4"),
        ("pathway_Z", "module_5"),
    ])
    for i, name in enumerate(var_names):
        parent = ["module_1", "module_2", "module_3", "module_4", "module_5"][i % 5]
        var_tree.add_edge(parent, name)

    # Create TreeData (label=None prevents adding "tree" column to obs/var)
    tdata = TreeData(
        X=np.random.randn(n_obs, n_vars).astype(np.float32),
        obs=pd.DataFrame(
            {"cell_type": pd.Categorical(["T cell", "B cell"] * (n_obs // 2))},
            index=obs_names,
        ),
        var=pd.DataFrame({"gene_name": var_names}, index=var_names),
        obst={"phylogeny": obs_tree},
        vart={"gene_ontology": var_tree},
        label=None,  # Don't add "tree" column
    )

    # Add standard annotations
    tdata.uns["cell_type_colors"] = ["#e41a1c", "#377eb8"]
    tdata.obsm["X_pca"] = np.random.randn(n_obs, 10).astype(np.float32)
    tdata.layers["raw"] = np.random.randn(n_obs, n_vars).astype(np.float32)

    return tdata


def create_test_anndata() -> AnnData:
    """Create a comprehensive test AnnData with all features."""
    n_obs, n_vars = 100, 50

    adata = AnnData(
        sp.random(n_obs, n_vars, density=0.1, format="csr", dtype=np.float32),
        obs=pd.DataFrame({
            "batch": pd.Categorical(["A", "B"] * (n_obs // 2)),
            "n_counts": np.random.randint(1000, 10000, n_obs),
            "cell_type": pd.Categorical(
                ["T cell", "B cell", "NK cell"] * (n_obs // 3) + ["T cell"] * (n_obs % 3)
            ),
            "string_col": ["sample_" + str(i % 5) for i in range(n_obs)],  # Will warn
            "donor_id": [f"donor_{i}" for i in range(n_obs)],  # All unique, no warn
        }),
        var=pd.DataFrame({
            "gene_name": [f"gene_{i}" for i in range(n_vars)],
            "highly_variable": np.random.choice([True, False], n_vars),
            "mean_expression": np.random.randn(n_vars).astype(np.float32),
        }),
    )

    # Add color annotations (matching)
    adata.uns["batch_colors"] = ["#FF6B6B", "#4ECDC4"]
    adata.uns["cell_type_colors"] = ["#FF6B6B", "#4ECDC4", "#45B7D1"]

    # Add color annotations (mismatched - should warn)
    adata.uns["mismatched_colors"] = ["#FF0000"]  # Not matching any categorical

    # Add nested AnnData
    inner_adata = AnnData(np.zeros((10, 5)))
    inner_adata.obs["inner_col"] = list(range(10))
    adata.uns["nested_adata"] = inner_adata

    # Add various uns types
    adata.uns["neighbors"] = {"params": {"n_neighbors": 15, "method": "umap"}}
    adata.uns["pca"] = {
        "variance": np.random.randn(50).astype(np.float32),
        "variance_ratio": np.random.randn(50).astype(np.float32),
    }
    adata.uns["string_value"] = "This is a test string"
    adata.uns["int_value"] = 42
    adata.uns["float_value"] = 3.14159
    adata.uns["list_value"] = [1, 2, 3, "mixed", {"nested": True}]

    # Add unserializable type (should warn)
    class CustomClass:
        def __repr__(self):
            return "CustomClass()"
    adata.uns["unserializable"] = CustomClass()

    # Add obsm/varm
    adata.obsm["X_pca"] = np.random.randn(n_obs, 50).astype(np.float32)
    adata.obsm["X_umap"] = np.random.randn(n_obs, 2).astype(np.float32)
    adata.obsm["X_tsne"] = np.random.randn(n_obs, 2).astype(np.float32)
    adata.obsm["cell_metadata"] = pd.DataFrame({
        "spatial_x": np.random.randn(n_obs),
        "spatial_y": np.random.randn(n_obs),
        "area": np.random.rand(n_obs) * 100,
    }, index=adata.obs_names)
    adata.varm["PCs"] = np.random.randn(n_vars, 50).astype(np.float32)

    # Add layers
    adata.layers["raw"] = sp.random(n_obs, n_vars, density=0.1, format="csr")
    adata.layers["normalized"] = np.random.randn(n_obs, n_vars).astype(np.float32)

    # Add obsp/varp
    adata.obsp["distances"] = sp.random(n_obs, n_obs, density=0.01, format="csr")
    adata.obsp["connectivities"] = sp.random(n_obs, n_obs, density=0.01, format="csr")
    adata.varp["gene_corr"] = sp.random(n_vars, n_vars, density=0.1, format="csr")

    return adata


def create_html_page(sections: list[tuple[str, str, str | None]]) -> str:
    """Create a full HTML page with multiple test cases.

    Parameters
    ----------
    sections
        List of (title, html_content, description) tuples.
        Description can be None for no description box.
    """
    html_parts = [
        """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="Content-Security-Policy" content="default-src 'self' 'unsafe-inline' 'unsafe-eval' data: https:; style-src 'self' 'unsafe-inline';">
    <title>AnnData _repr_html_ Visual Test</title>
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background: #f5f5f5;
        }
        h1 {
            color: #333;
            border-bottom: 2px solid #0d6efd;
            padding-bottom: 10px;
        }
        h2 {
            color: #555;
            margin-top: 40px;
        }
        .test-case {
            background: white;
            border-radius: 8px;
            padding: 20px;
            margin: 20px 0;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .test-case h3 {
            margin-top: 0;
            color: #0d6efd;
        }
        .description {
            color: #666;
            font-size: 0.9em;
            margin-bottom: 15px;
            padding: 10px;
            background: #f8f9fa;
            border-radius: 4px;
        }
        /* Dark mode toggle */
        .dark-mode-toggle {
            position: fixed;
            top: 20px;
            right: 20px;
            padding: 10px 20px;
            background: #333;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
        }
        body.dark-mode {
            background: #1a1a1a;
            color: #e0e0e0;
        }
        body.dark-mode h1, body.dark-mode h2 {
            color: #e0e0e0;
        }
        body.dark-mode .test-case {
            background: #2d2d2d;
        }
        body.dark-mode .description {
            background: #333;
            color: #aaa;
        }
    </style>
</head>
<body>
    <button class="dark-mode-toggle" onclick="document.body.classList.toggle('dark-mode')">
        Toggle Dark Mode
    </button>

    <h1>AnnData _repr_html_ Visual Test</h1>
    <p>This page displays various AnnData configurations to visually verify the HTML representation.</p>
"""
    ]

    for item in sections:
        title = item[0]
        html_content = item[1]
        description = item[2] if len(item) > 2 else None

        desc_html = ""
        if description:
            desc_html = f'<div class="description">{description}</div>'

        html_parts.append(f"""
    <div class="test-case">
        <h3>{title}</h3>
        {desc_html}
        <div class="repr-output">
            {html_content}
        </div>
    </div>
""")

    html_parts.append("""
</body>
</html>
""")

    return "".join(html_parts)


def strip_script_tags(html: str) -> str:
    """Remove <script>...</script> tags from HTML to simulate no-JS environment."""
    import re
    return re.sub(r'<script>.*?</script>', '', html, flags=re.DOTALL)


def main():
    """Generate visual test HTML file."""
    print("Generating visual test cases...")

    sections = []

    # Test 1: Full AnnData
    print("  1. Full AnnData with all features")
    adata_full = create_test_anndata()
    sections.append((
        "1. Full AnnData (all features)",
        adata_full._repr_html_(),
    ))

    # Test 2: Empty AnnData
    print("  2. Empty AnnData")
    adata_empty = AnnData()
    sections.append((
        "2. Empty AnnData",
        adata_empty._repr_html_(),
    ))

    # Test 3: Minimal AnnData
    print("  3. Minimal AnnData (just X)")
    adata_minimal = AnnData(np.zeros((10, 5)))
    sections.append((
        "3. Minimal AnnData (just X matrix)",
        adata_minimal._repr_html_(),
    ))

    # Test 4: View
    print("  4. AnnData View")
    view = adata_full[0:20, 0:10]
    sections.append((
        "4. AnnData View (subset)",
        view._repr_html_(),
    ))

    # Test 5: Dense matrix
    print("  5. Dense matrix")
    adata_dense = AnnData(np.random.randn(50, 30).astype(np.float32))
    adata_dense.obs["cluster"] = pd.Categorical(["A", "B", "C", "D", "E"] * 10)
    adata_dense.uns["cluster_colors"] = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
    sections.append((
        "5. Dense Matrix with Categories",
        adata_dense._repr_html_(),
    ))

    # Test 6: Many columns (collapsed sections)
    print("  6. Many columns (tests folding)")
    adata_many = AnnData(np.zeros((20, 10)))
    for i in range(15):
        adata_many.obs[f"column_{i}"] = list(range(20))
    for i in range(12):
        adata_many.obsm[f"X_embedding_{i}"] = np.random.randn(20, 2).astype(np.float32)
    sections.append((
        "6. Many Columns (tests auto-folding)",
        adata_many._repr_html_(),
    ))

    # Test 7: Special characters
    print("  7. Special characters in names")
    adata_special = AnnData(np.zeros((5, 3)))
    adata_special.obs["column<with>html"] = list(range(5))
    adata_special.obs["column&ampersand"] = list(range(5))
    adata_special.uns["key\"with'quotes"] = "value"
    adata_special.uns["unicode_日本語"] = "japanese"
    sections.append((
        "7. Special Characters (XSS/Unicode test)",
        adata_special._repr_html_(),
    ))

    # Test 8: Dask array (if available)
    if HAS_DASK:
        print("  8. Dask array")
        X_dask = da.zeros((100, 50), chunks=(10, 50))
        adata_dask = AnnData(X_dask)
        sections.append((
            "8. Dask Array (lazy/chunked)",
            adata_dask._repr_html_(),
        ))

    # Test 9: Backed AnnData (H5AD file)
    print("  9. Backed AnnData (H5AD file)")
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp_path = tmp.name
    adata_backed = None
    try:
        adata_to_save = AnnData(np.random.randn(50, 20).astype(np.float32))
        adata_to_save.obs["cluster"] = pd.Categorical(["A", "B", "C"] * 16 + ["A", "B"])
        adata_to_save.var["gene_name"] = [f"gene_{i}" for i in range(20)]
        adata_to_save.write_h5ad(tmp_path)
        adata_backed = ad.read_h5ad(tmp_path, backed="r")
        sections.append((
            "9. Backed AnnData (H5AD file)",
            adata_backed._repr_html_(),
            f"File path: {tmp_path}. Shows badge with format and status, plus file path.",
        ))
    finally:
        if adata_backed is not None:
            adata_backed.file.close()
        os.unlink(tmp_path)

    # Test 10: Nested AnnData at depth
    print("  10. Deeply nested AnnData")
    inner3 = AnnData(np.zeros((3, 2)))
    inner2 = AnnData(np.zeros((5, 3)))
    inner2.uns["level3"] = inner3
    inner1 = AnnData(np.zeros((10, 5)))
    inner1.uns["level2"] = inner2
    outer = AnnData(np.zeros((20, 10)))
    outer.uns["level1"] = inner1
    sections.append((
        "10. Deeply Nested AnnData (tests max depth)",
        outer._repr_html_(),
    ))

    # Test 11: Many categories (tests truncation and wrap button)
    # Default max_categories is 100, but we set it to 20 here to test truncation
    print("  11. Many categories (tests category truncation)")
    adata_many_cats = AnnData(np.zeros((100, 10)))
    # 30 categories - with max_categories=20 should show first 20 + '...+10'
    many_cat_values = [f"type_{i}" for i in range(30)] * (100 // 30) + [f"type_{i}" for i in range(100 % 30)]
    adata_many_cats.obs["cell_type"] = pd.Categorical(many_cat_values)
    # Add colors for the categories
    adata_many_cats.uns["cell_type_colors"] = [
        "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
        "#ffff33", "#a65628", "#f781bf", "#999999", "#66c2a5",
        "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f",
        "#e5c494", "#b3b3b3", "#1b9e77", "#d95f02", "#7570b3",
        "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666",
        "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
    ]
    # Also add a column with exactly 20 categories
    adata_many_cats.obs["batch"] = pd.Categorical([f"batch_{i}" for i in range(20)] * 5)
    adata_many_cats.uns["batch_colors"] = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
        "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
    ]
    # Use lower max_categories (default is 100) to demonstrate truncation
    original_max_cats = ad.settings.repr_html_max_categories
    ad.settings.repr_html_max_categories = 20
    sections.append((
        "11. Many Categories (tests truncation)",
        adata_many_cats._repr_html_(),
        "Default max_categories is 100; set to 20 here to test truncation. "
        "cell_type has 30 categories (shows 20 + '...+10'). Click ⋯ button to expand.",
    ))
    ad.settings.repr_html_max_categories = original_max_cats

    # Test 12: Uns value previews and custom TypeFormatter
    print("  12. Uns value previews and type hints")

    # Register a custom TypeFormatter for tagged data in uns
    @register_formatter
    class AnalysisHistoryFormatter(TypeFormatter):
        """Example TypeFormatter for analysis history data with embedded type hint."""
        priority = 100  # High priority to check before fallback

        def can_format(self, obj):
            hint, _ = extract_uns_type_hint(obj)
            return hint == "example.history"

        def format(self, obj, context):
            import json
            _hint, value = extract_uns_type_hint(obj)

            # Parse JSON if string, otherwise use as-is
            if isinstance(value, str):
                try:
                    data = json.loads(value)
                except json.JSONDecodeError:
                    data = {"raw": value}
            else:
                data = value if isinstance(value, dict) else {"data": value}

            # Build a rich HTML preview
            runs = data.get("runs", [])
            params = data.get("params", {})

            html_parts = ['<div style="font-size:11px;">']
            if runs:
                html_parts.append(f'<strong>{len(runs)} runs</strong>')
            if params:
                param_str = ", ".join(f"{k}={v}" for k, v in list(params.items())[:3])
                if len(params) > 3:
                    param_str += "..."
                html_parts.append(f' · params: {param_str}')
            html_parts.append('</div>')

            return FormattedOutput(
                type_name="analysis history",
                html_content="".join(html_parts),
            )

    adata_uns = AnnData(np.zeros((10, 5)))
    # Simple types with previews
    adata_uns.uns["string_param"] = "A short string value"
    adata_uns.uns["long_string"] = "This is a very long string that should be truncated in the preview because it exceeds the maximum length allowed for display in the meta column"
    adata_uns.uns["int_param"] = 42
    adata_uns.uns["float_param"] = 3.14159265359
    adata_uns.uns["bool_param"] = True
    adata_uns.uns["none_param"] = None
    adata_uns.uns["small_list"] = [1, 2, 3]
    adata_uns.uns["small_dict"] = {"a": 1, "b": 2}
    adata_uns.uns["larger_dict"] = {"key1": "val1", "key2": "val2", "key3": "val3", "key4": "val4", "key5": "val5"}

    # Type hint WITH registered renderer (shows custom HTML)
    adata_uns.uns["analysis_history"] = {
        "__anndata_repr__": "example.history",
        "runs": [{"id": 1}, {"id": 2}, {"id": 3}],
        "params": {"method": "umap", "n_neighbors": 15, "metric": "euclidean"},
    }

    # Type hint WITHOUT registered renderer (shows fallback with import hint)
    adata_uns.uns["unregistered_data"] = {
        "__anndata_repr__": "otherpackage.custom_type",
        "data": {"some": "data", "values": [1, 2, 3]},
    }
    # String format type hint (also unregistered)
    adata_uns.uns["string_hint"] = "__anndata_repr__:otherpackage.config::{'setting': 'value'}"

    sections.append((
        "12. Uns Value Previews and Type Hints",
        adata_uns._repr_html_(),
        "Shows: (1) preview values for simple types, (2) 'analysis_history' with registered TypeFormatter "
        "(shows '3 runs · params: method=umap...'), (3) unregistered type hints show 'import X to enable' message.",
    ))

    # Test 13: No JavaScript (graceful degradation)
    print("  13. No JavaScript (graceful degradation)")
    adata_nojs = AnnData(np.random.randn(30, 15).astype(np.float32))
    adata_nojs.obs["group"] = pd.Categorical(["X", "Y", "Z"] * 10)
    adata_nojs.uns["group_colors"] = ["#e41a1c", "#377eb8", "#4daf4a"]
    for i in range(8):
        adata_nojs.obs[f"metric_{i}"] = np.random.randn(30)
    adata_nojs.obsm["X_pca"] = np.random.randn(30, 10).astype(np.float32)
    adata_nojs.layers["raw"] = np.random.randn(30, 15).astype(np.float32)
    # Add a DataFrame with many columns to test column list wrapping without JS
    adata_nojs.obsm["cell_measurements"] = pd.DataFrame({
        "area": np.random.rand(30) * 500,
        "perimeter": np.random.rand(30) * 100,
        "circularity": np.random.rand(30),
        "eccentricity": np.random.rand(30),
        "solidity": np.random.rand(30),
        "extent": np.random.rand(30),
        "major_axis_length": np.random.rand(30) * 50,
        "minor_axis_length": np.random.rand(30) * 30,
        "orientation": np.random.rand(30) * 180,
        "mean_intensity": np.random.rand(30) * 255,
        "max_intensity": np.random.rand(30) * 255,
        "min_intensity": np.random.rand(30) * 50,
        "std_intensity": np.random.rand(30) * 30,
        "centroid_x": np.random.randn(30) * 100,
        "centroid_y": np.random.randn(30) * 100,
        "bbox_area": np.random.rand(30) * 600,
        "convex_area": np.random.rand(30) * 550,
        "euler_number": np.random.randint(-2, 3, 30),
        "equivalent_diameter": np.random.rand(30) * 25,
        "filled_area": np.random.rand(30) * 500,
    }, index=adata_nojs.obs_names)
    # Strip script tags to simulate no-JS environment
    nojs_html = strip_script_tags(adata_nojs._repr_html_())
    sections.append((
        "13. No JavaScript (graceful degradation)",
        nojs_html,
        "This example has script tags removed to simulate environments where JS is disabled. "
        "All content should be visible, sections should be expanded, category lists and "
        "DataFrame column lists should wrap naturally to multiple lines, and interactive buttons "
        "(fold icons, copy buttons, search, expand, wrap toggle) should be hidden. "
        "The obsm 'cell_measurements' DataFrame has 20 columns to test column list wrapping.",
    ))

    # Test 14: Custom sections example using TreeData (if available)
    if HAS_TREEDATA:
        print("  14. Custom Sections (TreeData example)")
        tdata = create_test_treedata()
        sections.append((
            "14. Custom Sections (TreeData example)",
            tdata._repr_html_(),
            "Demonstrates how to add custom sections using SectionFormatter. "
            "This example uses <a href='https://treedata.readthedocs.io/en/latest/' target='_blank'>TreeData</a> "
            "to show 'obst' and 'vart' sections "
            "(<a href='https://github.com/scverse/ecosystem-packages/pull/282' target='_blank'>scverse ecosystem PR</a>). "
            "The 'obst' section appears after 'obsm' and 'vart' after 'varm'. "
            "Click the expand button (▸) to see an SVG tree visualization. "
            "Trees with >30 leaves show a text message instead of the full preview.",
        ))

    # Test 15: Expandable DataFrame in obsm
    print("  15. Expandable DataFrame in obsm")
    # Enable DataFrame expansion for this test
    original_expand = ad.settings.repr_html_dataframe_expand
    ad.settings.repr_html_dataframe_expand = True
    try:
        adata_df = AnnData(np.random.randn(30, 10).astype(np.float32))
        adata_df.obs["group"] = pd.Categorical(["A", "B", "C"] * 10)
        # Add a wide DataFrame to obsm with many columns
        adata_df.obsm["spatial_metrics"] = pd.DataFrame({
            "x_centroid": np.random.randn(30) * 100,
            "y_centroid": np.random.randn(30) * 100,
            "area": np.random.rand(30) * 500,
            "perimeter": np.random.rand(30) * 100,
            "circularity": np.random.rand(30),
            "eccentricity": np.random.rand(30),
            "solidity": np.random.rand(30),
            "extent": np.random.rand(30),
            "major_axis": np.random.rand(30) * 50,
            "minor_axis": np.random.rand(30) * 30,
            "orientation": np.random.rand(30) * 180,
            "intensity_mean": np.random.rand(30) * 255,
        }, index=adata_df.obs_names)
        adata_df.obsm["X_pca"] = np.random.randn(30, 5).astype(np.float32)
        sections.append((
            "15. Expandable DataFrame in obsm",
            adata_df._repr_html_(),
            "When <code>anndata.settings.repr_html_dataframe_expand = True</code>, "
            "DataFrames in obsm/varm show an 'Expand' button. Click to see pandas <code>_repr_html_()</code> output "
            "(styled table with zebra striping and hover). Configure pandas display options: "
            "<code>pd.set_option('display.max_rows', 10)</code>. "
            "Column names are shown in the rightmost column (meta column).",
        ))
    finally:
        ad.settings.repr_html_dataframe_expand = original_expand

    # Test 16: Very long field names
    print("  16. Very long field names")
    adata_long = AnnData(np.random.randn(20, 10).astype(np.float32))
    # Add columns with very long names to test field name column width calculation
    adata_long.obs["short"] = np.random.randn(20)
    adata_long.obs["this_is_a_moderately_long_column_name"] = np.random.randn(20)
    adata_long.obs["this_is_an_extremely_long_column_name_that_should_test_the_max_width_setting"] = np.random.randn(20)
    adata_long.obs["cell_type_annotation_from_automated_classifier_v2"] = pd.Categorical(["A", "B"] * 10)
    adata_long.obsm["X_pca_computed_with_highly_variable_genes_batch_corrected"] = np.random.randn(20, 5).astype(np.float32)
    adata_long.uns["preprocessing_parameters_for_normalization_and_scaling"] = {"method": "log1p", "scale": True}
    adata_long.layers["raw_counts_before_any_preprocessing_steps"] = np.random.randn(20, 10).astype(np.float32)
    sections.append((
        "16. Very Long Field Names",
        adata_long._repr_html_(),
        "Tests the dynamic field name column width calculation. The longest field name is "
        "'this_is_an_extremely_long_column_name_that_should_test_the_max_width_setting' (77 chars). "
        "The name column width should expand to fit longer names but be capped by "
        "<code>repr_html_max_field_width</code> (default: 250px). Names exceeding the max width "
        "should show ellipsis (...) on hover.",
    ))

    # Generate HTML file
    output_path = Path(__file__).parent / "repr_html_visual_test.html"
    html_content = create_html_page(sections)
    output_path.write_text(html_content)

    print(f"\nVisual test file generated: {output_path}")
    print("Open this file in a browser to inspect the HTML representation.")


if __name__ == "__main__":
    main()
