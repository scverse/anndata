"""
Visual inspection script for AnnData HTML representation.

Run this script to generate an HTML file that can be opened in a browser
to visually inspect the _repr_html_ output.

Usage:
    python tests/visual_inspect_repr_html.py

Then open tests/repr_html_visual_test.html in your browser.
"""

from __future__ import annotations

import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

import anndata as ad

# Suppress anndata warning about string index transformation (not relevant for visual tests)
from anndata._warnings import ImplicitModificationWarning

warnings.filterwarnings(
    "ignore",
    message="Transforming to str index",
    category=ImplicitModificationWarning,
)
from anndata import AnnData  # noqa: E402
from anndata._repr import (  # noqa: E402
    FormattedOutput,
    TypeFormatter,
    extract_uns_type_hint,
    register_formatter,
)

# Check optional dependencies
try:
    import dask.array as da

    HAS_DASK = True
except ImportError:
    HAS_DASK = False

try:
    import xarray  # noqa: F401

    from anndata.experimental import read_lazy

    HAS_XARRAY = True
except ImportError:
    HAS_XARRAY = False

try:
    import networkx as nx
    from treedata import TreeData

    from anndata._repr import (
        FormattedEntry,
        FormattedOutput,
        FormatterContext,
        SectionFormatter,
        register_formatter,
    )

    HAS_TREEDATA = True

    def _render_tree_svg(
        tree: nx.DiGraph, max_leaves: int = 30, width: int = 300, height: int = 150
    ) -> str:
        """Render a tree as an SVG visualization.

        Uses a simple top-down layout similar to pycea's approach but generates SVG.
        """
        # Find root and leaves
        roots = [n for n in tree.nodes() if tree.in_degree(n) == 0]
        if not roots:
            return (
                "<span style='color:#888;font-size:11px;'>Invalid tree (no root)</span>"
            )
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
            svg_parts.append(
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}"/>'
            )

        svg_parts.append("</svg>")
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

except (ImportError, AttributeError):
    # AttributeError can occur on Python 3.14+ with incompatible networkx versions
    HAS_TREEDATA = False

# Check for MuData
try:
    from mudata import MuData

    from anndata._repr import (
        FormattedEntry,
        FormattedOutput,
        FormatterContext,
        SectionFormatter,
        register_formatter,
    )
    from anndata._repr.html import generate_repr_html
    from anndata._repr.utils import format_number

    HAS_MUDATA = True

    # Suppress MuData's internal mapping attributes using a SectionFormatter
    # that handles multiple sections and returns empty (suppresses them)
    @register_formatter
    class MuDataInternalSectionsFormatter(SectionFormatter):
        """Suppress MuData's internal mapping attributes."""

        section_names = ("obsmap", "varmap", "axis")

        @property
        def section_name(self) -> str:
            return self.section_names[0]  # Primary name for compatibility

        def should_show(self, obj) -> bool:
            return False  # Never show these sections

        def get_entries(self, obj, context):
            return []  # No entries

    # Register a SectionFormatter for MuData's .mod section
    # This allows generate_repr_html() to work directly on MuData objects
    @register_formatter
    class ModSectionFormatter(SectionFormatter):
        """
        SectionFormatter for MuData's .mod attribute.

        This demonstrates how external packages (like mudata) can extend
        anndata's HTML repr to add new sections. The .mod section contains
        AnnData objects for each modality, similar to how .uns can contain
        nested AnnData objects.
        """

        section_name = "mod"
        priority = 200  # High priority to show before other sections

        @property
        def after_section(self) -> str:
            return "X"  # Show right after X (before obs)

        @property
        def doc_url(self) -> str:
            return "https://mudata.readthedocs.io/en/latest/api/generated/mudata.MuData.html"

        @property
        def tooltip(self) -> str:
            return "Modalities (MuData)"

        def should_show(self, obj) -> bool:
            return hasattr(obj, "mod") and len(obj.mod) > 0

        def get_entries(self, obj, context: FormatterContext) -> list[FormattedEntry]:
            entries = []
            for mod_name, adata in obj.mod.items():
                shape_str = (
                    f"{format_number(adata.n_obs)} × {format_number(adata.n_vars)}"
                )
                # Generate nested HTML for expandable content
                can_expand = context.depth < context.max_depth
                nested_html = None
                if can_expand:
                    nested_html = generate_repr_html(
                        adata,
                        depth=context.depth + 1,
                        max_depth=context.max_depth,
                        show_header=True,
                        show_search=False,
                    )
                output = FormattedOutput(
                    type_name=f"AnnData ({shape_str})",
                    css_class="dtype-anndata",
                    tooltip=f"Modality: {mod_name}",
                    details={
                        "n_obs": adata.n_obs,
                        "n_vars": adata.n_vars,
                    },
                    html_content=nested_html,
                    is_expandable=can_expand,
                    is_serializable=True,
                )
                entries.append(FormattedEntry(key=mod_name, output=output))
            return entries

except ImportError:
    HAS_MUDATA = False
    MuData = None  # type: ignore[assignment,misc]


# =============================================================================
# SpatialData Example: Building custom _repr_html_ using anndata's building blocks
# =============================================================================
# This demonstrates how packages like SpatialData can create their own _repr_html_
# while reusing anndata's CSS, JavaScript, and rendering helpers.
#
# KEY BUILDING BLOCKS USED:
#   - get_css()                  : Reuse anndata's CSS (dark mode, styling)
#   - get_javascript(id)         : Reuse anndata's JS (fold, search, copy)
#   - render_section()           : Render a collapsible section
#   - render_formatted_entry()   : Render a table row
#   - FormattedEntry/Output      : Data classes for entry configuration
#   - generate_repr_html()       : Embed nested AnnData objects
#   - FormatterRegistry          : (Optional) Allow third-party extensions

try:
    import uuid

    from anndata._repr import (
        FormattedEntry,
        FormattedOutput,
        FormatterContext,
        FormatterRegistry,
        SectionFormatter,
        TypeFormatter,
        escape_html,
        format_number,
        get_css,
        get_javascript,
        render_badge,
        render_formatted_entry,
        render_search_box,
        render_section,
    )
    from anndata._repr.html import generate_repr_html

    HAS_SPATIALDATA_EXAMPLE = True

    # =========================================================================
    # MockSpatialData: Minimal example of custom _repr_html_
    # =========================================================================

    class MockSpatialData:
        """
        Mock SpatialData demonstrating custom _repr_html_ with anndata's building blocks.

        This is a simplified example showing the essential pattern. A real
        implementation would have more complex data structures.
        """

        def __init__(
            self,
            *,
            images: dict | None = None,
            labels: dict | None = None,
            points: dict | None = None,
            shapes: dict | None = None,
            tables: dict | None = None,  # Contains AnnData objects
            coordinate_systems: list | None = None,
            path: str | None = None,
        ):
            self.images = images or {}
            self.labels = labels or {}
            self.points = points or {}
            self.shapes = shapes or {}
            self.tables = tables or {}
            self.coordinate_systems = coordinate_systems or []
            self.path = path

        def _repr_html_(self) -> str:
            """
            Build HTML using anndata's building blocks.

            Pattern:
                1. get_css() - include styling
                2. Container div with unique ID
                3. Custom header (optional)
                4. Coordinate systems preview (like obs_names/var_names in AnnData)
                5. Sections using render_section() + render_formatted_entry()
                6. Custom sections via FormatterRegistry (optional)
                7. get_javascript(id) - include interactivity
            """
            container_id = f"spatialdata-{uuid.uuid4().hex[:8]}"
            parts = []

            # --- STEP 1: Include anndata's CSS ---
            parts.append(get_css())

            # --- STEP 2: Container with anndata-repr class ---
            parts.append(
                f'<div class="anndata-repr" id="{container_id}" data-depth="0" '
                f'style="--anndata-name-col-width: 150px; --anndata-type-col-width: 300px;">'
            )

            # --- STEP 3: Custom header (SpatialData has no shape) ---
            parts.append(self._build_header(container_id))

            # --- STEP 4: Coordinate systems preview (alternative to obs_names/var_names) ---
            parts.append(self._build_coordinate_systems_preview())

            # --- STEP 5: Sections using render_section() ---
            parts.append('<div class="adata-sections">')
            parts.append(self._build_images_section())
            parts.append(self._build_labels_section())
            parts.append(self._build_points_section())
            parts.append(self._build_shapes_section())
            parts.append(self._build_tables_section())  # Nested AnnData
            # --- STEP 6: Custom sections from FormatterRegistry ---
            parts.append(self._build_custom_sections())
            parts.append("</div>")

            parts.append("</div>")

            # --- STEP 7: Include anndata's JavaScript ---
            parts.append(get_javascript(container_id))

            return "\n".join(parts)

        def _build_header(self, container_id: str) -> str:
            """Custom header - shows 'SpatialData' with Zarr badge and file path."""
            parts = ['<div class="anndata-hdr">']
            parts.append('<span class="adata-type">SpatialData</span>')

            # Zarr badge using render_badge() helper
            if self.path:
                parts.append(
                    render_badge("Zarr", "adata-badge-backed", "Backed by Zarr storage")
                )
                parts.append(
                    f'<span class="adata-file-path" style="font-family:ui-monospace,monospace;'
                    f'font-size:11px;color:var(--anndata-text-secondary, #6c757d);">'
                    f"{escape_html(self.path)}</span>"
                )

            # Search box using render_search_box() helper
            parts.append('<span style="flex-grow:1;"></span>')
            parts.append(render_search_box(container_id))
            parts.append("</div>")
            return "\n".join(parts)

        def _build_coordinate_systems_preview(self) -> str:
            """
            Build coordinate systems preview - SpatialData's equivalent to obs_names/var_names.

            Simple list of coordinate system names with element details in tooltips.
            """
            if not self.coordinate_systems:
                return ""

            # Collect element names for tooltips
            all_elements = []
            if self.images:
                all_elements.extend([f"{k} (Images)" for k in self.images])
            if self.labels:
                all_elements.extend([f"{k} (Labels)" for k in self.labels])
            if self.points:
                all_elements.extend([f"{k} (Points)" for k in self.points])
            if self.shapes:
                all_elements.extend([f"{k} (Shapes)" for k in self.shapes])

            elements_str = ", ".join(all_elements) if all_elements else "no elements"

            # Build simple inline list
            parts = ['<div class="adata-index-preview" style="padding:2px 8px;">']
            parts.append(
                '<span style="color:var(--anndata-text-secondary, #6c757d);'
                'font-size:12px;">coordinate_systems: </span>'
            )

            # Render coordinate systems as simple badges with tooltips
            cs_parts = []
            for cs_name in self.coordinate_systems:
                tooltip = f"Elements: {elements_str}"
                cs_parts.append(
                    f'<span title="{escape_html(tooltip)}" style="'
                    f"font-family:ui-monospace,monospace;font-size:11px;"
                    f'color:var(--anndata-accent, #0d6efd);cursor:help;">'
                    f"'{escape_html(cs_name)}'</span>"
                )

            parts.append(", ".join(cs_parts))
            parts.append("</div>")
            return "".join(parts)

        def _build_images_section(self) -> str:
            """
            Build images section using render_section() + render_formatted_entry().

            This is the core pattern: create FormattedEntry objects and render them.
            """
            rows = []
            for name, info in self.images.items():
                # Build meta content (dimensions info) for the META column
                dims_str = ", ".join(info.get("dims", ["y", "x"]))
                meta = f'<span class="adata-meta-info">[{dims_str}]</span>'

                # Create a FormattedEntry with FormattedOutput
                entry = FormattedEntry(
                    key=name,
                    output=FormattedOutput(
                        type_name=f"DataArray {info['shape']} {info['dtype']}",
                        css_class="dtype-ndarray",
                        meta_content=meta,  # Content in meta column (rightmost)
                    ),
                )
                # render_formatted_entry() creates the table row HTML
                rows.append(render_formatted_entry(entry))

            # render_section() wraps rows in a collapsible section
            return render_section(
                "images",
                "\n".join(rows),
                n_items=len(self.images),
                tooltip="Image data (xarray.DataArray)",
            )

        def _build_labels_section(self) -> str:
            """Build labels section - same pattern as images."""
            rows = []
            for name, info in self.labels.items():
                dims_str = ", ".join(info.get("dims", ["y", "x"]))
                meta = f'<span class="adata-meta-info">[{dims_str}]</span>'

                entry = FormattedEntry(
                    key=name,
                    output=FormattedOutput(
                        type_name=f"Labels {info['shape']} {info['dtype']}",
                        css_class="dtype-ndarray",
                        meta_content=meta,
                    ),
                )
                rows.append(render_formatted_entry(entry))

            return render_section(
                "labels",
                "\n".join(rows),
                n_items=len(self.labels),
                tooltip="Segmentation masks (xarray.DataArray)",
            )

        def _build_points_section(self) -> str:
            """Build points section."""
            rows = []
            for name, info in self.points.items():
                meta = f'<span class="adata-meta-info">{info["n_dims"]}D coordinates</span>'

                entry = FormattedEntry(
                    key=name,
                    output=FormattedOutput(
                        type_name=f"dask.DataFrame ({format_number(info['n_points'])} × {info['n_dims']})",
                        css_class="dtype-dataframe",
                        meta_content=meta,
                    ),
                )
                rows.append(render_formatted_entry(entry))

            return render_section(
                "points",
                "\n".join(rows),
                n_items=len(self.points),
                tooltip="Point annotations (dask.DataFrame)",
            )

        def _build_shapes_section(self) -> str:
            """Build shapes section."""
            rows = []
            for name, info in self.shapes.items():
                meta = f'<span class="adata-meta-info">{info["geometry_type"]}</span>'

                entry = FormattedEntry(
                    key=name,
                    output=FormattedOutput(
                        type_name=f"GeoDataFrame ({format_number(info['n_shapes'])} shapes)",
                        css_class="dtype-dataframe",
                        meta_content=meta,
                    ),
                )
                rows.append(render_formatted_entry(entry))

            return render_section(
                "shapes",
                "\n".join(rows),
                n_items=len(self.shapes),
                tooltip="Vector shapes (geopandas.GeoDataFrame)",
            )

        def _build_tables_section(self) -> str:
            """
            Build tables section with NESTED AnnData objects.

            Uses generate_repr_html() to embed full AnnData representations
            that are expandable with all standard features.
            """
            rows = []
            for name, adata in self.tables.items():
                # generate_repr_html() creates nested AnnData HTML
                nested_html = generate_repr_html(
                    adata,
                    depth=1,  # Nested level
                    max_depth=3,
                    show_header=True,
                    show_search=False,
                )

                # FormattedOutput with is_expandable=True makes it collapsible
                entry = FormattedEntry(
                    key=name,
                    output=FormattedOutput(
                        type_name=f"AnnData ({adata.n_obs} × {adata.n_vars})",
                        css_class="dtype-anndata",
                        html_content=nested_html,
                        is_expandable=True,  # Makes the nested content collapsible
                    ),
                )
                rows.append(render_formatted_entry(entry))

            return render_section(
                "tables",
                "\n".join(rows),
                n_items=len(self.tables),
                tooltip="Annotation tables (AnnData)",
            )

        def _build_custom_sections(self) -> str:
            """
            Render custom sections from FormatterRegistry.

            This demonstrates how third-party packages can add new sections
            by registering SectionFormatters with spatialdata_formatter_registry.
            """
            parts = []
            context = FormatterContext()

            for (
                section_name
            ) in spatialdata_formatter_registry.get_registered_sections():
                formatter = spatialdata_formatter_registry.get_section_formatter(
                    section_name
                )
                if formatter is None or not formatter.should_show(self):
                    continue

                entries = formatter.get_entries(self, context)
                if not entries:
                    continue

                rows = [render_formatted_entry(entry) for entry in entries]
                section_html = render_section(
                    formatter.section_name,
                    "\n".join(rows),
                    n_items=len(entries),
                    tooltip=getattr(formatter, "tooltip", ""),
                )
                parts.append(section_html)

            return "\n".join(parts)

    # =========================================================================
    # OPTIONAL: FormatterRegistry for third-party extensibility
    # =========================================================================
    # SpatialData can create its own registry to allow plugins to add
    # custom type formatters or new sections. This mirrors anndata's pattern.

    # Create SpatialData's own formatter registry
    spatialdata_formatter_registry = FormatterRegistry()

    # Example: TypeFormatter for custom value rendering
    class DataTreeFormatter(TypeFormatter):
        """Example: format xarray DataTree objects."""

        priority = 100

        def can_format(self, obj) -> bool:
            return isinstance(obj, dict) and "shape" in obj and "dtype" in obj

        def format(self, obj, context: FormatterContext) -> FormattedOutput:
            return FormattedOutput(
                type_name=f"DataTree {obj['shape']} {obj['dtype']}",
                css_class="dtype-ndarray",
            )

    spatialdata_formatter_registry.register_type_formatter(DataTreeFormatter())

    # Example: SectionFormatter to add new sections
    class TransformsSectionFormatter(SectionFormatter):
        """Example: add a 'transforms' section."""

        section_name = "transforms"

        def should_show(self, obj) -> bool:
            return (
                hasattr(obj, "coordinate_systems") and len(obj.coordinate_systems) > 1
            )

        def get_entries(self, obj, context: FormatterContext) -> list[FormattedEntry]:
            cs = list(obj.coordinate_systems)
            return [
                FormattedEntry(
                    key=f"{cs[i]} → {cs[i + 1]}",
                    output=FormattedOutput(type_name="Affine (3×3)"),
                )
                for i in range(len(cs) - 1)
            ]

    spatialdata_formatter_registry.register_section_formatter(
        TransformsSectionFormatter()
    )

    # =========================================================================
    # Test data factory
    # =========================================================================

    def create_test_spatialdata():
        """Create a mock SpatialData object for testing."""
        # Create nested AnnData tables
        cell_table = AnnData(
            np.random.randn(150, 30).astype(np.float32),
            obs=pd.DataFrame({
                "cell_type": pd.Categorical(["Tumor", "Immune", "Stromal"] * 50),
                "area": np.random.uniform(100, 1000, 150),
            }),
        )
        cell_table.obsm["spatial"] = np.random.randn(150, 2).astype(np.float32)

        transcript_table = AnnData(
            np.random.randn(80, 10).astype(np.float32),
            obs=pd.DataFrame({
                "gene": pd.Categorical(
                    np.random.choice([f"gene_{i}" for i in range(10)], 80)
                ),
            }),
        )

        return MockSpatialData(
            images={
                "raw_image": {
                    "shape": (3, 2048, 2048),
                    "dims": ("c", "y", "x"),
                    "dtype": "uint16",
                },
                "processed": {
                    "shape": (3, 1024, 1024),
                    "dims": ("c", "y", "x"),
                    "dtype": "float32",
                },
            },
            labels={
                "cell_segmentation": {
                    "shape": (2048, 2048),
                    "dims": ("y", "x"),
                    "dtype": "int32",
                },
                "nucleus_segmentation": {
                    "shape": (2048, 2048),
                    "dims": ("y", "x"),
                    "dtype": "int32",
                },
            },
            points={
                "transcripts": {"n_points": 50000, "n_dims": 3},
            },
            shapes={
                "cell_boundaries": {"n_shapes": 150, "geometry_type": "Polygon"},
                "roi_annotations": {"n_shapes": 5, "geometry_type": "Polygon"},
            },
            tables={
                "cell_annotations": cell_table,
                "transcript_counts": transcript_table,
            },
            coordinate_systems=["global", "aligned", "microscope"],
            path="/data/experiment_001.zarr",
        )

except (ImportError, AttributeError):
    HAS_SPATIALDATA_EXAMPLE = False


def create_test_mudata():
    """Create a comprehensive test MuData with multiple modalities."""
    if not HAS_MUDATA:
        return None

    np.random.seed(42)

    # RNA modality
    n_cells = 100
    n_genes = 50
    rna = AnnData(
        np.random.randn(n_cells, n_genes).astype(np.float32),
        obs=pd.DataFrame({
            "cell_type": pd.Categorical(
                ["T cell", "B cell", "NK cell"] * 33 + ["T cell"]
            ),
            "n_counts": np.random.randint(1000, 10000, n_cells),
        }),
        var=pd.DataFrame({
            "gene_name": [f"gene_{i}" for i in range(n_genes)],
            "highly_variable": np.random.choice([True, False], n_genes),
        }),
    )
    rna.uns["cell_type_colors"] = ["#e41a1c", "#377eb8", "#4daf4a"]
    rna.obsm["X_pca"] = np.random.randn(n_cells, 10).astype(np.float32)
    rna.obsm["X_umap"] = np.random.randn(n_cells, 2).astype(np.float32)
    rna.layers["raw"] = np.random.randn(n_cells, n_genes).astype(np.float32)

    # ATAC modality (same cells, different features)
    n_peaks = 30
    atac = AnnData(
        np.random.randn(n_cells, n_peaks).astype(np.float32),
        obs=pd.DataFrame({
            "peak_count": np.random.randint(500, 5000, n_cells),
            "tss_enrichment": np.random.uniform(2, 10, n_cells),
        }),
        var=pd.DataFrame({
            "peak_name": [f"peak_{i}" for i in range(n_peaks)],
            "chr": [f"chr{i % 22 + 1}" for i in range(n_peaks)],
        }),
    )
    atac.obsm["X_lsi"] = np.random.randn(n_cells, 15).astype(np.float32)

    # Protein modality (subset of cells)
    n_prot_cells = 80
    n_proteins = 20
    prot = AnnData(
        np.random.randn(n_prot_cells, n_proteins).astype(np.float32),
        obs=pd.DataFrame({
            "protein_count": np.random.randint(100, 1000, n_prot_cells),
        }),
        var=pd.DataFrame({
            "protein_name": [f"CD{i}" for i in range(n_proteins)],
            "isotype_control": [i < 3 for i in range(n_proteins)],
        }),
    )

    # Create MuData
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mdata = MuData({"rna": rna, "atac": atac, "prot": prot})

    # Add shared annotations
    mdata.uns["experiment"] = "multiome_sample_001"
    mdata.uns["processing_date"] = "2024-03-15"

    return mdata


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
        ("root", "clade_A"),
        ("root", "clade_B"),
        ("clade_A", "subA1"),
        ("clade_A", "subA2"),
        ("clade_B", "subB1"),
        ("clade_B", "subB2"),
    ])
    for i, name in enumerate(obs_names):
        parent = ["subA1", "subA2", "subB1", "subB2"][i % 4]
        obs_tree.add_edge(parent, name)

    # Create variable tree (gene ontology-like structure, >30 leaves)
    var_tree = nx.DiGraph()
    var_tree.add_edges_from([
        ("all_genes", "pathway_X"),
        ("all_genes", "pathway_Y"),
        ("all_genes", "pathway_Z"),
        ("pathway_X", "module_1"),
        ("pathway_X", "module_2"),
        ("pathway_Y", "module_3"),
        ("pathway_Y", "module_4"),
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
                ["T cell", "B cell", "NK cell"] * (n_obs // 3)
                + ["T cell"] * (n_obs % 3)
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

    # Add unserializable type in uns (should warn)
    class CustomClass:
        def __repr__(self):
            return "CustomClass()"

    adata.uns["unserializable"] = CustomClass()

    # Add obsm/varm
    adata.obsm["X_pca"] = np.random.randn(n_obs, 50).astype(np.float32)
    adata.obsm["X_umap"] = np.random.randn(n_obs, 2).astype(np.float32)
    adata.obsm["X_tsne"] = np.random.randn(n_obs, 2).astype(np.float32)
    adata.obsm["cell_metadata"] = pd.DataFrame(
        {
            "spatial_x": np.random.randn(n_obs),
            "spatial_y": np.random.randn(n_obs),
            "area": np.random.rand(n_obs) * 100,
        },
        index=adata.obs_names,
    )
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
    # Generate TOC entries
    toc_items = []
    for item in sections:
        title = item[0]
        # Create anchor ID from title
        anchor_id = title.lower().replace(" ", "-").replace("(", "").replace(")", "")
        anchor_id = "".join(c for c in anchor_id if c.isalnum() or c == "-")
        toc_items.append(f'<a href="#{anchor_id}">{title}</a>')

    toc_html = "\n            ".join(toc_items)

    html_parts = [
        f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="Content-Security-Policy" content="default-src 'self' 'unsafe-inline' 'unsafe-eval' data: https:; style-src 'self' 'unsafe-inline';">
    <title>AnnData _repr_html_ Visual Test</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            max-width: 1000px;
            margin-left: 230px;
            margin-right: 20px;
            padding: 20px;
            background: #f5f5f5;
        }}
        h1 {{
            color: #333;
            border-bottom: 2px solid #0d6efd;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #555;
            margin-top: 40px;
        }}
        .test-case {{
            background: white;
            border-radius: 8px;
            padding: 20px;
            margin: 20px 0;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .test-case h3 {{
            margin-top: 0;
            color: #0d6efd;
        }}
        .description {{
            color: #666;
            font-size: 0.9em;
            margin-bottom: 15px;
            padding: 10px;
            background: #f8f9fa;
            border-radius: 4px;
        }}
        /* TOC sidebar */
        .toc {{
            position: fixed;
            top: 20px;
            left: 20px;
            width: 180px;
            max-height: calc(100vh - 40px);
            overflow-y: auto;
            background: white;
            border-radius: 8px;
            padding: 15px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.15);
            font-size: 12px;
            z-index: 100;
        }}
        .toc a {{
            display: block;
            padding: 4px 0;
            color: #555;
            text-decoration: none;
            border-bottom: 1px solid #eee;
        }}
        .toc a:last-child {{
            border-bottom: none;
        }}
        .toc a:hover {{
            color: #0d6efd;
        }}
        /* Dark mode toggle */
        .dark-mode-toggle {{
            position: fixed;
            top: 20px;
            right: 20px;
            padding: 10px 20px;
            background: #333;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
        }}
        body.dark-mode {{
            background: #1a1a1a;
            color: #e0e0e0;
        }}
        body.dark-mode h1, body.dark-mode h2 {{
            color: #e0e0e0;
        }}
        body.dark-mode .test-case {{
            background: #2d2d2d;
        }}
        body.dark-mode .description {{
            background: #333;
            color: #aaa;
        }}
        body.dark-mode .toc {{
            background: #2d2d2d;
        }}
        body.dark-mode .toc a {{
            color: #aaa;
            border-bottom-color: #444;
        }}
        body.dark-mode .toc a:hover {{
            color: #6ea8fe;
        }}
    </style>
</head>
<body>
    <nav class="toc">
        {toc_html}
    </nav>

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

        # Create anchor ID from title (same logic as TOC generation)
        anchor_id = title.lower().replace(" ", "-").replace("(", "").replace(")", "")
        anchor_id = "".join(c for c in anchor_id if c.isalnum() or c == "-")

        desc_html = ""
        if description:
            desc_html = f'<div class="description">{description}</div>'

        html_parts.append(f"""
    <div id="{anchor_id}" class="test-case">
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

    return re.sub(r"<script>.*?</script>", "", html, flags=re.DOTALL)


def main():  # noqa: PLR0915, PLR0912
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
    adata_dense.uns["cluster_colors"] = [
        "#e41a1c",
        "#377eb8",
        "#4daf4a",
        "#984ea3",
        "#ff7f00",
    ]
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

    # Test 8: Dask array (if available) - demonstrates lazy loading safety
    if HAS_DASK:
        print("  8. Dask array (lazy loading safety)")
        # Create Dask arrays in multiple sections to show lazy handling
        X_dask = da.random.random((1000, 500), chunks=(100, 100))
        adata_dask = AnnData(X_dask)
        adata_dask.obs["cluster"] = pd.Categorical(["A", "B", "C"] * 333 + ["A"])
        adata_dask.var["gene_name"] = [f"gene_{i}" for i in range(500)]
        # Dask arrays in layers and obsm
        adata_dask.layers["counts"] = da.random.randint(
            0, 100, (1000, 500), chunks=(100, 100)
        )
        adata_dask.obsm["X_pca"] = da.random.random((1000, 50), chunks=(100, 50))
        adata_dask.varm["loadings"] = da.random.random((500, 50), chunks=(100, 50))
        sections.append((
            "8. Dask Arrays (Lazy Loading Safety)",
            adata_dask._repr_html_(),
            "<strong>No data is loaded from disk!</strong> The repr only reads metadata:<br>"
            "<ul>"
            "<li><code>X</code>: Shows shape (1,000 × 500), dtype, chunks — no <code>.compute()</code></li>"
            "<li><code>layers['counts']</code>: Same — reads <code>.shape</code>, <code>.dtype</code>, <code>.chunks</code></li>"
            "<li><code>obsm['X_pca']</code>: Shows (1,000 × 50) from <code>.shape</code> attribute</li>"
            "<li><code>varm['loadings']</code>: Shows (500 × 50) from <code>.shape</code> attribute</li>"
            "</ul>"
            "This is safe for terabyte-scale datasets — displaying the repr never triggers computation.",
        ))

    # Test 8b: Lazy AnnData (experimental) - fully lazy obs/var
    # Tests the lazy category loading behavior:
    # - Categorical columns with few categories: load and display categories
    # - Categorical columns with too many categories: show "(lazy)" to avoid loading
    # - Categorical columns with colors in uns: display color swatches
    # - Non-categorical columns: show "(lazy)" for all
    if HAS_XARRAY:
        print("  8b. Lazy AnnData (experimental read_lazy)")
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
            tmp_path = tmp.name
        adata_lazy = None
        h5_file = None
        try:
            import h5py

            # Create a comprehensive test file for lazy loading behavior
            adata_to_save = AnnData(
                sp.random(1000, 500, density=0.1, format="csr", dtype=np.float32)
            )

            # --- Categorical columns ---
            # 1. Small categorical WITH colors (should show categories + color dots)
            adata_to_save.obs["cell_type"] = pd.Categorical(
                np.random.choice(["T cell", "B cell", "Monocyte", "NK cell"], 1000)
            )
            adata_to_save.uns["cell_type_colors"] = [
                "#e41a1c",  # T cell - red
                "#377eb8",  # B cell - blue
                "#4daf4a",  # Monocyte - green
                "#984ea3",  # NK cell - purple
            ]

            # 2. Small categorical WITHOUT colors (should show categories only)
            adata_to_save.obs["cluster"] = pd.Categorical(
                np.random.choice(["C0", "C1", "C2", "C3", "C4"], 1000)
            )

            # 3. Medium categorical (50 cats) - will show truncation with max_lazy_categories=30
            medium_categories = [f"sample_{i}" for i in range(50)]
            adata_to_save.obs["sample_id"] = pd.Categorical(
                np.random.choice(medium_categories, 1000),
                categories=medium_categories,  # Ensure all 50 categories exist
            )

            # --- Non-categorical columns (all should show "(lazy)") ---
            adata_to_save.obs["n_genes"] = np.random.randint(500, 5000, 1000)
            adata_to_save.obs["total_counts"] = np.random.randint(1000, 50000, 1000)

            # --- var columns ---
            adata_to_save.var["gene_symbol"] = [f"GENE{i}" for i in range(500)]
            adata_to_save.var["highly_variable"] = np.random.choice([True, False], 500)
            adata_to_save.var["mean_expression"] = np.random.uniform(0, 10, 500)

            # --- obsm/varm ---
            adata_to_save.obsm["X_pca"] = np.random.randn(1000, 50).astype(np.float32)
            adata_to_save.obsm["X_umap"] = np.random.randn(1000, 2).astype(np.float32)
            adata_to_save.varm["PCs"] = np.random.randn(500, 50).astype(np.float32)

            # --- uns with array (to show dask array WITH size in uns) ---
            adata_to_save.uns["neighbors"] = {
                "connectivities_key": "connectivities",
                "distances_key": "distances",
            }
            adata_to_save.uns["pca_variance"] = np.random.rand(50).astype(np.float32)

            adata_to_save.write_h5ad(tmp_path)

            # Read with experimental lazy loading
            h5_file = h5py.File(tmp_path, "r")
            adata_lazy = read_lazy(h5_file)

            # Use setting to demonstrate truncation behavior (default is 100)
            # - cell_type (4 cats): all shown
            # - cluster (5 cats): all shown
            # - sample_id (50 cats): first 30 shown + "...+20"
            original_max_lazy_cats = ad.settings.repr_html_max_lazy_categories
            ad.settings.repr_html_max_lazy_categories = 30
            custom_lazy_html = adata_lazy._repr_html_()
            ad.settings.repr_html_max_lazy_categories = original_max_lazy_cats

            sections.append((
                "8b. Lazy AnnData (Experimental)",
                custom_lazy_html,
                "<code>anndata.experimental.read_lazy()</code><br>"
                "<p style='margin: 5px 0;'><b>Design: Truly lazy category loading</b></p>"
                "<p style='margin: 5px 0; font-size: 0.9em;'>"
                "Categories are read directly from storage without loading the full column data. "
                "This example sets <code>ad.settings.repr_html_max_lazy_categories = 30</code> "
                "(default: 100) to illustrate truncation:</p>"
                "<ul style='margin: 5px 0; padding-left: 20px; font-size: 0.9em;'>"
                "<li><b>cell_type</b> (4 cats): all shown with colors</li>"
                "<li><b>cluster</b> (5 cats): all shown</li>"
                '<li><b>sample_id</b> (50 cats): first 30 shown + "...+20"</li>'
                "</ul>"
                "<p style='margin: 5px 0;'><b>What the repr avoids loading:</b></p>"
                "<ul style='margin: 5px 0; padding-left: 20px;'>"
                "<li>Numeric data (dask arrays not computed)</li>"
                "<li>Category codes (only category labels are read)</li>"
                "<li>Categories beyond <code>max_lazy_categories</code> limit</li>"
                "<li>Colors for truncated categories (only loads colors for displayed cats)</li>"
                "</ul>",
            ))

            # Test 8c: Lazy AnnData with max_lazy_categories=0 (metadata-only mode)
            print("  8c. Lazy AnnData (metadata-only mode)")

            # Use setting to disable category loading (instead of parameter)
            original_max_lazy_cats = ad.settings.repr_html_max_lazy_categories
            ad.settings.repr_html_max_lazy_categories = 0
            metadata_only_html = adata_lazy._repr_html_()
            ad.settings.repr_html_max_lazy_categories = original_max_lazy_cats

            sections.append((
                "8c. Lazy AnnData (Metadata-Only Mode)",
                metadata_only_html,
                "<code>ad.settings.repr_html_max_lazy_categories = 0</code><br>"
                "<p style='margin: 5px 0;'><b>Zero disk I/O mode:</b> Set <code>repr_html_max_lazy_categories=0</code> to "
                "skip all data loading from disk.</p>"
                "<p style='margin: 5px 0;'><b>What's shown:</b></p>"
                "<ul style='margin: 5px 0; padding-left: 20px;'>"
                "<li>All categorical columns show <code>(N categories)</code> - count from dtype metadata (already in memory)</li>"
                "<li>Category labels are NOT loaded from disk</li>"
                "<li>Colors from <code>uns</code> are NOT loaded (skips <code>.compute()</code> on dask arrays)</li>"
                "<li>Type column still shows <code>category (N, lazy)</code> with the count</li>"
                "</ul>"
                "<p style='margin: 5px 0;'><b>Use case:</b> Quick inspection of lazy AnnData when you want to avoid "
                "all data loading from disk, for the fastest possible repr.</p>",
            ))

        except (OSError, ImportError, TypeError) as e:
            print(f"    Warning: Failed to create lazy example: {e}")
        finally:
            if h5_file is not None:
                h5_file.close()
            Path(tmp_path).unlink()
    else:
        print("  8b. Lazy AnnData (skipped - xarray not installed)")

    # Test 9: Backed AnnData (H5AD file) - demonstrates on-disk safety
    print("  9. Backed AnnData (H5AD file)")
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp_path = tmp.name
    adata_backed = None
    try:
        adata_to_save = AnnData(
            sp.random(500, 200, density=0.1, format="csr", dtype=np.float32)
        )
        adata_to_save.obs["cluster"] = pd.Categorical(
            ["A", "B", "C"] * 166 + ["A", "B"]
        )
        adata_to_save.obs["n_counts"] = np.random.randint(1000, 10000, 500)
        adata_to_save.var["gene_name"] = [f"gene_{i}" for i in range(200)]
        adata_to_save.var["highly_variable"] = np.random.choice([True, False], 200)
        adata_to_save.obsm["X_pca"] = np.random.randn(500, 50).astype(np.float32)
        adata_to_save.write_h5ad(tmp_path)
        adata_backed = ad.read_h5ad(tmp_path, backed="r")
        sections.append((
            "9. Backed AnnData (H5AD File)",
            adata_backed._repr_html_(),
            "<strong>File-backed mode — X matrix stays on disk!</strong><br>"
            f"<code>{tmp_path}</code><br><br>"
            "What the repr reads (safe, small metadata):<br>"
            "<ul>"
            "<li><code>X.shape</code>, <code>X.dtype</code>, <code>X.nnz</code> — from HDF5 attributes</li>"
            "<li><code>obs</code>/<code>var</code> DataFrames — loaded in memory (typically small)</li>"
            "<li><code>obsm</code>/<code>varm</code> shapes — from HDF5 dataset attributes</li>"
            "</ul>"
            "What it does NOT load:<br>"
            "<ul>"
            "<li>The actual X matrix data (stays memory-mapped on disk)</li>"
            "<li>No slicing or indexing into X</li>"
            "</ul>",
        ))
    finally:
        if adata_backed is not None:
            adata_backed.file.close()
        Path(tmp_path).unlink()

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
    many_cat_values = [f"type_{i}" for i in range(30)] * (100 // 30) + [
        f"type_{i}" for i in range(100 % 30)
    ]
    adata_many_cats.obs["cell_type"] = pd.Categorical(many_cat_values)
    # Add colors for the categories
    adata_many_cats.uns["cell_type_colors"] = [
        "#e41a1c",
        "#377eb8",
        "#4daf4a",
        "#984ea3",
        "#ff7f00",
        "#ffff33",
        "#a65628",
        "#f781bf",
        "#999999",
        "#66c2a5",
        "#fc8d62",
        "#8da0cb",
        "#e78ac3",
        "#a6d854",
        "#ffd92f",
        "#e5c494",
        "#b3b3b3",
        "#1b9e77",
        "#d95f02",
        "#7570b3",
        "#e7298a",
        "#66a61e",
        "#e6ab02",
        "#a6761d",
        "#666666",
        "#8dd3c7",
        "#ffffb3",
        "#bebada",
        "#fb8072",
        "#80b1d3",
    ]
    # Also add a column with exactly 20 categories
    adata_many_cats.obs["batch"] = pd.Categorical([f"batch_{i}" for i in range(20)] * 5)
    adata_many_cats.uns["batch_colors"] = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
        "#aec7e8",
        "#ffbb78",
        "#98df8a",
        "#ff9896",
        "#c5b0d5",
        "#c49c94",
        "#f7b6d2",
        "#c7c7c7",
        "#dbdb8d",
        "#9edae5",
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
                html_parts.append(f"<strong>{len(runs)} runs</strong>")
            if params:
                param_str = ", ".join(f"{k}={v}" for k, v in list(params.items())[:3])
                if len(params) > 3:
                    param_str += "..."
                html_parts.append(f" · params: {param_str}")
            html_parts.append("</div>")

            return FormattedOutput(
                type_name="analysis history",
                html_content="".join(html_parts),
            )

    adata_uns = AnnData(np.zeros((10, 5)))
    # Simple types with previews
    adata_uns.uns["string_param"] = "A short string value"
    adata_uns.uns["long_string"] = (
        "This is a very long string that should be truncated in the preview because it exceeds the maximum length allowed for display in the meta column"
    )
    adata_uns.uns["int_param"] = 42
    adata_uns.uns["float_param"] = 3.14159265359
    adata_uns.uns["bool_param"] = True
    adata_uns.uns["none_param"] = None
    adata_uns.uns["small_list"] = [1, 2, 3]
    adata_uns.uns["small_dict"] = {"a": 1, "b": 2}
    adata_uns.uns["larger_dict"] = {
        "key1": "val1",
        "key2": "val2",
        "key3": "val3",
        "key4": "val4",
        "key5": "val5",
    }

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
    adata_uns.uns["string_hint"] = (
        "__anndata_repr__:otherpackage.config::{'setting': 'value'}"
    )

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
    adata_nojs.obsm["cell_measurements"] = pd.DataFrame(
        {
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
        },
        index=adata_nojs.obs_names,
    )
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
        adata_df.obsm["spatial_metrics"] = pd.DataFrame(
            {
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
            },
            index=adata_df.obs_names,
        )
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
    adata_long.obs[
        "this_is_an_extremely_long_column_name_that_should_test_the_max_width_setting"
    ] = np.random.randn(20)
    adata_long.obs["cell_type_annotation_from_automated_classifier_v2"] = (
        pd.Categorical(["A", "B"] * 10)
    )
    adata_long.obsm["X_pca_computed_with_highly_variable_genes_batch_corrected"] = (
        np.random.randn(20, 5).astype(np.float32)
    )
    adata_long.uns["preprocessing_parameters_for_normalization_and_scaling"] = {
        "method": "log1p",
        "scale": True,
    }
    adata_long.layers["raw_counts_before_any_preprocessing_steps"] = np.random.randn(
        20, 10
    ).astype(np.float32)
    sections.append((
        "16. Very Long Field Names",
        adata_long._repr_html_(),
        "Tests the dynamic field name column width calculation. The longest field name is "
        "'this_is_an_extremely_long_column_name_that_should_test_the_max_width_setting' (77 chars). "
        "The name column width should expand to fit longer names but be capped by "
        "<code>repr_html_max_field_width</code> (default: 400px). Names exceeding the max width "
        "show an ellipsis (...) via CSS text-overflow; hover over truncated names to see the "
        "full name in a tooltip.",
    ))

    # Test 17: README icon
    print("  17. README icon")
    adata_readme = AnnData(np.random.randn(50, 20).astype(np.float32))
    adata_readme.obs["cluster"] = pd.Categorical(["A", "B", "C", "D", "E"] * 10)
    adata_readme.uns["cluster_colors"] = [
        "#e41a1c",
        "#377eb8",
        "#4daf4a",
        "#984ea3",
        "#ff7f00",
    ]
    adata_readme.obsm["X_pca"] = np.random.randn(50, 10).astype(np.float32)
    adata_readme.uns["README"] = """# Human Lung Adenocarcinoma - Patient LU-A047

Single-cell RNA sequencing of a *primary* lung adenocarcinoma tumor sample. This dataset was generated as part of a study investigating **tumor heterogeneity** and ***immune cell infiltration*** patterns in early-stage lung cancer.

## Sample Information
- **Tissue**: Primary lung tumor, right upper lobe
- **Diagnosis**: Adenocarcinoma, stage IIA (T2aN0M0)
- **Collection date**: 2024-03-15
- **Patient ID**: LU-A047 (IRB #2023-0892)

## Wet Lab Protocol

### Tissue Dissociation
1. Fresh tissue dissociation using Miltenyi Tumor Dissociation Kit
2. Dead cell removal via MACS Dead Cell Removal Kit
3. Red blood cell lysis (ACK buffer, 2 min)

### Quality Control
4. Filtered through 40 µm cell strainer
5. Viability assessment: **92%** (Trypan Blue)

#### Technical Notes
Cell viability was *above threshold* for sequencing. Data stored in `adata.obs['viability']`.

## 10x Genomics Processing
- **Chemistry**: Chromium Next GEM 3' v3.1
- **Target cells**: 10,000
- **Cells recovered**: 8,247
- **Sequencing**: NovaSeq 6000, 28×90 bp

## Quality Metrics
| Metric | Value |
|--------|-------|
| Median genes/cell | 2,847 |
| Median UMIs/cell | 8,392 |
| Sequencing saturation | 78.2% |

## Usage Example
```python
import scanpy as sc

adata = sc.read_h5ad("LU-A047.h5ad")
sc.pp.filter_cells(adata, min_genes=200)
sc.pl.umap(adata, color=["cell_type", "cluster"])
```

## References
- [10x Genomics Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression)
- [Scanpy documentation](https://scanpy.readthedocs.io/)
- GEO accession: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE_example_data

## Contact
For questions about this dataset: `genome-lab@example-hospital.org`

> **Note**: All patient identifiers have been de-identified per HIPAA guidelines.
"""
    sections.append((
        "17. README Icon",
        adata_readme._repr_html_(),
        "When <code>uns['README']</code> contains a string, a small ⓘ icon appears in the header. "
        "Click the icon to open a modal with the README content rendered as markdown. "
        "The modal supports: headers (h1-h4), bold/italic, code blocks, inline code, "
        "ordered lists (numbered), unordered lists (bulleted), links, tables, and blockquotes. "
        "Press Escape or click outside to close.",
    ))

    # Test 18: README icon in No-JS mode
    print("  18. README icon in No-JS mode")
    adata_readme_nojs = AnnData(np.random.randn(20, 10).astype(np.float32))
    adata_readme_nojs.obs["batch"] = pd.Categorical(["batch1", "batch2"] * 10)
    adata_readme_nojs.uns["README"] = """# Dataset Information

This dataset contains processed single-cell data.

## Key Features
- 20 cells, 10 genes
- 2 batches

For more details, see the full documentation.
"""
    nojs_readme_html = strip_script_tags(adata_readme_nojs._repr_html_())
    sections.append((
        "18. README Icon in No-JS Mode",
        nojs_readme_html,
        "In no-JS mode, the README icon still appears but clicking it won't open a modal. "
        "Instead, **hover over the icon** to see the README content as a tooltip (browser's "
        "native title attribute). The tooltip shows the first 500 characters of the README.",
    ))

    # Test 19: MuData (multimodal data)
    # This demonstrates how MuData can reuse anndata's repr by:
    # 1. Registering a SectionFormatter for the .mod attribute (done at import time above)
    # 2. Calling generate_repr_html() directly on the MuData object
    if HAS_MUDATA:
        print("  19. MuData (multimodal data)")
        from anndata._repr.html import generate_repr_html

        mdata = create_test_mudata()
        if mdata is not None:
            sections.append((
                "19. MuData (Multimodal Data)",
                generate_repr_html(mdata),
                "Demonstrates how MuData can reuse anndata's HTML repr machinery by simply "
                "registering a <code>SectionFormatter</code> for the <code>.mod</code> attribute. "
                "The <code>mod</code> section shows each modality as an expandable nested AnnData "
                "(click the arrow to expand). All standard sections (obs, var, obsm, varm, uns, etc.) "
                "work automatically. This example has 3 modalities: RNA (100×50), ATAC (100×30), "
                "and Protein (80×20).",
            ))
    else:
        print("  19. MuData (skipped - mudata not installed)")

    # Test 20: SpatialData (custom _repr_html_ using anndata's building blocks)
    # This demonstrates how packages with completely different structures can build
    # their own _repr_html_ while reusing anndata's CSS, JavaScript, and formatters.
    if HAS_SPATIALDATA_EXAMPLE:
        print("  20. SpatialData (custom _repr_html_ using anndata's building blocks)")
        sdata = create_test_spatialdata()
        sections.append((
            "20. SpatialData (Custom _repr_html_)",
            sdata._repr_html_(),
            "Demonstrates how packages like <a href='https://spatialdata.scverse.org/' "
            "target='_blank'>SpatialData</a> can build custom <code>_repr_html_</code> "
            "using anndata's building blocks: "
            "<ul>"
            "<li><code>get_css()</code> / <code>get_javascript()</code> - reuse styling and interactivity</li>"
            "<li><code>render_section()</code> - create collapsible sections (images, labels, points, shapes, tables)</li>"
            "<li><code>render_formatted_entry()</code> with <code>meta_content</code> - table rows with meta column</li>"
            "<li><code>generate_repr_html()</code> - embed nested AnnData (see 'tables' section)</li>"
            "<li><code>FormatterRegistry</code> - custom 'transforms' section added via SectionFormatter</li>"
            "</ul>"
            "Note the meta column shows dimension info like <code>[c, y, x]</code>. "
            "The nested AnnData objects in <code>tables</code> are fully interactive (click Expand). "
            "Hover over coordinate system names to see associated elements.",
        ))
    else:
        print("  20. SpatialData (skipped - example failed to load)")

    # Test 21: Raw section with detailed info
    print("  21. Raw section (unprocessed data)")
    # Create an AnnData that simulates a typical workflow:
    # 1. Start with more genes (raw)
    # 2. Filter to fewer genes (current)
    n_obs, n_vars_raw = 100, 2000
    n_vars_filtered = 500
    adata_raw = AnnData(
        # Current filtered data
        np.random.randn(n_obs, n_vars_filtered).astype(np.float32),
        obs=pd.DataFrame({
            "cell_type": pd.Categorical(
                ["T cell", "B cell", "NK cell"] * 33 + ["T cell"]
            ),
            "n_counts": np.random.randint(1000, 10000, n_obs),
        }),
        var=pd.DataFrame({
            "gene_name": [f"HVG_{i}" for i in range(n_vars_filtered)],
            "highly_variable": [True] * n_vars_filtered,
            "mean_expression": np.random.randn(n_vars_filtered).astype(np.float32),
        }),
    )
    # Set raw to have more genes (simulating pre-filtering state)
    raw_X = np.random.randn(n_obs, n_vars_raw).astype(np.float32)
    raw_var = pd.DataFrame(
        {
            "gene_name": [f"gene_{i}" for i in range(n_vars_raw)],
            "highly_variable": [i < n_vars_filtered for i in range(n_vars_raw)],
        },
        index=[f"gene_{i}" for i in range(n_vars_raw)],
    )
    adata_raw.raw = AnnData(raw_X, var=raw_var)
    # Add varm to raw
    adata_raw.raw.varm["PCs"] = np.random.randn(n_vars_raw, 50).astype(np.float32)
    sections.append((
        "21a. Raw Section - Dense Matrix with var and varm",
        adata_raw._repr_html_(),
        "Shows the <code>.raw</code> attribute with dense matrix, var columns, and varm. "
        "Click the raw row to expand and see the full Raw repr with X, var, and varm sections.",
    ))

    # Test 21b: Raw with sparse matrix
    print("  21b. Raw section (sparse matrix)")
    adata_sparse_raw = AnnData(
        sp.random(n_obs, n_vars_filtered, density=0.1, format="csr", dtype=np.float32),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars_filtered)]),
    )
    sparse_raw_X = sp.random(
        n_obs, n_vars_raw, density=0.05, format="csr", dtype=np.float32
    )
    adata_sparse_raw.raw = AnnData(sparse_raw_X, var=raw_var)
    sections.append((
        "21b. Raw Section - Sparse Matrix",
        adata_sparse_raw._repr_html_(),
        "Raw with sparse CSR matrix. The X section should show sparsity info.",
    ))

    # Test 21c: Raw with no varm (minimal raw)
    print("  21c. Raw section (minimal - no varm)")
    adata_minimal_raw = AnnData(
        np.random.randn(50, 100).astype(np.float32),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(100)]),
    )
    minimal_raw_var = pd.DataFrame(
        {"gene_symbol": [f"GENE{i}" for i in range(200)]},
        index=[f"gene_{i}" for i in range(200)],
    )
    adata_minimal_raw.raw = AnnData(
        np.random.randn(50, 200).astype(np.float32),
        var=minimal_raw_var,
    )
    sections.append((
        "21c. Raw Section - Minimal (no varm)",
        adata_minimal_raw._repr_html_(),
        "Raw with only X and var (no varm). Shows graceful handling of optional attributes.",
    ))

    # Test 21d: Raw with empty var columns
    print("  21d. Raw section (empty var columns)")
    adata_empty_var_raw = AnnData(
        np.random.randn(30, 50).astype(np.float32),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(50)]),
    )
    empty_raw_var = pd.DataFrame(index=[f"gene_{i}" for i in range(80)])  # No columns
    adata_empty_var_raw.raw = AnnData(
        np.random.randn(30, 80).astype(np.float32),
        var=empty_raw_var,
    )
    sections.append((
        "21d. Raw Section - Empty var columns",
        adata_empty_var_raw._repr_html_(),
        "Raw where var has no columns (only index). The meta info should not show 'var: 0 cols'.",
    ))

    # Test 22: Unknown sections and error handling
    print("  22. Unknown sections and error handling")

    class ExtendedAnnData(AnnData):
        """AnnData subclass with custom mapping attributes."""

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._custom_mappings = {}

        @property
        def custom_data(self):
            """Custom mapping-like attribute."""
            return self._custom_mappings

        @property
        def failing_data(self):
            """Property that raises an error when accessed."""
            msg = "This property intentionally fails for testing"
            raise RuntimeError(msg)

    adata_extended = ExtendedAnnData(
        np.random.randn(50, 100).astype(np.float32),
        obs=pd.DataFrame({"cluster": pd.Categorical(["A", "B"] * 25)}),
    )
    # Add some custom data
    adata_extended._custom_mappings = {
        "embedding": np.random.randn(50, 2),
        "config": {"param1": 1, "param2": "value"},
    }
    sections.append((
        "22a. Unknown Sections",
        adata_extended._repr_html_(),
        "Demonstrates two important features for scientific accuracy:<br>"
        "<ol>"
        "<li><strong>Unknown sections:</strong> The <code>custom_data</code> mapping attribute "
        "appears in an 'other' section at the bottom, ensuring no data is silently hidden.</li>"
        "<li><strong>Error handling:</strong> The <code>failing_data</code> property raises an "
        "error when accessed. Instead of silently hiding it, the repr shows it as 'inaccessible' "
        "in the 'other' section.</li>"
        "</ol>"
        "This ensures researchers always know what data exists, even if it can't be rendered.",
    ))

    # Test 22b: Failing section rendering (using real rendering mechanism)
    # Use unittest.mock to safely patch properties
    print("  22b. Failing section rendering (real errors)")
    from unittest.mock import PropertyMock, patch

    # Create a mapping-like object that raises an error when accessed
    class FailingMapping:
        """A mapping that raises an error when iterated."""

        def __init__(self, error_msg: str):
            self._error_msg = error_msg

        def keys(self):
            raise RuntimeError(self._error_msg)

        def __len__(self):
            return 1  # Report as non-empty so it tries to render

        def __iter__(self):
            raise RuntimeError(self._error_msg)

    # Create a real AnnData with data
    rng = np.random.default_rng(42)
    adata_failing = ad.AnnData(
        X=rng.random((100, 50)),
        obs=pd.DataFrame(
            {"cell_type": ["A", "B", "C"] * 33 + ["A"]},
            index=[f"cell_{i}" for i in range(100)],
        ),
        var=pd.DataFrame(
            {"gene_name": [f"gene_{i}" for i in range(50)]},
            index=[f"gene_{i}" for i in range(50)],
        ),
        obsm={"X_pca": rng.random((100, 10)), "X_umap": rng.random((100, 2))},
        varm={"loadings": rng.random((50, 10))},
        layers={"counts": rng.integers(0, 100, (100, 50))},
        obsp={"distances": sp.csr_matrix(rng.random((100, 100)))},
        uns={"method": "test", "params": {"k": 10}},
    )

    # Create failing mappings
    failing_varm = FailingMapping(
        "Failed to decompress data block (corrupted zarr chunk)"
    )
    failing_layers = FailingMapping(
        "IOError: [Errno 5] Input/output error reading '/data/counts.h5'"
    )

    # Use unittest.mock.patch to safely patch properties (auto-restores on exit)
    with (
        patch.object(
            type(adata_failing),
            "varm",
            new_callable=PropertyMock,
            return_value=failing_varm,
        ),
        patch.object(
            type(adata_failing),
            "layers",
            new_callable=PropertyMock,
            return_value=failing_layers,
        ),
    ):
        # Use the real repr mechanism - this will trigger errors on varm and layers
        failing_html = adata_failing._repr_html_()

    sections.append((
        "22b. Failing Section Rendering (Real Errors)",
        failing_html,
        "Demonstrates the <strong>actual error handling mechanism</strong> in the repr. "
        "This test uses <code>unittest.mock.patch</code> to make <code>varm</code> and "
        "<code>layers</code> properties raise real exceptions when accessed.<br>"
        "<ul>"
        "<li><strong>Normal sections:</strong> X, obs, var, uns, obsm, obsp, varp render correctly</li>"
        "<li><strong>varm (error):</strong> Real RuntimeError from corrupted data simulation</li>"
        "<li><strong>layers (error):</strong> Real IOError from I/O failure simulation</li>"
        "</ul>"
        "This tests the actual <code>_render_section</code> try/except error handling "
        "using the real <code>generate_repr_html</code> pipeline.",
    ))

    # Test 23: Serialization warnings and edge cases
    print("  23. Serialization warnings and edge cases")

    class CustomObject:
        def __repr__(self):
            return "CustomObject()"

    adata_serial = ad.AnnData(X=np.eye(5))

    # Cases that SHOULD warn (and fail to serialize)
    adata_serial.obs["list_values"] = [["a", "b"], ["c"], ["d"], ["e"], ["f"]]
    adata_serial.obs["dict_values"] = [{"k": 1}, {"k": 2}, {"k": 3}, {"k": 4}, {"k": 5}]
    adata_serial.obs["custom_obj"] = [CustomObject() for _ in range(5)]
    adata_serial.obs["datetime_col"] = pd.to_datetime([
        "2024-01-01",
        "2024-01-02",
        "2024-01-03",
        "2024-01-04",
        "2024-01-05",
    ])
    adata_serial.obs["timedelta_col"] = pd.to_timedelta([
        "1 days",
        "2 days",
        "3 days",
        "4 days",
        "5 days",
    ])
    adata_serial.obs["path/slash"] = ["a", "b", "c", "d", "e"]  # Slash in name
    adata_serial.obs[("tuple", "name")] = [1, 2, 3, 4, 5]  # Non-string name

    # Cases that should NOT warn (and serialize fine)
    adata_serial.var["gène_名前"] = ["a", "b", "c", "d", "e"]  # Non-ASCII is OK
    adata_serial.var["normal_col"] = [1.0, 2.0, 3.0, 4.0, 5.0]  # Normal floats
    adata_serial.var["string_col"] = ["a", "b", "c", "d", "e"]  # Strings OK

    # Mapping sections (layers, obsm, etc.) - same rules apply
    adata_serial.layers[("tuple", "key")] = np.eye(5)  # Non-string - fails NOW
    adata_serial.obsm["path/embed"] = np.random.randn(5, 2)  # Slash - deprecated

    sections.append((
        "23. Serialization Warnings",
        adata_serial._repr_html_(),
        "<strong>Fails NOW (red):</strong>"
        "<ul>"
        "<li><code>obs.list_values</code> - Contains list</li>"
        "<li><code>obs.dict_values</code> - Contains dict</li>"
        "<li><code>obs.custom_obj</code> - Contains CustomObject</li>"
        "<li><code>obs.datetime_col</code> - datetime64 not serializable</li>"
        "<li><code>obs.timedelta_col</code> - timedelta64 not serializable</li>"
        "<li><code>obs.('tuple', 'name')</code> - Non-string column name</li>"
        "<li><code>layers.('tuple', 'key')</code> - Non-string key</li>"
        "</ul>"
        "<strong>Will fail in future (yellow):</strong>"
        "<ul>"
        "<li><code>obs.path/slash</code> - Slash in name</li>"
        "<li><code>obsm.path/embed</code> - Slash in key</li>"
        "</ul>"
        "<strong>Serializes fine (no warning):</strong>"
        "<ul>"
        "<li><code>var.gène_名前</code> - Non-ASCII is valid UTF-8</li>"
        "<li><code>var.normal_col</code> - Normal float values</li>"
        "<li><code>var.string_col</code> - String values</li>"
        "</ul>",
    ))

    # Generate HTML file
    output_path = Path(__file__).parent / "repr_html_visual_test.html"
    html_content = create_html_page(sections)
    output_path.write_text(html_content)

    print(f"\nVisual test file generated: {output_path}")
    print("Open this file in a browser to inspect the HTML representation.")


if __name__ == "__main__":
    main()
