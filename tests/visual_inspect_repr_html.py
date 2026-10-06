"""
Visual inspection script for the AnnData HTML representation (``_repr_html_``).

Generates one self-contained HTML page with every scenario below, for manual
review in a browser (light/dark, with/without JS/CSS, embedded themes, ...).

Usage::

    python tests/visual_inspect_repr_html.py              # all cases
    python tests/visual_inspect_repr_html.py --list       # list cases, don't render
    python tests/visual_inspect_repr_html.py --only lazy --only 5.   # filter
    python tests/visual_inspect_repr_html.py -o /tmp/repr.html --strict
    python tests/visual_inspect_repr_html.py --only theme --browser-check  # headless Chrome

Then open ``tests/repr_html_visual_test.html`` in a browser.

Each case is a small function registered with ``@case(category, title, ...)``.
It returns the HTML to show and declares what the reviewer should look for
(``expect``). Cases are numbered ``<category>.<n>`` in registration order;
anchors use the explicit ``slug`` so links stay stable when cases are added.
A case that raises shows its traceback in the page instead of aborting the run;
missing optional dependencies show as "skipped". ``[[slug]]`` in ``expect``
or ``notes`` renders as a link to that case.

Theme cases embed the repr in iframes that mimic each host page; OS dark mode
is simulated per iframe, and each pane self-reports PASS/FAIL. The summary at
the top collects them, and ``--browser-check`` prints them from headless Chrome.

Categories:

1. Core structure: full/empty/minimal objects, ``X`` vs ``layers[None]``, raw,
   nested AnnData, README.
2. Data types: dense/sparse, pandas extension dtypes, categoricals & colors,
   dask, awkward, array-API devices, ``uns`` value types.
3. Storage states: views, backed h5ad, ``read_lazy`` (h5ad and zarr).
4. Scale & truncation: folding, many columns/categories/entries, wide DataFrames,
   huge shapes, long names, README size limit.
5. Environments & theming: no JS, no CSS, HTML repr disabled, Jupyter, VS Code,
   Furo, pydata-sphinx-theme, theme-less pages in OS light/dark mode.
6. Robustness & security: XSS, broken properties, failing formatters,
   serialization warnings, the "evil" AnnData.
7. Extensibility / ecosystem: uns type hints, TreeData/MuData/SpatialData,
   ecosystem TypeFormatter, AnnData subclasses, custom array types.

See also:
- src/anndata/_repr/registry.py: TypeFormatter and SectionFormatter APIs
- Reviewer's Guide gist for architecture overview
"""

# ruff: noqa: EM101
# EM101: Exception string literals are used intentionally in evil test objects
# RUF001/RUF003: Unicode lookalike characters are intentional - testing confusable chars

from __future__ import annotations

import argparse
import html as html_mod
import importlib.util
import platform
import re
import shutil
import sys
import tempfile
import time
import traceback
import warnings
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path
from typing import TYPE_CHECKING, Protocol

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
    FormattedEntry,
    FormattedOutput,
    FormatterContext,
    SectionFormatter,
    TypeFormatter,
    escape_html,
    extract_uns_type_hint,
    formatter_registry,
    register_formatter,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator, Sequence
    from typing import Any, Literal

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
    import networkx as nx  # type: ignore[import-untyped]

    HAS_NETWORKX = True
except (ImportError, AttributeError):
    # AttributeError can occur on Python 3.14+ with incompatible networkx versions
    HAS_NETWORKX = False

try:
    from treedata import TreeData  # type: ignore[import-not-found]

    HAS_TREEDATA = HAS_NETWORKX
except (ImportError, AttributeError):
    HAS_TREEDATA = False

if HAS_NETWORKX:

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
            (
                f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg" '
                f'style="background:#fafafa;border-radius:4px;border:1px solid #e0e0e0;">'
            )
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

    # Register TreeData section formatters (what treedata would do at import time)
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
                    css_class="anndata-dtype--tree",
                    tooltip=f"Phylogenetic tree with {n_nodes} total nodes",
                    expanded_html=svg_html,
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
                    css_class="anndata-dtype--tree",
                    tooltip=f"Phylogenetic tree with {n_nodes} total nodes",
                    expanded_html=svg_html,
                )
                entries.append(FormattedEntry(key=key, output=output))
            return entries

    @register_formatter
    class TreeMetadataSectionFormatter(SectionFormatter):
        """Section formatter for TreeData metadata as a compact inline row.

        This demonstrates a fully custom section representation using
        ``render_html()`` instead of the standard ``get_entries()`` path.
        When ``render_html()`` is defined, it takes precedence and the
        returned HTML is inserted directly — no ``<details>`` wrapping,
        no entry grid. This is useful for compact metadata that doesn't
        fit the "list of entries" pattern.

        Renders like the X entry — a single non-foldable line showing
        key=value pairs for label, alignment, and allow_overlap.
        All values are escaped via ``escape_html(repr(val))``.
        """

        @property
        def section_name(self) -> str:
            return "tree_metadata"

        @property
        def display_name(self) -> str:
            return "tree"

        @property
        def after_section(self) -> str:
            return "X"

        @property
        def doc_url(self) -> str:
            return TREEDATA_DOCS

        @property
        def tooltip(self) -> str:
            return "Tree configuration parameters (TreeData)"

        def should_show(self, obj) -> bool:
            return hasattr(obj, "_tree_label")

        def get_entries(self, obj, context: FormatterContext) -> list[FormattedEntry]:
            """Fallback if render_html fails (e.g., missing optional dependency)."""
            entries = []
            for attr, label in [
                ("_tree_label", "label"),
                ("_alignment", "alignment"),
                ("_allow_overlap", "allow_overlap"),
            ]:
                val = getattr(obj, attr, None)
                if val is not None:
                    output = FormattedOutput(
                        type_name=type(val).__name__,
                        preview=repr(val),
                    )
                    entries.append(FormattedEntry(key=label, output=output))
            return entries

        def render_html(self, obj, context: FormatterContext) -> str:
            """Render as a compact line instead of a foldable section."""
            pairs = []
            for attr, label in [
                ("_tree_label", "label"),
                ("_alignment", "alignment"),
                ("_allow_overlap", "allow_overlap"),
            ]:
                val = getattr(obj, attr, None)
                if val is not None:
                    pairs.append(
                        f'<span style="color:var(--anndata-text-secondary,#6c757d);">{label}=</span>'
                        f"{escape_html(repr(val))}"
                    )
            summary = " &nbsp; ".join(pairs)
            return (
                '<div class="anndata-x__entry">'
                f"<span>tree</span>"
                f"<span>{summary}</span>"
                "</div>"
            )

    class TreeDataStandIn(AnnData):
        """Minimal stand-in exposing TreeData's attributes (``obst``, ``vart``, tree metadata).

        Used when the real ``treedata`` package is missing or incompatible with the
        installed anndata, so the SectionFormatter demo still renders.
        """

        def __init__(
            self,
            *args,
            obst: dict[str, nx.DiGraph],
            vart: dict[str, nx.DiGraph],
            label: str,
            alignment: str,
            allow_overlap: bool,
            **kwargs,
        ) -> None:
            super().__init__(*args, **kwargs)
            self._obst = obst
            self._vart = vart
            self._tree_label = label
            self._alignment = alignment
            self._allow_overlap = allow_overlap

        @property
        def obst(self) -> dict[str, nx.DiGraph]:
            return self._obst

        @property
        def vart(self) -> dict[str, nx.DiGraph]:
            return self._vart


# Check for MuData
try:
    from mudata import MuData

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
                    css_class="anndata-dtype--anndata",
                    tooltip=f"Modality: {mod_name}",
                    expanded_html=nested_html if can_expand else None,
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
        FormatterRegistry,
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
            parts.append('<div class="anndata-repr__sections">')
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
            parts = ['<div class="anndata-header">']
            parts.append('<span class="anndata-header__type">SpatialData</span>')

            # Zarr badge using render_badge() helper
            if self.path:
                parts.append(
                    render_badge(
                        "Zarr", "anndata-badge--backed", "Backed by Zarr storage"
                    )
                )
                parts.append(
                    f'<span class="anndata-file-path" style="font-family:ui-monospace,monospace;'
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
            parts = ['<div class="anndata-index-preview" style="padding:2px 8px;">']
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
                meta = f'<span class="anndata-meta-info">[{dims_str}]</span>'

                # Create a FormattedEntry with FormattedOutput
                entry = FormattedEntry(
                    key=name,
                    output=FormattedOutput(
                        type_name=f"DataArray {info['shape']} {info['dtype']}",
                        css_class="anndata-dtype--ndarray",
                        preview_html=meta,  # Content in preview column (rightmost)
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
                meta = f'<span class="anndata-meta-info">[{dims_str}]</span>'

                entry = FormattedEntry(
                    key=name,
                    output=FormattedOutput(
                        type_name=f"Labels {info['shape']} {info['dtype']}",
                        css_class="anndata-dtype--ndarray",
                        preview_html=meta,
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
                meta = f'<span class="anndata-meta-info">{info["n_dims"]}D coordinates</span>'

                entry = FormattedEntry(
                    key=name,
                    output=FormattedOutput(
                        type_name=f"dask.DataFrame ({format_number(info['n_points'])} × {info['n_dims']})",
                        css_class="anndata-dtype--dataframe",
                        preview_html=meta,
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
                meta = f'<span class="anndata-meta-info">{info["geometry_type"]}</span>'

                entry = FormattedEntry(
                    key=name,
                    output=FormattedOutput(
                        type_name=f"GeoDataFrame ({format_number(info['n_shapes'])} shapes)",
                        css_class="anndata-dtype--dataframe",
                        preview_html=meta,
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

                # FormattedOutput with expanded_html makes it collapsible
                entry = FormattedEntry(
                    key=name,
                    output=FormattedOutput(
                        type_name=f"AnnData ({adata.n_obs} × {adata.n_vars})",
                        css_class="anndata-dtype--anndata",
                        expanded_html=nested_html,  # Makes the nested content collapsible
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

        def can_format(self, obj, context):
            return isinstance(obj, dict) and "shape" in obj and "dtype" in obj

        def format(self, obj, context: FormatterContext) -> FormattedOutput:
            return FormattedOutput(
                type_name=f"DataTree {obj['shape']} {obj['dtype']}",
                css_class="anndata-dtype--ndarray",
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


# =============================================================================
# Case registry
# =============================================================================

CaseOutput = str | tuple[str, str]
"""What a case returns: the HTML to show, optionally with a runtime note."""


class SkipCase(Exception):
    """Raise inside a case to skip it; the message is shown in the page."""


@dataclass(frozen=True)
class Category:
    key: str
    title: str
    blurb: str


CATEGORIES: tuple[Category, ...] = (
    Category(
        "core",
        "Core structure",
        "Section layout and ordering for typical objects: X, obs/var, *m, *p, layers, uns, raw.",
    ),
    Category(
        "dtypes",
        "Data types",
        "How individual values are typed and previewed in each section.",
    ),
    Category(
        "storage",
        "Storage states (views, backed, lazy)",
        "Badges, file paths, and that nothing is loaded or computed just to render.",
    ),
    Category(
        "scale",
        "Scale & truncation",
        "Folding, truncation indicators, number formatting and column widths for big objects.",
    ),
    Category(
        "env",
        "Environments & theming",
        "Graceful degradation and following the host page's light/dark theme.",
    ),
    Category(
        "robust",
        "Robustness & security",
        "Escaping, broken objects, failing formatters: the repr must never crash or execute input.",
    ),
    Category(
        "ext",
        "Extensibility / ecosystem",
        "Public extension points used by downstream packages.",
    ),
)
_CATEGORY_INDEX = {c.key: i for i, c in enumerate(CATEGORIES, start=1)}


@dataclass
class Case:
    slug: str
    category: str
    title: str
    func: Callable[[], CaseOutput]
    expect: tuple[str, ...]
    notes: str = ""
    requires: tuple[str, ...] = ()
    tags: tuple[str, ...] = ()
    number: str = ""


CASES: list[Case] = []


def case(
    category: str,
    title: str,
    *,
    slug: str,
    expect: Sequence[str],
    notes: str = "",
    requires: Sequence[str] = (),
    tags: Sequence[str] = (),
) -> Callable[[Callable[[], CaseOutput]], Callable[[], CaseOutput]]:
    """Register a visual test case.

    Parameters
    ----------
    category
        Key from :data:`CATEGORIES`.
    title
        Short human-readable title.
    slug
        Stable anchor id (kebab-case). Do not change once published.
    expect
        Checklist of what the reviewer should see (HTML allowed, ``[[slug]]`` links).
    notes
        Optional longer background (HTML allowed), shown collapsed.
    requires
        Importable module names; the case is skipped if any is missing.
    tags
        Coverage tags (sections, dtypes, states, ...) for the coverage index.
    """
    if category not in _CATEGORY_INDEX:
        msg = f"Unknown category {category!r}"
        raise ValueError(msg)
    if not re.fullmatch(r"[a-z0-9]+(-[a-z0-9]+)*", slug):
        msg = f"Slug must be kebab-case: {slug!r}"
        raise ValueError(msg)

    def decorator(func: Callable[[], CaseOutput]) -> Callable[[], CaseOutput]:
        if any(c.slug == slug for c in CASES):
            msg = f"Duplicate case slug {slug!r}"
            raise ValueError(msg)
        CASES.append(
            Case(
                slug=slug,
                category=category,
                title=title,
                func=func,
                expect=tuple(expect),
                notes=notes,
                requires=tuple(requires),
                tags=tuple(tags),
            )
        )
        return func

    return decorator


def numbered_cases() -> list[Case]:
    """All cases ordered by category, numbered ``<category>.<n>``."""
    ordered = sorted(CASES, key=lambda c: _CATEGORY_INDEX[c.category])
    counters: dict[str, int] = {}
    for c in ordered:
        counters[c.category] = counters.get(c.category, 0) + 1
        c.number = f"{_CATEGORY_INDEX[c.category]}.{counters[c.category]}"
    return ordered


@dataclass
class CaseResult:
    case: Case
    status: Literal["ok", "skipped", "failed"]
    html: str = ""
    note: str = ""
    reason: str = ""
    seconds: float = 0.0
    warnings: list[str] = field(default_factory=list)


def run_case(c: Case) -> CaseResult:
    """Run one case, capturing skips, exceptions and emitted warnings."""
    missing = [m for m in c.requires if importlib.util.find_spec(m) is None]
    if missing:
        return CaseResult(
            c, "skipped", reason=f"skipped: {', '.join(missing)} not installed"
        )
    start = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        warnings.filterwarnings("ignore", message="Transforming to str index")
        try:
            out = c.func()
        except SkipCase as e:
            result = CaseResult(c, "skipped", reason=f"skipped: {e}")
        except Exception:  # noqa: BLE001
            result = CaseResult(c, "failed", reason=traceback.format_exc())
        else:
            html_out, note = out if isinstance(out, tuple) else (out, "")
            result = CaseResult(c, "ok", html=html_out, note=note)
    result.seconds = time.perf_counter() - start
    counts: dict[str, int] = {}
    for w in caught:
        msg = f"{w.category.__name__}: {str(w.message)[:300]}"
        counts[msg] = counts.get(msg, 0) + 1
    result.warnings = [f"{m} (×{n})" if n > 1 else m for m, n in counts.items()]
    return result


# =============================================================================
# Helpers
# =============================================================================


class _HasReprHtml(Protocol):
    def _repr_html_(self) -> str | None: ...


def render(obj: _HasReprHtml) -> str:
    """Return ``obj._repr_html_()``, failing loudly if it returned ``None``."""
    out = obj._repr_html_()
    if out is None:
        msg = f"{type(obj).__name__}._repr_html_() returned None"
        raise RuntimeError(msg)
    return out


def strip_script_tags(html: str) -> str:
    """Remove <script>...</script> tags from HTML to simulate no-JS environment."""
    return re.sub(r"<script>.*?</script>", "", html, flags=re.DOTALL)


def strip_style_and_script_tags(html: str) -> str:
    """Remove <style> and <script> tags to simulate GitHub/untrusted notebook rendering."""
    html = re.sub(r"<style[^>]*>.*?</style>", "", html, flags=re.DOTALL)
    return strip_script_tags(html)


def iframe(doc: str, *, title: str, style: str = "") -> str:
    """Embed a full HTML document in an auto-sized, CSS-isolated ``<iframe srcdoc>``.

    The harness page script resizes ``iframe.vt-frame`` to its content.
    """
    return (
        f'<iframe class="vt-frame" title="{html_mod.escape(title)}" '
        f'srcdoc="{html_mod.escape(doc, quote=True)}" style="{style}"></iframe>'
    )


@contextmanager
def temporarily_registered(
    formatter: TypeFormatter | SectionFormatter,
) -> Iterator[None]:
    """Register a formatter instance for the duration of one case.

    Real packages register at import time with ``@register_formatter``; scoping
    it here keeps one case's formatter from leaking into later cases.
    """
    register_formatter(formatter)
    try:
        yield
    finally:
        if isinstance(formatter, TypeFormatter):
            formatter_registry.unregister_type_formatter(formatter)
        else:
            unregister = getattr(
                formatter_registry, "unregister_section_formatter", None
            )
            for name in formatter.section_names:
                if unregister is not None:
                    unregister(name)
                else:  # older registry without the public method
                    formatter_registry._section_formatters.pop(name, None)


@contextmanager
def tmp_path(suffix: str) -> Iterator[Path]:
    """Yield a path in a fresh temporary directory, removed afterwards."""
    d = Path(tempfile.mkdtemp(prefix="anndata-repr-visual-"))
    try:
        yield d / f"data{suffix}"
    finally:
        shutil.rmtree(d, ignore_errors=True)


def palette(n: int) -> list[str]:
    """``n`` distinct hex colors (tab20 + set1-ish, cycled)."""
    base = [
        "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628",
        "#f781bf", "#999999", "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854",
        "#ffd92f", "#e5c494", "#b3b3b3", "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
        "#66a61e", "#e6ab02", "#a6761d", "#666666", "#8dd3c7", "#ffffb3", "#bebada",
        "#fb8072", "#80b1d3",
    ]  # fmt: skip
    return [base[i % len(base)] for i in range(n)]


# =============================================================================
# Shared test data factories
# =============================================================================


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
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mdata = MuData({"rna": rna, "atac": atac, "prot": prot})

    # Add shared annotations
    mdata.uns["experiment"] = "multiome_sample_001"
    mdata.uns["processing_date"] = "2024-03-15"

    return mdata


def create_test_treedata() -> tuple[AnnData, str]:
    """Create a TreeData object with observation and variable trees.

    Falls back to :class:`TreeDataStandIn` if ``treedata`` is missing or
    incompatible; the second return value then explains why.
    """
    if not HAS_NETWORKX:
        raise SkipCase("networkx not installed")

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

    kwargs: dict[str, Any] = dict(
        X=np.random.randn(n_obs, n_vars).astype(np.float32),
        obs=pd.DataFrame(
            {"cell_type": pd.Categorical(["T cell", "B cell"] * (n_obs // 2))},
            index=obs_names,
        ),
        var=pd.DataFrame({"gene_name": var_names}, index=var_names),
        obst={"phylogeny": obs_tree},
        vart={"gene_ontology": var_tree},
        label="phylogeny",
        alignment="leaves",
        allow_overlap=False,
    )
    note = ""
    tdata: AnnData
    if HAS_TREEDATA:
        try:
            tdata = TreeData(**kwargs)
        except Exception as e:  # noqa: BLE001
            note = (
                f"Real <code>treedata</code> failed with this anndata "
                f"(<code>{escape_html(f'{type(e).__name__}: {e}')}</code>); "
                "rendering <code>TreeDataStandIn</code> instead."
            )
            tdata = TreeDataStandIn(**kwargs)
    else:
        note = "treedata not installed; rendering <code>TreeDataStandIn</code>."
        tdata = TreeDataStandIn(**kwargs)

    # Add standard annotations
    tdata.uns["cell_type_colors"] = ["#e41a1c", "#377eb8"]
    tdata.obsm["X_pca"] = np.random.randn(n_obs, 10).astype(np.float32)
    tdata.layers["raw"] = np.random.randn(n_obs, n_vars).astype(np.float32)

    return tdata, note


def create_test_anndata() -> AnnData:
    """Create a comprehensive test AnnData with all features.

    This showcases common patterns from real single-cell analysis workflows:
    - Sparse X matrix with typical density
    - Categorical columns with color annotations
    - Numeric QC metrics
    - String columns (some will trigger serialization warnings)
    - Datetime columns (will trigger serialization warnings)
    - Boolean columns
    - Cluster assignments (louvain, leiden)
    - Dimensionality reductions (PCA, UMAP, t-SNE)
    - Neighbor graphs
    - Layers (raw counts, normalized)
    - Raw attribute (unprocessed data)
    - Various uns types (dicts, arrays, nested AnnData)
    """
    n_obs, n_vars = 100, 50

    # Main AnnData with sparse X
    # obs: 5 columns (stays expanded below fold_threshold)
    # var: more columns to demonstrate folding and various types
    adata = AnnData(
        sp.random(n_obs, n_vars, density=0.1, format="csr", dtype=np.float32),
        obs=pd.DataFrame({
            # Categorical with colors (5 categories)
            "cell_type": pd.Categorical(
                ["T cell", "B cell", "NK cell", "Monocyte", "DC"] * (n_obs // 5)
            ),
            # Categorical with colors (8 clusters)
            "louvain": pd.Categorical([
                f"cluster_{i}" for i in (np.random.randint(0, 8, n_obs))
            ]),
            # Numeric QC metric
            "n_counts": np.random.randint(1000, 50000, n_obs),
            # Float QC metric
            "percent_mito": np.random.uniform(0, 15, n_obs).astype(np.float32),
            # Boolean column
            "is_doublet": np.random.choice([True, False], n_obs, p=[0.1, 0.9]),
        }),
        var=pd.DataFrame({
            # Basic gene info
            "gene_symbol": [f"GN{i}" for i in range(n_vars)],
            "highly_variable": np.random.choice([True, False], n_vars, p=[0.2, 0.8]),
            "means": np.random.exponential(1, n_vars).astype(np.float32),
            "dispersions": np.random.exponential(0.5, n_vars).astype(np.float32),
            # Categorical column
            "chromosome": pd.Categorical([f"chr{i % 22 + 1}" for i in range(n_vars)]),
            # String column (will trigger categorical conversion warning)
            "gene_biotype": ["protein_coding"] * (n_vars - 5) + ["lncRNA"] * 5,
            # Datetime column (will trigger serialization warning)
            "annotation_date": pd.to_datetime(["2024-01-15"] * n_vars),
            # All unique strings (no warning - too many unique values)
            "ensembl_id": [f"ENSG{i:011d}" for i in range(n_vars)],
        }),
    )

    # === Color annotations ===
    adata.uns["cell_type_colors"] = [
        "#FF6B6B",
        "#4ECDC4",
        "#45B7D1",
        "#96CEB4",
        "#FFEAA7",
    ]
    adata.uns["louvain_colors"] = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
    ]

    # === Uns: Analysis results (typical scanpy output) ===
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": {"n_neighbors": 15, "method": "umap", "metric": "euclidean"},
    }
    adata.uns["pca"] = {
        "variance": np.random.exponential(10, 50).astype(np.float32),
        "variance_ratio": np.sort(np.random.uniform(0, 0.1, 50))[::-1].astype(
            np.float32
        ),
    }
    adata.uns["umap"] = {"params": {"min_dist": 0.5, "spread": 1.0}}
    adata.uns["louvain"] = {"params": {"resolution": 1.0, "random_state": 0}}

    # === Uns: Simple values ===
    adata.uns["experiment_id"] = "EXP_2024_001"
    adata.uns["n_highly_variable"] = int(adata.var["highly_variable"].sum())
    adata.uns["total_counts"] = float(adata.obs["n_counts"].sum())
    adata.uns["processing_steps"] = [
        "filtering",
        "normalization",
        "hvg",
        "pca",
        "neighbors",
        "umap",
        "clustering",
    ]

    # === Uns: Nested AnnData ===
    inner_adata = AnnData(np.zeros((10, 5)))
    inner_adata.obs["inner_cluster"] = pd.Categorical(["A", "B"] * 5)
    inner_adata.var["gene"] = [f"gene_{i}" for i in range(5)]
    adata.uns["subset_adata"] = inner_adata

    # === Uns: Unserializable type (will warn) ===
    class CustomAnalysisResult:
        def __repr__(self):
            return "CustomAnalysisResult(n_clusters=8)"

    adata.uns["custom_result"] = CustomAnalysisResult()

    # === Obsm: Embeddings and metadata ===
    adata.obsm["X_pca"] = np.random.randn(n_obs, 50).astype(np.float32)
    adata.obsm["X_umap"] = np.random.randn(n_obs, 2).astype(np.float32)
    adata.obsm["X_tsne"] = np.random.randn(n_obs, 2).astype(np.float32)
    # DataFrame in obsm (spatial coordinates)
    adata.obsm["spatial"] = pd.DataFrame(
        {
            "x": np.random.uniform(0, 1000, n_obs),
            "y": np.random.uniform(0, 1000, n_obs),
            "z": np.random.uniform(0, 100, n_obs),
            "area": np.random.uniform(50, 500, n_obs),
            "perimeter": np.random.uniform(20, 100, n_obs),
        },
        index=adata.obs_names,
    )

    # === Varm: Gene loadings ===
    adata.varm["PCs"] = np.random.randn(n_vars, 50).astype(np.float32)

    # === Layers: Different normalizations ===
    adata.layers["counts"] = sp.random(
        n_obs, n_vars, density=0.1, format="csr", dtype=np.float32
    )
    adata.layers["normalized"] = np.random.randn(n_obs, n_vars).astype(np.float32)
    adata.layers["log1p"] = np.log1p(np.abs(np.random.randn(n_obs, n_vars))).astype(
        np.float32
    )

    # === Obsp/Varp: Graphs ===
    adata.obsp["distances"] = sp.random(
        n_obs, n_obs, density=0.05, format="csr", dtype=np.float32
    )
    adata.obsp["connectivities"] = sp.random(
        n_obs, n_obs, density=0.05, format="csr", dtype=np.float32
    )
    adata.varp["gene_correlation"] = sp.random(
        n_vars, n_vars, density=0.1, format="csr", dtype=np.float32
    )

    # === Raw: Unprocessed data (common in scanpy workflows) ===
    raw_adata = AnnData(
        sp.random(n_obs, n_vars + 20, density=0.1, format="csr", dtype=np.float32),
        var=pd.DataFrame({
            "gene_name": [f"Gene_{i}" for i in range(n_vars + 20)],
            "n_cells": np.random.randint(1, n_obs, n_vars + 20),
        }),
    )
    adata.raw = raw_adata

    return adata


def create_theme_demo_anndata() -> AnnData:
    """Compact AnnData exercising most colored elements (for theme panes)."""
    rng = np.random.default_rng(0)
    X = sp.random(40, 12, density=0.2, format="csr", dtype=np.float32, rng=rng)
    adata = AnnData(
        X,
        obs=pd.DataFrame({
            "cell_type": pd.Categorical(["T", "B", "NK", "Mono"] * 10),
            "n_counts": rng.integers(100, 1000, 40),
            "is_doublet": rng.choice([True, False], 40),
        }),
        var=pd.DataFrame({
            "gene_symbol": [f"G{i}" for i in range(12)],
            "date": pd.to_datetime(["2024-01-01"] * 12),  # serialization warning
        }),
    )
    adata.uns["cell_type_colors"] = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"]
    adata.uns["params"] = {"k": 15, "metric": "cosine"}
    adata.uns["nested"] = AnnData(np.zeros((3, 2)))
    adata.uns["README"] = "# Theme demo\n\nREADME icon and modal colors."
    adata.obsm["X_umap"] = rng.standard_normal((40, 2)).astype(np.float32)
    adata.layers["counts"] = X.copy()
    return adata


def write_lazy_demo_h5ad(path: Path) -> None:
    """Write the file used by the lazy h5ad cases."""
    adata = AnnData(sp.random(1000, 500, density=0.1, format="csr", dtype=np.float32))

    # --- Categorical columns ---
    # 1. Small categorical WITH colors (should show categories + color dots)
    adata.obs["cell_type"] = pd.Categorical(
        np.random.choice(["T cell", "B cell", "Monocyte", "NK cell"], 1000)
    )
    adata.uns["cell_type_colors"] = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"]
    # 2. Small categorical WITHOUT colors (should show categories only)
    adata.obs["cluster"] = pd.Categorical(
        np.random.choice(["C0", "C1", "C2", "C3", "C4"], 1000)
    )
    # 3. Medium categorical (50 cats) - will show truncation with max_lazy_categories=30
    medium_categories = [f"sample_{i}" for i in range(50)]
    adata.obs["sample_id"] = pd.Categorical(
        np.random.choice(medium_categories, 1000),
        categories=medium_categories,  # Ensure all 50 categories exist
    )
    # --- Non-categorical columns (all should show "(lazy)") ---
    adata.obs["n_genes"] = np.random.randint(500, 5000, 1000)
    adata.obs["total_counts"] = np.random.randint(1000, 50000, 1000)
    # --- var columns ---
    adata.var["gene_symbol"] = [f"GENE{i}" for i in range(500)]
    adata.var["highly_variable"] = np.random.choice([True, False], 500)
    adata.var["mean_expression"] = np.random.uniform(0, 10, 500)
    # --- obsm/varm ---
    adata.obsm["X_pca"] = np.random.randn(1000, 50).astype(np.float32)
    adata.obsm["X_umap"] = np.random.randn(1000, 2).astype(np.float32)
    adata.varm["PCs"] = np.random.randn(500, 50).astype(np.float32)
    # --- uns with array (to show dask array WITH size in uns) ---
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
    }
    adata.uns["pca_variance"] = np.random.rand(50).astype(np.float32)
    adata.write_h5ad(path)


def write_lazy_demo_zarr(path: Path) -> None:
    """Write the file used by the lazy zarr case."""
    adata = AnnData(sp.random(800, 400, density=0.1, format="csr", dtype=np.float32))
    adata.obs["tissue"] = pd.Categorical(
        np.random.choice(["Brain", "Heart", "Liver", "Lung", "Kidney"], 800)
    )
    adata.uns["tissue_colors"] = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
    adata.obs["donor"] = pd.Categorical(
        np.random.choice([f"D{i}" for i in range(10)], 800)
    )
    adata.obs["n_counts"] = np.random.randint(1000, 50000, 800)
    adata.var["gene_name"] = [f"GENE{i}" for i in range(400)]
    adata.obsm["X_umap"] = np.random.randn(800, 2).astype(np.float32)
    adata.write_zarr(path)


def write_backed_demo_h5ad(path: Path, *, dense: bool = False) -> None:
    """Write the file used by the backed cases."""
    X = (
        np.random.randn(500, 200).astype(np.float32)
        if dense
        else sp.random(500, 200, density=0.1, format="csr", dtype=np.float32)
    )
    adata = AnnData(X)
    adata.obs["cluster"] = pd.Categorical(["A", "B", "C"] * 166 + ["A", "B"])
    adata.uns["cluster_colors"] = ["#e41a1c", "#377eb8", "#4daf4a"]
    adata.obs["n_counts"] = np.random.randint(1000, 10000, 500)
    adata.var["gene_name"] = [f"gene_{i}" for i in range(200)]
    adata.var["highly_variable"] = np.random.choice([True, False], 200)
    adata.obsm["X_pca"] = np.random.randn(500, 50).astype(np.float32)
    adata.layers["counts"] = sp.random(
        500, 200, density=0.1, format="csr", dtype=np.float32
    )
    adata.write_h5ad(path)


README_LUNG = """# Human Lung Adenocarcinoma - Patient LU-A047

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


# =============================================================================
# 1. Core structure
# =============================================================================


@case(
    "core",
    "Full AnnData (all features)",
    slug="full-anndata",
    tags=(
        "X",
        "obs",
        "var",
        "obsm",
        "varm",
        "obsp",
        "varp",
        "layers",
        "uns",
        "raw",
        "sparse",
        "colors",
        "nested-anndata",
    ),
    expect=(
        "Section order: X first, then obs, var, obsm, varm, obsp, varp, layers, raw, uns.",
        "X appears exactly once; <code>layers</code> lists counts/normalized/log1p only (no <code>None</code> entry).",
        "Categorical columns show color dots from <code>uns['*_colors']</code>.",
        "<code>var.annotation_date</code> (datetime) and <code>uns.custom_result</code> are flagged as not serializable.",
        "Each section header has a <b>?</b> link to the docs; hovering a section name shows a tooltip.",
    ),
    notes="Baseline reference for a typical annotated dataset; several other cases reuse this object.",
)
def _full_anndata() -> CaseOutput:
    return render(create_test_anndata())


@case(
    "core",
    "Empty AnnData",
    slug="empty-anndata",
    tags=("empty",),
    expect=(
        "Header shows 0 × 0 and no crash.",
        "Empty sections are hidden or show 'No entries'; no stray X row with a bogus shape.",
    ),
)
def _empty_anndata() -> CaseOutput:
    return render(AnnData())


@case(
    "core",
    "Minimal AnnData (just X)",
    slug="minimal-x-only",
    tags=("X", "dense"),
    expect=(
        "Only the X row carries information; obs/var have default integer names and no columns.",
        "No <code>layers</code> section entry for X (X lives in <code>layers[None]</code> internally).",
    ),
)
def _minimal() -> CaseOutput:
    return render(AnnData(np.zeros((10, 5))))


@case(
    "core",
    "X stored as layers[None], with named layers",
    slug="x-layers-none-with-layers",
    tags=("X", "layers"),
    expect=(
        "X row shown once, at the top, as a dense float32 ndarray.",
        "<code>layers</code> lists exactly <code>counts</code> and <code>log1p</code> (2 items, not 3).",
        "No entry named <code>None</code> anywhere; the section count says 2.",
    ),
    notes="Since upstream #1707 <code>adata.X is adata.layers[None]</code>; the repr must hide that alias.",
)
def _x_layers_none_with_layers() -> CaseOutput:
    X = np.random.randn(30, 8).astype(np.float32)
    adata = AnnData(X)
    adata.layers["counts"] = sp.random(30, 8, density=0.3, format="csr")
    adata.layers["log1p"] = np.log1p(np.abs(X))
    assert adata.layers[None] is adata.X
    return render(adata)


@case(
    "core",
    "X is None, only named layers",
    slug="x-none-with-layers",
    tags=("X", "layers", "empty"),
    expect=(
        "X row says <code>None</code> (or is absent); shape in header is still 30 × 8 (from obs/var).",
        "<code>layers</code> shows <code>spliced</code> and <code>unspliced</code> only.",
    ),
)
def _x_none_with_layers() -> CaseOutput:
    adata = AnnData(
        obs=pd.DataFrame(index=[f"c{i}" for i in range(30)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(8)]),
        layers={
            "spliced": sp.random(30, 8, density=0.3, format="csr"),
            "unspliced": sp.random(30, 8, density=0.3, format="csr"),
        },
    )
    return render(adata)


@case(
    "core",
    "Dense matrix with categories",
    slug="dense-with-categories",
    tags=("X", "dense", "categorical", "colors"),
    expect=(
        "X shows <code>ndarray</code> (not CSR/CSC).",
        "<code>obs.cluster</code> shows 5 categories with color dots.",
    ),
)
def _dense_categories() -> CaseOutput:
    adata = AnnData(np.random.randn(50, 30).astype(np.float32))
    adata.obs["cluster"] = pd.Categorical(["A", "B", "C", "D", "E"] * 10)
    adata.uns["cluster_colors"] = palette(5)
    return render(adata)


@case(
    "core",
    "Raw: dense matrix with var and varm",
    slug="raw-dense",
    tags=("raw", "dense"),
    expect=(
        "Current shape 100 × 500; the <code>raw</code> row reports 100 × 2,000.",
        "Expanding raw shows X (dense), var (2 columns) and varm (PCs).",
    ),
    notes="Typical workflow: filtered AnnData with <code>.raw</code> preserving all genes.",
)
def _raw_dense() -> CaseOutput:
    n_obs, n_vars_raw, n_vars_filtered = 100, 2000, 500
    adata = AnnData(
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
    raw_var = pd.DataFrame(
        {
            "gene_name": [f"gene_{i}" for i in range(n_vars_raw)],
            "highly_variable": [i < n_vars_filtered for i in range(n_vars_raw)],
        },
        index=[f"gene_{i}" for i in range(n_vars_raw)],
    )
    raw = AnnData(np.random.randn(n_obs, n_vars_raw).astype(np.float32), var=raw_var)
    raw.varm["PCs"] = np.random.randn(n_vars_raw, 50).astype(np.float32)
    adata.raw = raw
    return render(adata)


@case(
    "core",
    "Raw: sparse matrix",
    slug="raw-sparse",
    tags=("raw", "sparse"),
    expect=(
        "Both X (10% density) and raw.X (5% density) are CSR; the expanded raw shows CSR type and sparsity.",
        "Compare with [[raw-dense]].",
    ),
)
def _raw_sparse() -> CaseOutput:
    n_obs, n_vars_raw, n_vars_filtered = 100, 2000, 500
    adata = AnnData(
        sp.random(n_obs, n_vars_filtered, density=0.1, format="csr", dtype=np.float32),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars_filtered)]),
    )
    raw_var = pd.DataFrame(
        {"gene_name": [f"gene_{i}" for i in range(n_vars_raw)]},
        index=[f"gene_{i}" for i in range(n_vars_raw)],
    )
    adata.raw = AnnData(
        sp.random(n_obs, n_vars_raw, density=0.05, format="csr", dtype=np.float32),
        var=raw_var,
    )
    return render(adata)


@case(
    "core",
    "Raw: minimal (no varm)",
    slug="raw-minimal",
    tags=("raw",),
    expect=("Expanded raw shows only X and var; no empty varm section.",),
)
def _raw_minimal() -> CaseOutput:
    adata = AnnData(
        np.random.randn(50, 100).astype(np.float32),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(100)]),
    )
    adata.raw = AnnData(
        np.random.randn(50, 200).astype(np.float32),
        var=pd.DataFrame(
            {"gene_symbol": [f"GENE{i}" for i in range(200)]},
            index=[f"gene_{i}" for i in range(200)],
        ),
    )
    return render(adata)


@case(
    "core",
    "Raw: var without columns",
    slug="raw-empty-var",
    tags=("raw", "empty"),
    expect=(
        "Raw row shows the shape (30 × 80) but no 'var: 0 cols' text and no empty var section.",
    ),
)
def _raw_empty_var() -> CaseOutput:
    adata = AnnData(
        np.random.randn(30, 50).astype(np.float32),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(50)]),
    )
    adata.raw = AnnData(
        np.random.randn(30, 80).astype(np.float32),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(80)]),
    )
    return render(adata)


@case(
    "core",
    "Deeply nested AnnData (max depth)",
    slug="nested-max-depth",
    tags=("uns", "nested-anndata"),
    expect=(
        "outer → level1 → level2 expand as nested reprs.",
        "level3 is beyond <code>repr_html_max_depth</code> (default 3) and shows as a non-expandable entry.",
    ),
)
def _nested_depth() -> CaseOutput:
    inner3 = AnnData(np.zeros((3, 2)))
    inner2 = AnnData(np.zeros((5, 3)))
    inner2.uns["level3"] = inner3
    inner1 = AnnData(np.zeros((10, 5)))
    inner1.uns["level2"] = inner2
    outer = AnnData(np.zeros((20, 10)))
    outer.uns["level1"] = inner1
    return render(outer)


@case(
    "core",
    "README icon",
    slug="readme-icon",
    tags=("uns", "readme"),
    expect=(
        "A small ⓘ icon in the header; clicking opens a modal with the raw markdown as plain text.",
        "Escape or clicking outside closes the modal.",
    ),
)
def _readme_icon() -> CaseOutput:
    adata = AnnData(np.random.randn(50, 20).astype(np.float32))
    adata.obs["cluster"] = pd.Categorical(["A", "B", "C", "D", "E"] * 10)
    adata.uns["cluster_colors"] = palette(5)
    adata.obsm["X_pca"] = np.random.randn(50, 10).astype(np.float32)
    adata.uns["README"] = README_LUNG
    return render(adata)


# =============================================================================
# 2. Data types
# =============================================================================


@case(
    "dtypes",
    "Sparse formats and numeric dtypes",
    slug="sparse-formats",
    tags=("X", "layers", "obsm", "obsp", "sparse", "dense"),
    expect=(
        "X: <code>csr_matrix</code>; layers distinguish <code>csc_matrix</code>, <code>csr_array</code>, <code>csc_array</code>.",
        "Sparse entries show density / '% sparse' and stored count; dense layers show their dtype (int8, uint16, bool, float16, complex64).",
        "obsm sparse embedding and obsp <code>csr_array</code> render like layers.",
    ),
)
def _sparse_formats() -> CaseOutput:
    n_obs, n_vars = 60, 40
    adata = AnnData(sp.random(n_obs, n_vars, density=0.05, format="csr"))
    adata.layers["csc_matrix"] = sp.random(n_obs, n_vars, density=0.2, format="csc")
    adata.layers["csr_array"] = sp.csr_array(sp.random(n_obs, n_vars, density=0.5))
    adata.layers["csc_array_int"] = sp.csc_array(
        sp.random(n_obs, n_vars, density=0.01, format="csc", dtype=np.float32)
    ).astype(np.int32)
    adata.layers["int8"] = np.zeros((n_obs, n_vars), dtype=np.int8)
    adata.layers["uint16"] = np.ones((n_obs, n_vars), dtype=np.uint16)
    adata.layers["bool"] = np.zeros((n_obs, n_vars), dtype=bool)
    adata.layers["float16"] = np.zeros((n_obs, n_vars), dtype=np.float16)
    adata.layers["complex64"] = np.zeros((n_obs, n_vars), dtype=np.complex64)
    adata.obsm["sparse_embedding"] = sp.random(n_obs, 100, density=0.01, format="csr")
    adata.obsp["knn"] = sp.csr_array(sp.random(n_obs, n_obs, density=0.1))
    return render(adata)


@case(
    "dtypes",
    "pandas extension and nullable dtypes in obs",
    slug="pandas-dtypes",
    tags=("obs", "nullable", "extension-dtype", "string", "categorical"),
    expect=(
        "Nullable <code>Int64</code>/<code>Float64</code>/<code>boolean</code> and <code>string</code> columns show their pandas dtype names.",
        "Ordered categorical is distinguishable (or at least renders), unused categories are listed.",
        "datetime with tz, period, interval and object-mixed columns are flagged as not serializable where writing would fail.",
        "No crash on an all-NaN categorical.",
    ),
)
def _pandas_dtypes() -> CaseOutput:
    n = 12
    adata = AnnData(np.zeros((n, 3)))
    adata.obs["Int64"] = pd.array([1, None, *range(n - 2)], dtype="Int64")
    adata.obs["Float64"] = pd.array([0.5, None, *np.arange(n - 2.0)], dtype="Float64")
    adata.obs["boolean"] = pd.array([True, None] + [False] * (n - 2), dtype="boolean")
    adata.obs["string"] = pd.array([f"s{i}" for i in range(n)], dtype="string")
    if importlib.util.find_spec("pyarrow") is not None:
        adata.obs["string_pyarrow"] = pd.array(
            [f"s{i % 3}" for i in range(n)], dtype="string[pyarrow]"
        )
    adata.obs["ordered_cat"] = pd.Categorical(
        ["low", "mid", "high"] * (n // 3),
        categories=["low", "mid", "high"],
        ordered=True,
    )
    adata.obs["unused_cats"] = pd.Categorical(["a"] * n, categories=["a", "b", "c"])
    adata.obs["int_cats"] = pd.Categorical([1, 2, 3] * (n // 3))
    adata.obs["bool_cats"] = pd.Categorical([True, False] * (n // 2))
    adata.obs["all_nan_cat"] = pd.Categorical([np.nan] * n, categories=["x"])
    adata.obs["datetime_tz"] = pd.date_range(
        "2024-01-01", periods=n, tz="Europe/Berlin"
    )
    adata.obs["period"] = pd.period_range("2024-01", periods=n, freq="M")
    adata.obs["interval"] = pd.interval_range(0, n)
    adata.obs["mixed_object"] = pd.Series([1, "a", 2.0, None] * (n // 4), dtype=object)
    adata.obs["uint8"] = np.arange(n, dtype=np.uint8)
    return render(adata)


@case(
    "dtypes",
    "Categorical colors: formats and edge cases",
    slug="categorical-colors",
    tags=("obs", "var", "categorical", "colors"),
    expect=(
        "Hex (#rgb, #rrggbb, #rrggbbaa), named CSS colors, and numpy-array colors all render as dots.",
        "Colors for a var categorical (<code>var.chrom</code>) are looked up in <code>uns['chrom_colors']</code> too.",
        "Category with NaN values: NaN is not a category and must not shift colors.",
    ),
    notes="Malformed / hostile color arrays are covered in [[evil-anndata]].",
)
def _categorical_colors() -> CaseOutput:
    adata = AnnData(np.zeros((12, 6)))
    adata.obs["hex_short"] = pd.Categorical(["a", "b", "c"] * 4)
    adata.uns["hex_short_colors"] = ["#f00", "#0f0", "#00f"]
    adata.obs["hex_alpha"] = pd.Categorical(["a", "b"] * 6)
    adata.uns["hex_alpha_colors"] = ["#ff000080", "#0000ff80"]
    adata.obs["named"] = pd.Categorical(["x", "y", "z", "w"] * 3)
    adata.uns["named_colors"] = ["tomato", "steelblue", "gold", "black"]
    adata.obs["np_array_colors"] = pd.Categorical(["p", "q"] * 6)
    adata.uns["np_array_colors_colors"] = np.array(["#1b9e77", "#d95f02"])
    adata.obs["with_nan"] = pd.Categorical(["u", np.nan, "v"] * 4)
    adata.uns["with_nan_colors"] = ["#e41a1c", "#377eb8"]
    adata.var["chrom"] = pd.Categorical(["chr1", "chr2", "chrX"] * 2)
    adata.uns["chrom_colors"] = ["#66c2a5", "#fc8d62", "#8da0cb"]
    return render(adata)


@case(
    "dtypes",
    "Dask arrays (no compute)",
    slug="dask-arrays",
    requires=("dask",),
    tags=("X", "layers", "obsm", "varm", "dask", "sparse"),
    expect=(
        "X, layers, obsm and varm show <code>dask.array</code> with shape, dtype and chunks.",
        "Rendering is fast: nothing is computed (no <code>.compute()</code>).",
        "<code>layers['sparse_chunks']</code> has CSR chunks: ideally hinted as sparse.",
    ),
    notes=(
        "Regular in-memory AnnData whose arrays are dask arrays. obs/var are plain pandas; "
        "compare with the file-backed lazy cases in category 3."
    ),
)
def _dask() -> CaseOutput:
    X_dask = da.random.random((1000, 500), chunks=(100, 100))
    adata = AnnData(X_dask)
    adata.obs["cluster"] = pd.Categorical(["A", "B", "C"] * 333 + ["A"])
    adata.var["gene_name"] = [f"gene_{i}" for i in range(500)]
    adata.layers["counts"] = da.random.randint(0, 100, (1000, 500), chunks=(100, 100))
    adata.layers["sparse_chunks"] = da.from_array(
        sp.random(1000, 500, density=0.01, format="csr"),
        chunks=(250, 500),
        asarray=False,
    )
    adata.obsm["X_pca"] = da.random.random((1000, 50), chunks=(100, 50))
    adata.varm["loadings"] = da.random.random((500, 50), chunks=(100, 50))
    return render(adata)


@case(
    "dtypes",
    "Awkward arrays in obsm",
    slug="awkward-arrays",
    requires=("awkward",),
    tags=("obsm", "awkward"),
    expect=(
        "Ragged list and record arrays show <code>awkward.Array</code> with record count (ideally the type, e.g. <code>var * int64</code>).",
        "No crash; awkward styling distinct from numpy arrays.",
    ),
)
def _awkward() -> CaseOutput:
    import awkward as ak

    n = 20
    adata = AnnData(np.zeros((n, 4)))
    adata.obsm["ragged"] = ak.Array([list(range(i % 4)) for i in range(n)])
    adata.obsm["records"] = ak.Array([
        {"x": float(i), "tags": ["a"] * (i % 3)} for i in range(n)
    ])
    return render(adata)


@case(
    "dtypes",
    "Array-API arrays with device info",
    slug="array-api-devices",
    tags=("obsm", "uns", "array-api"),
    expect=(
        "Device appears inline as <code>dtype · device</code> in the type column (no hover needed).",
        "<code>X_jax_gpu</code> cuda:0, <code>X_jax_tpu</code> tpu:0, <code>X_jax_cpu</code> cpu, <code>X_cupy_gpu</code> GPU:0 (GPU-green), <code>uns['gpu_embedding']</code> cuda:1.",
    ),
    notes="Uses mock objects satisfying the <code>SupportsArrayApi</code> protocol; no GPU or JAX needed.",
)
def _array_api() -> CaseOutput:
    def make_mock(module, *, shape, dtype, device="cpu"):
        """Create a mock array satisfying the SupportsArrayApi protocol."""
        ns_module = type("Namespace", (), {"__name__": module.split(".")[0]})()
        cls = type(
            "MockArrayAPI",
            (),
            {
                "shape": shape,
                "dtype": dtype,
                "ndim": len(shape),
                "size": int(np.prod(shape)),
                "device": device,
                "__array_namespace__": lambda self, **kw: ns_module,
                "to_device": lambda self, dev, /, **kw: self,
                "__dlpack__": lambda self, **kw: None,
                "__dlpack_device__": lambda self: (1, 0),
                "__getitem__": lambda self, k, /: self,
            },
        )
        cls.__module__ = module
        return cls()

    n_obs, n_vars = 100, 50
    adata = AnnData(
        np.random.randn(n_obs, n_vars).astype(np.float32),
        obs=pd.DataFrame(
            {"cell_type": pd.Categorical(["T cell", "B cell"] * (n_obs // 2))},
            index=[f"cell_{i}" for i in range(n_obs)],
        ),
        var=pd.DataFrame(
            {"gene_name": [f"gene_{i}" for i in range(n_vars)]},
            index=[f"gene_{i}" for i in range(n_vars)],
        ),
    )
    adata.obsm["X_jax_gpu"] = make_mock(
        "jax.numpy", shape=(n_obs, 30), dtype=np.dtype("float32"), device="cuda:0"
    )
    adata.obsm["X_jax_tpu"] = make_mock(
        "jax.numpy", shape=(n_obs, 10), dtype=np.dtype("float16"), device="tpu:0"
    )
    adata.obsm["X_jax_cpu"] = make_mock(
        "jax.numpy", shape=(n_obs, 50), dtype=np.dtype("float64"), device="cpu"
    )

    class _MockGPUDevice:
        id = 0

    adata.obsm["X_cupy_gpu"] = make_mock(
        "cupy._core.core",
        shape=(n_obs, 20),
        dtype=np.dtype("float32"),
        device=_MockGPUDevice(),
    )
    adata.uns["gpu_embedding"] = make_mock(
        "jax.numpy", shape=(20, 5), dtype=np.dtype("float32"), device="cuda:1"
    )
    return render(adata)


@case(
    "dtypes",
    "uns value types",
    slug="uns-value-types",
    tags=("uns", "dense", "sparse", "string"),
    expect=(
        "Every entry has a sensible type label and, where cheap, a value preview.",
        "numpy scalar vs Python scalar, 0-d array, bytes, string array, structured array, DataFrame, sparse matrix, nested list, tuple, set, empty containers.",
        "Non-serializable values (tuple? set, complex) flagged consistently with what <code>write_h5ad</code> does.",
    ),
)
def _uns_value_types() -> CaseOutput:
    adata = AnnData(np.zeros((5, 3)))
    adata.uns["np_float32"] = np.float32(1.5)
    adata.uns["np_int64"] = np.int64(7)
    adata.uns["np_bool"] = np.True_
    adata.uns["py_complex"] = 1 + 2j
    adata.uns["zero_d_array"] = np.array(3)
    adata.uns["bytes"] = b"\x00\x01binary"
    adata.uns["str_array"] = np.array(["a", "bb", "ccc"])
    adata.uns["object_array"] = np.array(["a", 1, None], dtype=object)
    adata.uns["structured"] = np.array(
        [(1, 2.0), (3, 4.0)], dtype=[("a", "i4"), ("b", "f8")]
    )
    adata.uns["dataframe"] = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})
    adata.uns["sparse"] = sp.csr_matrix(np.eye(4))
    adata.uns["nested_list"] = [[1, 2], [3, [4, 5]]]
    adata.uns["tuple"] = (1, "two", 3.0)
    adata.uns["set"] = {1, 2, 3}
    adata.uns["empty_dict"] = {}
    adata.uns["empty_list"] = []
    adata.uns["empty_string"] = ""
    adata.uns["multiline_string"] = "line one\nline two\n\ttabbed"
    adata.uns["nan"] = float("nan")
    adata.uns["inf"] = float("-inf")
    return render(adata)


# =============================================================================
# 3. Storage states (views, backed, lazy)
# =============================================================================


@case(
    "storage",
    "AnnData view (subset)",
    slug="view",
    tags=("view", "sparse", "raw"),
    expect=(
        "Header shows a 'View' badge and the subset shape 20 × 10.",
        "All sections of [[full-anndata]] are present with subset shapes; raw keeps its own var count.",
    ),
)
def _view() -> CaseOutput:
    return render(create_test_anndata()[0:20, 0:10])


@case(
    "storage",
    "View of a view (boolean mask)",
    slug="view-of-view",
    tags=("view", "categorical"),
    expect=(
        "Still a single 'View' badge; shape reflects both subsets.",
        "Categoricals keep all categories (unused ones included) and colors stay aligned.",
    ),
)
def _view_of_view() -> CaseOutput:
    adata = create_test_anndata()
    first = adata[adata.obs["cell_type"].isin(["T cell", "B cell"])]
    return render(first[:, 5:25])


@case(
    "storage",
    "Backed AnnData (h5ad, sparse X)",
    slug="backed-h5ad",
    requires=("h5py",),
    tags=("backed", "sparse", "X", "layers"),
    expect=(
        "Header shows a backed/H5AD badge and the file path.",
        "X is a backed sparse dataset (shape, dtype, nnz from HDF5 metadata; data not loaded).",
        "obs/var are in memory (regular columns, categories with colors).",
    ),
    notes="Backed mode loads obs/var fully, while lazy mode ([[lazy-h5ad]]) keeps them as dask-backed xarray.",
)
def _backed_h5ad() -> CaseOutput:
    with tmp_path(".h5ad") as path:
        write_backed_demo_h5ad(path)
        adata = ad.read_h5ad(path, backed="r")
        try:
            return render(adata)
        finally:
            adata.file.close()


@case(
    "storage",
    "Backed AnnData (h5ad, dense X)",
    slug="backed-h5ad-dense",
    requires=("h5py",),
    tags=("backed", "dense", "X"),
    expect=("X shows an h5py Dataset (dense) with shape/dtype, not loaded.",),
)
def _backed_h5ad_dense() -> CaseOutput:
    with tmp_path(".h5ad") as path:
        write_backed_demo_h5ad(path, dense=True)
        adata = ad.read_h5ad(path, backed="r")
        try:
            return render(adata)
        finally:
            adata.file.close()


@case(
    "storage",
    "View of a backed AnnData",
    slug="backed-view",
    requires=("h5py",),
    tags=("backed", "view"),
    expect=("Both the backed badge/path and the View badge appear; shape is 50 × 20.",),
)
def _backed_view() -> CaseOutput:
    with tmp_path(".h5ad") as path:
        write_backed_demo_h5ad(path)
        adata = ad.read_h5ad(path, backed="r")
        try:
            return render(adata[:50, :20])
        finally:
            adata.file.close()


@case(
    "storage",
    "Lazy AnnData (read_lazy, h5ad)",
    slug="lazy-h5ad",
    requires=("xarray", "h5py"),
    tags=("lazy", "categorical", "colors", "lazy-categorical"),
    expect=(
        "Header: <b>Lazy (H5AD)</b> badge and file path.",
        "<code>cell_type</code>: 4 labels + colors; <code>cluster</code>: 5 labels; <code>sample_id</code>: first 30 of 50 (<code>max_lazy_categories=30</code>).",
        "Non-categorical columns show '(lazy)'; arrays show shape/dtype only.",
    ),
    notes=(
        "Only category labels (and colors from uns) are read from disk; codes, numeric values and "
        "categories beyond the limit are not. Compare with [[lazy-metadata-only]]."
    ),
)
def _lazy_h5ad() -> CaseOutput:
    import h5py

    with tmp_path(".h5ad") as path:
        write_lazy_demo_h5ad(path)
        with (
            h5py.File(path, "r") as f,
            ad.settings.override(repr_html_max_lazy_categories=30),
        ):
            return render(read_lazy(f))


@case(
    "storage",
    "Lazy AnnData, metadata-only (max_lazy_categories=0)",
    slug="lazy-metadata-only",
    requires=("xarray", "h5py"),
    tags=("lazy", "lazy-categorical"),
    expect=(
        "Same object as [[lazy-h5ad]], but zero disk I/O for the repr.",
        "Categoricals only show '(N categories)' from dtype metadata; no labels, no color dots.",
    ),
)
def _lazy_metadata_only() -> CaseOutput:
    import h5py

    with tmp_path(".h5ad") as path:
        write_lazy_demo_h5ad(path)
        with (
            h5py.File(path, "r") as f,
            ad.settings.override(repr_html_max_lazy_categories=0),
        ):
            return render(read_lazy(f))


@case(
    "storage",
    "Lazy AnnData (read_lazy, zarr)",
    slug="lazy-zarr",
    requires=("xarray", "zarr"),
    tags=("lazy", "zarr", "lazy-categorical", "colors"),
    expect=(
        "Header: <b>Lazy (Zarr)</b> badge and the zarr directory path.",
        "Same lazy behavior as [[lazy-h5ad]]: labels on demand, numeric columns '(lazy)'.",
    ),
)
def _lazy_zarr() -> CaseOutput:
    import zarr

    with tmp_path(".zarr") as path:
        write_lazy_demo_zarr(path)
        return render(read_lazy(zarr.open_group(path, mode="r")))


@case(
    "storage",
    "Lazy AnnData subset (view)",
    slug="lazy-view",
    requires=("xarray", "h5py"),
    tags=("lazy", "view"),
    expect=(
        "Subset shape 100 × 50 with lazy badges; no data loaded beyond category labels.",
    ),
)
def _lazy_view() -> CaseOutput:
    import h5py

    with tmp_path(".h5ad") as path:
        write_lazy_demo_h5ad(path)
        with h5py.File(path, "r") as f:
            return render(read_lazy(f)[:100, :50])


# =============================================================================
# 4. Scale & truncation
# =============================================================================


@case(
    "scale",
    "Auto-folding sections",
    slug="auto-folding",
    tags=("obs", "obsm"),
    expect=(
        "Sections with more than <code>repr_html_fold_threshold</code> (default 5) entries start collapsed: obs (15) and obsm (12).",
        "Clicking the header or fold icon expands/collapses.",
    ),
)
def _auto_folding() -> CaseOutput:
    adata = AnnData(np.zeros((20, 10)))
    for i in range(15):
        adata.obs[f"column_{i}"] = list(range(20))
    for i in range(12):
        adata.obsm[f"X_embedding_{i}"] = np.random.randn(20, 2).astype(np.float32)
    return render(adata)


@case(
    "scale",
    "Many obs columns (beyond max_items)",
    slug="many-obs-columns",
    tags=("obs", "var"),
    expect=(
        "obs has 250 columns: only the first <code>repr_html_max_items</code> (200) are listed, then a '… +50 more' indicator.",
        "Search still finds columns; the section header count says 250.",
    ),
)
def _many_obs_columns() -> CaseOutput:
    n = 30
    obs = pd.DataFrame(
        {f"qc_metric_{i:03d}": np.arange(n, dtype=float) for i in range(250)},
        index=[f"c{i}" for i in range(n)],
    )
    return render(AnnData(np.zeros((n, 4)), obs=obs))


@case(
    "scale",
    "max_items setting (30 layers, max 10)",
    slug="max-items-setting",
    tags=("layers",),
    expect=(
        "layers shows 10 entries and a truncation indicator for the remaining 20.",
    ),
)
def _max_items() -> CaseOutput:
    adata = AnnData(np.zeros((10, 5)))
    for i in range(30):
        adata.layers[f"layer_{i:02d}"] = np.zeros((10, 5), dtype=np.float32)
    with ad.settings.override(repr_html_max_items=10):
        return render(adata)


@case(
    "scale",
    "Many categories (truncation)",
    slug="many-categories",
    tags=("obs", "categorical", "colors"),
    expect=(
        "With <code>max_categories=20</code>: <code>cell_type</code> (30) shows 20 + '…+10'; <code>batch</code> (exactly 20) shows all.",
        "The ▼ button appears only when truncated and expands the full list.",
    ),
)
def _many_categories() -> CaseOutput:
    adata = AnnData(np.zeros((100, 10)))
    values = [f"type_{i}" for i in range(30)] * (100 // 30) + [
        f"type_{i}" for i in range(100 % 30)
    ]
    adata.obs["cell_type"] = pd.Categorical(values)
    adata.uns["cell_type_colors"] = palette(30)
    adata.obs["batch"] = pd.Categorical([f"batch_{i}" for i in range(20)] * 5)
    adata.uns["batch_colors"] = palette(20)[::-1]
    with ad.settings.override(repr_html_max_categories=20):
        return render(adata)


@case(
    "scale",
    "Wide DataFrame in obsm (150 columns)",
    slug="wide-obsm-dataframe",
    tags=("obsm", "dataframe"),
    expect=(
        "The type column shows <code>DataFrame (40 × 150)</code>.",
        "The column-name preview lists at most 100 names then '…+50', and is clipped to one line (the wrap toggle expands it).",
    ),
)
def _wide_obsm_df() -> CaseOutput:
    n = 40
    adata = AnnData(np.zeros((n, 5)))
    adata.obsm["morphology"] = pd.DataFrame(
        np.random.rand(n, 150),
        columns=[f"feature_{i:03d}_intensity" for i in range(150)],
        index=adata.obs_names,
    )
    adata.obsm["X_pca"] = np.random.rand(n, 10)
    return render(adata)


@case(
    "scale",
    "Expandable DataFrame in obsm",
    slug="expandable-obsm-dataframe",
    tags=("obsm", "dataframe"),
    expect=(
        "With <code>repr_html_dataframe_expand=True</code> the DataFrame row has an 'Expand' control.",
        "Expanding shows pandas' <code>_repr_html_()</code> table (respects <code>pd.options.display.max_rows</code>).",
    ),
)
def _expandable_df() -> CaseOutput:
    adata = AnnData(np.random.randn(30, 10).astype(np.float32))
    adata.obs["group"] = pd.Categorical(["A", "B", "C"] * 10)
    adata.obsm["spatial_metrics"] = pd.DataFrame(
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
        index=adata.obs_names,
    )
    adata.obsm["X_pca"] = np.random.randn(30, 5).astype(np.float32)
    with ad.settings.override(repr_html_dataframe_expand=True):
        return render(adata)


@case(
    "scale",
    "Large uns containers",
    slug="large-uns",
    tags=("uns",),
    expect=(
        "A 100,000-element list, a 2,000-key dict and a 1e6-element array render instantly with counts, not full contents.",
        "Long string preview is truncated with an ellipsis.",
    ),
)
def _large_uns() -> CaseOutput:
    adata = AnnData(np.zeros((5, 3)))
    adata.uns["big_list"] = list(range(100_000))
    adata.uns["big_dict"] = {f"key_{i:04d}": i for i in range(2_000)}
    adata.uns["big_array"] = np.zeros(1_000_000, dtype=np.float32)
    adata.uns["list_of_strings"] = [f"sample_{i}" for i in range(5_000)]
    adata.uns["long_string"] = "lorem ipsum " * 2_000
    return render(adata)


@case(
    "scale",
    "Huge shape (millions of cells)",
    slug="huge-shape",
    tags=("X", "sparse", "obs"),
    expect=(
        "Header shows 1,234,567 × 33,538 with thousands separators.",
        "Sparse X density/nnz and memory estimate are formatted readably; render stays fast.",
    ),
)
def _huge_shape() -> CaseOutput:
    n_obs, n_vars = 1_234_567, 33_538
    X = sp.csr_matrix((n_obs, n_vars), dtype=np.float32)
    obs = pd.DataFrame(index=pd.RangeIndex(n_obs).astype(str))
    obs["batch"] = pd.Categorical(
        np.repeat(["a", "b"], [n_obs // 2, n_obs - n_obs // 2])
    )
    var = pd.DataFrame(index=[f"g{i}" for i in range(n_vars)])
    return render(AnnData(X, obs=obs, var=var))


@case(
    "scale",
    "Unique-count limit",
    slug="unique-limit",
    tags=("obs", "string"),
    expect=(
        "With <code>repr_html_unique_limit=100</code> and 1,000 rows, string/numeric columns do not show '(N unique)'.",
        "Categoricals still show their categories.",
    ),
)
def _unique_limit() -> CaseOutput:
    n = 1000
    adata = AnnData(np.zeros((n, 2)))
    adata.obs["donor"] = [f"d{i % 7}" for i in range(n)]
    adata.obs["n_counts"] = np.arange(n)
    adata.obs["cat"] = pd.Categorical(["x", "y"] * (n // 2))
    with ad.settings.override(repr_html_unique_limit=100):
        return render(adata)


@case(
    "scale",
    "Very long field names",
    slug="long-field-names",
    tags=("obs", "obsm", "uns", "layers"),
    expect=(
        "Name column widens for long names but is capped by <code>repr_html_max_field_width</code> (400px).",
        "Overlong names end in an ellipsis; hover shows the full name; copy button copies the full name.",
    ),
)
def _long_names() -> CaseOutput:
    adata = AnnData(np.random.randn(20, 10).astype(np.float32))
    adata.obs["short"] = np.random.randn(20)
    adata.obs["this_is_a_moderately_long_column_name"] = np.random.randn(20)
    adata.obs[
        "this_is_an_extremely_long_column_name_that_should_test_the_max_width_setting"
    ] = np.random.randn(20)
    adata.obs["cell_type_annotation_from_automated_classifier_v2"] = pd.Categorical(
        ["A", "B"] * 10
    )
    adata.obsm["X_pca_computed_with_highly_variable_genes_batch_corrected"] = (
        np.random.randn(20, 5).astype(np.float32)
    )
    adata.uns["preprocessing_parameters_for_normalization_and_scaling"] = {
        "method": "log1p",
        "scale": True,
    }
    adata.layers["raw_counts_before_any_preprocessing_steps"] = np.random.randn(
        20, 10
    ).astype(np.float32)
    return render(adata)


@case(
    "scale",
    "README truncation (large README)",
    slug="readme-truncation",
    tags=("uns", "readme"),
    expect=(
        "The README is ~140,000 characters; the modal shows the first 100,000 (<code>repr_html_max_readme_size</code>) plus a truncation note.",
    ),
)
def _readme_truncation() -> CaseOutput:
    adata = AnnData(np.zeros((5, 5)))
    adata.uns["README"] = "# Large README\n\n" + "This is a very long README. " * 5000
    return render(adata)


# =============================================================================
# 5. Environments & theming
# =============================================================================


@case(
    "env",
    "No JavaScript (graceful degradation)",
    slug="no-javascript",
    tags=("no-js", "raw", "nested-anndata", "dataframe"),
    expect=(
        "All content visible; sections expanded; interactive buttons (fold, copy, search, wrap) hidden.",
        "Category lists and the 20-column <code>cell_measurements</code> DataFrame column list wrap naturally.",
        "Nested AnnData in uns and raw still expand via native <code>&lt;details&gt;</code>.",
        "A small 'interactive features require JavaScript' hint is shown.",
    ),
)
def _no_js() -> CaseOutput:
    adata = AnnData(np.random.randn(30, 15).astype(np.float32))
    adata.obs["group"] = pd.Categorical(["X", "Y", "Z"] * 10)
    adata.uns["group_colors"] = ["#e41a1c", "#377eb8", "#4daf4a"]
    for i in range(8):
        adata.obs[f"metric_{i}"] = np.random.randn(30)
    adata.obsm["X_pca"] = np.random.randn(30, 10).astype(np.float32)
    adata.layers["raw"] = np.random.randn(30, 15).astype(np.float32)
    adata.uns["nested_adata"] = AnnData(
        np.zeros((5, 3)),
        obs=pd.DataFrame({"label": ["A", "B", "C", "D", "E"]}),
    )
    adata.raw = adata.copy()
    measurements = [
        "area", "perimeter", "circularity", "eccentricity", "solidity", "extent",
        "major_axis_length", "minor_axis_length", "orientation", "mean_intensity",
        "max_intensity", "min_intensity", "std_intensity", "centroid_x", "centroid_y",
        "bbox_area", "convex_area", "euler_number", "equivalent_diameter", "filled_area",
    ]  # fmt: skip
    adata.obsm["cell_measurements"] = pd.DataFrame(
        np.random.rand(30, len(measurements)),
        columns=measurements,
        index=adata.obs_names,
    )
    return strip_script_tags(render(adata))


@case(
    "env",
    "No CSS (GitHub / untrusted notebook)",
    slug="no-css",
    tags=("no-css", "no-js"),
    expect=(
        "Rendered in an isolated iframe with all &lt;style&gt; and &lt;script&gt; removed.",
        "Still readable: one entry per line, monospace, comma-separated categories, a 'styled representation available' hint.",
        "Sections fold/unfold via native <code>&lt;details&gt;</code>/<code>&lt;summary&gt;</code>.",
    ),
)
def _no_css() -> CaseOutput:
    nocss = strip_style_and_script_tags(render(create_test_anndata()))
    return iframe(
        nocss,
        title="No-CSS repr",
        style="width:100%;border:1px solid #ccc;border-radius:4px;background:white;",
    )


@case(
    "env",
    "README icon without JavaScript",
    slug="readme-no-js",
    tags=("no-js", "readme"),
    expect=(
        "The ⓘ icon is present but does not open a modal.",
        "Hovering it shows the first ~500 characters of the README as a native tooltip.",
    ),
)
def _readme_no_js() -> CaseOutput:
    adata = AnnData(np.random.randn(20, 10).astype(np.float32))
    adata.obs["batch"] = pd.Categorical(["batch1", "batch2"] * 10)
    adata.uns["README"] = (
        "# Dataset Information\n\nThis dataset contains processed single-cell data.\n\n"
        "## Key Features\n- 20 cells, 10 genes\n- 2 batches\n\n"
        "For more details, see the full documentation.\n"
    )
    return strip_script_tags(render(adata))


@case(
    "env",
    "HTML repr disabled",
    slug="html-disabled",
    tags=("settings",),
    expect=(
        "<code>_repr_html_()</code> returns <code>None</code> with <code>repr_html_enabled=False</code>, so Jupyter falls back to the text repr shown below.",
    ),
)
def _html_disabled() -> CaseOutput:
    adata = create_test_anndata()
    with ad.settings.override(repr_html_enabled=False):
        result = adata._repr_html_()
    if result is not None:
        return f"<p style='color:#cf222e'>Expected None, got {len(result)} chars of HTML.</p>{result}"
    return f"<pre>_repr_html_() -> None\n\n{escape_html(repr(adata))}</pre>"


@dataclass(frozen=True)
class ThemeEnv:
    """A simulated host page for one theme pane."""

    label: str
    expect: Literal["light", "dark"]
    os: Literal["light", "dark"] = "light"
    html_attrs: str = ""
    body_attrs: str = ""
    page_css: str = ""


_THEME_PROBE_JS = """
(() => {
  const probe = document.getElementById("vt-probe");
  const os = matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  const repr = document.querySelector(".anndata-repr");
  let got = "?";
  if (repr) {
    const t = document.createElement("span");
    t.style.color = "var(--anndata-text-primary)";
    repr.appendChild(t);
    const m = getComputedStyle(t).color.match(/[\\d.]+/g);
    if (m) got = (+m[0] * 299 + +m[1] * 587 + +m[2] * 114) / 1000 > 128 ? "dark" : "light";
    t.remove();
  }
  const ok = got === EXPECTED;
  const sim = os === SIM_OS ? "" : " (OS simulation unsupported, toggle your OS theme)";
  probe.innerHTML = `OS ${os}${sim} · repr ${got} · expected ${EXPECTED} ` +
    `<b style="color:${ok ? "#1a7f37" : "#cf222e"}">${ok ? "PASS" : "FAIL"}</b>`;
})();
"""


def render_theme_panes(envs: Sequence[ThemeEnv], repr_html: str) -> str:
    """Render ``repr_html`` inside one iframe per simulated host environment.

    OS dark mode is simulated by setting ``color-scheme: dark`` on the iframe
    element: per CSS Color Adjust, the embedded document's
    ``prefers-color-scheme`` follows the embedding element's used color scheme.
    Each pane self-checks whether the repr picked the expected scheme.
    """
    panes = []
    for env in envs:
        probe = _THEME_PROBE_JS.replace("EXPECTED", repr(env.expect)).replace(
            "SIM_OS", repr(env.os)
        )
        doc = (
            f"<!doctype html><html {env.html_attrs}><head><meta charset='utf-8'>"
            "<style>body{margin:0;padding:10px;font:13px system-ui,sans-serif}"
            "#vt-probe{font:11px ui-monospace,monospace;margin-bottom:6px;opacity:.85}"
            f"{env.page_css}</style></head><body {env.body_attrs}>"
            f"<div id='vt-probe'>probe needs JS</div>{repr_html}"
            f"<script>{probe}</script></body></html>"
        )
        panes.append(
            '<figure class="vt-pane">'
            f"<figcaption>{escape_html(env.label)} · OS {env.os} · expect {env.expect}</figcaption>"
            + iframe(doc, title=env.label, style=f"color-scheme:{env.os};")
            + "</figure>"
        )
    return f'<div class="vt-panes">{"".join(panes)}</div>'


_THEME_EXPECT_COMMON = (
    "Each pane prints <b>PASS</b>/<b>FAIL</b>: whether the repr's text color matches the expected scheme.",
    "Visually: no light repr box on a dark page (or vice versa); category dots, warnings and badges stay legible.",
)


@case(
    "env",
    "Jupyter (JupyterLab / Notebook 7)",
    slug="theme-jupyter",
    tags=("theme", "dark-mode"),
    expect=(
        *_THEME_EXPECT_COMMON,
        "An explicit light Jupyter theme wins over OS dark mode.",
    ),
    notes="JupyterLab sets <code>data-jp-theme-light</code> and <code>jp-Theme-*</code> on <code>&lt;body&gt;</code>.",
)
def _theme_jupyter() -> CaseOutput:
    light = 'data-jp-theme-light="true" data-jp-theme-name="JupyterLab Light" class="jp-Theme-Light"'
    dark = 'data-jp-theme-light="false" data-jp-theme-name="JupyterLab Dark" class="jp-Theme-Dark"'
    envs = [
        ThemeEnv("JupyterLab Light", "light", body_attrs=light, page_css="body{background:#fff;color:#000}"),
        ThemeEnv("JupyterLab Light", "light", os="dark", body_attrs=light, page_css="body{background:#fff;color:#000}"),
        ThemeEnv("JupyterLab Dark", "dark", body_attrs=dark, page_css="body{background:#111;color:#ddd}"),
        ThemeEnv("JupyterLab Dark", "dark", os="dark", body_attrs=dark, page_css="body{background:#111;color:#ddd}"),
    ]  # fmt: skip
    return render_theme_panes(envs, render(create_theme_demo_anndata()))


@case(
    "env",
    "VS Code notebooks",
    slug="theme-vscode",
    tags=("theme", "dark-mode"),
    expect=(
        *_THEME_EXPECT_COMMON,
        "High-contrast themes: HC dark should render dark, HC light should render light.",
    ),
    notes=(
        "VS Code webviews set <code>body.vscode-{light,dark,high-contrast,high-contrast-light}</code> "
        "and <code>data-vscode-theme-kind</code>."
    ),
)
def _theme_vscode() -> CaseOutput:
    def attrs(kind: str) -> str:
        return f'class="{kind}" data-vscode-theme-kind="{kind}"'

    envs = [
        ThemeEnv("VS Code Light+", "light", os="dark", body_attrs=attrs("vscode-light"), page_css="body{background:#fff;color:#3b3b3b}"),
        ThemeEnv("VS Code Dark+", "dark", body_attrs=attrs("vscode-dark"), page_css="body{background:#1f1f1f;color:#ccc}"),
        ThemeEnv("VS Code High Contrast (dark)", "dark", body_attrs=attrs("vscode-high-contrast"), page_css="body{background:#000;color:#fff}"),
        ThemeEnv("VS Code High Contrast Light", "light", os="dark", body_attrs=attrs("vscode-high-contrast-light"), page_css="body{background:#fff;color:#292929}"),
    ]  # fmt: skip
    return render_theme_panes(envs, render(create_theme_demo_anndata()))


@case(
    "env",
    "Sphinx docs: Furo",
    slug="theme-furo",
    tags=("theme", "dark-mode"),
    expect=(
        *_THEME_EXPECT_COMMON,
        "Furo 'auto' follows the OS: auto + OS dark is a dark page, so the repr should be dark.",
    ),
    notes="Furo sets <code>body[data-theme=light|dark|auto]</code>; 'auto' uses a <code>prefers-color-scheme</code> media query.",
)
def _theme_furo() -> CaseOutput:
    auto_css = (
        "body{background:#fff;color:#000}"
        "@media (prefers-color-scheme: dark){body[data-theme=auto]{background:#131416;color:#cfd0d0}}"
    )
    envs = [
        ThemeEnv("Furo light", "light", os="dark", body_attrs='data-theme="light"', page_css="body{background:#fff;color:#000}"),
        ThemeEnv("Furo dark", "dark", body_attrs='data-theme="dark"', page_css="body{background:#131416;color:#cfd0d0}"),
        ThemeEnv("Furo auto", "light", body_attrs='data-theme="auto"', page_css=auto_css),
        ThemeEnv("Furo auto", "dark", os="dark", body_attrs='data-theme="auto"', page_css=auto_css),
    ]  # fmt: skip
    return render_theme_panes(envs, render(create_theme_demo_anndata()))


@case(
    "env",
    "Sphinx docs: pydata-sphinx-theme / sphinx-book-theme",
    slug="theme-pydata",
    tags=("theme", "dark-mode"),
    expect=(*_THEME_EXPECT_COMMON,),
    notes=(
        "pydata-sphinx-theme (and sphinx-book-theme, used by scverse docs) resolve 'auto' in JS "
        "and set <code>html[data-theme=light|dark]</code> plus <code>data-mode</code>."
    ),
)
def _theme_pydata() -> CaseOutput:
    envs = [
        ThemeEnv("pydata light", "light", os="dark", html_attrs='data-theme="light" data-mode="light"', page_css="body{background:#fff;color:#222832}"),
        ThemeEnv("pydata dark", "dark", html_attrs='data-theme="dark" data-mode="dark"', page_css="body{background:#14181e;color:#ced6dd}"),
    ]  # fmt: skip
    return render_theme_panes(envs, render(create_theme_demo_anndata()))


@case(
    "env",
    "Theme-less pages (OS light vs dark)",
    slug="theme-none",
    tags=("theme", "dark-mode"),
    expect=(
        *_THEME_EXPECT_COMMON,
        "A page that declares no theme and no <code>color-scheme</code> stays light in OS dark mode (white canvas), so the repr must stay light too.",
        "A page that opts into <code>color-scheme: light dark</code> turns dark in OS dark mode, so the repr should follow.",
    ),
    notes=(
        "Covers classic Notebook, nbconvert/static HTML, GitHub-like pages and the top of this harness. "
        "Use the 'Toggle dark mode' button (adds <code>body.dark-mode</code>) to check the harness itself."
    ),
)
def _theme_none() -> CaseOutput:
    scheme_css = (
        ":root{color-scheme:light dark}body{background:Canvas;color:CanvasText}"
    )
    envs = [
        ThemeEnv("plain page", "light"),
        ThemeEnv("plain page", "light", os="dark"),
        ThemeEnv("page with color-scheme: light dark", "light", page_css=scheme_css),
        ThemeEnv(
            "page with color-scheme: light dark", "dark", os="dark", page_css=scheme_css
        ),
    ]
    return render_theme_panes(envs, render(create_theme_demo_anndata()))


# =============================================================================
# 6. Robustness & security
# =============================================================================


@case(
    "robust",
    "Special characters in names",
    slug="special-characters",
    tags=("obs", "uns", "xss", "unicode"),
    expect=(
        "<code>column&lt;with&gt;html</code>, ampersands, quotes and Japanese characters render literally.",
        "No broken layout and no HTML interpretation.",
    ),
)
def _special_chars() -> CaseOutput:
    adata = AnnData(np.zeros((5, 3)))
    adata.obs["column<with>html"] = list(range(5))
    adata.obs["column&ampersand"] = list(range(5))
    adata.uns["key\"with'quotes"] = "value"
    adata.uns["unicode_日本語"] = "japanese"
    return render(adata)


@case(
    "robust",
    "Serialization warnings",
    slug="serialization-warnings",
    tags=("obs", "var", "layers", "obsm", "uns", "serialization"),
    expect=(
        "<b>Red (fails now):</b> obs list/dict/custom-object/tuple-named columns; var datetime/timedelta; <code>layers[('tuple','key')]</code>; uns custom object, lambda, nested bad value.",
        "<b>Yellow (fails in future):</b> <code>obs['path/slash']</code>, <code>obsm['path/embed']</code>.",
        "<b>No warning:</b> <code>var['gène_名前']</code>, normal floats/strings, <code>uns.valid_dict</code>.",
    ),
)
def _serialization() -> CaseOutput:
    class CustomObject:
        def __repr__(self):
            return "CustomObject()"

    adata = AnnData(X=np.eye(5))

    def obj_col(values: list[object]) -> pd.Series:
        return pd.Series(values, index=adata.obs_names, dtype=object)

    adata.obs["list_values"] = obj_col([["a", "b"], ["c"], ["d"], ["e"], ["f"]])
    adata.obs["dict_values"] = obj_col([{"k": i} for i in range(1, 6)])
    adata.obs["custom_obj"] = obj_col([CustomObject() for _ in range(5)])
    adata.obs["path/slash"] = ["a", "b", "c", "d", "e"]
    adata.obs[("tuple", "name")] = [1, 2, 3, 4, 5]
    adata.var["datetime_col"] = pd.to_datetime([f"2024-01-0{i}" for i in range(1, 6)])
    adata.var["timedelta_col"] = pd.to_timedelta([f"{i} days" for i in range(1, 6)])
    adata.var["gène_名前"] = ["a", "b", "c", "d", "e"]
    adata.var["normal_col"] = [1.0, 2.0, 3.0, 4.0, 5.0]
    adata.var["string_col"] = ["a", "b", "c", "d", "e"]
    adata.layers[("tuple", "key")] = np.eye(5)  # type: ignore[index]  # intentionally invalid key
    adata.obsm["path/embed"] = np.random.randn(5, 2)
    adata.uns["custom_obj"] = CustomObject()
    adata.uns["lambda_func"] = lambda x: x
    adata.uns["nested_bad"] = {"ok": 1, "bad": CustomObject()}
    adata.uns["valid_dict"] = {"a": 1, "b": [1, 2, 3]}
    return render(adata)


@case(
    "robust",
    "Subclass with unknown and failing attributes",
    slug="unknown-sections",
    tags=("custom", "subclass", "errors"),
    expect=(
        "The <code>custom_data</code> mapping appears in an 'other' section at the bottom (nothing silently hidden).",
        "The <code>failing_data</code> property raises; it is shown as inaccessible in 'other', not dropped.",
    ),
)
def _unknown_sections() -> CaseOutput:
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

    adata = ExtendedAnnData(
        np.random.randn(50, 100).astype(np.float32),
        obs=pd.DataFrame({"cluster": pd.Categorical(["A", "B"] * 25)}),
    )
    adata._custom_mappings = {
        "embedding": np.random.randn(50, 2),
        "config": {"param1": 1, "param2": "value"},
    }
    return render(adata)


@case(
    "robust",
    "Failing section accessors (real errors)",
    slug="failing-sections",
    tags=("varm", "layers", "errors"),
    expect=(
        "obs, var, uns, obsm, obsp, varp render normally.",
        "<code>varm</code> and <code>layers</code> show an error row with the exception message instead of crashing.",
        "X lives in <code>layers[None]</code>, so it shows the same I/O error (with its message), not a bare exception type.",
    ),
    notes="Uses <code>unittest.mock.patch</code> so the real <code>generate_repr_html</code> pipeline hits the exceptions.",
)
def _failing_sections() -> CaseOutput:
    from unittest.mock import PropertyMock, patch

    class FailingMapping:
        """A mapping that raises an error when iterated."""

        isbacked = False  # AnnData.isbacked / .X consult layers since X is layers[None]

        def __init__(self, error_msg: str):
            self._error_msg = error_msg

        def keys(self):
            raise RuntimeError(self._error_msg)

        def __len__(self):
            return 1  # Report as non-empty so it tries to render

        def __iter__(self):
            raise RuntimeError(self._error_msg)

        def __getitem__(self, key):
            raise RuntimeError(self._error_msg)

        def get(self, key, default=None):
            raise RuntimeError(self._error_msg)

        def __contains__(self, key):
            raise RuntimeError(self._error_msg)

    rng = np.random.default_rng(42)
    adata = AnnData(
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
    failing_varm = FailingMapping(
        "Failed to decompress data block (corrupted zarr chunk)"
    )
    failing_layers = FailingMapping(
        "IOError: [Errno 5] Input/output error reading '/data/counts.h5'"
    )
    with (
        patch.object(
            type(adata), "varm", new_callable=PropertyMock, return_value=failing_varm
        ),
        patch.object(
            type(adata),
            "layers",
            new_callable=PropertyMock,
            return_value=failing_layers,
        ),
    ):
        return render(adata)


class _Exploding:
    """Marker type routed to the failing formatters below."""

    shape = (3, 3)
    dtype = np.dtype("float32")

    def __repr__(self) -> str:
        return "_Exploding()"


class _FailingTypeFormatter(TypeFormatter):
    """A third-party TypeFormatter with a bug in ``format()``."""

    priority = 10_000

    def can_format(self, obj, context):
        return isinstance(obj, _Exploding)

    def format(self, obj, context) -> FormattedOutput:
        raise ZeroDivisionError("third-party formatter bug <script>alert(1)</script>")


class _FailingCanFormatFormatter(TypeFormatter):
    """A third-party TypeFormatter whose ``can_format()`` raises for everything."""

    priority = 10_001

    def can_format(self, obj, context):
        if isinstance(obj, str) and obj == "trigger-can-format-bug":
            raise KeyError("can_format bug")
        return False

    def format(self, obj, context) -> FormattedOutput:  # pragma: no cover
        raise NotImplementedError


class _FailingSectionFormatter(SectionFormatter):
    """A third-party SectionFormatter whose ``get_entries()`` raises."""

    section_name = "broken_plugin"

    @property
    def after_section(self) -> str:
        return "obsm"

    def should_show(self, obj) -> bool:
        return "broken_plugin_marker" in getattr(obj, "uns", {})

    def get_entries(self, obj, context) -> list[FormattedEntry]:
        raise RuntimeError("plugin section exploded")


@case(
    "robust",
    "Failing third-party formatters",
    slug="failing-formatters",
    tags=("obsm", "uns", "layers", "custom", "errors", "xss"),
    expect=(
        "Entries handled by the buggy TypeFormatter fall back to an error row (or the default formatter), the rest renders.",
        "The exception message containing <code>&lt;script&gt;</code> is escaped.",
        "A <code>can_format()</code> that raises does not break other entries.",
        "The <code>broken_plugin</code> section shows an error instead of aborting the whole repr.",
    ),
)
def _failing_formatters() -> CaseOutput:
    adata = AnnData(np.zeros((3, 3)))
    adata.obsm._data["exploding"] = _Exploding()  # type: ignore[union-attr]  # bypass validation
    adata.layers["fine"] = np.ones((3, 3))
    adata.uns["exploding"] = _Exploding()
    adata.uns["triggers_can_format"] = "trigger-can-format-bug"
    adata.uns["broken_plugin_marker"] = True
    with (
        temporarily_registered(_FailingTypeFormatter()),
        temporarily_registered(_FailingCanFormatFormatter()),
        temporarily_registered(_FailingSectionFormatter()),
    ):
        return render(adata)


@case(
    "robust",
    "Evil AnnData (adversarial robustness)",
    slug="evil-anndata",
    tags=(
        "obs",
        "var",
        "obsm",
        "varm",
        "varp",
        "layers",
        "uns",
        "xss",
        "unicode",
        "errors",
        "colors",
        "readme",
        "nested-anndata",
        "circular",
    ),
    expect=(
        "No crash, no script execution, no layout breakout (the card below stays intact).",
        "Errors in <span style='color:#dc3545'>red</span>, warnings in <span style='color:#d29922'>orange</span>, rows tinted accordingly.",
        "All XSS payloads in column names, category values, DataFrame columns, uns keys/values, type and exception names show as literal text.",
        "Bad colors (too many/few, invalid, CSS/url injection, 1000-char strings) never produce dots with injected styles.",
        "varp shows 200 of 300 entries plus a truncation indicator; 10k-category column, 50 KB string and 500-item dict are truncated.",
        "Circular references (dict, self-referencing AnnData) terminate.",
    ),
    notes=(
        "<b>Crashing objects in uns:</b> exploding_repr/len/str, lying_object, infinite_len (10^18), "
        "exploding_shape/dtype, xss_via_exception, xss_via_type_name, long_error_object_uns, unknown_anndata_type (orange).<br>"
        "<b>Evil README:</b> displayed via textContent, so nothing can fire; includes script/style tags, RTL override, null bytes, template injection and a 50 KB bomb.<br>"
        "<b>varm:</b> long_error_object (error should truncate).<br>"
        "<b>Circular:</b> circular_dict, self_reference, child_with_parent_ref.<br>"
        "<b>Nesting:</b> deeply_nested_15_levels; nested_adata_with_errors shows tinted rows inside.<br>"
        "<b>obs names:</b> &lt;script&gt;, &lt;img onerror&gt;, onclick=, &lt;svg onload&gt;, javascript:, emoji, CJK, RTL override, null bytes.<br>"
        "<b>var names:</b> &lt;/style&gt;&lt;script&gt;, &lt;/div&gt; breakout.<br>"
        "<b>uns strings:</b> SVG XSS, mutation XSS, UTF-7, BOM prefix."
    ),
)
def _evil() -> CaseOutput:  # noqa: PLR0915
    class ExplodingRepr:
        """Object whose __repr__ crashes."""

        def __repr__(self):
            raise RuntimeError("BOOM! __repr__ exploded")

    class ExplodingLen:
        """Object whose __len__ crashes."""

        def __len__(self):
            raise MemoryError("BOOM! __len__ exploded")

    class ExplodingStr:
        """Object whose __str__ crashes."""

        def __str__(self):
            raise ValueError("BOOM! __str__ exploded")

        def __repr__(self):
            return "ExplodingStr(str crashes)"

    class LyingObject:
        """Object that lies about all its properties."""

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

    class InfiniteLen:
        """Object claiming impossibly large length."""

        def __len__(self):
            return 10**18  # 1 quintillion items

        def __repr__(self):
            return "InfiniteLen(10^18 items)"

    class ExplodingShape:
        """Object whose shape property explodes."""

        @property
        def shape(self):
            raise TypeError("BOOM! shape exploded")

        def __repr__(self):
            return "ExplodingShape(.shape crashes)"

    class ExplodingDtype:
        """Object whose dtype property explodes."""

        shape = (10, 10)  # Normal shape

        @property
        def dtype(self):
            raise TypeError("BOOM! dtype exploded")

        def __repr__(self):
            return "ExplodingDtype(.dtype crashes)"

    adata_evil = AnnData(X=np.random.rand(50, 30).astype(np.float32))

    # XSS injection attempts (6 variants)
    adata_evil.obs["normal_column"] = np.random.choice(["A", "B", "C"], size=50)
    adata_evil.obs['<script>alert("XSS")</script>'] = np.random.randint(0, 10, size=50)
    adata_evil.obs["<img onerror=alert(1)>"] = np.random.rand(50)
    adata_evil.obs['onclick="evil()"'] = np.random.rand(50)
    adata_evil.obs["<svg onload=alert(1)>"] = np.random.rand(50)
    adata_evil.obs["javascript:alert(1)"] = np.random.rand(50)

    # Unicode bombs
    adata_evil.obs["emoji_\U0001f4a9_poop"] = np.random.rand(50)
    adata_evil.obs["chinese_\u4e2d\u6587"] = pd.Categorical(
        np.random.choice(["cat", "dog", "bird"], size=50)
    )
    adata_evil.obs["rtl_\u202eEVIL\u202c_override"] = np.random.rand(50)
    adata_evil.obs["null\x00byte\x00col"] = np.random.rand(50)

    # XSS in CATEGORY VALUES (not just column names)
    xss_categories = [
        '<script>alert("cat")</script>',
        "<img onerror=alert(1)>",
        '<svg onload="evil()">',
        "normal_category",
        "<div onclick=bad()>",
        "javascript:void(0)",
    ]
    adata_evil.obs["xss_category_values"] = pd.Categorical(
        np.random.choice(xss_categories, size=50), categories=xss_categories
    )

    # HTML/CSS breakout
    adata_evil.var["gene_normal"] = [f"gene_{i}" for i in range(30)]
    adata_evil.var["</style><script>bad()</script>"] = np.random.rand(30)
    adata_evil.var["</div></div></div>breakout"] = np.random.rand(30)

    # CRASHING OBJECTS in uns
    adata_evil.uns["normal"] = {"key": "value", "nested": {"a": 1, "b": 2}}
    adata_evil.uns["exploding_repr"] = ExplodingRepr()
    adata_evil.uns["exploding_len"] = ExplodingLen()
    adata_evil.uns["exploding_str"] = ExplodingStr()
    adata_evil.uns["lying_object"] = LyingObject()
    adata_evil.uns["infinite_len"] = InfiniteLen()
    adata_evil.uns["exploding_shape"] = ExplodingShape()
    adata_evil.uns["exploding_dtype"] = ExplodingDtype()

    # XSS VIA EXCEPTION CLASS NAME - tests that error messages escape __name__
    class XSSException(Exception):
        pass

    XSSException.__name__ = "<img src=x onerror=alert('exception')>"

    class XSSViaException:
        """Object whose .shape raises exception with XSS payload in __name__."""

        @property
        def shape(self):
            raise XSSException("gotcha")

    adata_evil.uns["xss_via_exception"] = XSSViaException()

    # XSS VIA TYPE NAME - tests that type names are escaped
    class XSSViaTypeName:
        pass

    XSSViaTypeName.__name__ = "<script>alert('type')</script>"
    adata_evil.uns["xss_via_type_name"] = XSSViaTypeName()

    # UNKNOWN TYPE WARNING (orange text) - object pretending to be from anndata package
    class FakeAnndataType:
        """Unknown type from anndata package triggers warning (not error)."""

        __module__ = "anndata.experimental.fake"

        def __repr__(self):
            return "FakeAnndataType()"

    adata_evil.uns["unknown_anndata_type"] = FakeAnndataType()

    # EVIL README - displayed as plain text via textContent (not innerHTML),
    # so none of these vectors can fire. This verifies the data-readme
    # attribute handles edge-case content without breaking HTML structure.
    adata_evil.uns["README"] = (
        """Evil README - XSS and Injection Test

<script>alert('XSS in readme!')</script>
<img src=x onerror="alert('img onerror')">
</div></div></table></section>
<style>body { display: none !important; }</style>

Unicode: \u202eSIHT DAER\u202c
Null bytes: before\x00after
Emoji bomb: \U0001f480\U0001f480\U0001f480\U0001f480\U0001f480

{{constructor.constructor('alert(1)')()}}
${alert('template_literal')}

Size bomb below (50KB):
"""
        + "A" * 50000
    )

    # CIRCULAR REFERENCES
    circular_dict: dict = {"level1": {"level2": {}}}
    circular_dict["level1"]["level2"]["back_to_start"] = circular_dict
    adata_evil.uns["circular_dict"] = circular_dict
    adata_evil.uns["self_reference"] = adata_evil
    child_adata = AnnData(np.zeros((5, 5)))
    child_adata.uns["parent_ref"] = adata_evil
    adata_evil.uns["child_with_parent_ref"] = child_adata

    # DEEPLY NESTED - 15 levels
    deeply_nested: dict = {}
    current = deeply_nested
    for i in range(15):
        current[f"level_{i}"] = {}
        current = current[f"level_{i}"]
    current["bottom"] = "reached the bottom!"
    adata_evil.uns["deeply_nested_15_levels"] = deeply_nested

    # XSS in uns keys, size bombs
    adata_evil.uns["<script>evil()</script>"] = "XSS key"
    adata_evil.uns["giant_string_50kb"] = "X" * 50_000
    adata_evil.uns["many_items_500"] = {f"item_{i:04d}": i for i in range(500)}
    huge_cats = [f"category_{i:05d}" for i in range(10000)]
    adata_evil.obs["huge_categorical_10k"] = pd.Categorical(
        np.random.choice(huge_cats[:50], size=50), categories=huge_cats
    )

    # MANY ENTRIES IN ONE SECTION (default max is 200). Populate the internal
    # store directly after one validated entry to skip slow validation.
    tiny_sparse = sp.csr_matrix(([1.0], ([0], [0])), shape=(30, 30))
    adata_evil.varp["varp_000"] = tiny_sparse
    for i in range(1, 300):
        adata_evil.varp._data[f"varp_{i:03d}"] = tiny_sparse  # type: ignore[union-attr]

    # BAD COLORS - various malformed color arrays
    bad_colors: dict[str, tuple[list[str], list[str]]] = {
        "cat_too_many_colors": (
            ["A", "B", "C"],
            ["red", "green", "blue", "yellow", "purple", "orange"],
        ),
        "cat_too_few_colors": (["X", "Y", "Z", "W"], ["red"]),
        "cat_bad_colors": (["alpha", "beta"], ["not_a_color", "also_invalid"]),
        "cat_strange_colors": (
            ["one", "two", "three"],
            ["#FF0000", "rgb(0,255,0)", "rgba(0,0,255,0.5)"],
        ),
        "cat_empty_colors": (["p", "q"], []),
        # CSS injection attempts (must be blocked by whitelist)
        "cat_css_injection": (
            ["x", "y", "z"],
            [
                "#ff0000",
                "blue; } .adata-table { display:none } .x {",
                "red; background-image: url(https://evil.com/steal)",
            ],
        ),
        # URL/expression injection (must be blocked)
        "cat_url_injection": (
            ["a", "b"],
            ["url(https://evil.com/track)", "expression(alert(1))"],
        ),
        # Very long color strings (DoS protection)
        "cat_long_colors": (["m", "n"], ["red" + "x" * 1000, "blue"]),
    }
    for col, (cats, colors) in bad_colors.items():
        adata_evil.obs[col] = pd.Categorical(np.random.choice(cats, size=50))
        adata_evil.uns[f"{col}_colors"] = colors

    # NESTED ANNDATA WITH ERRORS (should show yellow/red rows in nested content)
    nested_with_errors = AnnData(np.zeros((10, 10)))
    nested_with_errors.uns["bad_obj_in_nested"] = ExplodingRepr()
    nested_with_errors.uns["another_bad"] = LyingObject()
    adata_evil.uns["nested_adata_with_errors"] = nested_with_errors

    class VeryLongErrorObject:
        """Object that produces a very long error message."""

        @property
        def shape(self):
            raise TypeError(
                "This is a VERY LONG ERROR MESSAGE that should be properly truncated. "
                * 10
                + "It contains lots of details about what went wrong: "
                + "ValueError: The input array has shape (100, 200, 300) but expected (50, 100). "
                + "Additional context: This error occurred while processing the data matrix. "
                + "Stack trace would go here with many lines of debugging information. "
                * 5
            )

        def __repr__(self):
            return "VeryLongErrorObject(produces long error)"

    adata_evil.uns["long_error_object_uns"] = VeryLongErrorObject()

    # SVG XSS, mutation XSS, encoding attacks
    adata_evil.uns["svg_script"] = "<svg><script>alert(1)</script></svg>"
    adata_evil.uns["svg_onload"] = '<svg onload="alert(1)">'
    adata_evil.uns["mxss_unclosed"] = "<img src=x onerror=alert(1)//"
    adata_evil.uns["mxss_nested"] = "<div<script>alert(1)</script>>"
    adata_evil.uns["utf7_script"] = "+ADw-script+AD4-alert(1)+ADw-/script+AD4-"
    adata_evil.uns["bom_prefix"] = "\ufeffmalicious_content"

    # Long error object in varm (bypass validation via internal store)
    adata_evil.varm["gene_scores"] = np.random.rand(30, 5)
    adata_evil.varm._data["long_error_object"] = VeryLongErrorObject()  # type: ignore[union-attr]

    # XSS in DataFrame COLUMN NAMES (shown in obsm preview)
    adata_evil.obsm["X_evil_df_cols"] = pd.DataFrame(
        {
            '<script>alert("col")</script>': np.random.rand(50),
            "<img onerror=alert(1)>": np.random.rand(50),
            "normal_col": np.random.rand(50),
        },
        index=adata_evil.obs_names,
    )

    # Standard sections to show they still work
    adata_evil.obsm["X_pca"] = np.random.rand(50, 10)
    adata_evil.obsm["X_umap"] = np.random.rand(50, 2)
    adata_evil.layers["raw_counts"] = np.random.randint(0, 100, (50, 30))
    adata_evil.obsp["connectivities"] = sp.random(50, 50, density=0.1, format="csr")

    # The warnings from crashing objects are expected; don't list them in the page.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return render(adata_evil)


# =============================================================================
# 7. Extensibility / ecosystem
# =============================================================================


class AnalysisHistoryFormatter(TypeFormatter):
    """Example TypeFormatter for analysis history data with embedded type hint.

    A package would decorate this with ``@register_formatter``; the case below
    registers it only while rendering.
    """

    priority = 100  # High priority to check before fallback

    def can_format(self, obj, context):
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

        runs = data.get("runs", [])
        params = data.get("params", {})

        html_parts = ['<div style="font-size:11px;">']
        if runs:
            html_parts.append(f"<strong>{len(runs)} runs</strong>")
        if params:
            param_str = ", ".join(f"{k}={v}" for k, v in list(params.items())[:3])
            if len(params) > 3:
                param_str += "..."
            html_parts.append(f" · params: {escape_html(param_str)}")
        html_parts.append("</div>")

        return FormattedOutput(
            type_name="analysis history",
            preview_html="".join(html_parts),  # Use preview_html for inline preview
        )


@case(
    "ext",
    "uns value previews and type hints",
    slug="uns-type-hints",
    tags=("uns", "type-formatter"),
    expect=(
        "Simple types (str, int, float, bool, None) show inline previews; <code>long_string</code> is truncated.",
        "<code>small_list</code>/<code>small_dict</code> show content; <code>larger_dict</code> shows a key count.",
        "<code>analysis_history</code>: custom TypeFormatter renders '3 runs · params: …'.",
        "<code>unregistered_data</code> and <code>string_hint</code>: hint without formatter → 'import otherpackage to enable'.",
    ),
    notes="The <code>__anndata_repr__</code> type hint lets packages register renderers for their data stored in uns.",
)
def _uns_type_hints() -> CaseOutput:
    adata = AnnData(np.zeros((10, 5)))
    adata.uns["string_param"] = "A short string value"
    adata.uns["long_string"] = (
        "This is a very long string that should be truncated in the preview because "
        "it exceeds the maximum length allowed for display in the meta column"
    )
    adata.uns["int_param"] = 42
    adata.uns["float_param"] = 3.14159265359
    adata.uns["bool_param"] = True
    adata.uns["none_param"] = None
    adata.uns["small_list"] = [1, 2, 3]
    adata.uns["small_dict"] = {"a": 1, "b": 2}
    adata.uns["larger_dict"] = {f"key{i}": f"val{i}" for i in range(1, 6)}
    # Type hint WITH registered renderer (shows custom HTML)
    adata.uns["analysis_history"] = {
        "__anndata_repr__": "example.history",
        "runs": [{"id": 1}, {"id": 2}, {"id": 3}],
        "params": {"method": "umap", "n_neighbors": 15, "metric": "euclidean"},
    }
    # Type hint WITHOUT registered renderer (shows fallback with import hint)
    adata.uns["unregistered_data"] = {
        "__anndata_repr__": "otherpackage.custom_type",
        "data": {"some": "data", "values": [1, 2, 3]},
    }
    adata.uns["string_hint"] = (
        "__anndata_repr__:otherpackage.config::{'setting': 'value'}"
    )
    with temporarily_registered(AnalysisHistoryFormatter()):
        return render(adata)


@case(
    "ext",
    "Custom sections (TreeData)",
    slug="treedata-sections",
    tags=("custom", "section-formatter"),
    expect=(
        "<code>obst</code> (after obsm) and <code>vart</code> (after varm) are foldable sections with SVG tree previews; trees with &gt;30 leaves show a text message.",
        "<code>tree</code> (right after X) is a compact non-foldable line rendered via <code>render_html()</code>.",
        "If a yellow note says a stand-in is used, the real treedata package is incompatible with this anndata.",
    ),
    notes=(
        "<a href='https://treedata.readthedocs.io/en/latest/' target='_blank'>TreeData</a> registers three "
        "SectionFormatters (<a href='https://github.com/scverse/ecosystem-packages/pull/282' "
        "target='_blank'>scverse ecosystem PR</a>)."
    ),
)
def _treedata() -> CaseOutput:
    tdata, note = create_test_treedata()
    return render(tdata), note


@case(
    "ext",
    "MuData (multimodal, SectionFormatter for .mod)",
    slug="mudata",
    requires=("mudata",),
    tags=("custom", "section-formatter", "nested-anndata"),
    expect=(
        "A <code>mod</code> section right after X lists rna (100×50), atac (100×30) and prot (80×20), each expandable.",
        "MuData's internal <code>obsmap</code>/<code>varmap</code>/<code>axis</code> are suppressed (not in 'other').",
    ),
    notes="MuData reuses anndata's repr by registering a SectionFormatter and calling <code>generate_repr_html()</code> on itself.",
)
def _mudata() -> CaseOutput:
    if not HAS_MUDATA:
        raise SkipCase("mudata failed to import")
    return generate_repr_html(create_test_mudata())


@case(
    "ext",
    "SpatialData (custom _repr_html_ from building blocks)",
    slug="spatialdata",
    tags=("custom", "building-blocks", "nested-anndata"),
    expect=(
        "Custom header 'SpatialData' with Zarr badge and path; coordinate systems list with tooltips.",
        "images/labels/points/shapes sections with <code>[c, y, x]</code>-style previews; tables embed expandable AnnData.",
        "A <code>transforms</code> section contributed via a separate <code>FormatterRegistry</code>.",
    ),
    notes=(
        "Uses <code>get_css()</code>, <code>get_javascript()</code>, <code>render_section()</code>, "
        "<code>render_formatted_entry()</code>, <code>render_badge()</code>, <code>render_search_box()</code>, "
        "<code>generate_repr_html()</code> and <code>FormatterRegistry</code>."
    ),
)
def _spatialdata() -> CaseOutput:
    if not HAS_SPATIALDATA_EXAMPLE:
        raise SkipCase("building blocks failed to import")
    return render(create_test_spatialdata())


ONTOLOGY_METADATA_KEY = "__ontology_annotations__"


class OntologyAnnotatedCategoricalFormatter(TypeFormatter):
    """
    Example TypeFormatter for columns annotated with ontology metadata.

    This demonstrates how ecosystem packages (bionty, lamindb, cellxgene, ...)
    can enhance the HTML repr for categorical columns in obs/var by:
    1. Checking if the column has ontology metadata in uns
    2. Rendering enhanced type info (registry name, validation status)
    3. Adding tooltips with ontology IDs

    To use this pattern in your package:
    1. Define a metadata convention (e.g., uns["__mypackage_annotations__"])
    2. Create a TypeFormatter with sections=("obs", "var")
    3. Use context.adata_ref and context.key to look up metadata
    4. Decorate it with ``@register_formatter``

    See: src/anndata/_repr/registry.py for TypeFormatter API
    """

    priority = 115  # Higher than CategoricalFormatter (110)
    sections = ("obs", "var")  # Only apply to obs/var columns

    def can_format(self, obj, context):
        """Check if this column has ontology metadata.

        The context parameter provides access to:
        - context.adata_ref: reference to root AnnData for uns lookups
        - context.key: current entry key being formatted
        - context.section: current section ("obs", "var", etc.)
        """
        if not (isinstance(obj, pd.Series) and hasattr(obj, "cat")):
            return False
        if context.adata_ref is None or context.key is None:
            return False
        annotations = context.adata_ref.uns.get(ONTOLOGY_METADATA_KEY, {})
        section_annotations = annotations.get(context.section, {})
        return context.key in section_annotations

    def format(self, obj, context):
        """Render the categorical with ontology information."""
        annotations = context.adata_ref.uns[ONTOLOGY_METADATA_KEY]
        col_info = annotations[context.section][context.key]

        registry = col_info.get("registry", "unknown")
        ontology_id = col_info.get("ontology_id", "")
        validated = col_info.get("validated", True)
        unmapped_count = col_info.get("unmapped_count", 0)

        n_cats = len(obj.cat.categories)
        type_name = f"category[{registry}] ({n_cats})"

        categories = list(obj.cat.categories[:5])
        cat_html = ", ".join(
            f'<span style="color: var(--anndata-category-color, #666);">{escape_html(str(c))}</span>'
            for c in categories
        )
        if n_cats > 5:
            cat_html += f' <span style="color: #888;">...+{n_cats - 5}</span>'
        if validated:
            cat_html += (
                ' <span style="color: #28a745;" title="All values validated">✓</span>'
            )
        else:
            cat_html += f' <span style="color: #fd7e14;" title="{unmapped_count} unmapped values">⚠ {unmapped_count} unmapped</span>'

        tooltip_parts = [f"Registry: {registry}"]
        if ontology_id:
            tooltip_parts.append(f"Ontology: {ontology_id}")
        tooltip_parts.append(f"Validated: {'Yes' if validated else 'No'}")
        if not validated:
            tooltip_parts.append(f"Unmapped: {unmapped_count} values")

        return FormattedOutput(
            type_name=type_name,
            css_class="anndata-dtype--category",
            tooltip="\n".join(tooltip_parts),
            preview_html=cat_html,
            warnings=[]
            if validated
            else [f"{unmapped_count} values not mapped to ontology"],
        )


@case(
    "ext",
    "Ecosystem TypeFormatter for obs/var columns",
    slug="ecosystem-type-formatter",
    tags=("obs", "var", "type-formatter", "categorical"),
    expect=(
        "<code>cell_type</code> → <code>category[bionty.CellType] (4)</code> ✓; <code>tissue</code> → ⚠ 2 unmapped (warning row); <code>assay</code> ✓; var <code>gene_symbol</code> ✓.",
        "<code>batch</code>, <code>n_counts</code>, <code>mean_expression</code> keep default rendering.",
        "Hovering annotated columns shows registry/ontology/validation tooltip.",
    ),
    notes=(
        "Pattern used by packages like <a href='https://lamin.ai/docs/bionty' target='_blank'>bionty</a> / "
        "<a href='https://docs.lamin.ai/' target='_blank'>lamindb</a>: store metadata in uns, register a "
        "TypeFormatter with <code>sections=('obs', 'var')</code> and <code>priority=115</code> (above the "
        "built-in CategoricalFormatter at 110), and look up metadata via <code>context.adata_ref</code> / "
        "<code>context.key</code> / <code>context.section</code>. See also [[uns-type-hints]] and [[treedata-sections]]."
    ),
)
def _ecosystem_formatter() -> CaseOutput:
    adata = AnnData(
        np.random.randn(100, 50).astype(np.float32),
        obs=pd.DataFrame({
            "cell_type": pd.Categorical(
                np.random.choice(["T cell", "B cell", "NK cell", "Monocyte"], 100)
            ),
            "tissue": pd.Categorical(
                np.random.choice(["blood", "spleen", "lymph_node", "bone_marrow"], 100)
            ),
            "assay": pd.Categorical(np.random.choice(["10x 3' v3", "Smart-seq2"], 100)),
            "batch": pd.Categorical(
                np.random.choice(["batch_1", "batch_2", "batch_3"], 100)
            ),
            "n_counts": np.random.randint(1000, 10000, 100),
        }),
        var=pd.DataFrame({
            "gene_symbol": pd.Categorical([f"GENE{i}" for i in range(50)]),
            "mean_expression": np.random.randn(50).astype(np.float32),
        }),
    )
    adata.uns[ONTOLOGY_METADATA_KEY] = {
        "obs": {
            "cell_type": {
                "registry": "bionty.CellType",
                "ontology_id": "cl",
                "validated": True,
                "unmapped_count": 0,
            },
            "tissue": {
                "registry": "bionty.Tissue",
                "ontology_id": "uberon",
                "validated": False,
                "unmapped_count": 2,
            },
            "assay": {
                "registry": "bionty.ExperimentalFactor",
                "ontology_id": "efo",
                "validated": True,
                "unmapped_count": 0,
            },
        },
        "var": {
            "gene_symbol": {
                "registry": "bionty.Gene",
                "ontology_id": "ensembl",
                "validated": True,
                "unmapped_count": 0,
            },
        },
    }
    adata.obsm["X_pca"] = np.random.randn(100, 10).astype(np.float32)
    adata.obsm["X_umap"] = np.random.randn(100, 2).astype(np.float32)
    adata.uns["cell_type_colors"] = palette(4)
    with temporarily_registered(OntologyAnnotatedCategoricalFormatter()):
        return render(adata)


class SpatialExperiment(AnnData):
    """A plain AnnData subclass, as downstream packages often define."""


@case(
    "ext",
    "Plain AnnData subclass",
    slug="anndata-subclass",
    tags=("subclass", "nested-anndata", "view"),
    expect=(
        "Header shows the subclass name <code>SpatialExperiment</code> (not 'AnnData'), once.",
        "All standard sections render as for AnnData; no spurious 'other' section.",
        "A view of the subclass and a subclass instance nested in uns render the same way.",
    ),
)
def _subclass() -> CaseOutput:
    adata = SpatialExperiment(
        np.random.randn(20, 8).astype(np.float32),
        obs=pd.DataFrame({"region": pd.Categorical(["cortex", "hippocampus"] * 10)}),
    )
    adata.obsm["spatial"] = np.random.rand(20, 2)
    adata.uns["region_colors"] = ["#1b9e77", "#d95f02"]
    adata.uns["sub_sample"] = SpatialExperiment(np.zeros((4, 2)))
    return render(adata) + "<hr>" + render(adata[:5])


class SparseCodes:
    """A custom duck array (e.g. from a downstream package) stored in obsm/layers."""

    def __init__(self, shape: tuple[int, int], n_codes: int) -> None:
        self.shape = shape
        self.dtype = np.dtype("uint8")
        self.ndim = 2
        self.n_codes = n_codes

    def __getitem__(self, idx: object) -> SparseCodes:
        return self


class SparseCodesFormatter(TypeFormatter):
    """TypeFormatter teaching the repr about :class:`SparseCodes`."""

    def can_format(self, obj, context):
        return isinstance(obj, SparseCodes)

    def format(self, obj, context) -> FormattedOutput:
        n, m = obj.shape
        return FormattedOutput(
            type_name=f"SparseCodes ({n} × {m}) · {obj.n_codes} codes",
            css_class="anndata-dtype--extension",
            tooltip="Custom array type from a downstream package",
            preview=f"codebook of {obj.n_codes}",
        )


@case(
    "ext",
    "Custom array type with TypeFormatter",
    slug="custom-array-type",
    tags=("obsm", "type-formatter"),
    expect=(
        "<code>obsm['codes']</code> shows 'SparseCodes (40 × 16) · 256 codes' with the extension color and a preview.",
        "Without the formatter (second repr) it falls back to a generic type entry, flagged as unknown/non-serializable rather than crashing.",
    ),
)
def _custom_array_type() -> CaseOutput:
    adata = AnnData(np.zeros((40, 6)))
    adata.obsm._data["codes"] = SparseCodes((40, 16), n_codes=256)  # type: ignore[union-attr]
    with temporarily_registered(SparseCodesFormatter()):
        with_formatter = render(adata)
    return with_formatter + "<hr><p><b>Without formatter:</b></p>" + render(adata)


# =============================================================================
# Page rendering
# =============================================================================

_PAGE_CSS = """
:root { --vt-bg:#f5f5f5; --vt-card:#fff; --vt-text:#222; --vt-muted:#666; --vt-border:#e3e3e3;
        --vt-accent:#0d6efd; --vt-desc:#f8f9fa; --vt-ok:#1a7f37; --vt-skip:#9a6700; --vt-fail:#cf222e; }
body.dark-mode { --vt-bg:#1a1a1a; --vt-card:#2a2a2a; --vt-text:#e0e0e0; --vt-muted:#aaa;
        --vt-border:#3a3a3a; --vt-accent:#6ea8fe; --vt-desc:#333; }
body { font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif; margin:0 20px 0 300px;
       padding:20px; background:var(--vt-bg); color:var(--vt-text); }
h1 { border-bottom:2px solid var(--vt-accent); padding-bottom:10px; margin-top:0; }
h2.vt-cat { margin:48px 0 4px; }
.vt-cat-blurb { color:var(--vt-muted); margin:0 0 12px; }
a { color:var(--vt-accent); }
.vt-toc { position:fixed; top:0; left:0; bottom:0; width:270px; overflow-y:auto; background:var(--vt-card);
          border-right:1px solid var(--vt-border); padding:14px; box-sizing:border-box; font-size:12px; }
.vt-toc h4 { margin:12px 0 4px; font-size:12px; text-transform:uppercase; letter-spacing:.04em; }
.vt-toc a { display:block; padding:2px 0; color:var(--vt-text); text-decoration:none; }
.vt-toc a:hover { color:var(--vt-accent); }
.vt-dot { display:inline-block; width:7px; height:7px; border-radius:50%; margin-right:5px; }
.vt-dot.ok { background:var(--vt-ok); } .vt-dot.skipped { background:var(--vt-skip); }
.vt-dot.failed { background:var(--vt-fail); }
.vt-toolbar { position:fixed; top:12px; right:20px; z-index:10; }
.vt-toolbar button { padding:8px 14px; border:0; border-radius:4px; background:#333; color:#fff; cursor:pointer; }
.vt-summary { background:var(--vt-card); border-radius:8px; padding:12px 16px; font-size:13px; }
.vt-case { background:var(--vt-card); border-radius:8px; padding:18px 20px; margin:16px 0;
           box-shadow:0 1px 3px rgba(0,0,0,.12); scroll-margin-top:10px; }
.vt-case h3 { margin:0 0 6px; }
.vt-case h3 .vt-num { color:var(--vt-accent); margin-right:6px; }
.vt-case h3 a.vt-anchor { color:var(--vt-muted); text-decoration:none; font-weight:normal; font-size:12px; margin-left:6px; }
.vt-meta { font-size:11px; color:var(--vt-muted); margin-bottom:8px; }
.vt-tag { display:inline-block; padding:0 6px; margin:0 3px 3px 0; border:1px solid var(--vt-border); border-radius:9px; }
.vt-expect { background:var(--vt-desc); border-left:3px solid var(--vt-accent); border-radius:4px;
             padding:8px 12px; margin-bottom:10px; font-size:13px; }
.vt-expect ul { margin:4px 0 0; padding-left:20px; }
.vt-notes, .vt-warnings { font-size:12px; color:var(--vt-muted); margin-bottom:10px; }
.vt-note { background:#fff3cd; color:#664d03; border-radius:4px; padding:6px 10px; font-size:12px; margin-bottom:10px; }
.vt-status { border-radius:4px; padding:8px 12px; font-size:13px; }
.vt-status.skipped { background:#fff8c5; color:#4d2d00; }
.vt-status.failed { background:#ffebe9; color:#82071e; }
.vt-status pre { white-space:pre-wrap; font-size:11px; margin:6px 0 0; }
.vt-panes { display:grid; grid-template-columns:repeat(auto-fit,minmax(380px,1fr)); gap:12px; }
.vt-pane { margin:0; }
.vt-pane figcaption { font-size:12px; color:var(--vt-muted); margin-bottom:4px; }
iframe.vt-frame { width:100%; min-height:120px; border:1px solid var(--vt-border); border-radius:4px; display:block; }
.vt-coverage { columns:3 260px; font-size:12px; }
.vt-coverage div { break-inside:avoid; margin-bottom:2px; }
@media (max-width: 900px) { body { margin-left:0; } .vt-toc { position:static; width:auto; border:0; } }
"""

_PAGE_JS = """
function vtFit(f) {
  try {
    const d = f.contentDocument;
    if (!d || !d.body) return;
    const fit = () => { f.style.height = (d.documentElement.scrollHeight + 2) + "px"; };
    fit();
    new ResizeObserver(fit).observe(d.body);
  } catch (e) { /* cross-origin: keep min-height */ }
}
document.querySelectorAll("iframe.vt-frame").forEach((f) => {
  f.addEventListener("load", () => vtFit(f));
  if (f.contentDocument && f.contentDocument.readyState === "complete") vtFit(f);
});
// Collect the theme panes' self-checks into the summary (read by --browser-check)
window.addEventListener("load", () => setTimeout(() => {
  const el = document.getElementById("vt-probe-summary");
  const rows = [...document.querySelectorAll(".vt-pane")].map((fg) => {
    let probe = "unavailable";
    try { probe = fg.querySelector("iframe").contentDocument.getElementById("vt-probe").textContent; } catch (e) {}
    return { slug: fg.closest(".vt-case").id, pane: fg.querySelector("figcaption").textContent, probe };
  });
  if (!el || !rows.length) return;
  const fails = rows.filter((r) => !r.probe.endsWith(" PASS"));
  el.dataset.results = JSON.stringify(rows);
  el.textContent = `Theme probes: ${rows.length - fails.length}/${rows.length} PASS`;
  const ul = document.createElement("ul");
  for (const r of fails) {
    const li = document.createElement("li");
    const a = document.createElement("a");
    a.href = "#" + r.slug;
    a.textContent = r.pane;
    li.append(a, ": " + r.probe);
    ul.append(li);
  }
  if (fails.length) el.append(ul);
}, 1000));
"""


def _link_refs(text: str, by_slug: dict[str, Case]) -> str:
    """Replace ``[[slug]]`` with a link to that case."""

    def repl(m: re.Match[str]) -> str:
        c = by_slug.get(m.group(1))
        if c is None:
            return f"<code>[[{m.group(1)}]] (unknown case)</code>"
        return f'<a href="#{c.slug}">{c.number} {escape_html(c.title)}</a>'

    return re.sub(r"\[\[([a-z0-9-]+)\]\]", repl, text)


def _render_case(r: CaseResult, by_slug: dict[str, Case]) -> str:
    c = r.case
    expect = "".join(f"<li>{_link_refs(e, by_slug)}</li>" for e in c.expect)
    tags = "".join(f'<span class="vt-tag">{escape_html(t)}</span>' for t in c.tags)
    timing = f" · {r.seconds:.2f}s" if r.status != "skipped" else ""
    size = f" · {len(r.html) / 1024:.0f} KB HTML" if r.html else ""
    parts = [
        f'<section id="{c.slug}" class="vt-case">',
        (
            f'<h3><span class="vt-num">{c.number}</span>{escape_html(c.title)}'
            f'<a class="vt-anchor" href="#{c.slug}">#{c.slug}</a></h3>'
        ),
        f'<div class="vt-meta">{tags}{timing}{size}</div>',
        f'<div class="vt-expect"><b>What to check</b><ul>{expect}</ul></div>',
    ]
    if c.notes:
        parts.append(
            f'<details class="vt-notes"><summary>Background</summary>'
            f"{_link_refs(c.notes, by_slug)}</details>"
        )
    if r.note:
        parts.append(f'<div class="vt-note">{r.note}</div>')
    if r.warnings:
        items = "".join(f"<li>{escape_html(w)}</li>" for w in r.warnings)
        parts.append(
            f'<details class="vt-warnings"><summary>{len(r.warnings)} Python warning(s) '
            f"emitted while building/rendering</summary><ul>{items}</ul></details>"
        )
    if r.status == "ok":
        parts.append(f'<div class="vt-output">{r.html}</div>')
    elif r.status == "skipped":
        parts.append(f'<div class="vt-status skipped">{escape_html(r.reason)}</div>')
    else:
        parts.append(
            '<div class="vt-status failed"><b>Case raised an exception</b>'
            f"<pre>{escape_html(r.reason)}</pre></div>"
        )
    parts.append("</section>")
    return "\n".join(parts)


def render_page(results: Sequence[CaseResult]) -> str:
    """Assemble the full HTML page: TOC, summary, coverage index, cases."""
    by_slug = {r.case.slug: r.case for r in results}
    status_of = {r.case.slug: r.status for r in results}

    toc = []
    body = []
    for cat in CATEGORIES:
        cat_results = [r for r in results if r.case.category == cat.key]
        if not cat_results:
            continue
        n = _CATEGORY_INDEX[cat.key]
        toc.append(
            f'<h4><a href="#cat-{cat.key}">{n}. {escape_html(cat.title)}</a></h4>'
        )
        toc.extend(
            f'<a href="#{r.case.slug}"><span class="vt-dot {r.status}"></span>'
            f"{r.case.number} {escape_html(r.case.title)}</a>"
            for r in cat_results
        )
        body.append(
            f'<h2 id="cat-{cat.key}" class="vt-cat">{n}. {escape_html(cat.title)}</h2>'
            f'<p class="vt-cat-blurb">{escape_html(cat.blurb)}</p>'
        )
        body.extend(_render_case(r, by_slug) for r in cat_results)

    counts = {
        s: sum(r.status == s for r in results) for s in ("ok", "skipped", "failed")
    }
    problems = [r for r in results if r.status != "ok"]
    problem_links = ", ".join(
        f'<a href="#{r.case.slug}">{r.case.number} {escape_html(r.case.title)}</a> ({r.status})'
        for r in problems
    )

    tag_index: dict[str, list[Case]] = {}
    for r in results:
        for t in r.case.tags:
            tag_index.setdefault(t, []).append(r.case)
    coverage = "".join(
        f"<div><b>{escape_html(t)}</b>: "
        + ", ".join(
            f'<a href="#{c.slug}" title="{escape_html(c.title)}">'
            f'<span class="vt-dot {status_of[c.slug]}"></span>{c.number}</a>'
            for c in cases
        )
        + "</div>"
        for t, cases in sorted(tag_index.items())
    )

    generated = datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<meta http-equiv="Content-Security-Policy" content="default-src 'self' 'unsafe-inline' 'unsafe-eval' data: https:; style-src 'self' 'unsafe-inline';">
<title>AnnData repr visual test</title>
<style>{_PAGE_CSS}</style>
</head>
<body>
<nav class="vt-toc"><b>AnnData _repr_html_</b>{"".join(toc)}</nav>
<div class="vt-toolbar"><button onclick="document.body.classList.toggle('dark-mode')">Toggle dark mode</button></div>
<h1>AnnData <code>_repr_html_</code> visual test</h1>
<div class="vt-summary">
<p style="margin-top:0">anndata {escape_html(version("anndata"))} · Python {platform.python_version()} ·
pandas {pd.__version__} · numpy {np.__version__} · generated {generated}</p>
<p><span class="vt-dot ok"></span>{counts["ok"]} rendered ·
<span class="vt-dot skipped"></span>{counts["skipped"]} skipped ·
<span class="vt-dot failed"></span>{counts["failed"]} failed
{f"<br>Not rendered: {problem_links}" if problems else ""}</p>
<p>Each case lists <b>what to check</b>. The dark-mode button adds <code>body.dark-mode</code>
(which the repr honors); embedded host themes are simulated in category
<a href="#cat-env">5</a>. Regenerate with
<code>python tests/visual_inspect_repr_html.py</code> (<code>--help</code> for filters).</p>
<p id="vt-probe-summary"></p>
<details><summary>Coverage index (tag → cases)</summary><div class="vt-coverage">{coverage}</div></details>
</div>
{"".join(body)}
<script>{_PAGE_JS}</script>
</body>
</html>
"""


# =============================================================================
# CLI
# =============================================================================


def find_chrome() -> str | None:
    """Locate a Chrome/Chromium binary (``$CHROME`` overrides)."""
    import os

    candidates = [
        os.environ.get("CHROME"),
        *(
            shutil.which(n)
            for n in (
                "google-chrome",
                "google-chrome-stable",
                "chromium",
                "chromium-browser",
                "chrome",
            )
        ),
        "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
        "/Applications/Chromium.app/Contents/MacOS/Chromium",
    ]
    return next((c for c in candidates if c and Path(c).exists()), None)


def browser_check(page: Path) -> int | None:
    """Load ``page`` in headless Chrome and report the theme panes' self-checks.

    Returns the number of failing panes, or ``None`` if the check could not run.
    """
    import json
    import subprocess

    chrome = find_chrome()
    if chrome is None:
        print(
            "--browser-check: no Chrome/Chromium found (set $CHROME).", file=sys.stderr
        )
        return None
    proc = subprocess.run(
        [
            chrome,
            "--headless=new",
            "--disable-gpu",
            "--virtual-time-budget=10000",
            "--dump-dom",
            page.resolve().as_uri(),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=300,
    )
    m = re.search(r'id="vt-probe-summary"[^>]*data-results="([^"]*)"', proc.stdout)
    if m is None:
        print(
            "--browser-check: no probe results in the rendered page.", file=sys.stderr
        )
        return None
    rows = json.loads(html_mod.unescape(m.group(1)))
    n_fail = 0
    for row in rows:
        ok = row["probe"].endswith(" PASS")
        n_fail += not ok
        print(
            f"  {'PASS' if ok else 'FAIL'}  #{row['slug']}  {row['pane']}: {row['probe']}"
        )
    print(f"Theme probes: {len(rows) - n_fail}/{len(rows)} PASS")
    return n_fail


def _matches(c: Case, patterns: Sequence[str]) -> bool:
    haystack = f"{c.number} {c.slug} {c.title} {c.category}".lower()
    return any(
        c.number == p or c.number.startswith(p if p.endswith(".") else f"{p}.")
        if re.fullmatch(r"\d+(\.\d*)?", p)
        else p.lower() in haystack
        for p in patterns
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Generate the visual test HTML file."""
    parser = argparse.ArgumentParser(
        description="Render AnnData _repr_html_ scenarios into one HTML page.",
    )
    parser.add_argument(
        "--only",
        action="append",
        default=[],
        metavar="PATTERN",
        help=(
            "Only run matching cases (repeatable). A number like '3' or '3.2' selects "
            "by category/case number; anything else is a case-insensitive substring "
            "of number, slug, title or category key."
        ),
    )
    parser.add_argument(
        "--list", action="store_true", help="List cases and exit without rendering."
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path(__file__).parent / "repr_html_visual_test.html",
        help="Output HTML path (default: %(default)s).",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit with status 1 if any case raised an exception (or a theme probe failed).",
    )
    parser.add_argument(
        "--browser-check",
        action="store_true",
        help="Afterwards, load the page in headless Chrome and print the theme probes.",
    )
    args = parser.parse_args(argv)

    cases = numbered_cases()
    if args.only:
        cases = [c for c in cases if _matches(c, args.only)]
        if not cases:
            print(f"No cases match {args.only}", file=sys.stderr)
            return 2

    if args.list:
        for c in cases:
            req = f"  [requires {', '.join(c.requires)}]" if c.requires else ""
            print(f"{c.number:>5}  {c.slug:<32} {c.title}{req}")
        return 0

    print(f"Rendering {len(cases)} visual test cases...")
    results = []
    for c in cases:
        r = run_case(c)
        detail = (
            f"{r.seconds:.2f}s"
            if r.status == "ok"
            else r.reason.strip().splitlines()[-1]
        )
        print(f"  {c.number:>5} {c.title} ... {r.status} ({detail})")
        results.append(r)

    args.output.write_text(render_page(results), encoding="utf-8")
    n_failed = sum(r.status == "failed" for r in results)
    print(f"\nVisual test file generated: {args.output}")
    if n_failed:
        print(f"{n_failed} case(s) failed; see the red boxes in the page.")
    if args.browser_check:
        n_failed += browser_check(args.output) or 0
    return 1 if args.strict and n_failed else 0


if __name__ == "__main__":
    sys.exit(main())
