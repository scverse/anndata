"""
CSS styles for AnnData HTML representation.

Provides inline CSS with:
- Light and dark mode support
- Jupyter notebook theme detection
- Responsive design
- Accessible color contrasts
"""

from __future__ import annotations


def get_css() -> str:
    """Get the complete CSS for the HTML representation."""
    return f"""<style>
{_CSS_CONTENT}
</style>"""


_CSS_CONTENT = """
/* AnnData HTML Representation Styles */
/* Scoped to .anndata-repr to avoid conflicts */

.anndata-repr {
    /* CSS Variables - Light mode (default) */
    --anndata-bg-primary: #ffffff;
    --anndata-bg-secondary: #f8f9fa;
    --anndata-bg-tertiary: #e9ecef;
    --anndata-highlight: #e7f1ff;
    --anndata-text-primary: #212529;
    --anndata-text-secondary: #6c757d;
    --anndata-text-muted: #adb5bd;
    --anndata-border-color: #dee2e6;
    --anndata-border-light: #e9ecef;
    --anndata-accent-color: #0d6efd;
    --anndata-accent-hover: #0b5ed7;
    --anndata-warning-color: #ffc107;
    --anndata-warning-bg: #fff3cd;
    --anndata-error-color: #dc3545;
    --anndata-error-bg: #f8d7da;
    --anndata-success-color: #198754;
    --anndata-info-color: #0dcaf0;
    --anndata-link-color: #0d6efd;
    --anndata-code-bg: #f8f9fa;
    --anndata-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
    --anndata-radius: 4px;
    --anndata-font-mono: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace;
    --anndata-font-size: 13px;
    --anndata-line-height: 1.4;

    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    font-size: var(--anndata-font-size);
    line-height: var(--anndata-line-height);
    color: var(--anndata-text-primary);
    background: var(--anndata-bg-primary);
    border: 1px solid var(--anndata-border-color);
    border-radius: var(--anndata-radius);
    padding: 0;
    margin: 8px 0;
    max-width: 100%;
    overflow: hidden;
}

/* Dark mode - via media query */
@media (prefers-color-scheme: dark) {
    .anndata-repr {
        --anndata-bg-primary: #1e1e1e;
        --anndata-bg-secondary: #252526;
        --anndata-bg-tertiary: #2d2d2d;
        --anndata-highlight: #264f78;
        --anndata-text-primary: #e0e0e0;
        --anndata-text-secondary: #a0a0a0;
        --anndata-text-muted: #707070;
        --anndata-border-color: #404040;
        --anndata-border-light: #333333;
        --anndata-accent-color: #58a6ff;
        --anndata-accent-hover: #79b8ff;
        --anndata-warning-color: #d29922;
        --anndata-warning-bg: #3d3200;
        --anndata-error-color: #f85149;
        --anndata-error-bg: #3d1a1a;
        --anndata-success-color: #3fb950;
        --anndata-info-color: #58a6ff;
        --anndata-link-color: #58a6ff;
        --anndata-code-bg: #2d2d2d;
        --anndata-shadow: 0 1px 3px rgba(0, 0, 0, 0.3);
    }
}

/* Jupyter dark theme detection + manual dark mode class */
[data-jp-theme-light="false"] .anndata-repr,
.jp-Theme-Dark .anndata-repr,
body.vscode-dark .anndata-repr,
body[data-vscode-theme-kind="vscode-dark"] .anndata-repr,
body.dark-mode .anndata-repr {
    --anndata-bg-primary: #1e1e1e;
    --anndata-bg-secondary: #252526;
    --anndata-bg-tertiary: #2d2d2d;
    --anndata-highlight: #264f78;
    --anndata-text-primary: #e0e0e0;
    --anndata-text-secondary: #a0a0a0;
    --anndata-text-muted: #707070;
    --anndata-border-color: #404040;
    --anndata-border-light: #333333;
    --anndata-accent-color: #58a6ff;
    --anndata-accent-hover: #79b8ff;
    --anndata-warning-color: #d29922;
    --anndata-warning-bg: #3d3200;
    --anndata-error-color: #f85149;
    --anndata-error-bg: #3d1a1a;
    --anndata-success-color: #3fb950;
    --anndata-info-color: #58a6ff;
    --anndata-link-color: #58a6ff;
    --anndata-code-bg: #2d2d2d;
    --anndata-shadow: 0 1px 3px rgba(0, 0, 0, 0.3);
}

/* Header */
.anndata-repr .anndata-hdr {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 8px;
    padding: 10px 12px;
    background: #f8f9fa; /* Fallback */
    background: var(--anndata-bg-secondary);
    border-bottom: 1px solid #dee2e6; /* Fallback */
    border-bottom: 1px solid var(--anndata-border-color);
}

.anndata-repr .adata-type {
    font-weight: 600;
    font-size: 14px;
    color: #212529; /* Fallback */
    color: var(--anndata-text-primary);
}

.anndata-repr .adata-shape {
    font-family: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace; /* Fallback */
    font-family: var(--anndata-font-mono);
    font-size: 12px;
    color: #6c757d; /* Fallback */
    color: var(--anndata-text-secondary);
}

.anndata-repr .adata-badge {
    display: inline-flex;
    align-items: center;
    gap: 4px;
    padding: 2px 8px;
    font-size: 11px;
    font-weight: 500;
    border-radius: 10px;
    white-space: nowrap;
}

.anndata-repr .adata-badge-view {
    background: var(--anndata-info-color);
    color: white;
}

.anndata-repr .adata-badge-backed {
    background: var(--anndata-success-color);
    color: white;
}

.anndata-repr .adata-badge-extension {
    background: var(--anndata-accent-color);
    color: white;
}

/* Search box */
.anndata-repr .adata-search {
    padding: 8px 12px;
    background: var(--anndata-bg-primary);
    border-bottom: 1px solid var(--anndata-border-light);
}

.anndata-repr .adata-search-input {
    width: 100%;
    max-width: 300px;
    padding: 6px 10px;
    font-size: 12px;
    border: 1px solid var(--anndata-border-color);
    border-radius: var(--anndata-radius);
    background: var(--anndata-bg-primary);
    color: var(--anndata-text-primary);
    outline: none;
    transition: border-color 0.15s;
}

.anndata-repr .adata-search-input:focus {
    border-color: var(--anndata-accent-color);
}

.anndata-repr .adata-search-input::placeholder {
    color: var(--anndata-text-muted);
}

.anndata-repr .adata-filter-indicator {
    display: none;
    margin-left: 8px;
    font-size: 11px;
    color: var(--anndata-accent-color);
}

.anndata-repr .adata-filter-indicator.active {
    display: inline;
}

/* Metadata bar */
.anndata-repr .adata-metadata {
    display: flex;
    flex-wrap: wrap;
    gap: 12px;
    padding: 6px 12px;
    font-size: 11px;
    color: var(--anndata-text-secondary);
    background: var(--anndata-bg-tertiary);
    border-bottom: 1px solid var(--anndata-border-light);
}

.anndata-repr .adata-metadata span {
    white-space: nowrap;
}

/* Index preview */
.anndata-repr .adata-index-preview {
    padding: 8px 12px;
    font-size: 11px;
    font-family: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace; /* Fallback */
    font-family: var(--anndata-font-mono);
    color: #6c757d; /* Fallback */
    color: var(--anndata-text-secondary);
    background: #ffffff; /* Fallback */
    background: var(--anndata-bg-primary);
    border-bottom: 1px solid #e9ecef; /* Fallback */
    border-bottom: 1px solid var(--anndata-border-light);
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
}

.anndata-repr .adata-index-preview strong {
    color: #212529; /* Fallback */
    color: var(--anndata-text-primary);
    font-weight: 500;
}

/* Sections container */
.anndata-repr .anndata-secs {
    padding: 0;
}

/* Individual section */
.anndata-repr .anndata-sec {
    border-bottom: 1px solid #e9ecef; /* Fallback */
    border-bottom: 1px solid var(--anndata-border-light);
}

.anndata-repr .anndata-sec:last-child {
    border-bottom: none;
}

.anndata-repr .anndata-sechdr {
    display: flex;
    align-items: center;
    gap: 8px;
    padding: 8px 12px;
    cursor: pointer;
    user-select: none;
    background: #ffffff; /* Fallback */
    background: var(--anndata-bg-primary);
    transition: background-color 0.15s;
}

.anndata-repr .anndata-sechdr:hover {
    background: #f8f9fa; /* Fallback */
    background: var(--anndata-bg-secondary);
}

.anndata-repr .adata-fold-icon {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 16px;
    height: 16px;
    font-size: 10px;
    color: #adb5bd; /* Fallback */
    color: var(--anndata-text-muted);
    transition: transform 0.15s;
    transform-origin: center;
    flex-shrink: 0;
}

.anndata-repr .anndata-sec.collapsed .adata-fold-icon {
    transform: rotate(-90deg);
}

.anndata-repr .anndata-sec-name {
    font-weight: 600;
    color: #212529; /* Fallback */
    color: var(--anndata-text-primary);
}

.anndata-repr .anndata-sec-count {
    font-size: 11px;
    color: #6c757d; /* Fallback */
    color: var(--anndata-text-secondary);
}

.anndata-repr .adata-help-link {
    margin-left: auto;
    padding: 2px 6px;
    font-size: 11px;
    color: #adb5bd; /* Fallback */
    color: var(--anndata-text-muted);
    text-decoration: none;
    border-radius: 4px; /* Fallback */
    border-radius: var(--anndata-radius);
    transition: color 0.15s, background-color 0.15s;
}

.anndata-repr .adata-help-link:hover {
    color: #0d6efd; /* Fallback */
    color: var(--anndata-accent-color);
    background: #e9ecef; /* Fallback */
    background: var(--anndata-bg-tertiary);
}

/* Section content */
.anndata-repr .anndata-seccontent {
    padding: 0;
    overflow: hidden;
    transition: max-height 0.2s ease-out;
}

.anndata-repr .anndata-sec.collapsed .anndata-seccontent {
    max-height: 0 !important;
    padding: 0;
}

/* Entries table */
.anndata-repr .adata-table {
    width: 100%;
    border-collapse: collapse;
    font-size: 12px;
}

.anndata-repr .adata-entry {
    border-bottom: 1px solid #e9ecef; /* Fallback */
    border-bottom: 1px solid var(--anndata-border-light);
    transition: background-color 0.1s;
}

.anndata-repr .adata-entry:last-child {
    border-bottom: none;
}

.anndata-repr .adata-entry:hover {
    background: #f8f9fa; /* Fallback */
    background: var(--anndata-bg-secondary);
}

.anndata-repr .adata-entry.hidden {
    display: none;
}

.anndata-repr .adata-entry.warning {
    background: var(--anndata-warning-bg);
}

.anndata-repr .adata-entry.error {
    background: var(--anndata-error-bg);
}

.anndata-repr .adata-entry td {
    padding: 6px 12px;
    vertical-align: middle;
}

.anndata-repr .adata-entry-name {
    font-family: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace; /* Fallback */
    font-family: var(--anndata-font-mono);
    font-weight: 500;
    color: #212529; /* Fallback */
    color: var(--anndata-text-primary);
    white-space: nowrap;
}

.anndata-repr .adata-entry-type {
    font-family: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace; /* Fallback */
    font-family: var(--anndata-font-mono);
    font-size: 11px;
    color: #6c757d; /* Fallback */
    color: var(--anndata-text-secondary);
}

.anndata-repr .adata-entry-meta {
    font-size: 11px;
    color: #adb5bd; /* Fallback */
    color: var(--anndata-text-muted);
    text-align: right;
}

/* Copy button */
.anndata-repr .adata-copy-btn {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 20px;
    height: 20px;
    margin-left: 4px;
    padding: 0;
    font-size: 11px;
    color: var(--anndata-text-muted);
    background: transparent;
    border: none;
    border-radius: 3px;
    cursor: pointer;
    opacity: 0;
    transition: opacity 0.15s, color 0.15s, background-color 0.15s;
}

.anndata-repr .adata-entry:hover .adata-copy-btn {
    opacity: 1;
}

.anndata-repr .adata-copy-btn:hover {
    color: var(--anndata-accent-color);
    background: var(--anndata-bg-tertiary);
}

.anndata-repr .adata-copy-btn.copied {
    color: var(--anndata-success-color);
}

/* Type-specific styling */
.anndata-repr .dtype-category { color: #8250df; }
.anndata-repr .dtype-int { color: #0550ae; }
.anndata-repr .dtype-float { color: #0550ae; }
.anndata-repr .dtype-bool { color: #cf222e; }
.anndata-repr .dtype-string { color: #0a3069; }
.anndata-repr .dtype-object { color: #6e7781; }
.anndata-repr .dtype-sparse { color: #1a7f37; }
.anndata-repr .dtype-array { color: #0550ae; }
.anndata-repr .dtype-dataframe { color: #8250df; }
.anndata-repr .dtype-anndata { color: #cf222e; font-weight: 600; }
.anndata-repr .dtype-unknown { color: #6e7781; font-style: italic; }
.anndata-repr .dtype-extension { color: #8250df; }
.anndata-repr .dtype-warning { color: var(--anndata-warning-color); }
.anndata-repr .dtype-dask { color: #fb8500; }
.anndata-repr .dtype-gpu { color: #76b900; }
.anndata-repr .dtype-awkward { color: #e85d04; }
.anndata-repr .dtype-array-api { color: #9a6700; } /* Array-API arrays (JAX, PyTorch, etc.) - PR #2063 */

/* Dark mode type colors */
@media (prefers-color-scheme: dark) {
    .anndata-repr .dtype-category { color: #d2a8ff; }
    .anndata-repr .dtype-int { color: #79c0ff; }
    .anndata-repr .dtype-float { color: #79c0ff; }
    .anndata-repr .dtype-bool { color: #ff7b72; }
    .anndata-repr .dtype-string { color: #a5d6ff; }
    .anndata-repr .dtype-sparse { color: #7ee787; }
    .anndata-repr .dtype-array { color: #79c0ff; }
    .anndata-repr .dtype-dataframe { color: #d2a8ff; }
    .anndata-repr .dtype-anndata { color: #ff7b72; }
    .anndata-repr .dtype-dask { color: #ffc168; }
    .anndata-repr .dtype-gpu { color: #a0db63; }
    .anndata-repr .dtype-awkward { color: #ff9d76; }
    .anndata-repr .dtype-array-api { color: #e6c400; } /* Array-API arrays - dark mode */
}

[data-jp-theme-light="false"] .anndata-repr .dtype-category,
.jp-Theme-Dark .anndata-repr .dtype-category,
body.dark-mode .anndata-repr .dtype-category { color: #d2a8ff; }
[data-jp-theme-light="false"] .anndata-repr .dtype-int,
.jp-Theme-Dark .anndata-repr .dtype-int,
body.dark-mode .anndata-repr .dtype-int { color: #79c0ff; }
[data-jp-theme-light="false"] .anndata-repr .dtype-float,
.jp-Theme-Dark .anndata-repr .dtype-float,
body.dark-mode .anndata-repr .dtype-float { color: #79c0ff; }
[data-jp-theme-light="false"] .anndata-repr .dtype-bool,
.jp-Theme-Dark .anndata-repr .dtype-bool,
body.dark-mode .anndata-repr .dtype-bool { color: #ff7b72; }
[data-jp-theme-light="false"] .anndata-repr .dtype-string,
.jp-Theme-Dark .anndata-repr .dtype-string,
body.dark-mode .anndata-repr .dtype-string { color: #a5d6ff; }
[data-jp-theme-light="false"] .anndata-repr .dtype-sparse,
.jp-Theme-Dark .anndata-repr .dtype-sparse,
body.dark-mode .anndata-repr .dtype-sparse { color: #7ee787; }
[data-jp-theme-light="false"] .anndata-repr .dtype-array,
.jp-Theme-Dark .anndata-repr .dtype-array,
body.dark-mode .anndata-repr .dtype-array { color: #79c0ff; }
[data-jp-theme-light="false"] .anndata-repr .dtype-dataframe,
.jp-Theme-Dark .anndata-repr .dtype-dataframe,
body.dark-mode .anndata-repr .dtype-dataframe { color: #d2a8ff; }
[data-jp-theme-light="false"] .anndata-repr .dtype-anndata,
.jp-Theme-Dark .anndata-repr .dtype-anndata,
body.dark-mode .anndata-repr .dtype-anndata { color: #ff7b72; }

/* Color swatches */
.anndata-repr .adata-color-swatches {
    display: inline-flex;
    gap: 2px;
    margin-left: 6px;
    vertical-align: middle;
}

.anndata-repr .adata-color-swatch {
    display: inline-block;
    width: 12px;
    height: 12px;
    border-radius: 2px;
    border: 1px solid var(--anndata-border-color);
}

/* Warning indicator */
.anndata-repr .adata-warning-icon {
    color: var(--anndata-warning-color);
    margin-left: 4px;
}

/* Expandable nested content */
.anndata-repr .adata-expand-btn {
    padding: 2px 8px;
    font-size: 11px;
    color: var(--anndata-accent-color);
    background: transparent;
    border: 1px solid var(--anndata-accent-color);
    border-radius: var(--anndata-radius);
    cursor: pointer;
    transition: background-color 0.15s, color 0.15s;
}

.anndata-repr .adata-expand-btn:hover {
    background: var(--anndata-accent-color);
    color: white;
}

/* Nested row (contains nested AnnData) */
.anndata-repr .adata-nested-row {
    display: none;
}

.anndata-repr .adata-nested-row.expanded {
    display: table-row;
}

.anndata-repr .adata-nested-content {
    padding: 8px 12px 8px 24px;
    background: var(--anndata-bg-secondary);
    border-top: 1px solid var(--anndata-border-light);
    /* Constrain width to prevent wide tables from expanding the layout */
    max-width: 1px;
    overflow: hidden;
}

/* Container for embedded content (e.g., pandas DataFrame HTML output) */
.anndata-repr .adata-custom-expanded {
    display: block;
    overflow-x: auto;
    overflow-y: hidden;
    /* Use width:0 + min-width:100% trick to prevent table from expanding container */
    width: 0;
    min-width: 100%;
    padding: 8px 12px;
    -webkit-overflow-scrolling: touch;
}

/* Jupyter-like table styling for embedded DataFrames */
.anndata-repr .adata-custom-expanded > div {
    display: inline-block;
}

.anndata-repr .adata-custom-expanded table {
    border-collapse: collapse;
    border-spacing: 0;
    border: none;
    font-size: 12px;
    font-family: var(--anndata-font-mono);
    table-layout: auto;
    margin: 0;
    white-space: nowrap;
}

.anndata-repr .adata-custom-expanded thead {
    border-bottom: 1px solid var(--anndata-border-color);
    vertical-align: bottom;
}

.anndata-repr .adata-custom-expanded th,
.anndata-repr .adata-custom-expanded td {
    vertical-align: middle;
    padding: 6px 10px;
    line-height: normal;
    border: none;
    text-align: right;
}

.anndata-repr .adata-custom-expanded th {
    font-weight: 600;
    color: var(--anndata-text-primary);
    background: var(--anndata-bg-secondary);
}

.anndata-repr .adata-custom-expanded td {
    color: var(--anndata-text-primary);
}

.anndata-repr .adata-custom-expanded tbody tr:nth-child(odd) {
    background: var(--anndata-bg-primary);
}

.anndata-repr .adata-custom-expanded tbody tr:nth-child(even) {
    background: var(--anndata-bg-secondary);
}

.anndata-repr .adata-custom-expanded tbody tr:hover {
    background: var(--anndata-highlight);
}

/* Nested AnnData */
.anndata-repr .adata-nested-anndata {
    margin: 8px 0;
    border: 1px solid var(--anndata-border-color);
    border-radius: var(--anndata-radius);
}

/* X section (special) */
.anndata-repr .adata-x-section {
    padding: 8px 12px;
    background: var(--anndata-bg-secondary);
    border-bottom: 1px solid var(--anndata-border-light);
}

.anndata-repr .adata-x-info {
    display: grid;
    grid-template-columns: auto 1fr;
    gap: 4px 12px;
    font-size: 12px;
}

.anndata-repr .adata-x-info dt {
    color: var(--anndata-text-secondary);
    font-weight: 500;
}

.anndata-repr .adata-x-info dd {
    margin: 0;
    font-family: var(--anndata-font-mono);
    color: var(--anndata-text-primary);
}

/* Truncation indicator */
.anndata-repr .adata-truncated {
    padding: 8px 12px;
    font-size: 11px;
    color: var(--anndata-text-muted);
    text-align: center;
    font-style: italic;
}

/* Max depth indicator */
.anndata-repr .adata-max-depth {
    padding: 8px 12px;
    font-size: 11px;
    color: var(--anndata-text-muted);
    background: var(--anndata-bg-tertiary);
    border-radius: var(--anndata-radius);
    text-align: center;
}

/* Empty section indicator */
.anndata-repr .adata-empty {
    padding: 8px 12px;
    font-size: 11px;
    color: var(--anndata-text-muted);
    font-style: italic;
}

/* Footer */
.anndata-repr .anndata-ftr {
    color: var(--anndata-text-muted);
    border-top: 1px solid var(--anndata-border-light);
}

/* Muted text helper */
.anndata-repr .adata-text-muted {
    color: var(--anndata-text-muted);
}

/* X entry row */
.anndata-repr .adata-x-entry {
    border-bottom: 1px solid #e9ecef; /* Fallback */
    border-bottom: 1px solid var(--anndata-border-light);
    color: #6c757d; /* Fallback */
    color: var(--anndata-text-secondary);
}

/* Tooltip */
.anndata-repr [data-tooltip] {
    position: relative;
}

.anndata-repr [data-tooltip]:hover::after {
    content: attr(data-tooltip);
    position: absolute;
    bottom: 100%;
    left: 50%;
    transform: translateX(-50%);
    padding: 4px 8px;
    font-size: 11px;
    font-weight: normal;
    color: white;
    background: #333;
    border-radius: var(--anndata-radius);
    white-space: nowrap;
    z-index: 1000;
    pointer-events: none;
}
"""
