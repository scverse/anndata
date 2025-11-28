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
    --ad-bg-primary: #ffffff;
    --ad-bg-secondary: #f8f9fa;
    --ad-bg-tertiary: #e9ecef;
    --ad-text-primary: #212529;
    --ad-text-secondary: #6c757d;
    --ad-text-muted: #adb5bd;
    --ad-border-color: #dee2e6;
    --ad-border-light: #e9ecef;
    --ad-accent-color: #0d6efd;
    --ad-accent-hover: #0b5ed7;
    --ad-warning-color: #ffc107;
    --ad-warning-bg: #fff3cd;
    --ad-error-color: #dc3545;
    --ad-error-bg: #f8d7da;
    --ad-success-color: #198754;
    --ad-info-color: #0dcaf0;
    --ad-link-color: #0d6efd;
    --ad-code-bg: #f8f9fa;
    --ad-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
    --ad-radius: 4px;
    --ad-font-mono: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace;
    --ad-font-size: 13px;
    --ad-line-height: 1.4;

    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    font-size: var(--ad-font-size);
    line-height: var(--ad-line-height);
    color: var(--ad-text-primary);
    background: var(--ad-bg-primary);
    border: 1px solid var(--ad-border-color);
    border-radius: var(--ad-radius);
    padding: 0;
    margin: 8px 0;
    max-width: 100%;
    overflow: hidden;
}

/* Dark mode - via media query */
@media (prefers-color-scheme: dark) {
    .anndata-repr {
        --ad-bg-primary: #1e1e1e;
        --ad-bg-secondary: #252526;
        --ad-bg-tertiary: #2d2d2d;
        --ad-text-primary: #e0e0e0;
        --ad-text-secondary: #a0a0a0;
        --ad-text-muted: #707070;
        --ad-border-color: #404040;
        --ad-border-light: #333333;
        --ad-accent-color: #58a6ff;
        --ad-accent-hover: #79b8ff;
        --ad-warning-color: #d29922;
        --ad-warning-bg: #3d3200;
        --ad-error-color: #f85149;
        --ad-error-bg: #3d1a1a;
        --ad-success-color: #3fb950;
        --ad-info-color: #58a6ff;
        --ad-link-color: #58a6ff;
        --ad-code-bg: #2d2d2d;
        --ad-shadow: 0 1px 3px rgba(0, 0, 0, 0.3);
    }
}

/* Jupyter dark theme detection + manual dark mode class */
[data-jp-theme-light="false"] .anndata-repr,
.jp-Theme-Dark .anndata-repr,
body.vscode-dark .anndata-repr,
body[data-vscode-theme-kind="vscode-dark"] .anndata-repr,
body.dark-mode .anndata-repr {
    --ad-bg-primary: #1e1e1e;
    --ad-bg-secondary: #252526;
    --ad-bg-tertiary: #2d2d2d;
    --ad-text-primary: #e0e0e0;
    --ad-text-secondary: #a0a0a0;
    --ad-text-muted: #707070;
    --ad-border-color: #404040;
    --ad-border-light: #333333;
    --ad-accent-color: #58a6ff;
    --ad-accent-hover: #79b8ff;
    --ad-warning-color: #d29922;
    --ad-warning-bg: #3d3200;
    --ad-error-color: #f85149;
    --ad-error-bg: #3d1a1a;
    --ad-success-color: #3fb950;
    --ad-info-color: #58a6ff;
    --ad-link-color: #58a6ff;
    --ad-code-bg: #2d2d2d;
    --ad-shadow: 0 1px 3px rgba(0, 0, 0, 0.3);
}

/* Header */
.anndata-repr .ad-header {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 8px;
    padding: 10px 12px;
    background: var(--ad-bg-secondary);
    border-bottom: 1px solid var(--ad-border-color);
}

.anndata-repr .ad-type {
    font-weight: 600;
    font-size: 14px;
    color: var(--ad-accent-color);
}

.anndata-repr .ad-shape {
    font-family: var(--ad-font-mono);
    font-size: 12px;
    color: var(--ad-text-secondary);
}

.anndata-repr .ad-badge {
    display: inline-flex;
    align-items: center;
    gap: 4px;
    padding: 2px 8px;
    font-size: 11px;
    font-weight: 500;
    border-radius: 10px;
    white-space: nowrap;
}

.anndata-repr .ad-badge-view {
    background: var(--ad-info-color);
    color: white;
}

.anndata-repr .ad-badge-backed {
    background: var(--ad-success-color);
    color: white;
}

.anndata-repr .ad-badge-extension {
    background: var(--ad-accent-color);
    color: white;
}

/* Search box */
.anndata-repr .ad-search {
    padding: 8px 12px;
    background: var(--ad-bg-primary);
    border-bottom: 1px solid var(--ad-border-light);
}

.anndata-repr .ad-search-input {
    width: 100%;
    max-width: 300px;
    padding: 6px 10px;
    font-size: 12px;
    border: 1px solid var(--ad-border-color);
    border-radius: var(--ad-radius);
    background: var(--ad-bg-primary);
    color: var(--ad-text-primary);
    outline: none;
    transition: border-color 0.15s;
}

.anndata-repr .ad-search-input:focus {
    border-color: var(--ad-accent-color);
}

.anndata-repr .ad-search-input::placeholder {
    color: var(--ad-text-muted);
}

.anndata-repr .ad-filter-indicator {
    display: none;
    margin-left: 8px;
    font-size: 11px;
    color: var(--ad-accent-color);
}

.anndata-repr .ad-filter-indicator.active {
    display: inline;
}

/* Metadata bar */
.anndata-repr .ad-metadata {
    display: flex;
    flex-wrap: wrap;
    gap: 12px;
    padding: 6px 12px;
    font-size: 11px;
    color: var(--ad-text-secondary);
    background: var(--ad-bg-tertiary);
    border-bottom: 1px solid var(--ad-border-light);
}

.anndata-repr .ad-metadata span {
    white-space: nowrap;
}

/* Index preview */
.anndata-repr .ad-index-preview {
    padding: 8px 12px;
    font-size: 11px;
    font-family: var(--ad-font-mono);
    color: var(--ad-text-secondary);
    background: var(--ad-bg-primary);
    border-bottom: 1px solid var(--ad-border-light);
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
}

.anndata-repr .ad-index-preview strong {
    color: var(--ad-text-primary);
    font-weight: 500;
}

/* Sections container */
.anndata-repr .ad-sections {
    padding: 0;
}

/* Individual section */
.anndata-repr .ad-section {
    border-bottom: 1px solid var(--ad-border-light);
}

.anndata-repr .ad-section:last-child {
    border-bottom: none;
}

.anndata-repr .ad-section-header {
    display: flex;
    align-items: center;
    gap: 8px;
    padding: 8px 12px;
    cursor: pointer;
    user-select: none;
    background: var(--ad-bg-primary);
    transition: background-color 0.15s;
}

.anndata-repr .ad-section-header:hover {
    background: var(--ad-bg-secondary);
}

.anndata-repr .ad-fold-icon {
    width: 16px;
    font-size: 10px;
    color: var(--ad-text-muted);
    transition: transform 0.15s;
}

.anndata-repr .ad-section.collapsed .ad-fold-icon {
    transform: rotate(-90deg);
}

.anndata-repr .ad-section-name {
    font-weight: 600;
    color: var(--ad-text-primary);
}

.anndata-repr .ad-section-count {
    font-size: 11px;
    color: var(--ad-text-secondary);
}

.anndata-repr .ad-help-link {
    margin-left: auto;
    padding: 2px 6px;
    font-size: 11px;
    color: var(--ad-text-muted);
    text-decoration: none;
    border-radius: var(--ad-radius);
    transition: color 0.15s, background-color 0.15s;
}

.anndata-repr .ad-help-link:hover {
    color: var(--ad-accent-color);
    background: var(--ad-bg-tertiary);
}

/* Section content */
.anndata-repr .ad-section-content {
    padding: 0;
    overflow: hidden;
    transition: max-height 0.2s ease-out;
}

.anndata-repr .ad-section.collapsed .ad-section-content {
    max-height: 0 !important;
    padding: 0;
}

/* Entries table */
.anndata-repr .ad-table {
    width: 100%;
    border-collapse: collapse;
    font-size: 12px;
}

.anndata-repr .ad-entry {
    border-bottom: 1px solid var(--ad-border-light);
    transition: background-color 0.1s;
}

.anndata-repr .ad-entry:last-child {
    border-bottom: none;
}

.anndata-repr .ad-entry:hover {
    background: var(--ad-bg-secondary);
}

.anndata-repr .ad-entry.hidden {
    display: none;
}

.anndata-repr .ad-entry.warning {
    background: var(--ad-warning-bg);
}

.anndata-repr .ad-entry.error {
    background: var(--ad-error-bg);
}

.anndata-repr .ad-entry td {
    padding: 6px 12px;
    vertical-align: middle;
}

.anndata-repr .ad-entry-name {
    font-family: var(--ad-font-mono);
    font-weight: 500;
    color: var(--ad-text-primary);
    white-space: nowrap;
}

.anndata-repr .ad-entry-type {
    font-family: var(--ad-font-mono);
    font-size: 11px;
    color: var(--ad-text-secondary);
}

.anndata-repr .ad-entry-meta {
    font-size: 11px;
    color: var(--ad-text-muted);
    text-align: right;
}

/* Copy button */
.anndata-repr .ad-copy-btn {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 20px;
    height: 20px;
    margin-left: 4px;
    padding: 0;
    font-size: 11px;
    color: var(--ad-text-muted);
    background: transparent;
    border: none;
    border-radius: 3px;
    cursor: pointer;
    opacity: 0;
    transition: opacity 0.15s, color 0.15s, background-color 0.15s;
}

.anndata-repr .ad-entry:hover .ad-copy-btn {
    opacity: 1;
}

.anndata-repr .ad-copy-btn:hover {
    color: var(--ad-accent-color);
    background: var(--ad-bg-tertiary);
}

.anndata-repr .ad-copy-btn.copied {
    color: var(--ad-success-color);
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
.anndata-repr .dtype-warning { color: var(--ad-warning-color); }
.anndata-repr .dtype-dask { color: #fb8500; }
.anndata-repr .dtype-gpu { color: #76b900; }
.anndata-repr .dtype-awkward { color: #e85d04; }

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
.anndata-repr .ad-color-swatches {
    display: inline-flex;
    gap: 2px;
    margin-left: 6px;
    vertical-align: middle;
}

.anndata-repr .ad-color-swatch {
    display: inline-block;
    width: 12px;
    height: 12px;
    border-radius: 2px;
    border: 1px solid var(--ad-border-color);
}

/* Warning indicator */
.anndata-repr .ad-warning-icon {
    color: var(--ad-warning-color);
    margin-left: 4px;
}

/* Expandable nested content */
.anndata-repr .ad-expand-btn {
    padding: 2px 8px;
    font-size: 11px;
    color: var(--ad-accent-color);
    background: transparent;
    border: 1px solid var(--ad-accent-color);
    border-radius: var(--ad-radius);
    cursor: pointer;
    transition: background-color 0.15s, color 0.15s;
}

.anndata-repr .ad-expand-btn:hover {
    background: var(--ad-accent-color);
    color: white;
}

/* Nested row (contains nested AnnData) */
.anndata-repr .ad-nested-row {
    display: none;
}

.anndata-repr .ad-nested-row.expanded {
    display: table-row;
}

.anndata-repr .ad-nested-content {
    padding: 8px 12px 8px 24px;
    background: var(--ad-bg-secondary);
    border-top: 1px solid var(--ad-border-light);
}

/* Nested AnnData */
.anndata-repr .ad-nested-anndata {
    margin: 8px 0;
    border: 1px solid var(--ad-border-color);
    border-radius: var(--ad-radius);
}

/* X section (special) */
.anndata-repr .ad-x-section {
    padding: 8px 12px;
    background: var(--ad-bg-secondary);
    border-bottom: 1px solid var(--ad-border-light);
}

.anndata-repr .ad-x-info {
    display: grid;
    grid-template-columns: auto 1fr;
    gap: 4px 12px;
    font-size: 12px;
}

.anndata-repr .ad-x-info dt {
    color: var(--ad-text-secondary);
    font-weight: 500;
}

.anndata-repr .ad-x-info dd {
    margin: 0;
    font-family: var(--ad-font-mono);
    color: var(--ad-text-primary);
}

/* Truncation indicator */
.anndata-repr .ad-truncated {
    padding: 8px 12px;
    font-size: 11px;
    color: var(--ad-text-muted);
    text-align: center;
    font-style: italic;
}

/* Max depth indicator */
.anndata-repr .ad-max-depth {
    padding: 8px 12px;
    font-size: 11px;
    color: var(--ad-text-muted);
    background: var(--ad-bg-tertiary);
    border-radius: var(--ad-radius);
    text-align: center;
}

/* Empty section indicator */
.anndata-repr .ad-empty {
    padding: 8px 12px;
    font-size: 11px;
    color: var(--ad-text-muted);
    font-style: italic;
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
    border-radius: var(--ad-radius);
    white-space: nowrap;
    z-index: 1000;
    pointer-events: none;
}
"""
