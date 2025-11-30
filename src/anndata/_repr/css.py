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
    /* Column widths are set dynamically via inline style on container */

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

/* README icon in header */
.anndata-repr .adata-readme-icon {
    cursor: pointer;
    font-size: 14px;
    opacity: 0.7;
    transition: opacity 0.15s;
    margin-left: 4px;
}

.anndata-repr .adata-readme-icon:hover {
    opacity: 1;
}

/* README modal overlay */
.anndata-repr .adata-readme-overlay {
    position: fixed;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    background: rgba(0, 0, 0, 0.5);
    display: flex;
    align-items: center;
    justify-content: center;
    z-index: 10000;
    padding: 20px;
}

/* README modal */
.anndata-repr .adata-readme-modal {
    background: var(--anndata-bg-primary);
    border: 1px solid var(--anndata-border-color);
    border-radius: 8px;
    box-shadow: 0 4px 20px rgba(0, 0, 0, 0.2);
    max-width: 700px;
    max-height: 80vh;
    width: 100%;
    display: flex;
    flex-direction: column;
    overflow: hidden;
}

.anndata-repr .adata-readme-header {
    display: flex;
    align-items: center;
    justify-content: space-between;
    padding: 12px 16px;
    border-bottom: 1px solid var(--anndata-border-color);
    background: var(--anndata-bg-secondary);
}

.anndata-repr .adata-readme-header h3 {
    margin: 0;
    font-size: 14px;
    font-weight: 600;
    color: var(--anndata-text-primary);
}

.anndata-repr .adata-readme-close {
    background: none;
    border: none;
    font-size: 20px;
    cursor: pointer;
    color: var(--anndata-text-secondary);
    padding: 0 4px;
    line-height: 1;
}

.anndata-repr .adata-readme-close:hover {
    color: var(--anndata-text-primary);
}

.anndata-repr .adata-readme-content {
    padding: 16px;
    overflow-y: auto;
    font-size: 13px;
    line-height: 1.6;
    color: var(--anndata-text-primary);
}

/* Markdown rendering in README */
.anndata-repr .adata-readme-content h1,
.anndata-repr .adata-readme-content h2,
.anndata-repr .adata-readme-content h3,
.anndata-repr .adata-readme-content h4 {
    margin-top: 16px;
    margin-bottom: 8px;
    font-weight: 600;
    color: var(--anndata-text-primary);
}

.anndata-repr .adata-readme-content h1 {
    font-size: 1.5em;
    padding-bottom: 8px;
    border-bottom: 1px solid var(--anndata-border-color);
}
.anndata-repr .adata-readme-content h2 { font-size: 1.3em; }
.anndata-repr .adata-readme-content h3 { font-size: 1.1em; }

.anndata-repr .adata-readme-content p {
    margin: 0 0 12px 0;
}

.anndata-repr .adata-readme-content code {
    background: var(--anndata-code-bg);
    padding: 2px 6px;
    border-radius: 3px;
    font-family: var(--anndata-font-mono);
    font-size: 0.9em;
}

.anndata-repr .adata-readme-content pre {
    background: var(--anndata-code-bg);
    padding: 12px;
    border-radius: var(--anndata-radius);
    overflow-x: auto;
    margin: 12px 0;
}

.anndata-repr .adata-readme-content pre code {
    background: none;
    padding: 0;
}

.anndata-repr .adata-readme-content ul,
.anndata-repr .adata-readme-content ol {
    margin: 0 0 12px 24px;
    padding-left: 0;
    list-style-position: outside;
}

.anndata-repr .adata-readme-content li {
    margin: 4px 0;
    padding-left: 4px;
}

.anndata-repr .adata-readme-content a {
    color: var(--anndata-link-color);
    text-decoration: none;
}

.anndata-repr .adata-readme-content a:hover {
    text-decoration: underline;
}

.anndata-repr .adata-readme-content blockquote {
    margin: 12px 0;
    padding: 8px 16px;
    border-left: 3px solid var(--anndata-accent-color);
    background: var(--anndata-bg-secondary);
    color: var(--anndata-text-secondary);
}

.anndata-repr .adata-readme-content table {
    border-collapse: collapse;
    margin: 12px 0;
    font-size: 0.9em;
    width: auto;
}

.anndata-repr .adata-readme-content th,
.anndata-repr .adata-readme-content td {
    border: 1px solid var(--anndata-border-color);
    padding: 6px 12px;
    text-align: left;
}

.anndata-repr .adata-readme-content th {
    background: var(--anndata-bg-secondary);
    font-weight: 600;
}

.anndata-repr .adata-readme-content tbody tr:hover {
    background: var(--anndata-bg-secondary);
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

/* When JS is enabled, use fixed layout so meta column fills remaining space */
.anndata-repr.js-enabled .adata-table {
    table-layout: fixed;
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
    font-size: 13px; /* Fallback */
    font-size: var(--anndata-font-size);
    font-weight: 500;
    color: #212529; /* Fallback */
    color: var(--anndata-text-primary);
    white-space: nowrap;
    text-align: left;
    /* Fixed width for name column */
    width: var(--anndata-name-col-width, 150px);
}

/* Name cell inner container: text truncates, copy button stays visible */
.anndata-repr .adata-entry-name-inner {
    display: flex;
    align-items: center;
    gap: 4px;
    min-width: 0;  /* Allow flex child to shrink below content size */
}

.anndata-repr .adata-name-text {
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
    flex: 1;
    min-width: 0;  /* Allow text to shrink and show ellipsis */
    font-size: inherit;  /* Prevent mobile browsers from auto-sizing truncated text */
}

.anndata-repr .adata-entry-type {
    font-family: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace; /* Fallback */
    font-family: var(--anndata-font-mono);
    font-size: 11px;
    color: #6c757d; /* Fallback */
    color: var(--anndata-text-secondary);
    text-align: left;
    white-space: nowrap;
    /* Fixed width for consistent alignment (dtype names are predictable) */
    width: var(--anndata-type-col-width, 180px);
    min-width: var(--anndata-type-col-width, 180px);
}

.anndata-repr .adata-entry-meta {
    font-size: 11px;
    color: #adb5bd; /* Fallback */
    color: var(--anndata-text-muted);
    text-align: right;
    /* Default: allow wrapping for graceful no-JS degradation */
    white-space: normal;
    word-break: break-word;
    /* Takes remaining space */
}

/* When JS is enabled, meta cell truncates with ellipsis */
.anndata-repr.js-enabled .adata-entry-meta {
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
}

.anndata-repr .adata-entry-meta.expanded {
    white-space: normal;
    overflow: visible;
}

/* Copy button - uses CSS to draw two overlapping squares icon */
.anndata-repr .adata-copy-btn {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 16px;
    height: 16px;
    padding: 0;
    background: transparent;
    border: none;
    border-radius: 2px;
    cursor: pointer;
    opacity: 0;
    transition: opacity 0.15s, background-color 0.15s;
    flex-shrink: 0;  /* Don't shrink the button */
    position: relative;
}

/* Two overlapping squares icon using pseudo-elements */
.anndata-repr .adata-copy-btn::before,
.anndata-repr .adata-copy-btn::after {
    content: '';
    position: absolute;
    width: 7px;
    height: 8px;
    border: 1.5px solid var(--anndata-text-muted);
    border-radius: 1px;
    background: var(--anndata-bg-primary);
}

.anndata-repr .adata-copy-btn::before {
    top: 2px;
    left: 2px;
}

.anndata-repr .adata-copy-btn::after {
    top: 5px;
    left: 5px;
}

.anndata-repr .adata-entry:hover .adata-copy-btn {
    opacity: 1;
}

.anndata-repr .adata-copy-btn:hover::before,
.anndata-repr .adata-copy-btn:hover::after {
    border-color: var(--anndata-accent-color);
}

/* When copied, hide squares and show checkmark */
.anndata-repr .adata-copy-btn.copied::before,
.anndata-repr .adata-copy-btn.copied::after {
    display: none;
}

.anndata-repr .adata-copy-btn.copied::before {
    display: block;
    content: '✓';
    width: auto;
    height: auto;
    border: none;
    background: none;
    color: var(--anndata-success-color);
    font-size: 12px;
    font-weight: bold;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%);
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
/* Warning class kept for backwards compatibility but no longer used for coloring */
.anndata-repr .dtype-warning {
    /* Type name keeps its original color; only the icon is styled */
}
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

/* Category list with wrap toggle */
/* Default: multi-line for no-JS graceful degradation */
.anndata-repr .adata-cats-list {
    display: inline;
    white-space: normal;
    word-break: break-word;
}

/* When JS is enabled, categories inherit truncation from parent cell */
.anndata-repr.js-enabled .adata-cats-list {
    display: inline;
    white-space: nowrap;
    word-break: normal;
}

/* Wrapped state (toggled by JS) */
.anndata-repr.js-enabled .adata-cats-list.wrapped {
    display: inline;
    max-width: none;
    white-space: normal;
    word-break: break-word;
    overflow: visible;
    text-overflow: clip;
}

/* Wrap buttons - hidden by default for no-JS graceful degradation */
.anndata-repr .adata-cats-wrap-btn,
.anndata-repr .adata-cols-wrap-btn {
    display: none;
    background: transparent;
    border: none;
    color: var(--anndata-text-muted);
    cursor: pointer;
    font-size: 11px;
    padding: 0 4px;
    margin-left: 4px;
    transition: color 0.15s;
    vertical-align: middle;
}

/* Show wrap buttons when JS is enabled */
.anndata-repr.js-enabled .adata-cats-wrap-btn,
.anndata-repr.js-enabled .adata-cols-wrap-btn {
    display: inline-block;
}

.anndata-repr .adata-cats-wrap-btn:hover,
.anndata-repr .adata-cols-wrap-btn:hover {
    color: var(--anndata-accent-color);
}

/* DataFrame columns list with wrap toggle */
/* Default: multi-line for no-JS graceful degradation */
.anndata-repr .adata-cols-list {
    display: inline;
    white-space: normal;
    word-break: break-word;
}

/* When JS is enabled, columns inherit truncation from parent cell */
.anndata-repr.js-enabled .adata-cols-list {
    display: inline;
    white-space: nowrap;
    word-break: normal;
}

/* Wrapped state (toggled by JS) */
.anndata-repr.js-enabled .adata-cols-list.wrapped {
    display: inline;
    max-width: none;
    white-space: normal;
    word-break: break-word;
    overflow: visible;
    text-overflow: clip;
}

/* Warning indicator */
.anndata-repr .adata-warning-icon {
    font-size: 10px;
    font-weight: 600;
    color: var(--anndata-warning-color);
    cursor: help;
    margin-left: 2px;
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
    display: flex;
    justify-content: space-between;
    padding: 4px 12px;
    font-size: 10px;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    color: var(--anndata-text-muted);
    border-top: 1px solid var(--anndata-border-light);
}

/* Muted text helper */
.anndata-repr .adata-text-muted {
    color: var(--anndata-text-muted);
}

/* X entry row */
.anndata-repr .adata-x-entry {
    display: flex;
    align-items: center;
    gap: 12px;
    padding: 6px 12px;
    border-bottom: 1px solid #e9ecef; /* Fallback */
    border-bottom: 1px solid var(--anndata-border-light);
    color: #6c757d; /* Fallback */
    color: var(--anndata-text-secondary);
}

.anndata-repr .adata-x-entry > span:first-child {
    font-family: var(--anndata-font-mono);
    font-weight: 600;
    min-width: 60px;
}

.anndata-repr .adata-x-entry > span:last-child {
    font-family: var(--anndata-font-mono);
    font-size: 11px;
}

/* Category item in list */
.anndata-repr .adata-cat-item {
    display: inline-flex;
    align-items: center;
    gap: 3px;
    margin-right: 8px;
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
