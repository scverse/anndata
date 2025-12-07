"""
JavaScript for AnnData HTML representation interactivity.

Provides:
- Section folding/unfolding
- Search/filter functionality across all levels
- Copy to clipboard
- Nested content expansion
- README modal with markdown rendering
"""

from __future__ import annotations

from anndata._repr.markdown import get_markdown_parser_js


def get_javascript(container_id: str) -> str:
    """
    Get the JavaScript code for a specific container.

    Parameters
    ----------
    container_id
        Unique ID for the container element

    Returns
    -------
    JavaScript code wrapped in script tags
    """
    markdown_parser = get_markdown_parser_js()
    return f"""<script>
(function() {{
    // Scoped to avoid conflicts with multiple repr instances
    const container = document.getElementById('{container_id}');
    if (!container) return;

    {_JS_CONTENT}

    {markdown_parser}
}})();
</script>"""


_JS_CONTENT = """
    // Mark container as JS-enabled (shows interactive elements)
    container.classList.add('js-enabled');

    // Show interactive elements (hidden by default for no-JS graceful degradation)
    container.querySelectorAll('.adata-fold-icon').forEach(icon => {
        icon.style.display = 'inline-flex';
    });
    container.querySelectorAll('.adata-copy-btn').forEach(btn => {
        btn.style.display = 'inline-flex';
    });
    container.querySelectorAll('.adata-search-input').forEach(input => {
        input.style.display = 'inline-block';
    });
    container.querySelectorAll('.adata-expand-btn').forEach(btn => {
        btn.style.display = 'inline-block';
    });
    // Filter indicator is shown via CSS .active class, no need to set display here

    // Apply initial collapse state from data attributes
    container.querySelectorAll('.anndata-sec[data-should-collapse="true"]').forEach(section => {
        section.classList.add('collapsed');
        const content = section.querySelector('.anndata-seccontent');
        if (content) {
            content.style.maxHeight = '0';
            content.style.overflow = 'hidden';
        }
    });

    // Toggle section fold/unfold
    function toggleSection(header) {
        const section = header.closest('.anndata-sec');
        if (!section) return;

        const isCollapsed = section.classList.toggle('collapsed');
        const content = section.querySelector('.anndata-seccontent');
        if (content) {
            content.style.maxHeight = isCollapsed ? '0' : '';
            content.style.overflow = isCollapsed ? 'hidden' : '';
            content.setAttribute('aria-hidden', isCollapsed);
        }

        // Rotate fold icon
        const icon = section.querySelector('.adata-fold-icon');
        if (icon) {
            icon.style.transform = isCollapsed ? 'rotate(-90deg)' : '';
        }
    }

    // Attach click handlers to section headers
    container.querySelectorAll('.anndata-sechdr').forEach(header => {
        header.addEventListener('click', (e) => {
            // Don't toggle if clicking on help link
            if (e.target.closest('.adata-help-link')) return;
            toggleSection(header);
        });

        // Keyboard accessibility
        header.setAttribute('tabindex', '0');
        header.setAttribute('role', 'button');
        header.addEventListener('keydown', (e) => {
            if (e.key === 'Enter' || e.key === ' ') {
                e.preventDefault();
                toggleSection(header);
            }
        });
    });

    // Search/filter functionality
    const searchInput = container.querySelector('.adata-search-input');
    const filterIndicator = container.querySelector('.adata-filter-indicator');

    if (searchInput) {
        let debounceTimer;

        searchInput.addEventListener('input', (e) => {
            clearTimeout(debounceTimer);
            debounceTimer = setTimeout(() => {
                filterEntries(e.target.value.toLowerCase().trim());
            }, 150);
        });

        // Clear on Escape
        searchInput.addEventListener('keydown', (e) => {
            if (e.key === 'Escape') {
                searchInput.value = '';
                filterEntries('');
            }
        });
    }

    function filterEntries(query) {
        let totalMatches = 0;
        let totalEntries = 0;

        // Get all entries
        const allEntries = container.querySelectorAll('.adata-entry');

        // First pass: determine which entries match directly
        const matchingEntries = new Set();
        allEntries.forEach(entry => {
            if (!query) {
                matchingEntries.add(entry);
                return;
            }

            const key = (entry.dataset.key || '').toLowerCase();
            const dtype = (entry.dataset.dtype || '').toLowerCase();
            const text = entry.textContent.toLowerCase();

            if (key.includes(query) || dtype.includes(query) || text.includes(query)) {
                matchingEntries.add(entry);
            }
        });

        // Second pass: for matching entries inside nested content, also show their parent entries
        // This ensures hierarchical search works (e.g., searching for a field in nested AnnData
        // keeps the parent entry visible)
        const entriesToAdd = [];
        matchingEntries.forEach(entry => {
            // Walk up the DOM to find if this entry is inside nested content
            let element = entry.parentElement;
            let iterations = 0;
            const maxIterations = 100; // Guard against infinite loops
            while (element && element !== container && iterations < maxIterations) {
                iterations++;
                if (element.classList.contains('adata-nested-content')) {
                    // Found nested content - find the parent entry (previous sibling of nested-row)
                    const nestedRow = element.closest('.adata-nested-row');
                    if (nestedRow && nestedRow.previousElementSibling &&
                        nestedRow.previousElementSibling.classList.contains('adata-entry')) {
                        entriesToAdd.push(nestedRow.previousElementSibling);
                    }
                    break; // Found the nesting level, stop here
                }
                element = element.parentElement;
            }
        });
        entriesToAdd.forEach(e => matchingEntries.add(e));

        // Apply visibility to entries
        allEntries.forEach(entry => {
            totalEntries++;

            if (matchingEntries.has(entry)) {
                entry.classList.remove('hidden');
                totalMatches++;

                // Expand parent sections to show match
                const section = entry.closest('.anndata-sec');
                if (section && section.classList.contains('collapsed')) {
                    section.classList.remove('collapsed');
                }

                // Expand nested content if match is inside nested area
                const nestedContent = entry.closest('.adata-nested-content');
                if (nestedContent) {
                    const nestedRow = nestedContent.closest('.adata-nested-row');
                    if (nestedRow && !nestedRow.classList.contains('expanded')) {
                        nestedRow.classList.add('expanded');
                    }
                }
            } else {
                entry.classList.add('hidden');
            }
        });

        // Also filter X entries in nested AnnData (they use adata-x-entry class, not adata-entry)
        // This prevents orphaned X rows from showing when their sibling entries are hidden
        if (query) {
            container.querySelectorAll('.adata-nested-content .adata-x-entry').forEach(xEntry => {
                // Check if the nested AnnData has any visible entries
                const nestedRepr = xEntry.closest('.anndata-repr');
                if (nestedRepr) {
                    const hasVisibleEntries = nestedRepr.querySelector('.adata-entry:not(.hidden)');
                    xEntry.style.display = hasVisibleEntries ? '' : 'none';
                }
            });
        } else {
            // Reset X entries when no query
            container.querySelectorAll('.adata-nested-content .adata-x-entry').forEach(xEntry => {
                xEntry.style.display = '';
            });
        }

        // Update filter indicator
        if (filterIndicator) {
            if (query) {
                filterIndicator.classList.add('active');
                filterIndicator.textContent = `Showing ${totalMatches} of ${totalEntries}`;
            } else {
                filterIndicator.classList.remove('active');
            }
        }

        // Hide sections with no visible entries
        container.querySelectorAll('.anndata-sec').forEach(section => {
            const visibleEntries = section.querySelectorAll('.adata-entry:not(.hidden)');

            if (query && visibleEntries.length === 0) {
                section.style.display = 'none';
            } else {
                section.style.display = '';
            }
        });
    }

    // Copy to clipboard
    container.querySelectorAll('.adata-copy-btn').forEach(btn => {
        btn.addEventListener('click', async (e) => {
            e.stopPropagation();

            const text = btn.dataset.copy;
            if (!text) return;

            try {
                await navigator.clipboard.writeText(text);

                // Visual feedback (icon turns green via CSS)
                btn.classList.add('copied');
                setTimeout(() => btn.classList.remove('copied'), 1500);
            } catch (err) {
                // Fallback for older browsers
                const textarea = document.createElement('textarea');
                textarea.value = text;
                textarea.style.position = 'fixed';
                textarea.style.opacity = '0';
                document.body.appendChild(textarea);
                textarea.select();

                try {
                    document.execCommand('copy');
                    btn.classList.add('copied');
                    setTimeout(() => btn.classList.remove('copied'), 1500);
                } catch (e) {
                    console.error('Copy failed:', e);
                }

                document.body.removeChild(textarea);
            }
        });
    });

    // Expand/collapse nested content
    container.querySelectorAll('.adata-expand-btn').forEach(btn => {
        btn.addEventListener('click', (e) => {
            e.stopPropagation();

            const entry = btn.closest('.adata-entry');
            if (!entry) return;

            // The nested content is in a sibling <tr class="adata-nested-row">
            // which contains <td class="adata-nested-content">
            const nestedRow = entry.nextElementSibling;
            if (!nestedRow || !nestedRow.classList.contains('adata-nested-row')) return;

            const nestedContent = nestedRow.querySelector('.adata-nested-content');
            if (!nestedContent) return;

            const isExpanded = nestedRow.classList.toggle('expanded');

            btn.textContent = isExpanded ? 'Collapse ▲' : 'Expand ▼';
            btn.setAttribute('aria-expanded', isExpanded);
            nestedRow.setAttribute('aria-hidden', !isExpanded);
        });
    });

    // Expand all / Collapse all (if buttons exist)
    const expandAllBtn = container.querySelector('.adata-expand-all');
    const collapseAllBtn = container.querySelector('.adata-collapse-all');

    if (expandAllBtn) {
        expandAllBtn.addEventListener('click', () => {
            container.querySelectorAll('.anndata-sec.collapsed').forEach(section => {
                section.classList.remove('collapsed');
            });
        });
    }

    if (collapseAllBtn) {
        collapseAllBtn.addEventListener('click', () => {
            container.querySelectorAll('.anndata-sec:not(.collapsed)').forEach(section => {
                section.classList.add('collapsed');
            });
        });
    }

    // Helper to check if element is overflowing
    function isOverflowing(el) {
        return el.scrollWidth > el.clientWidth;
    }

    // Helper to update wrap button visibility based on overflow
    function updateWrapButtonVisibility(btn, list, metaCell) {
        if (!list || !metaCell) {
            btn.style.display = 'none';
            return;
        }
        // Show button only if content is overflowing or currently wrapped
        const isWrapped = list.classList.contains('wrapped');
        const overflows = isOverflowing(metaCell);
        btn.style.display = (overflows || isWrapped) ? 'inline' : 'none';
    }

    // Toggle category list wrap mode (single line vs multi-line)
    container.querySelectorAll('.adata-cats-wrap-btn').forEach(btn => {
        const typeCell = btn.closest('.adata-entry-type');
        const metaCell = typeCell ? typeCell.nextElementSibling : null;
        const catsList = metaCell ? metaCell.querySelector('.adata-cats-list') : null;

        // Initial visibility check
        updateWrapButtonVisibility(btn, catsList, metaCell);

        btn.addEventListener('click', (e) => {
            e.stopPropagation();
            if (!catsList || !metaCell) return;

            const isWrapped = catsList.classList.toggle('wrapped');
            metaCell.classList.toggle('expanded', isWrapped);
            btn.textContent = isWrapped ? '▲' : '⋯';
            btn.title = isWrapped ? 'Collapse to single line' : 'Toggle multi-line view';
            // Always show button when wrapped
            btn.style.display = 'inline';
        });
    });

    // Toggle DataFrame columns list wrap mode (single line vs multi-line)
    container.querySelectorAll('.adata-cols-wrap-btn').forEach(btn => {
        const typeCell = btn.closest('.adata-entry-type');
        const metaCell = typeCell ? typeCell.nextElementSibling : null;
        const colsList = metaCell ? metaCell.querySelector('.adata-cols-list') : null;

        // Initial visibility check
        updateWrapButtonVisibility(btn, colsList, metaCell);

        btn.addEventListener('click', (e) => {
            e.stopPropagation();
            if (!colsList || !metaCell) return;

            const isWrapped = colsList.classList.toggle('wrapped');
            metaCell.classList.toggle('expanded', isWrapped);
            btn.textContent = isWrapped ? '▲' : '⋯';
            btn.title = isWrapped ? 'Collapse to single line' : 'Toggle multi-line view';
            // Always show button when wrapped
            btn.style.display = 'inline';
        });
    });

    // Update button visibility on container resize (works for JupyterLab panes too)
    function updateAllWrapButtons() {
        container.querySelectorAll('.adata-cats-wrap-btn').forEach(btn => {
            const typeCell = btn.closest('.adata-entry-type');
            const metaCell = typeCell ? typeCell.nextElementSibling : null;
            const catsList = metaCell ? metaCell.querySelector('.adata-cats-list') : null;
            updateWrapButtonVisibility(btn, catsList, metaCell);
        });
        container.querySelectorAll('.adata-cols-wrap-btn').forEach(btn => {
            const typeCell = btn.closest('.adata-entry-type');
            const metaCell = typeCell ? typeCell.nextElementSibling : null;
            const colsList = metaCell ? metaCell.querySelector('.adata-cols-list') : null;
            updateWrapButtonVisibility(btn, colsList, metaCell);
        });
    }

    // Use ResizeObserver for robust resize detection (pane resizes, not just window)
    if (typeof ResizeObserver !== 'undefined') {
        let resizeTimer;
        const resizeObserver = new ResizeObserver(() => {
            clearTimeout(resizeTimer);
            resizeTimer = setTimeout(updateAllWrapButtons, 100);
        });
        resizeObserver.observe(container);
    } else {
        // Fallback for older browsers
        let resizeTimer;
        window.addEventListener('resize', () => {
            clearTimeout(resizeTimer);
            resizeTimer = setTimeout(updateAllWrapButtons, 100);
        });
    }

    // README modal functionality
    const readmeIcon = container.querySelector('.adata-readme-icon');
    if (readmeIcon) {
        // Ensure accessibility attributes
        readmeIcon.setAttribute('role', 'button');
        readmeIcon.setAttribute('tabindex', '0');
        readmeIcon.setAttribute('aria-label', 'View README');

        readmeIcon.addEventListener('click', (e) => {
            e.stopPropagation();
            const readmeContent = readmeIcon.dataset.readme;
            if (!readmeContent) return;

            // Create modal overlay
            const overlay = document.createElement('div');
            overlay.className = 'adata-readme-overlay';

            // Create modal with accessibility attributes
            const modal = document.createElement('div');
            modal.className = 'adata-readme-modal';
            modal.setAttribute('role', 'dialog');
            modal.setAttribute('aria-modal', 'true');
            modal.setAttribute('aria-labelledby', 'readme-modal-title');

            // Header
            const header = document.createElement('div');
            header.className = 'adata-readme-header';
            header.innerHTML = '<h3 id="readme-modal-title">README</h3>';

            const closeBtn = document.createElement('button');
            closeBtn.className = 'adata-readme-close';
            closeBtn.textContent = '×';
            closeBtn.setAttribute('aria-label', 'Close');
            header.appendChild(closeBtn);

            // Content
            const content = document.createElement('div');
            content.className = 'adata-readme-content';

            // Parse markdown to HTML (simple conversion)
            content.innerHTML = parseMarkdown(readmeContent);

            modal.appendChild(header);
            modal.appendChild(content);
            overlay.appendChild(modal);

            // Add to container (scoped styles apply)
            container.appendChild(overlay);

            // Close handlers
            const closeModal = () => {
                overlay.remove();
            };

            closeBtn.addEventListener('click', closeModal);
            overlay.addEventListener('click', (e) => {
                if (e.target === overlay) closeModal();
            });

            // Escape key closes modal
            const escHandler = (e) => {
                if (e.key === 'Escape') {
                    closeModal();
                    document.removeEventListener('keydown', escHandler);
                }
            };
            document.addEventListener('keydown', escHandler);

            // Focus trap
            closeBtn.focus();
        });

        // Keyboard accessibility for the icon
        readmeIcon.addEventListener('keydown', (e) => {
            if (e.key === 'Enter' || e.key === ' ') {
                e.preventDefault();
                readmeIcon.click();
            }
        });
    }
"""
