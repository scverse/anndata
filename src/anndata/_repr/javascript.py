"""
JavaScript for AnnData HTML representation interactivity.

Provides:
- Section folding/unfolding
- Search/filter functionality across all levels
- Copy to clipboard
- Nested content expansion
"""

from __future__ import annotations


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
    return f"""<script>
(function() {{
    // Scoped to avoid conflicts with multiple repr instances
    const container = document.getElementById('{container_id}');
    if (!container) return;

    {_JS_CONTENT}
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

        // Filter all entries
        container.querySelectorAll('.adata-entry').forEach(entry => {
            totalEntries++;

            if (!query) {
                entry.classList.remove('hidden');
                totalMatches++;
                return;
            }

            const key = (entry.dataset.key || '').toLowerCase();
            const dtype = (entry.dataset.dtype || '').toLowerCase();
            const text = entry.textContent.toLowerCase();

            const matches = key.includes(query) || dtype.includes(query) || text.includes(query);

            if (matches) {
                entry.classList.remove('hidden');
                totalMatches++;

                // Expand parent sections to show match
                const section = entry.closest('.anndata-sec');
                if (section && section.classList.contains('collapsed')) {
                    section.classList.remove('collapsed');
                }

                // Expand nested content if match is inside
                const nestedContent = entry.closest('.adata-nested-content');
                if (nestedContent && !nestedContent.classList.contains('expanded')) {
                    nestedContent.classList.add('expanded');
                }
            } else {
                entry.classList.add('hidden');
            }
        });

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
            const sectionHeader = section.querySelector('.anndata-sechdr');

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

                // Visual feedback
                btn.classList.add('copied');
                const originalText = btn.textContent;
                btn.textContent = '✓';

                setTimeout(() => {
                    btn.classList.remove('copied');
                    btn.textContent = originalText;
                }, 1500);
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
            container.querySelectorAll('.adata-section.collapsed').forEach(section => {
                section.classList.remove('collapsed');
            });
        });
    }

    if (collapseAllBtn) {
        collapseAllBtn.addEventListener('click', () => {
            container.querySelectorAll('.adata-section:not(.collapsed)').forEach(section => {
                section.classList.add('collapsed');
            });
        });
    }
"""
