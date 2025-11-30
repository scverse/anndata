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

    // Update button visibility on window resize
    let resizeTimer;
    window.addEventListener('resize', () => {
        clearTimeout(resizeTimer);
        resizeTimer = setTimeout(() => {
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
        }, 100);
    });

    // README modal functionality
    const readmeIcon = container.querySelector('.adata-readme-icon');
    if (readmeIcon) {
        readmeIcon.addEventListener('click', (e) => {
            e.stopPropagation();
            const readmeContent = readmeIcon.dataset.readme;
            if (!readmeContent) return;

            // Create modal overlay
            const overlay = document.createElement('div');
            overlay.className = 'adata-readme-overlay';

            // Create modal
            const modal = document.createElement('div');
            modal.className = 'adata-readme-modal';

            // Header
            const header = document.createElement('div');
            header.className = 'adata-readme-header';
            header.innerHTML = '<h3>README</h3>';

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

    // Simple markdown parser (handles common elements)
    function parseMarkdown(text) {
        // Decode HTML entities first (the content was escaped)
        const textarea = document.createElement('textarea');
        textarea.innerHTML = text;
        text = textarea.value;

        // Escape HTML to prevent XSS, then apply markdown
        text = text
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;');

        // Code blocks (``` ... ```)
        text = text.replace(/```(\\w*)\\n([\\s\\S]*?)```/g, '<pre><code>$2</code></pre>');

        // Inline code (`...`)
        text = text.replace(/`([^`]+)`/g, '<code>$1</code>');

        // Headers
        text = text.replace(/^#### (.+)$/gm, '<h4>$1</h4>');
        text = text.replace(/^### (.+)$/gm, '<h3>$1</h3>');
        text = text.replace(/^## (.+)$/gm, '<h2>$1</h2>');
        text = text.replace(/^# (.+)$/gm, '<h1>$1</h1>');

        // Bold and italic
        text = text.replace(/\\*\\*\\*(.+?)\\*\\*\\*/g, '<strong><em>$1</em></strong>');
        text = text.replace(/\\*\\*(.+?)\\*\\*/g, '<strong>$1</strong>');
        text = text.replace(/\\*(.+?)\\*/g, '<em>$1</em>');
        text = text.replace(/_(.+?)_/g, '<em>$1</em>');

        // Links [text](url)
        text = text.replace(/\\[([^\\]]+)\\]\\(([^)]+)\\)/g, '<a href="$2" target="_blank" rel="noopener">$1</a>');

        // Blockquotes
        text = text.replace(/^&gt; (.+)$/gm, '<blockquote>$1</blockquote>');

        // Unordered lists - use temporary marker to distinguish from ordered
        text = text.replace(/^[\\*\\-] (.+)$/gm, '<ul-li>$1</ul-li>');

        // Ordered lists - use temporary marker
        text = text.replace(/^\\d+\\. (.+)$/gm, '<ol-li>$1</ol-li>');

        // Wrap consecutive unordered items in <ul>
        text = text.replace(/(<ul-li>.*<\\/ul-li>\\n?)+/g, '<ul>$&</ul>');
        text = text.replace(/<ul-li>/g, '<li>');
        text = text.replace(/<\\/ul-li>/g, '</li>');

        // Wrap consecutive ordered items in <ol>
        text = text.replace(/(<ol-li>.*<\\/ol-li>\\n?)+/g, '<ol>$&</ol>');
        text = text.replace(/<ol-li>/g, '<li>');
        text = text.replace(/<\\/ol-li>/g, '</li>');

        // Paragraphs (double newlines)
        text = text.replace(/\\n\\n+/g, '</p><p>');
        text = '<p>' + text + '</p>';

        // Clean up empty paragraphs
        text = text.replace(/<p><\\/p>/g, '');
        text = text.replace(/<p>(<h[1-4]>)/g, '$1');
        text = text.replace(/(<\\/h[1-4]>)<\\/p>/g, '$1');
        text = text.replace(/<p>(<ul>)/g, '$1');
        text = text.replace(/(<\\/ul>)<\\/p>/g, '$1');
        text = text.replace(/<p>(<ol>)/g, '$1');
        text = text.replace(/(<\\/ol>)<\\/p>/g, '$1');
        text = text.replace(/<p>(<pre>)/g, '$1');
        text = text.replace(/(<\\/pre>)<\\/p>/g, '$1');
        text = text.replace(/<p>(<blockquote>)/g, '$1');
        text = text.replace(/(<\\/blockquote>)<\\/p>/g, '$1');

        // Single newlines to <br> within paragraphs
        text = text.replace(/([^>])\\n([^<])/g, '$1<br>$2');

        return text;
    }
"""
