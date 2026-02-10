// AnnData HTML Representation JavaScript
// This file provides interactivity for the HTML repr.
// The {container_id} placeholder is replaced at runtime.

// Mark container as JS-enabled (shows interactive elements)
container.classList.add("anndata-repr--js");

// Show interactive elements (hidden by default for no-JS graceful degradation)
container.querySelectorAll(".anndata-entry__copy").forEach((btn) => {
    btn.style.display = "inline-flex";
});
container.querySelectorAll(".anndata-search__box").forEach((box) => {
    box.style.display = "inline-flex";
});
container.querySelectorAll(".anndata-entry__expand").forEach((btn) => {
    btn.style.display = "inline-block";
});
container.querySelectorAll(".anndata-search__toggle").forEach((btn) => {
    btn.style.display = "inline-flex";
});
// Filter indicator is shown via CSS .active class, no need to set display here

// Section collapse is handled natively by <details>/<summary> elements.
// No JS needed for section toggle — the browser handles open/close state.

// Search/filter functionality
const searchBox = container.querySelector(".anndata-search__box");
const searchInput = container.querySelector(".anndata-search__input");
const filterIndicator = container.querySelector(".anndata-search__indicator");
const caseToggle = container.querySelector(".anndata-search__toggle--case");
const regexToggle = container.querySelector(".anndata-search__toggle--regex");

// Search state
let caseSensitive = false;
let useRegex = false;

if (searchInput) {
    let debounceTimer;

    const triggerFilter = () => {
        clearTimeout(debounceTimer);
        debounceTimer = setTimeout(() => {
            filterEntries(searchInput.value.trim());
        }, 150);
    };

    searchInput.addEventListener("input", triggerFilter);

    // Clear on Escape
    searchInput.addEventListener("keydown", (e) => {
        if (e.key === "Escape") {
            searchInput.value = "";
            filterEntries("");
        }
    });

    // Toggle button handlers
    if (caseToggle) {
        caseToggle.addEventListener("click", (e) => {
            e.stopPropagation();
            caseSensitive = !caseSensitive;
            caseToggle.classList.toggle("anndata--active", caseSensitive);
            caseToggle.setAttribute("aria-pressed", caseSensitive);
            triggerFilter();
        });
    }

    if (regexToggle) {
        regexToggle.addEventListener("click", (e) => {
            e.stopPropagation();
            useRegex = !useRegex;
            regexToggle.classList.toggle("anndata--active", useRegex);
            regexToggle.setAttribute("aria-pressed", useRegex);
            triggerFilter();
        });
    }
}

// Helper: test if text matches query (respects case sensitivity and regex mode)
function matchesQuery(text, query) {
    if (!query) return true;
    if (useRegex) {
        try {
            const flags = caseSensitive ? "" : "i";
            const regex = new RegExp(query, flags);
            if (searchBox)
                searchBox.classList.remove("anndata-search__box--error");
            return regex.test(text);
        } catch (e) {
            // Invalid regex - show error state but don't crash
            if (searchBox)
                searchBox.classList.add("anndata-search__box--error");
            return false;
        }
    } else {
        if (searchBox) searchBox.classList.remove("anndata-search__box--error");
        if (caseSensitive) {
            return text.includes(query);
        } else {
            return text.toLowerCase().includes(query.toLowerCase());
        }
    }
}

function filterEntries(query) {
    let totalMatches = 0;
    let totalEntries = 0;

    // First pass: mark all entries as hidden or not based on direct match
    const entries = container.querySelectorAll(".anndata-entry");
    const directMatches = new Set();

    entries.forEach((entry) => {
        totalEntries++;

        const key = entry.dataset.key || "";
        const dtype = entry.dataset.dtype || "";
        const text = entry.textContent;

        const matches =
            !query ||
            matchesQuery(key, query) ||
            matchesQuery(dtype, query) ||
            matchesQuery(text, query);

        if (matches) {
            directMatches.add(entry);
            entry.classList.remove("anndata-entry--hidden");
            totalMatches++;

            // Expand parent sections to show match
            const section = entry.closest(".anndata-section");
            if (section && !section.open) {
                section.open = true;
            }

            // Expand nested content if match is inside nested area
            const nestedContent = entry.closest(
                ".anndata-entry__nested-content",
            );
            if (nestedContent) {
                const nestedRow = nestedContent.closest(
                    ".anndata-entry--nested",
                );
                if (
                    nestedRow &&
                    !nestedRow.classList.contains("anndata-entry--expanded")
                ) {
                    nestedRow.classList.add("anndata-entry--expanded");
                }
            }
        } else {
            entry.classList.add("anndata-entry--hidden");
        }
    });

    // Second pass: if a nested entry matches, show all ancestor entry rows
    // This ensures that when searching for something inside a nested AnnData,
    // all parent rows remain visible so the user can expand them to see the match
    if (query) {
        directMatches.forEach((matchedEntry) => {
            // Walk up the DOM tree to find and show all parent entry rows
            // Safety limit prevents infinite loops (max nesting depth is typically 3)
            let element = matchedEntry;
            let iterations = 0;
            const maxIterations = 20;

            while (
                element &&
                element !== container &&
                iterations < maxIterations
            ) {
                iterations++;
                // Check if we're inside a nested content container
                const nestedContainer = element.closest(
                    ".anndata-entry__nested-content",
                );
                if (!nestedContainer) break;

                // Find the parent row that contains this nested content
                // Structure: tr.anndata-entry > tr.anndata-entry--nested > td.anndata-entry__nested-content
                const nestedRow = nestedContainer.closest(
                    ".anndata-entry--nested",
                );
                if (!nestedRow) break;

                // The parent entry is the previous sibling row
                const parentEntry = nestedRow.previousElementSibling;
                if (
                    parentEntry &&
                    parentEntry.classList.contains("anndata-entry")
                ) {
                    if (
                        parentEntry.classList.contains("anndata-entry--hidden")
                    ) {
                        parentEntry.classList.remove("anndata-entry--hidden");
                        totalMatches++;
                    }
                }
                // Continue searching from the parent's container
                element = nestedRow.parentElement;
            }
        });
    }

    // Also filter X entries in nested AnnData (they use anndata-x__entry class, not anndata-entry)
    // This prevents orphaned X rows from showing when their sibling entries are hidden
    if (query) {
        container
            .querySelectorAll(
                ".anndata-entry__nested-content .anndata-x__entry",
            )
            .forEach((xEntry) => {
                // Check if the nested AnnData has any visible entries
                const nestedRepr = xEntry.closest(".anndata-repr");
                if (nestedRepr) {
                    const hasVisibleEntries = nestedRepr.querySelector(
                        ".anndata-entry:not(.anndata-entry--hidden)",
                    );
                    xEntry.style.display = hasVisibleEntries ? "" : "none";
                }
            });
    } else {
        // Reset X entries when no query
        container
            .querySelectorAll(
                ".anndata-entry__nested-content .anndata-x__entry",
            )
            .forEach((xEntry) => {
                xEntry.style.display = "";
            });
    }

    // Update filter indicator
    if (filterIndicator) {
        if (query) {
            filterIndicator.classList.add("anndata--active");
            filterIndicator.textContent = `Showing ${totalMatches} of ${totalEntries}`;
        } else {
            filterIndicator.classList.remove("anndata--active");
        }
    }

    // Hide sections with no visible entries
    container.querySelectorAll(".anndata-section").forEach((section) => {
        const visibleEntries = section.querySelectorAll(
            ".anndata-entry:not(.anndata-entry--hidden)",
        );

        if (query && visibleEntries.length === 0) {
            section.style.display = "none";
        } else {
            section.style.display = "";
        }
    });
}

// Copy to clipboard
container.querySelectorAll(".anndata-entry__copy").forEach((btn) => {
    btn.addEventListener("click", async (e) => {
        e.stopPropagation();

        const text = btn.dataset.copy;
        if (!text) return;

        try {
            await navigator.clipboard.writeText(text);

            // Visual feedback (icon turns green via CSS)
            btn.classList.add("anndata-entry__copy--copied");
            setTimeout(
                () => btn.classList.remove("anndata-entry__copy--copied"),
                1500,
            );
        } catch (err) {
            // Fallback for older browsers
            const textarea = document.createElement("textarea");
            textarea.value = text;
            textarea.style.position = "fixed";
            textarea.style.opacity = "0";
            document.body.appendChild(textarea);
            textarea.select();

            try {
                document.execCommand("copy");
                btn.classList.add("anndata-entry__copy--copied");
                setTimeout(
                    () => btn.classList.remove("anndata-entry__copy--copied"),
                    1500,
                );
            } catch (e) {
                console.error("Copy failed:", e);
            }

            document.body.removeChild(textarea);
        }
    });
});

// Expand/collapse nested content
container.querySelectorAll(".anndata-entry__expand").forEach((btn) => {
    btn.addEventListener("click", (e) => {
        e.stopPropagation();

        const entry = btn.closest(".anndata-entry");
        if (!entry) return;

        // The nested content is in a sibling <tr class="anndata-entry--nested">
        // which contains <td class="anndata-entry__nested-content">
        const nestedRow = entry.nextElementSibling;
        if (
            !nestedRow ||
            !nestedRow.classList.contains("anndata-entry--nested")
        )
            return;

        const nestedContent = nestedRow.querySelector(
            ".anndata-entry__nested-content",
        );
        if (!nestedContent) return;

        const isExpanded = nestedRow.classList.toggle(
            "anndata-entry--expanded",
        );

        btn.textContent = isExpanded ? "Collapse ▲" : "Expand ▼";
        btn.setAttribute("aria-expanded", isExpanded);
        nestedRow.setAttribute("aria-hidden", !isExpanded);
    });
});

// Helper to check if element is overflowing
function isOverflowing(el) {
    return el.scrollWidth > el.clientWidth;
}

// Helper to update wrap button visibility based on overflow
function updateWrapButtonVisibility(btn, list, metaCell, wrappedClass) {
    if (!list || !metaCell) {
        btn.style.display = "none";
        return;
    }
    // Show button only if content is overflowing or currently wrapped
    const isWrapped = list.classList.contains(wrappedClass);
    const overflows = isOverflowing(metaCell);
    btn.style.display = overflows || isWrapped ? "inline" : "none";
}

// Factory function to set up wrap button handlers (DRY pattern for cats/cols buttons)
function setupWrapButtons(buttonSelector, listSelector, wrappedClass) {
    container.querySelectorAll(buttonSelector).forEach((btn) => {
        const typeCell = btn.closest(".anndata-entry__type");
        const metaCell = typeCell ? typeCell.nextElementSibling : null;
        const list = metaCell ? metaCell.querySelector(listSelector) : null;

        // Initial visibility check
        updateWrapButtonVisibility(btn, list, metaCell, wrappedClass);

        btn.addEventListener("click", (e) => {
            e.stopPropagation();
            if (!list || !metaCell) return;

            const isWrapped = list.classList.toggle(wrappedClass);
            metaCell.classList.toggle("anndata-entry--expanded", isWrapped);
            btn.textContent = isWrapped ? "▲" : "▼";
            btn.title = isWrapped
                ? "Collapse to single line"
                : "Expand to multi-line view";
            // Always show button when wrapped
            btn.style.display = "inline";
        });
    });
}

// Set up wrap buttons for categories and columns lists
setupWrapButtons(
    ".anndata-categories__wrap",
    ".anndata-categories",
    "anndata-categories--wrapped",
);
setupWrapButtons(
    ".anndata-columns__wrap",
    ".anndata-columns",
    "anndata-columns--wrapped",
);

// Update button visibility on container resize (works for JupyterLab panes too)
// Uses the same selector pairs as setupWrapButtons for consistency
function updateAllWrapButtons() {
    [
        [
            ".anndata-categories__wrap",
            ".anndata-categories",
            "anndata-categories--wrapped",
        ],
        [
            ".anndata-columns__wrap",
            ".anndata-columns",
            "anndata-columns--wrapped",
        ],
    ].forEach(([btnSel, listSel, wrappedClass]) => {
        container.querySelectorAll(btnSel).forEach((btn) => {
            const typeCell = btn.closest(".anndata-entry__type");
            const metaCell = typeCell ? typeCell.nextElementSibling : null;
            const list = metaCell ? metaCell.querySelector(listSel) : null;
            updateWrapButtonVisibility(btn, list, metaCell, wrappedClass);
        });
    });
}

// Use ResizeObserver for robust resize detection (pane resizes, not just window)
if (typeof ResizeObserver !== "undefined") {
    let resizeTimer;
    const resizeObserver = new ResizeObserver(() => {
        clearTimeout(resizeTimer);
        resizeTimer = setTimeout(updateAllWrapButtons, 100);
    });
    resizeObserver.observe(container);
} else {
    // Fallback for older browsers
    let resizeTimer;
    window.addEventListener("resize", () => {
        clearTimeout(resizeTimer);
        resizeTimer = setTimeout(updateAllWrapButtons, 100);
    });
}

// README modal functionality
const readmeIcon = container.querySelector(".anndata-readme__icon");
if (readmeIcon) {
    // Ensure accessibility attributes
    readmeIcon.setAttribute("role", "button");
    readmeIcon.setAttribute("tabindex", "0");
    readmeIcon.setAttribute("aria-label", "View README");

    readmeIcon.addEventListener("click", (e) => {
        e.stopPropagation();
        const readmeContent = readmeIcon.dataset.readme;
        if (!readmeContent) return;

        // Create modal overlay
        const overlay = document.createElement("div");
        overlay.className = "anndata-readme__overlay";

        // Create modal with accessibility attributes
        // Use container.id to make IDs unique across multiple cells
        const modalTitleId = container.id + "-readme-modal-title";
        const modal = document.createElement("div");
        modal.className = "anndata-readme__modal";
        modal.setAttribute("role", "dialog");
        modal.setAttribute("aria-modal", "true");
        modal.setAttribute("aria-labelledby", modalTitleId);

        // Header
        const header = document.createElement("div");
        header.className = "anndata-readme__header";
        header.innerHTML = '<h3 id="' + modalTitleId + '">README</h3>';

        const closeBtn = document.createElement("button");
        closeBtn.className = "anndata-readme__close";
        closeBtn.textContent = "×";
        closeBtn.setAttribute("aria-label", "Close");
        header.appendChild(closeBtn);

        // Content — plain text (no markdown parsing, XSS-safe via textContent)
        const content = document.createElement("div");
        content.className = "anndata-readme__content";
        const pre = document.createElement("pre");
        pre.textContent = readmeContent;
        content.appendChild(pre);

        modal.appendChild(header);
        modal.appendChild(content);
        overlay.appendChild(modal);

        // Add to container (scoped styles apply)
        container.appendChild(overlay);

        // Close handlers
        const closeModal = () => {
            overlay.remove();
        };

        closeBtn.addEventListener("click", closeModal);
        overlay.addEventListener("click", (e) => {
            if (e.target === overlay) closeModal();
        });

        // Escape key closes modal
        const escHandler = (e) => {
            if (e.key === "Escape") {
                closeModal();
                document.removeEventListener("keydown", escHandler);
            }
        };
        document.addEventListener("keydown", escHandler);

        // Focus trap
        closeBtn.focus();
    });

    // Keyboard accessibility for the icon
    readmeIcon.addEventListener("keydown", (e) => {
        if (e.key === "Enter" || e.key === " ") {
            e.preventDefault();
            readmeIcon.click();
        }
    });
}
