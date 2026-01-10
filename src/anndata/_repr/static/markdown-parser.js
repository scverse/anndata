// Simple markdown parser for README rendering.
// This is a minimal implementation to avoid external dependencies.
// It can be replaced with a full-featured library like marked.js or markdown-it.
// See src/anndata/_repr/markdown.py for replacement instructions.
function parseMarkdown(text) {
    // Decode HTML entities first (the content was escaped)
    const textarea = document.createElement("textarea");
    textarea.innerHTML = text;
    text = textarea.value;

    // Escape HTML to prevent XSS, then apply markdown
    // SECURITY: Must escape ALL characters that could break HTML context
    text = text
        .replace(/&/g, "&amp;")
        .replace(/</g, "&lt;")
        .replace(/>/g, "&gt;")
        .replace(/"/g, "&quot;")
        .replace(/'/g, "&#39;");

    // Tables: |col1|col2| format with |---|---| separator
    text = text.replace(/(^\|.+\|\n?)+/gm, function (tableBlock) {
        const rows = tableBlock
            .trim()
            .split("\n")
            .filter((r) => r.trim());
        if (rows.length < 2) return tableBlock;

        // Find separator row index (contains only |, -, :, and spaces)
        const sepIndex = rows.findIndex((r) => /^\|[\s:\-|]+\|$/.test(r));
        if (sepIndex < 1) return tableBlock; // No valid separator found

        let html = "<table>";

        // Header rows (before separator)
        html += "<thead>";
        for (let i = 0; i < sepIndex; i++) {
            const cells = rows[i].split("|").slice(1, -1);
            html +=
                "<tr>" +
                cells.map((c) => "<th>" + c.trim() + "</th>").join("") +
                "</tr>";
        }
        html += "</thead>";

        // Body rows (after separator)
        if (sepIndex < rows.length - 1) {
            html += "<tbody>";
            for (let i = sepIndex + 1; i < rows.length; i++) {
                const cells = rows[i].split("|").slice(1, -1);
                html +=
                    "<tr>" +
                    cells.map((c) => "<td>" + c.trim() + "</td>").join("") +
                    "</tr>";
            }
            html += "</tbody>";
        }

        html += "</table>";
        return html;
    });

    // Extract code blocks and inline code first to protect from formatting
    const codeBlocks = [];
    const inlineCodes = [];

    // Code blocks (``` ... ```) - extract and replace with placeholder
    text = text.replace(/```(\w*)\n([\s\S]*?)```/g, (match, lang, code) => {
        codeBlocks.push("<pre><code>" + code + "</code></pre>");
        return "%%CODEBLOCK" + (codeBlocks.length - 1) + "%%";
    });

    // Inline code (`...`) - extract and replace with placeholder
    text = text.replace(/`([^`]+)`/g, (match, code) => {
        inlineCodes.push("<code>" + code + "</code>");
        return "%%INLINECODE" + (inlineCodes.length - 1) + "%%";
    });

    // Links [text](url) - extract to protect URLs from formatting
    const links = [];
    text = text.replace(/\[([^\]]+)\]\(([^)]+)\)/g, (match, linkText, url) => {
        links.push({ text: linkText, url: url });
        return "%%LINK" + (links.length - 1) + "%%";
    });

    // Auto-link raw URLs (https://... or http://...)
    text = text.replace(/(https?:\/\/[^\s<>)"']+)/g, (match, url) => {
        // Don't match if already inside a placeholder
        if (url.includes("%%")) return match;
        links.push({ text: url, url: url });
        return "%%LINK" + (links.length - 1) + "%%";
    });

    // Headers
    text = text.replace(/^#### (.+)$/gm, "<h4>$1</h4>");
    text = text.replace(/^### (.+)$/gm, "<h3>$1</h3>");
    text = text.replace(/^## (.+)$/gm, "<h2>$1</h2>");
    text = text.replace(/^# (.+)$/gm, "<h1>$1</h1>");

    // Bold and italic
    text = text.replace(/\*\*\*(.+?)\*\*\*/g, "<strong><em>$1</em></strong>");
    text = text.replace(/\*\*(.+?)\*\*/g, "<strong>$1</strong>");
    text = text.replace(/\*(.+?)\*/g, "<em>$1</em>");
    text = text.replace(/_(.+?)_/g, "<em>$1</em>");

    // Blockquotes
    text = text.replace(/^&gt; (.+)$/gm, "<blockquote>$1</blockquote>");

    // Unordered lists - use temporary marker to distinguish from ordered
    text = text.replace(/^[\*\-] (.+)$/gm, "<ul-li>$1</ul-li>");

    // Ordered lists - use temporary marker
    text = text.replace(/^\d+\. (.+)$/gm, "<ol-li>$1</ol-li>");

    // Wrap consecutive unordered items in <ul>
    text = text.replace(/(<ul-li>.*<\/ul-li>\n?)+/g, "<ul>$&</ul>");
    text = text.replace(/<ul-li>/g, "<li>");
    text = text.replace(/<\/ul-li>/g, "</li>");

    // Wrap consecutive ordered items in <ol>
    text = text.replace(/(<ol-li>.*<\/ol-li>\n?)+/g, "<ol>$&</ol>");
    text = text.replace(/<ol-li>/g, "<li>");
    text = text.replace(/<\/ol-li>/g, "</li>");

    // Paragraphs (double newlines)
    text = text.replace(/\n\n+/g, "</p><p>");
    text = "<p>" + text + "</p>";

    // Clean up empty paragraphs
    text = text.replace(/<p><\/p>/g, "");
    text = text.replace(/<p>(<h[1-4]>)/g, "$1");
    text = text.replace(/(<\/h[1-4]>)<\/p>/g, "$1");
    text = text.replace(/<p>(<ul>)/g, "$1");
    text = text.replace(/(<\/ul>)<\/p>/g, "$1");
    text = text.replace(/<p>(<ol>)/g, "$1");
    text = text.replace(/(<\/ol>)<\/p>/g, "$1");
    text = text.replace(/<p>(<pre>)/g, "$1");
    text = text.replace(/(<\/pre>)<\/p>/g, "$1");
    text = text.replace(/<p>(<blockquote>)/g, "$1");
    text = text.replace(/(<\/blockquote>)<\/p>/g, "$1");
    text = text.replace(/<p>(<table>)/g, "$1");
    text = text.replace(/(<\/table>)<\/p>/g, "$1");

    // Single newlines become spaces (standard markdown behavior)
    // Text flows together unless there's a double newline (paragraph break)
    text = text.replace(/([^>])\n([^<])/g, "$1 $2");

    // Restore links (apply formatting to link text, but URL is protected)
    links.forEach((link, i) => {
        // Apply formatting to link text only
        let formattedText = link.text;
        formattedText = formattedText.replace(
            /\*\*(.+?)\*\*/g,
            "<strong>$1</strong>",
        );
        formattedText = formattedText.replace(/\*(.+?)\*/g, "<em>$1</em>");
        // Sanitize URL: only allow safe protocols (http, https, mailto, relative)
        let safeUrl = link.url;
        if (!/^(https?:|mailto:|\/|#)/i.test(safeUrl)) {
            safeUrl = "#"; // Block unsafe protocols
        }
        const linkHtml =
            '<a href="' +
            safeUrl +
            '" target="_blank" rel="noopener">' +
            formattedText +
            "</a>";
        text = text.replace("%%LINK" + i + "%%", linkHtml);
    });

    // Restore code blocks and inline code
    codeBlocks.forEach((block, i) => {
        text = text.replace("%%CODEBLOCK" + i + "%%", block);
    });
    inlineCodes.forEach((code, i) => {
        text = text.replace("%%INLINECODE" + i + "%%", code);
    });

    return text;
}
