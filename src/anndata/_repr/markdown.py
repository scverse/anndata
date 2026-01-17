"""
Markdown parsing for README rendering in HTML repr.

Current Implementation
----------------------
README markdown is currently parsed client-side using JavaScript
(see static/markdown-parser.js). This avoids adding Python dependencies
but has limitations (no nested lists, limited table support, etc.).

Server-Side Python Alternatives
-------------------------------
To replace the JavaScript parser with Python server-side rendering,
consider these libraries (in order of recommendation):

1. **markdown-it-py** (recommended)
   - Port of markdown-it, CommonMark compliant, extensible
   - Install: ``pip install markdown-it-py``
   - Usage::

       from markdown_it import MarkdownIt

       md = MarkdownIt("commonmark", {"html": False})  # Disable raw HTML for security
       html = md.render(readme_text)

   - Pros: Fast, CommonMark compliant, plugin ecosystem
   - Docs: https://markdown-it-py.readthedocs.io/

2. **mistune**
   - Fast, pure Python, customizable renderer
   - Install: ``pip install mistune``
   - Usage::

       import mistune

       html = mistune.html(readme_text)

   - Pros: Very fast, no dependencies, easy customization
   - Docs: https://mistune.lepture.com/

3. **markdown** (Python-Markdown)
   - Original Python markdown library, many extensions
   - Install: ``pip install markdown``
   - Usage::

       import markdown

       md = markdown.Markdown(extensions=["tables", "fenced_code"])
       html = md.convert(readme_text)

   - Pros: Mature, extensive extension library
   - Docs: https://python-markdown.github.io/

4. **cmarkgfm**
   - GitHub Flavored Markdown via C bindings (fastest)
   - Install: ``pip install cmarkgfm``
   - Usage::

       import cmarkgfm

       html = cmarkgfm.github_flavored_markdown_to_html(readme_text)

   - Pros: Fastest, GitHub-compatible
   - Cons: Requires C compiler for installation
   - Docs: https://github.com/theacodes/cmarkgfm

Integration Steps
-----------------
To switch to server-side Python parsing:

1. Add the chosen library as an optional dependency in pyproject.toml::

       [project.optional - dependencies]
       markdown = ["markdown-it-py>=3.0"]

2. Update ``_repr/sections.py`` or ``_repr/html.py`` to:
   - Parse README content server-side before embedding in HTML
   - Remove the parseMarkdown() JavaScript call from the README modal
   - Sanitize output HTML to prevent XSS (use bleach or html-sanitizer)

3. Update ``static/repr.js`` to display pre-rendered HTML directly
   instead of calling the client-side parser

4. Remove ``static/markdown-parser.js`` once migration is complete

Security Considerations
-----------------------
README content often comes from h5ad/zarr files shared by others, making
proper sanitization essential to prevent XSS attacks.

Required security measures:

1. **HTML sanitization with nh3** (recommended) or bleach::

       import nh3

       # After markdown parsing
       safe_html = nh3.clean(
           raw_html,
           tags={
               "p",
               "h1",
               "h2",
               "h3",
               "h4",
               "a",
               "code",
               "pre",
               "strong",
               "em",
               "ul",
               "ol",
               "li",
               "blockquote",
               "table",
               "thead",
               "tbody",
               "tr",
               "th",
               "td",
           },
           attributes={"a": {"href", "rel", "target"}},
           link_rel="noopener noreferrer",
       )

2. **Disable raw HTML passthrough** in the markdown parser::

       md = MarkdownIt("commonmark", {"html": False})

3. **Sanitize link URLs** to block javascript: and data: protocols::

       url_rel = "noopener noreferrer"  # nh3 handles this
       # Or manually check: if not url.startswith(("http://", "https://", "/")):

Add nh3 as a dependency in pyproject.toml::

    [project.optional - dependencies]
    markdown = ["markdown-it-py>=3.0", "nh3>=0.2"]
"""

from __future__ import annotations
