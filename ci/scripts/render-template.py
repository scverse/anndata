#!/usr/bin/env python3
# /// script
# dependencies = ["sphinx", "jinja2"]
# ///

from __future__ import annotations

import jinja2
from sphinx.ext.autosummary import generate

loader = jinja2.FileSystemLoader("docs/_templates/autosummary")
env = jinja2.Environment(loader=loader)
env.filters.update(underline=generate._underline)
t = env.get_template("class.rst")

print(t.render(fullname="a.Test", name="Test", attributes=["a", "b"], methods=["c"]))
