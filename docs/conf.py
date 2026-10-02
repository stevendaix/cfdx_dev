"""Sphinx configuration for the CFDX scientific documentation."""

from __future__ import annotations

project = "CFDX"
copyright = "CFDX contributors"
author = "CFDX contributors"
release = "0.7"

extensions = [
    "myst_parser",
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",
]

templates_path = []
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}

myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    "colon_fence",
    "deflist",
]

bibtex_bibfiles = ["references/bibliography.bib"]
bibtex_default_style = "unsrt"

html_theme = "pydata_sphinx_theme"
html_title = "CFDX Documentation"
html_theme_options = {
    "show_toc_level": 2,
}
