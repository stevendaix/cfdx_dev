"""Sphinx configuration for the CFDX executable scientific documentation."""

from __future__ import annotations

project = "CFDX"
copyright = "CFDX contributors"
author = "CFDX contributors"
release = "0.7"

extensions = [
    "myst_parser",
    "myst_nb",
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",
]

templates_path = []
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
    ".py": "jupyter_notebook",
}

myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    "colon_fence",
    "deflist",
]

nb_custom_formats = {
    ".py": ["jupytext.reads", {"fmt": "py:percent"}],
}
nb_execution_mode = "off"

bibtex_bibfiles = ["references/bibliography.bib"]
bibtex_default_style = "unsrt"

html_theme = "pydata_sphinx_theme"
html_title = "CFDX Documentation"
html_theme_options = {"show_toc_level": 2}
