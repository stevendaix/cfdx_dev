"""Sphinx configuration for the CFDX executable scientific documentation."""

from __future__ import annotations

project = "CFDX"
copyright = "CFDX contributors"
author = "CFDX contributors"
release = "0.7"

extensions = [
    "myst_nb",
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",
    "sphinx.ext.githubpages",
]

templates_path = []

# Publish only the maintained documentation model. Legacy migration material and
# executable support scripts remain in the repository but are not Sphinx sources.
include_patterns = [
    "index.md",
    "theory/README.md",
    "theory/*/README.md",
    "theory/*/chapter.py",
    "theory/04-gradients-reconstruction/*.md",
    "user/README.md",
    "developer/README.md",
    "developer/*/chapter.py",
    "vv/README.md",
    "vv/*/chapter.py",
    "validation/README.md",
    "validation/14-benchmarks/README.md",
    "validation/14-benchmarks/*/README.md",
    "validation/14-benchmarks/*/chapter.py",
    "validation/15-qualification/README.md",
    "references/README.md",
]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# MyST-NB parses each notebook Markdown cell independently. Section cells that
# intentionally begin at H2 therefore trigger the generic document-level
# heading warning even when the notebook has a valid H1 title cell.
suppress_warnings = ["myst.header"]

myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    "colon_fence",
    "deflist",
]

nb_custom_formats = {
    ".py": ["jupytext.reads", {"fmt": "py:percent"}],
}
nb_execution_mode = "auto"
nb_execution_timeout = 120

bibtex_bibfiles = ["references/bibliography.bib"]
bibtex_default_style = "unsrt"

html_theme = "pydata_sphinx_theme"
html_title = "CFDX Documentation"
html_theme_options = {"show_toc_level": 2}

# Pin MathJax v3 for deterministic rendering across the documentation toolchain.
mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"
