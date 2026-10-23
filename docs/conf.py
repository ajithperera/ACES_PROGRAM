# Sphinx configuration for the ACES program system manual.
import os
import subprocess
import sys

DOCS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(DOCS)

project = "ACES"
author = "Quantum Theory Project, University of Florida"
copyright = "1990-2026, the ACES authors"

extensions = ["myst_parser"]
myst_enable_extensions = ["colon_fence", "deflist", "dollarmath", "attrs_inline"]
myst_heading_anchors = 3

source_suffix = {".md": "markdown"}
exclude_patterns = ["_build", "reference/keyword_notes.md", "tools", "README.md"]

html_theme = "furo"
html_title = "ACES program system"

latex_documents = [("index", "aces_manual.tex", "The ACES Program System",
                    author, "manual")]


def generate_keywords(app):
    """Regenerate reference/keywords.md from aces2/joda/asv_ctl.c."""
    subprocess.run(
        [sys.executable, os.path.join(DOCS, "tools", "gen_keywords.py"),
         os.path.join(ROOT, "aces2", "joda", "asv_ctl.c"),
         os.path.join(DOCS, "reference", "keyword_notes.md"),
         os.path.join(DOCS, "reference", "keywords.md")],
        check=True)


def setup(app):
    app.connect("builder-inited", generate_keywords)
