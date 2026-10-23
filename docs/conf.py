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
exclude_patterns = ["_build", "reference/keyword_notes.md", "programmer/jobarc_notes.md",
                    "user/aces3_jobflows.md",
                    "tools", "README.md"]

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


def generate_jobarc(app):
    """Regenerate programmer/jobarc.md from the PUTREC/GETREC calls."""
    subprocess.run(
        [sys.executable, os.path.join(DOCS, "tools", "gen_jobarc.py"),
         os.path.join(ROOT, "aces2"),
         os.path.join(DOCS, "programmer", "jobarc_notes.md"),
         os.path.join(DOCS, "programmer", "jobarc.md")],
        check=True)


def generate_jobflows(app):
    """Regenerate user/aces3_jobflows.md from default_jobflows."""
    subprocess.run(
        [sys.executable, os.path.join(DOCS, "tools", "gen_jobflows.py"),
         os.path.join(ROOT, "aces3"),
         os.path.join(DOCS, "user", "aces3_jobflows.md")],
        check=True)


def setup(app):
    app.connect("builder-inited", generate_jobflows)
    app.connect("builder-inited", generate_keywords)
    app.connect("builder-inited", generate_jobarc)

# Published as a private web page: no copied sources, MathJax pinned to an exact version.
html_copy_source = False
html_show_sourcelink = False
mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-mml-chtml.js"
