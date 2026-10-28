# ACES manual sources

Sphinx + MyST Markdown. To build the HTML manual:

    python3 -m pip install --user -r requirements.txt
    ./build.sh            # output in _build/html/index.html

`reference/keywords.md` is generated at build time by `tools/gen_keywords.py`
from `aces2/joda/asv_ctl.c`; edit keyword descriptions in
`reference/keyword_notes.md` instead.
