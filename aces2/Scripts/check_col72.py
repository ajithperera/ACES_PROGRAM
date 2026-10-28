#!/usr/bin/env python3
"""List fixed-form source lines whose code runs past column 72.

The compiler silently drops everything after column 72 of a fixed-form
line, so a lost comma or parenthesis turns into a different program (two
arguments fused into one name, an argument list one short, ...). This
checks the preprocessed .f files that are actually compiled -- macros can
shorten or lengthen a .F line -- in every aces2 module directory (one
with a GNUmakefile). Comment lines, tab-format lines and text after an
inline ! comment are ignored. Exit status 1 when anything is found.

Usage: python3 aces2/Scripts/check_col72.py [ACES2_DIR]
"""
import glob
import os
import sys

root = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir)
found = 0
for d in sorted(os.listdir(root)):
    if d == "tests" or not os.path.exists(os.path.join(root, d, "GNUmakefile")):
        continue
    for src in sorted(glob.glob(os.path.join(root, d, "*.F"))):
        f = src[:-2] + ".f"
        if not os.path.exists(f):
            continue
        with open(f, errors="replace") as fh:
            for n, line in enumerate(fh.read().split("\n"), 1):
                if len(line) <= 72 or line[:1] in "cC*!dD" or "\t" in line[:6]:
                    continue
                if not line[72:].strip():
                    continue
                quote = None
                comment = False
                for ch in line[6:72]:
                    if quote is None and ch in "'\"":
                        quote = ch
                    elif ch == quote:
                        quote = None
                    elif ch == "!" and quote is None:
                        comment = True
                        break
                if comment:
                    continue
                print(f"{os.path.relpath(f, root)}:{n}: lost after column 72: "
                      f"{line[72:].strip()!r}")
                found += 1
print(f"check_col72: {found} line(s) with code past column 72")
sys.exit(1 if found else 0)
