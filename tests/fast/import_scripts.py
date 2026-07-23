#!/usr/bin/env python3
"""One-time importer for the archive's script.* cases -- these chain two
(or more) ACES2 jobs together in one shell script (e.g. run job A, save an
intermediate file like FCMINT outside the working dir, clean up, run job B
with !restart using the saved file). Unlike the zmat.* cases, these keep
their original shell-script logic verbatim (no need to reimplement the
chaining in Python) -- just strip the TEST.DAT trailer and record the
reference TOTENERG, same as import_archive.py does for ZMAT cases.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from import_archive import ARCHIVE, CASES_DIR, TEST_RESULTS, find_totenerg  # noqa: E402


def main():
    entries = []
    for fname in sorted(os.listdir(ARCHIVE)):
        if not fname.startswith("script."):
            continue
        path = os.path.join(ARCHIVE, fname)
        if not os.path.isfile(path):
            continue
        text = open(path, errors="replace").read()
        if "TEST.DAT" not in text:
            print(f"SKIP {fname}: no TEST.DAT block")
            continue
        lines = text.splitlines(keepends=True)
        split_idx = None
        for i, line in enumerate(lines):
            if line.strip() == "TEST.DAT":
                split_idx = i
                break
        if split_idx is None:
            continue
        script_body = "".join(lines[:split_idx]).rstrip() + "\n"
        totenerg = find_totenerg(lines[split_idx:])
        if totenerg is None:
            print(f"SKIP {fname}: no TOTENERG found in TEST.DAT block")
            continue
        case_name = "s_" + fname[len("script."):]
        entries.append((case_name, script_body, totenerg, fname))

    for case_name, _, _, _ in entries:
        case_dir = os.path.join(CASES_DIR, case_name)
        if os.path.isdir(case_dir):
            sys.exit(f"case name collision: {case_dir} already exists")

    for case_name, script_body, totenerg, src in entries:
        case_dir = os.path.join(CASES_DIR, case_name)
        os.makedirs(case_dir, exist_ok=True)
        with open(os.path.join(case_dir, "SCRIPT"), "w") as f:
            f.write(script_body)
        with open(TEST_RESULTS, "a") as tr:
            tr.write(f"{case_name} {src}\n")
            tr.write("TOTENERG 1\n")
            tr.write(f"{totenerg:.12E}\n")

    print(f"Imported {len(entries)} script cases into {CASES_DIR}")


if __name__ == "__main__":
    main()
