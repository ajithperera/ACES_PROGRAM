#!/usr/bin/env python3
"""One-time importer: pulls cases from the old A2_test_suite archive into
this suite's cases/<name>/ZMAT + test_results, using each case's own
embedded TEST.DAT TOTENERG value as the reference. Not part of the regular
test-running workflow -- run once, inspect the result, done.
"""
import os
import re
import sys

ARCHIVE = "/blue/bartlett/perera/A2_test_suite"
FAST_DIR = os.path.dirname(os.path.abspath(__file__))
CASES_DIR = os.path.join(FAST_DIR, "cases")
TEST_RESULTS = os.path.join(FAST_DIR, "test_results")

SKIP_DUPLICATES = {"zmat.073.mrcc"}  # byte-identical to zmat.072.mrcc

TOTENERG_RE = re.compile(r"^\s*d\s+TOTENERG\b", re.IGNORECASE)
VALUE_RE = re.compile(r"[-+]?\d*\.\d+[eEdD][-+]?\d+|[-+]?\d+\.\d+")


def find_totenerg(lines):
    for i, line in enumerate(lines):
        if TOTENERG_RE.match(line):
            for j in range(i + 1, len(lines)):
                if lines[j].strip():
                    m = VALUE_RE.search(lines[j])
                    if m:
                        return float(m.group(0).replace("D", "E").replace("d", "e"))
                    return None
    return None


def main():
    entries = []
    for fname in sorted(os.listdir(ARCHIVE)):
        if not fname.startswith("zmat."):
            continue
        if fname in SKIP_DUPLICATES:
            continue
        path = os.path.join(ARCHIVE, fname)
        if not os.path.isfile(path):
            continue
        text = open(path, errors="replace").read()
        if "*ACES2" not in text or "TEST.DAT" not in text or "*CRAPS" in text:
            continue
        lines = text.splitlines(keepends=True)
        split_idx = None
        for i, line in enumerate(lines):
            if line.strip() == "TEST.DAT":
                split_idx = i
                break
        if split_idx is None:
            continue
        zmat_body = "".join(lines[:split_idx]).rstrip() + "\n\n"
        totenerg = find_totenerg(lines[split_idx:])
        if totenerg is None:
            print(f"SKIP {fname}: no TOTENERG found in TEST.DAT block")
            continue
        case_name = fname[len("zmat."):]
        entries.append((case_name, zmat_body, totenerg, fname))

    os.makedirs(CASES_DIR, exist_ok=True)
    with open(TEST_RESULTS, "w") as tr:
        for case_name, zmat_body, totenerg, src in entries:
            case_dir = os.path.join(CASES_DIR, case_name)
            os.makedirs(case_dir, exist_ok=True)
            with open(os.path.join(case_dir, "ZMAT"), "w") as f:
                f.write(zmat_body)
            tr.write(f"{case_name} {src}\n")
            tr.write("TOTENERG 1\n")
            tr.write(f"{totenerg:.12E}\n")

    print(f"Imported {len(entries)} cases into {CASES_DIR}")
    print(f"Wrote references to {TEST_RESULTS}")


if __name__ == "__main__":
    main()
