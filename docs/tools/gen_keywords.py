#!/usr/bin/env python3
"""Generate the ACES II keyword reference from the source code.

The keyword names, types, defaults, units and allowed values come from
aces2/joda/asv_ctl.c (the ASV_nl[] table and asv_handle_proc), so the
reference cannot drift from what xjoda actually parses.  Descriptions come
from the hand-maintained docs/reference/keyword_notes.md, which groups
keywords by category under '# Category' headings and describes each one
under a '## KEYWORD' heading.

Usage: gen_keywords.py ASV_CTL_C NOTES_MD OUTPUT_MD
"""
import re
import sys

TYPE_NAMES = {
    "handle": "handle",
    "f_int": "integer",
    "string": "string",
    "double": "real",
    "f_int_array": "integer array",
    "double_array": "real array",
}

SWITCH_VALUES = ["OFF", "ON"]


def strip_c_comments(text):
    return re.sub(r"/\*.*?\*/", "", text, flags=re.S)


def parse_table(src):
    """Return the ASV_nl[] entries in table order (Fortran index = position)."""
    start = src.index("const ASV_nl_t ASV_nl[]")
    end = src.index("end ASV_nl[] definition", start)
    entry = re.compile(
        r'\{\s*/\*\s*(\d+)\s*\*/\s*"([^"]*)"\s*,\s*"([^"]*)"\s*,\s*'
        r'(h_ICHAR_\w+|0)\s*,\s*(-?\d+)\s*,\s*"([^"]*)"\s*\}')
    keys = []
    for m in entry.finditer(src[start:end]):
        index, alias, oppar, ichar, default, units = m.groups()
        if not oppar:
            continue
        keys.append(dict(index=int(index), alias=alias, oppar=oppar,
                         type=ichar.replace("h_ICHAR_", ""),
                         default=int(default), units=units))
    return keys


def parse_handles(src):
    """Map lower-case h_IOPPAR_ names to their allowed handle strings."""
    start = src.index("void asv_handle_proc(")
    body = src[start:]
    body = body[:body.index("end switch")]
    body = strip_c_comments(body)
    lists = {}
    pending = []
    for block in re.split(r"\bbreak\s*;", body):
        cases = re.findall(r"case\s+h_IOPPAR_(\w+)\s*:", block)
        m = re.search(r"handles\[\]\s*=\s*\{(.*?)\}\s*;", block, flags=re.S)
        if not cases:
            continue
        if m is None:
            pending.extend(cases)
            continue
        # adjacent literals ("A" "") are one string to the C compiler
        values = re.findall(r'"([^"]*)"', re.sub(r'"\s+"', "", m.group(1)))
        terminated = bool(values) and values[-1] == ""
        values = [v for v in values if v]
        for c in pending + cases:
            lists[c.lower()] = (values, terminated)
        pending = []
    return lists


def parse_flags(path):
    """Map Fortran ioppar index -> h_IOPPAR_ names (first name is primary)."""
    names = {}
    with open(path) as f:
        for line in f:
            m = re.match(r"#define\s+h_IOPPAR_(\w+)\s+(\d+)", line)
            if m:
                names.setdefault(int(m.group(2)), []).append(m.group(1).lower())
    return names


def parse_notes(path):
    """Return (ordered categories, {KEYWORD: (category, text)})."""
    categories, notes = [], {}
    category, key, buf = None, None, []

    def flush():
        if key is not None:
            notes[key] = (category, "\n".join(buf).strip())

    with open(path) as f:
        for line in f:
            if line.startswith("## "):
                flush()
                key, buf = line[3:].strip().upper(), []
            elif line.startswith("# "):
                flush()
                key, buf = None, []
                category = line[2:].strip()
                categories.append(category)
            elif key is not None:
                buf.append(line.rstrip("\n"))
    flush()
    return categories, notes


def keyword_name(oppar):
    return oppar.replace("#", "")


def min_abbrev(oppar):
    return oppar.split("#")[0]


def describe(key, handles, flag_names):
    names = flag_names.get(key["index"], [])
    values, terminated = None, True
    for n in names:
        if n in handles:
            values, terminated = handles[n]
            break
    if key["type"] == "handle" and values is None:
        values = SWITCH_VALUES
        kind = "switch"
    else:
        kind = TYPE_NAMES.get(key["type"], key["type"])
    default = str(key["default"])
    if values and 0 <= key["default"] < len(values):
        default = values[key["default"]]
    if key["units"] == "(tol)":
        kind = "tolerance"
        default = f"{key['default']} (10<sup>-{key['default']}</sup>)"
    elif key["units"]:
        default = f"{default} {key['units']}"
    return kind, default, values, terminated


def anchor(name):
    return "kw-" + re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")


def main(asv_ctl, notes_md, out_md):
    with open(asv_ctl) as f:
        src = f.read()
    flags_h = re.sub(r"joda/asv_ctl\.c$", "include/flags.h", asv_ctl)
    keys = parse_table(src)
    handles = parse_handles(src)
    flag_names = parse_flags(flags_h)
    categories, notes = parse_notes(notes_md)

    by_cat = {c: [] for c in categories}
    undocumented = []
    for k in keys:
        name = keyword_name(k["oppar"])
        if name in notes:
            by_cat[notes[name][0]].append(k)
        else:
            undocumented.append(k)
    code_names = {keyword_name(k["oppar"]) for k in keys}
    stale = sorted(n for n in notes if n not in code_names)

    out = []
    w = out.append
    w("<!-- GENERATED by docs/tools/gen_keywords.py -- do not edit. -->")
    w("<!-- Edit docs/reference/keyword_notes.md or aces2/joda/asv_ctl.c. -->\n")
    w("(keyword-reference)=")
    w("# ACES II keyword reference\n")
    w("Keywords go in the `*ACES2(...)` namelist of the ZMAT file "
      "(see {ref}`zmat-namelist`).")
    w("Names, types, defaults and allowed values on this page are generated "
      "from `aces2/joda/asv_ctl.c`, the table `xjoda` uses to parse the "
      "namelist, so they match the code. Descriptions are maintained in "
      "`docs/reference/keyword_notes.md`.\n")
    w("A keyword can be shortened to any prefix that is at least as long as "
      "its **minimum abbreviation** (for example `CALC` for `CALCLEVEL`). "
      "A *switch* given without a value is set to `ON`, and a leading `!` "
      "sets it to `OFF`. A *tolerance* N means 10<sup>-N</sup>. A *handle* "
      "takes one of the listed names (case-insensitive) or its position in "
      "the list counted from 0.\n")
    w(f"The code defines **{len(keys)}** keywords; "
      f"{len(keys) - len(undocumented)} are described below by category and "
      f"{len(undocumented)} are listed at the end without descriptions yet.\n")

    def entry(k):
        name = keyword_name(k["oppar"])
        kind, default, values, terminated = describe(k, handles, flag_names)
        w(f"({anchor(name)})=")
        w(f"### {name}\n")
        w(f"Minimum abbreviation `{min_abbrev(k['oppar'])}` · type *{kind}* · "
          f"default `{default}` · internal index {k['index']}\n")
        if values and kind != "switch":
            w("Values: " + ", ".join(f"`{v}`" for v in values) + "\n")
        if not terminated:
            w("```{warning}\nThe value list for this keyword in `asv_ctl.c` "
              "has no terminating empty string, so an unrecognized or numeric "
              "value reads past the end of the list.\n```\n")
        text = notes.get(name, (None, ""))[1]
        if text:
            w(text + "\n")

    for c in categories:
        if not by_cat[c]:
            continue
        w(f"## {c}\n")
        for k in by_cat[c]:
            entry(k)
    if undocumented:
        w("## Not yet described\n")
        w("These keywords exist in the code but have no description yet "
          "(most were added after the 2006 manual).\n")
        for k in undocumented:
            entry(k)
    if stale:
        w("## Keywords in the notes but not in the code\n")
        w("These appear in `keyword_notes.md` but `asv_ctl.c` no longer "
          "defines them (removed, renamed, or a separate namelist):\n")
        w(", ".join(f"`{n}`" for n in stale) + "\n")

    with open(out_md, "w") as f:
        f.write("\n".join(out))
    print(f"gen_keywords: {len(keys)} keywords, "
          f"{len(undocumented)} undocumented, {len(stale)} stale notes")


if __name__ == "__main__":
    main(*sys.argv[1:4])
