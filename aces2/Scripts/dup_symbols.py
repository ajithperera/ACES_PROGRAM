#!/usr/bin/env python3
"""Symbols defined by more than one member of the same aces2 archive.

GNU ld takes each archive member it needs in one pass, front to back. A
symbol that first becomes undefined while the pass is already past one
definition is resolved by the NEXT member that defines it, so with two
definitions in one archive the member that wins depends on archive order
and on which member first referenced the symbol. If the two copies differ
(arguments, COMMON blocks, behaviour) the binary can silently call the
wrong one -- see docs/programmer/duplicate_symbols.md.

For every duplicate this prints the archive members with their position,
source file and routine header, whether the bodies are identical, and the
files that call it. With --trace it also relinks every aces2 binary into a
scratch file with -Wl,--trace-symbol and reports which member each binary
actually gets (run it with the toolchain loaded, e.g. after
`source build_config.sh; aces_program_load_env`).

usage: aces2/Scripts/dup_symbols.py [--trace] [--machname NAME] [--json]
"""
import argparse, json, os, re, subprocess, tempfile
from collections import defaultdict

ACES2 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIB = os.path.join(ACES2, "lib")


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def archive_members(archive):
    return {m: i + 1 for i, m in enumerate(run(["ar", "t", archive]).stdout.split())}


def defined_symbols(archive):
    for line in run(["nm", "-A", "--defined-only", archive]).stdout.splitlines():
        m = re.match(r"[^:]+:([^:]+):\S*\s+([A-Za-z])\s+(\S+)$", line)
        if m and m.group(2) in "TDR":          # code and initialised data
            yield m.group(1), m.group(3)


def find_source(mdir, obj):
    for ext in (".F", ".f", ".c", ".cpp", ".F90", ".f90"):
        p = os.path.join(mdir, obj[:-2] + ext)
        if os.path.exists(p):
            return p
    return None


def routine(src, sym):
    """(header, normalised body) of Fortran routine `sym` in `src`."""
    name = None if sym == "MAIN__" else sym.rstrip("_").lower()
    lines = open(src, errors="replace").read().splitlines()
    hdr_re = re.compile(r"^\s{6,}.*?\b(subroutine|function|program)\s+(\w+)", re.I)
    for i, l in enumerate(lines):
        m = hdr_re.match(l) if l[:1] not in "cC*!" else None
        if m and (m.group(1).lower() == "program" if name is None
                  else m.group(2).lower() == name):
            break
    else:
        return None, None
    hdr, j = lines[i].strip(), i + 1
    while j < len(lines) and len(lines[j]) > 5 and lines[j][5] not in " 0" \
            and lines[j][:1] not in "cC*!":
        hdr += " " + lines[j][6:].strip()
        j += 1
    end = next((k + 1 for k in range(i + 1, len(lines))
                if re.match(r"^\s+end\s*$", lines[k], re.I)), len(lines))
    body = [re.sub(r"\s+", " ", l.strip()).lower() for l in lines[i:end]
            if l.strip() and l[:1] not in "cC*!"]
    return re.sub(r"\s+", " ", hdr), "\n".join(body)


def callers(sym):
    if sym == "MAIN__":
        return []
    name = sym.rstrip("_")
    files = run(["git", "grep", "-l", "-i", "-E",
                 r"call +%s *(\(|$)|= *%s *\(" % (name, name), "--",
                 "*.F", "*.f", "*.F90", "*.f90"], cwd=ACES2).stdout.split()
    # a tracked .f next to its .F is generated from it
    return sorted(f for f in files if not (
        f.endswith(".f") and os.path.exists(os.path.join(ACES2, f[:-2] + ".F"))))


def find_duplicates():
    dups = []
    for a in sorted(os.listdir(LIB)):
        if not a.endswith(".a"):
            continue
        path = os.path.join(LIB, a)
        order = archive_members(path)
        defs = defaultdict(set)
        for obj, sym in defined_symbols(path):
            defs[sym].add(obj)
        mdir = os.path.join(ACES2, a[3:-2])
        for sym, objs in sorted(defs.items()):
            if len(objs) < 2:
                continue
            members, bodies = [], set()
            for o in sorted(objs, key=lambda o: order.get(o, 0)):
                src = find_source(mdir, o) if os.path.isdir(mdir) else None
                hdr, body = routine(src, sym) if src else (None, None)
                if body is not None:
                    bodies.add(body)
                members.append({"member": o, "position": order.get(o),
                                "source": os.path.relpath(src, ACES2) if src else None,
                                "header": hdr})
            known = all(m["header"] for m in members)
            dups.append({
                "archive": a, "symbol": sym, "members": members,
                "headers_differ": len({m["header"] for m in members}) > 1,
                "bodies": "unknown" if not known else
                          "identical" if len(bodies) == 1 else "differ",
                "callers": callers(sym)})
    return dups


def trace(dups, machname):
    """Which member each aces2 binary links for each duplicated symbol."""
    wanted = defaultdict(set)
    for d in dups:
        wanted[d["archive"][3:-2]].add(d["symbol"])
    makeargs = ["FFLAGS_EXTRA=-fPIC", "F9XFLAGS_EXTRA=-fPIC",
                "CFLAGS_EXTRA=-fPIC", "CXXFLAGS_EXTRA=-fPIC"]
    if machname:
        makeargs.insert(0, "MACHNAME=" + machname)
    bins = []
    for line in open(os.path.join(ACES2, ".gitignore")):
        m = re.match(r"^([\w.-]+)/(x[\w.-]+)$", line.strip())
        if m and os.path.isdir(os.path.join(ACES2, m.group(1))):
            bins.append(m.groups())
    out = []
    with tempfile.TemporaryDirectory() as scratch:
        for mod, binary in bins:
            mdir = os.path.join(ACES2, mod)
            dry = run(["gmake", "-n", "-W", "lib%s.a" % mod] + makeargs + [binary],
                      cwd=mdir).stdout.splitlines()
            link = [l for l in dry if " -o %s " % binary in l]
            if not link:
                continue
            cmd = link[0].split(" 1>")[0]
            syms = sorted(set().union(*[wanted[l] for l in
                                        re.findall(r"(?<=\s)-l(\S+)", cmd) if l in wanted]))
            if not syms:
                continue
            cmd = cmd.replace(" -o %s " % binary,
                              " -o %s " % os.path.join(scratch, binary))
            cmd += "".join(" -Wl,--trace-symbol=%s" % s for s in syms)
            r = run(cmd, shell=True, cwd=mdir)
            got = {}
            for l in (r.stdout + r.stderr).splitlines():
                m = re.match(r"^(?:ld: )?(.*?): definition of (\S+)$", l)
                if m:
                    got[m.group(2)] = m.group(1).replace("ld: ", "")
            for s in syms:
                if s in got:
                    out.append({"binary": binary, "symbol": s, "member": got[s]})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--trace", action="store_true",
                    help="relink every binary and report the member each one gets")
    ap.add_argument("--machname", default=os.environ.get("ACES_PROGRAM_MACHNAME"),
                    help="MACHNAME for gmake (default: $ACES_PROGRAM_MACHNAME)")
    ap.add_argument("--json", action="store_true", help="JSON instead of text")
    args = ap.parse_args()

    dups = find_duplicates()
    links = trace(dups, args.machname) if args.trace else []
    if args.json:
        print(json.dumps({"duplicates": dups, "links": links}, indent=1))
        return
    for d in dups:
        flags = []
        if d["headers_differ"]:
            flags.append("HEADERS DIFFER")
        flags.append("bodies " + d["bodies"])
        print("%s  %s  (%s)" % (d["archive"], d["symbol"], ", ".join(flags)))
        for m in d["members"]:
            print("   %4s  %-32s %s" % (m["position"], m["member"], m["header"] or "-"))
        print("   callers: " + (" ".join(d["callers"]) or "none"))
        for l in links:
            if l["symbol"] == d["symbol"] and any(
                    l["member"].endswith(("(%s)" % m["member"], m["member"]))
                    for m in d["members"]):
                print("   linked:  %-16s <- %s" % (l["binary"], l["member"]))
    print("\n%d duplicated symbols" % len(dups))


if __name__ == "__main__":
    main()
