#!/usr/bin/env python3
"""Classify the ACES II fast test suite (~/Develop/ACESII_build/tests/fast) for
pyaces compatibility. pyaces currently exposes: Runints, Runscf, Runcc, Runhbar,
Runpccd, Runprops, Runee -- no gradient/Vdens path (rungrads/run1props excluded
from the pyaces build), so any case requiring an analytic/numerical gradient or
geometry optimization is out of scope for now. SCRIPT (chained multi-job) cases
are also out of scope for this first pass -- different driver shape entirely.
"""
import os
import re
import sys

SUITE = os.path.expanduser("~/Develop/ACESII_build/tests/fast")
CASES = os.path.join(SUITE, "cases")
TEST_RESULTS = os.path.join(SUITE, "test_results")

# CALC= values that route through Runcc (xvcc handles all of these in the
# traditional pipeline: MP2 is a VCC-internal special case too).
CC_METHODS = {
    "SCF", "HF", "MBPT(2)", "MP2", "MBPT(3)", "MP3", "MBPT(4)", "MP4",
    "SDQ-MBPT(4)", "CCD", "QCISD", "CCSD", "CCSD(T)", "CCSD+T(CCSD)",
    "CCSD(T)_L", "CC2", "CC3", "CCSDT-1", "CCSDT-1B", "CCSDT-2", "CCSDT-3",
    "CCSDT", "CI(SD)", "CISD", "QCISD(T)", "UCC(4)", "LCCD", "LCCSD",
    "CCD+ST(CCD)", "BD", "BD(T)", "CCSDTQ", "CCSDT-4",
}
PCCD_METHODS = {"PCCD", "PCCD(T)"}
EOM_MARKERS = {"EOM", "EOM-CCSD", "EOM-CCSD(T)", "EOMCC", "EOMEE", "EOMIP", "EOMEA"}


def parse_test_results(path):
    entries = {}
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f]
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        name, src = lines[i].split(None, 1)
        record_line = lines[i + 1].split()
        record = record_line[0]
        value = lines[i + 2].strip()
        entries[name] = {"src": src.strip(), "record": record, "value": value}
        i += 3
    return entries


def parse_zmat_keywords(zmat_path):
    """Return dict of ACES2 keyword=value pairs (uppercased keys), collapsed
    across the (possibly multi-line, comma-separated) *ACES2(...) block or the
    older bare '*ACES2' + following KEY=VAL lines form."""
    with open(zmat_path, errors="replace") as f:
        text = f.read()
    kw = {}
    start = text.find("*ACES2(")
    if start != -1:
        # Balanced-paren scan -- values like CALC=CCSD(T) contain their own
        # parens, so a non-greedy regex to the first ")" truncates the block.
        depth = 0
        i = start + len("*ACES2")
        body_start = i + 1
        body_end = None
        for j in range(i, len(text)):
            if text[j] == "(":
                depth += 1
            elif text[j] == ")":
                depth -= 1
                if depth == 0:
                    body_end = j
                    break
        body = text[body_start:body_end] if body_end else ""
    else:
        m2 = re.search(r"\*ACES2\s*\n(.*?)\n\s*\n", text, re.S)
        body = m2.group(1) if m2 else ""
    for tok in re.split(r"[,\n]", body):
        tok = tok.strip()
        if not tok or "=" not in tok:
            continue
        k, v = tok.split("=", 1)
        kw[k.strip().upper()] = v.strip().upper()

    # ACES II's classic convention: a '*' suffix on a Z-matrix variable name
    # (e.g. "H 1 R* 2 A*") marks it for optimization -- independent of, and
    # much more common in this archive than, an explicit GEOM_OPT/VIB
    # keyword. Only look at the geometry section (before *ACES2), and strip
    # comments first so a stray '#... *foo' note can't false-positive.
    geom_section = text.split("*ACES2", 1)[0]
    geom_section = "\n".join(l.split("#", 1)[0] for l in geom_section.splitlines())
    starred = bool(re.search(r"\b\w+\*", geom_section))

    return kw, starred


def classify(name, kw, starred):
    calc = kw.get("CALC", "SCF")
    ref = kw.get("REFERENCE", kw.get("REF", "RHF"))
    geom = kw.get("GEOM_OPT", kw.get("GEOM", ""))
    vib = kw.get("VIB", "")
    excite = kw.get("EXCITE", "")
    estate = kw.get("ESTATE_PROP", "")
    deriv_lvl = kw.get("DERIV_LEV", "")

    reasons = []
    if starred:
        reasons.append("starred Z-matrix variable(s) (implicit GEOM_OPT, needs gradients)")
    if vib not in ("", "OFF", "NO"):
        reasons.append(f"VIB={vib} (Hessian/gradient path)")
    if geom not in ("", "NOOPTIMIZE", "OFF"):
        reasons.append(f"GEOM={geom} (optimization needs gradients)")
    if deriv_lvl not in ("", "0"):
        reasons.append(f"DERIV_LEV={deriv_lvl}")
    if any(m in calc for m in ("MRCC",)):
        reasons.append("MRCC (already dropped from suite as unsupportable)")

    # QRHF (without MAKERHF) makes joda/gtflgs.F internally force the
    # reference to UHF/ROHF regardless of any explicit REFERENCE=/REF=
    # keyword (see gtflgs.F's own "Ajith 03/2000" comment) -- Runscf/Runcc
    # must be called with iuhf=1 for these cases even though the ZMAT
    # itself never says REFERENCE=UHF. MAKERHF (rare, "dangerous keyword",
    # unused anywhere in this suite) is the only thing that keeps it RHF.
    is_qrhf = any(k.startswith("QRHF_G") for k in kw) and "MAKERHF" not in kw
    iuhf = 0 if (ref == "RHF" and not is_qrhf) else 1

    if reasons:
        return "EXCLUDED", "; ".join(reasons), None, iuhf

    if excite not in ("", "NONE") or estate not in ("", "OFF") or any(
        m in calc for m in EOM_MARKERS
    ):
        return "EOM", f"CALC={calc} EXCITE={excite} ESTATE_PROP={estate}", "Runee", iuhf

    if calc in PCCD_METHODS:
        return "PCCD", f"CALC={calc}", "Runpccd", iuhf

    if calc in ("SCF", "HF"):
        return "SCF", f"CALC={calc}", "Runscf", iuhf

    # Default: any other correlated single-point method (CCSD variants, CI
    # variants, MBPT/MPn, ACCSD(T), CCSDTQ, BCCD, QCISD, ...) routes through
    # Vcc in the traditional pipeline too -- treat as CC unless known
    # otherwise, rather than maintaining an exhaustive whitelist.
    return "CC", f"CALC={calc}", "Runcc", iuhf


def main():
    entries = parse_test_results(TEST_RESULTS)
    rows = []
    for name, info in sorted(entries.items()):
        case_dir = os.path.join(CASES, name)
        zmat_path = os.path.join(case_dir, "ZMAT")
        if not os.path.isfile(zmat_path):
            rows.append((name, "SCRIPT", "chained multi-job case, out of scope for pass 1",
                         None, None, info["record"], info["value"]))
            continue
        kw, starred = parse_zmat_keywords(zmat_path)
        cat, reason, entrypoint, iuhf = classify(name, kw, starred)
        rows.append((name, cat, reason, entrypoint, iuhf, info["record"], info["value"]))

    counts = {}
    for r in rows:
        counts[r[1]] = counts.get(r[1], 0) + 1

    print("=== Category counts ===")
    for k, v in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {k}: {v}")
    print()

    out_path = os.path.join(os.path.dirname(SUITE), "..", "..", "acespy_build",
                             "suite_classification.tsv")
    out_path = os.path.abspath(os.path.join("/blue/ufhpc/perera/acespy_build",
                                             "suite_classification.tsv"))
    with open(out_path, "w") as f:
        f.write("name\tcategory\treason\tentrypoint\tiuhf\trecord\tref_value\n")
        for r in rows:
            f.write("\t".join("" if x is None else str(x) for x in r) + "\n")
    print(f"Wrote {len(rows)} rows to {out_path}")


if __name__ == "__main__":
    main()
