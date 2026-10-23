#!/usr/bin/env python3
"""Report where aces2py.so resolves a routine differently than the classic build.

aces2 is ~150 module libraries that were never meant to share a link unit: the
same routine name is defined in several of them, with genuinely different code
behind it, and the classic build gets away with it because each member
executable links only the handful of libraries named in that module's
ACES_LINK_LISTS entry -- so xvcc's DRAOLAD and xvee's DRAOLAD are simply
different routines. pyaces links everything into one shared object, where a
name can only mean one thing, and whichever archive the linker reaches first
silently wins for everybody.

That is the bug class this project has been paying for one failing test at a
time (FINISH, GSCHMIDT, MODF, EXPDEN, HBARXC, DT2INT2, Z2TOMO, W45T3DRV,
CMPENG, XGEMM ...). Some of those crash immediately; the expensive ones just
return a slightly wrong energy. This script finds the rest of them up front.

For every symbol defined by more than one module and referenced by a module
that does not define it, it works out
  * which definition the classic per-executable link would give that caller
    (walk its ACES_LINK_LISTS entry, first library that defines it wins), and
  * which definition actually won in aces2py.so (matched on symbol size),
and reports the disagreements whose two source files are not identical. A
disagreement between byte-identical copies is noise; one between different
implementations is a live wrong-answer risk on whatever path reaches it.

  python3 tools/collision_audit.py                 # report against the baseline
  python3 tools/collision_audit.py --update        # accept the current state
  python3 tools/collision_audit.py --calls         # fail if a pyaces call binds to data

WHAT THIS LIST IS NOT: a work queue. Every collision fixed so far started from a
test that failed and was verified one at a time against the tree's own reference
values. An attempt to go the other way -- take 42 entries the script judged
"safe" (single agreed-on owner, winner is not it) and rename the winners in one
pass -- dropped the pyaces suite from 45 passing to 31. The caller model here is
incomplete: ACES_LINK_LISTS only has entries for modules that build an
executable, so every reference from a pure library module (libr, libr2, librt3,
crust, tools, acescore ...) is invisible to it, and a "nobody else wants the
winner" verdict can simply be wrong. Modelling it properly means reproducing
what the linker does. Treat an entry as a lead worth checking when something
misbehaves on that path -- that is how CHECK_T2 was found before its test was
investigated -- not as a defect to fix on sight.

The baseline is the set of disagreements known and accepted at the time it was
written -- either already handled in build.sh's rename list, or on a path
nothing exercises yet. New entries appearing there mean a module (or a .skip
directory that has just been re-enabled) started shadowing something.
"""
import argparse, collections, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ACES2 = os.path.normpath(os.path.join(HERE, '..', '..', 'aces2'))
SO_GLOB = os.path.normpath(os.path.join(HERE, '..', 'Main', 'build'))
BASELINE = os.path.join(HERE, 'collision_baseline.txt')
COMMON_BASELINE = os.path.join(HERE, 'common_shadow_baseline.txt')


def find_so():
    for f in sorted(os.listdir(SO_GLOB)):
        if f.startswith('aces2py') and f.endswith('.so'):
            return os.path.join(SO_GLOB, f)
    sys.exit(f'no aces2py*.so under {SO_GLOB} -- build pyaces first')


def link_lists():
    mk = open(f'{ACES2}/Makefiles/GNUmakefile', errors='ignore').read()
    default = re.search(r'^LIST_DEF_LIBS\s*\?=\s*(\S+)', mk, re.M).group(1)
    block = re.search(r'^ACES_LINK_LISTS\s*=(.*?)(?=\n\n)', mk, re.M | re.S).group(1)
    out = {}
    for line in block.split('\n'):
        line = line.strip().rstrip('\\').strip().replace('${LIST_DEF_LIBS}', default)
        if not line or ':' not in line:
            continue
        parts = [p for p in line.split(':') if p]
        out[parts[0]] = parts
    return out


def nm(path, *flags):
    return subprocess.run(['nm', *flags, path], capture_output=True, text=True).stdout


def archives():
    return sorted(f for f in os.listdir(f'{ACES2}/lib') if f.endswith('.a'))


def symbol_tables():
    defs, refs, sizes = collections.defaultdict(set), collections.defaultdict(set), {}
    for f in archives():
        mod = f[3:-2]
        sizes[mod] = {}
        for line in nm(f'{ACES2}/lib/{f}', '-S', '--defined-only').split('\n'):
            p = line.split()
            if len(p) == 4 and p[2] == 'T':
                defs[p[3]].add(mod)
                sizes[mod][p[3]] = int(p[1], 16)
        for line in nm(f'{ACES2}/lib/{f}', '--undefined-only').split('\n'):
            p = line.split()
            if len(p) == 2 and p[0] == 'U':
                refs[p[1]].add(mod)
    return defs, refs, sizes


def linked_sizes(so):
    out = {}
    for line in nm(so, '-S').split('\n'):
        p = line.split()
        if len(p) == 4 and p[2] == 'T':
            out[p[3]] = int(p[1], 16)
    return out


def source_of(mod, sym):
    name = sym.rstrip('_')
    d = f'{ACES2}/{mod}'
    if not os.path.isdir(d):
        return None
    for ext in ('.F', '.f'):
        p = f'{d}/{name}{ext}'
        if os.path.isfile(p):
            return p
    hit = subprocess.run(
        ['grep', '-rlniE', rf'^ +(subroutine|.*function) +{name} *\(', d],
        capture_output=True, text=True).stdout.split('\n')
    for p in hit:
        if p.endswith(('.F', '.f')):
            return p
    return None


def code_of(path):
    if not path:
        return None
    return [l.rstrip().lower() for l in open(path, errors='ignore').read().split('\n')
            if l[:1] not in ('c', 'C', '*', '!') and l.strip()]



def renamed_pairs(defs):
    """(module, symbol) pairs that build.sh renames out of the way.

    Without this the audit reports its own fixes back as problems: build.sh
    renames lambda's DIISLST so vcc keeps its own, and an audit that only reads
    the archives then sees "the .so gives vcc's, but lambda wants lambda's" --
    true of the archives, false of the shared object lambda is actually linked
    into, where its calls go to pyaces_lambda_diislst_.
    """
    txt = open(os.path.join(HERE, '..', 'build.sh'), errors='ignore').read()
    pairs = set()
    m = re.search(r'EXTRA_LOCALIZE_SYMBOLS="([^"]*)"', txt, re.S)
    if m:
        for tok in m.group(1).split():
            if ':' in tok:
                mod, sym = tok.split(':', 1)
                pairs.add((mod, sym))
    m = re.search(r'NEWLY_LINKED_MODULES="([^"]*)"', txt)
    if m:
        for mod in m.group(1).split():
            for sym, owners in defs.items():
                if mod in owners and len(owners) > 1 and sym != 'MAIN__':
                    pairs.add((mod, sym))
    return pairs


def common_block_shadows(so):
    """Find subroutines shadowed by a same-named Fortran COMMON block.

    A DIFFERENT failure mode from the rest of this file. Everything above
    compares two subroutines that share a name and asks which implementation
    won. Here the competitor is not code at all: `common /abcd/` emits a
    global data symbol abcd_, the same name Fortran gives SUBROUTINE ABCD.

    A tentative (C) symbol loses to a strong (T) definition only if the
    linker has already pulled in the member holding it. When the commons
    come first they satisfy the reference outright, the archive member
    defining the subroutine is never extracted, and the name ends up in BSS.
    Calling it then jumps into zero-filled data -- an instruction-fetch
    fault, si_addr == $pc, not a wrong answer. There is no reading of this
    under which the program is correct, so anything reported here with a
    live caller is a real defect, not a judgement call.

    Found the expensive way: this is what blocked every KS case. vscf ->
    setnumint -> abcd faulted because seven modules (anti, eompccd, intprc,
    lcct, nmr, tdee, vee) declare `common /abcd/` and intgrt/abcd.F defines
    the 8-line SUBROUTINE ABCD that vscf actually wanted.

    Reported separately from the audit above because the fix differs. When
    the defining module is also the only caller, rename it (build.sh's
    EXTRA_LOCALIZE_SYMBOLS moves definition and caller together). When some
    OTHER module calls it, renaming the definition orphans that caller --
    the common blocks are what have to move instead.

    A name build.sh already renames is out of the collision whether or not
    the rename is visible in the .so. Most of these definers (mopac, nddo,
    vdint, pdens, ccd_light, aces2) contribute nothing to this link, so
    pyaces_<mod>_<sym> never appears there to be found -- reading only the
    .so would keep reporting a fix that is in place, and would have made
    the whole 20th-instance batch unverifiable. Read build.sh's list too.
    """
    text, common, refs = (collections.defaultdict(set) for _ in range(3))
    for f in archives():
        mod = f[3:-2]
        for line in nm(f'{ACES2}/lib/{f}').split('\n'):
            p = line.split()
            if len(p) == 2 and p[0] == 'U':
                refs[p[1]].add(mod)
            elif len(p) == 3:
                if p[1] == 'T':
                    text[p[2]].add(mod)
                elif p[1] == 'C':
                    common[p[2]].add(mod)

    won = {}
    for line in nm(so, '--defined-only').split('\n'):
        p = line.split()
        if len(p) == 3:
            won.setdefault(p[2], p[1])
    renamed = renamed_pairs(text)

    def moved(mod, sym):
        # A definer build.sh has already renamed is out of the collision: its
        # definition and its own callers both moved to pyaces_<mod>_<sym>, so
        # the plain name being data no longer costs anybody a call.
        return (mod, sym) in renamed or f'pyaces_{mod}_{sym}' in won

    out = []
    for sym in sorted(set(text) & set(common)):
        if won.get(sym, '') not in ('B', 'b', 'D', 'd'):
            continue
        definers = sorted(m for m in text[sym] if not moved(m, sym))
        if not definers:
            continue
        callers = sorted(m for m in refs[sym] if not moved(m, sym))
        if not callers:
            continue
        verdict = ('rename-definer' if set(callers) <= set(definers)
                   else 'rename-commons')
        out.append('%s %s %s %s' % (sym, ','.join(definers),
                                    ','.join(callers), verdict))
    return out


def audit():
    so = find_so()
    links = link_lists()
    defs, refs, sizes = symbol_tables()
    linked = linked_sizes(so)
    renamed = renamed_pairs(defs)
    found = set()
    for sym, owners in defs.items():
        if len(owners) < 2 or sym == 'MAIN__' or sym not in linked:
            continue
        live = {m for m in owners if (m, sym) not in renamed}
        if not live:
            continue
        winners = {m for m in live if sizes[m].get(sym) == linked[sym]}
        if len(winners) != 1:
            continue                     # cannot tell them apart by size
        win = next(iter(winners))
        for caller in refs.get(sym, ()):
            if caller not in links:
                continue
            if (caller, sym) in renamed:
                continue                 # keeps its own through the rename
            # A module that defines the routine AND calls it means its own copy;
            # the classic build guarantees that, because the linker satisfies the
            # reference from the archive it is already pulling the caller out of.
            # In one .so it can lose to another module's copy just like anyone
            # else -- that is how vcc ended up running lambda's INITIN.
            serial = caller if caller in owners else \
                next((m for m in links[caller] if m in owners), None)
            if serial is not None and (serial, sym) in renamed:
                # the classic answer is a copy this build renames away; the
                # caller now gets whatever else won, which is the point of it
                serial = next((m for m in links[caller] if m in live), None)
            if serial is None or serial == win:
                continue
            if code_of(source_of(win, sym)) == code_of(source_of(serial, sym)):
                continue                 # same code, different archive: harmless
            found.add(f'{sym} {win} {serial} {caller}')
    return sorted(found)


def calls_into_data(so):
    """Find calls from pyaces's own Fortran objects that bind to data.

    common_block_shadows() above needs a subroutine (T) somewhere to report a
    name lost to a COMMON block. A call to a name that is ONLY ever a common
    has no such subroutine, so it never shows up there -- and it is just as
    fatal. That is how Runhbar and Buildhbar shipped calling `Hbar`: aces2
    defines SUBROUTINE HBAR_ (hbar__), so hbar_ resolved to lambda's
    COMMON /HBAR/ and the call would have jumped into BSS.

    Fortran objects reference commons as C, never U, so every U in them is a
    call; one that the .so defines as data is a defect with no exceptions.
    The two C objects (the f2py module and fortranobject) are left out:
    they legitimately reference C data such as PyFortran_Type.
    """
    data = set()
    for line in nm(so, '--defined-only').split('\n'):
        p = line.split()
        if len(p) == 3 and p[1] in ('B', 'b', 'D', 'd'):
            data.add(p[2])
    out = []
    for f in sorted(os.listdir(SO_GLOB)):
        if not f.endswith('.o') or f in ('aces2pymodule.o', 'fortranobject.o'):
            continue
        for line in nm(os.path.join(SO_GLOB, f), '--undefined-only').split('\n'):
            p = line.split()
            if len(p) == 2 and p[0] == 'U' and p[1] in data:
                out.append(f'{p[1]} {f}')
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--update', action='store_true',
                    help='rewrite the baseline with the current findings')
    ap.add_argument('--calls', action='store_true',
                    help='only check that pyaces objects call code, not data; '
                         'exit 1 if any call binds to a COMMON block')
    args = ap.parse_args()

    if args.calls:
        bad = calls_into_data(find_so())
        print(f'@COLLISION_AUDIT: {len(bad)} calls from pyaces objects bound '
              f'to data instead of code')
        for l in bad:
            sym, obj = l.split()
            print(f'    {sym:22s} called from {obj}, defined as data in the .so')
        return 1 if bad else 0

    current = audit()
    shadows = common_block_shadows(find_so())
    if args.update:
        for path, header, rows in (
                (COMMON_BASELINE, '# symbol  defined-in  callers  verdict', shadows),
                (BASELINE, '# symbol  winner-in-so  what-the-classic-link-gives'
                           '  caller', current)):
            with open(path, 'w') as fh:
                fh.write('\n'.join([header, *rows]) + '\n')
        print(f'baseline updated: {len(current)} entries, '
              f'{len(shadows)} common-block shadows')
        return 0

    known_shadows = set()
    if os.path.isfile(COMMON_BASELINE):
        known_shadows = {l.strip() for l in open(COMMON_BASELINE)
                         if l.strip() and not l.startswith('#')}
    new_shadows = [l for l in shadows if l not in known_shadows]
    print(f'@COLLISION_AUDIT: {len(shadows)} subroutines shadowed by a '
          f'same-named COMMON block that build.sh does not rename out of the way')
    if new_shadows:
        print()
        print('@COLLISION_AUDIT: NEW common-block shadows -- a CALL now lands in '
              'BSS, which')
        print('                  faults on instruction fetch the moment it is '
              'reached.')
        for l in new_shadows:
            sym, definers, callers, verdict = l.split()
            print(f'    {sym:22s} subroutine in {definers}, called by {callers}'
                  f' -- {verdict}')
        print()
        print('                  rename-definer: add "<definer>:<symbol>" to '
              'EXTRA_LOCALIZE_SYMBOLS.')
        print('                  rename-commons: the callers are outside the '
              'defining module, so')
        print('                  renaming the definition would orphan them -- '
              'move the commons.')
        print()

    known = set()
    if os.path.isfile(BASELINE):
        known = {l.strip() for l in open(BASELINE)
                 if l.strip() and not l.startswith('#')}
    new = [l for l in current if l not in known]
    gone = [l for l in sorted(known) if l not in current]

    print(f'@COLLISION_AUDIT: {len(current)} resolution differences '
          f'({len(known)} in the baseline)')
    for l in gone:
        print(f'  resolved since the baseline: {l}')
    if new:
        print()
        print('@COLLISION_AUDIT: NEW resolution differences -- a module started '
              'shadowing a routine')
        print('                  another module owns. Each of these silently gives '
              'some caller the')
        print('                  wrong implementation; see build.sh for how the '
              'known ones are fixed.')
        for l in new:
            sym, win, serial, caller = l.split()
            print(f'    {sym:22s} .so gives {win}, {caller} expects {serial}')
        print()
        print('                  Fix by adding "<loser>:<symbol>" to '
              'EXTRA_LOCALIZE_SYMBOLS in build.sh,')
        print('                  or accept them with '
              'tools/collision_audit.py --update')
    return 0


if __name__ == '__main__':
    sys.exit(main())
