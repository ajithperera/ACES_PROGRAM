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


def audit():
    so = find_so()
    links = link_lists()
    defs, refs, sizes = symbol_tables()
    linked = linked_sizes(so)
    found = set()
    for sym, owners in defs.items():
        if len(owners) < 2 or sym == 'MAIN__' or sym not in linked:
            continue
        winners = {m for m in owners if sizes[m].get(sym) == linked[sym]}
        if len(winners) != 1:
            continue                     # cannot tell them apart by size
        win = next(iter(winners))
        for caller in refs.get(sym, ()):
            if caller not in links:
                continue
            # A module that defines the routine AND calls it means its own copy;
            # the classic build guarantees that, because the linker satisfies the
            # reference from the archive it is already pulling the caller out of.
            # In one .so it can lose to another module's copy just like anyone
            # else -- that is how vcc ended up running lambda's INITIN.
            serial = caller if caller in owners else \
                next((m for m in links[caller] if m in owners), None)
            if serial is None or serial == win:
                continue
            if code_of(source_of(win, sym)) == code_of(source_of(serial, sym)):
                continue                 # same code, different archive: harmless
            found.add(f'{sym} {win} {serial} {caller}')
    return sorted(found)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--update', action='store_true',
                    help='rewrite the baseline with the current findings')
    args = ap.parse_args()

    current = audit()
    if args.update:
        with open(BASELINE, 'w') as fh:
            fh.write('# symbol  winner-in-so  what-the-classic-link-gives  caller\n')
            fh.write('\n'.join(current) + '\n')
        print(f'baseline updated: {len(current)} entries')
        return 0

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
