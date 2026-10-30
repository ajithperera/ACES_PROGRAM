#!/bin/bash
set -e

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Toolchain, first match wins:
#  1. PYACES_MODULES, when set (the CMake superbuild; ../build.sh passes the
#     set ../configure recorded). "" uses the current environment (ifx,
#     MKLROOT, python3 with numpy already set up).
#  2. ACES_PROGRAM_MODULES from the environment, when non-empty.
#  3. Run on its own: ../build_config.sh's loader, i.e. ../configure's choice
#     of modules, a setvars script or PATH, and its Fortran and C compilers.
#  4. No build_config.sh: the HiPerGator module set.
source "$SCRIPT_DIR/../cmake/env.sh"
if [ -n "${PYACES_MODULES+set}" ]; then
    aces_load_modules "$PYACES_MODULES"
elif [ -n "${ACES_PROGRAM_MODULES:-}" ]; then
    aces_load_modules "$ACES_PROGRAM_MODULES"
elif [ -f "$SCRIPT_DIR/../build_config.sh" ]; then
    # shellcheck source=/dev/null
    source "$SCRIPT_DIR/../build_config.sh"
    if type aces_program_load_env >/dev/null 2>&1; then
        aces_program_load_env || exit 1
    else
        # Written by an older ./configure, which knew only modules.
        aces_load_modules "$ACES_PROGRAM_MODULES"
    fi
    PYACES_FC=${PYACES_FC:-${ACES_PROGRAM_FC:-}}
    PYACES_CC=${PYACES_CC:-${ACES_PROGRAM_CC:-}}
else
    aces_load_modules "intel/2025.1.0 openmpi/5.0.7 python/3.12"
fi
: "${MKLROOT:?is unset -- the loaded modules provide no Intel MKL}"
# Fortran compiler for every pyaces object and the link: PYACES_FC (../build.sh
# passes ./configure's --fortran choice), else ifx, else ifort.
if [ -z "${PYACES_FC:-}" ]; then
    if command -v ifx >/dev/null 2>&1; then PYACES_FC=ifx; else PYACES_FC=ifort; fi
fi
command -v "$PYACES_FC" >/dev/null 2>&1 || { echo "$PYACES_FC is not on PATH" >&2; exit 1; }
# C compiler for the f2py module and fortranobject.c (plain C, no MPI):
# PYACES_CC (../build.sh passes ./configure's choice), else icx, else gcc.
if [ -z "${PYACES_CC:-}" ]; then
    if command -v icx >/dev/null 2>&1; then PYACES_CC=icx; else PYACES_CC=gcc; fi
fi
command -v "$PYACES_CC" >/dev/null 2>&1 || { echo "$PYACES_CC is not on PATH" >&2; exit 1; }
echo "compilers: $(command -v $PYACES_FC), $(command -v $PYACES_CC)"
# The globs and sorts below set archive member and -l order, and with them which
# copy of a multiply-defined routine the link picks. Collate bytewise so that
# does not change with the caller's locale (PBS jobs often run with none).
export LC_ALL=C

WORK=$SCRIPT_DIR
MAIN=$WORK/Main
PYF=$WORK/a2_pyf
ACES=$SCRIPT_DIR/../aces2
NUMPY_INC=$(python3 -c "import numpy; print(numpy.get_include())")
PY_INC=$(python3 -c "import sysconfig; print(sysconfig.get_path('include'))")
# The extension suffix (cpython-312-x86_64-linux-gnu.so, ...) belongs to the
# interpreter doing the building; a hard-coded one only imports under 3.12.
EXT_SUFFIX=$(python3 -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")
# fortranobject.c/.h ship inside numpy itself, so ask the interpreter that is
# actually being built against rather than hardcoding one site's install path.
F2PY_SRC=$(python3 -c "import numpy, os; print(os.path.join(os.path.dirname(numpy.__file__), 'f2py', 'src'))")
MKL=$MKLROOT/lib/intel64
# oneAPI's own layout has the libraries directly under lib/ on some installs.
[ -d "$MKL" ] || MKL=$MKLROOT/lib

cd $MAIN
rm -rf gen build
mkdir -p build

# The Python interface is exactly the subroutines in these files: f2py scans
# them and wraps every routine it finds, taking names, argument types and
# intents from the source itself (Cf2py lines, e.g. "Cf2py intent(out) Length").
# There is no .pyf. The old Main/aces2.pyf was passed alongside these files,
# but f2py then built the module twice and the source scan overwrote the .pyf
# one, so its declarations never took effect. Passing EVERY a2_pyf source
# would wrap internal helpers too (and their all-caps F77 names hit a
# case-mangling mismatch), so only list files whose routines are meant to be
# called from Python.
PYF_SOURCES="$MAIN/init.F $PYF/dgetrecpy.F $PYF/igetrecpy.F $PYF/dputrecpy.F \
    $PYF/iputrecpy.F $PYF/cgetrecpy.F $PYF/cputrecpy.F $PYF/getlstpy.F \
    $PYF/putlstpy.F $PYF/getlstdimpy.F $PYF/getallpy.F $PYF/putallpy.F \
    $PYF/getreclenpy.F \
    $PYF/get1ehpy.F $PYF/get2ehpy.F $PYF/buildhbar.F $PYF/runints.F \
    $PYF/runscf.F $PYF/runcc.F $PYF/runhbar.F $PYF/runprops.F $PYF/runpccd.F \
    $PYF/runee.F $PYF/rundirmp2.F $PYF/runintgrt.F $PYF/acesjafin.F $PYF/acesjainit.F"

# Full source list for compilation+linking (everything the above call, plus
# their internal Fortran-level dependencies, actually need at link time).
SOURCES=""
for f in $PYF/*.F; do
    base=$(basename "$f" .F)
    case "$base" in
        rungrads|run1props|run) continue ;;
    esac
    SOURCES="$SOURCES $f"
done
SOURCES="$SOURCES $MAIN/init.F $MAIN/aces2_main.F $WORK/aux_pyaces/*.F"

echo "=== step 1: f2py generate wrappers ==="
python3 -m numpy.f2py $PYF_SOURCES -m aces2py --build-dir gen > $WORK/f2py_gen.log 2>&1

# f2py's F_FUNC(f,F) macro uses "f" (whatever case the Fortran source used
# for the routine name) as its lowercase form, assuming the compiler
# preserves source case. ifort/gfortran always fold external symbols to
# lowercase regardless of source case, so any routine not written in
# all-lowercase source (e.g. "Subroutine Init()") ends up with a mismatched
# undefined symbol at import time (e.g. "Init_" instead of "init_"). Force
# the first F_FUNC/F_FUNC_US argument to lowercase to match reality.
sed -i -E 's/(F_FUNC(_US)?\()([A-Za-z0-9_]+)(,)/\1\L\3\E\4/g' gen/aces2pymodule.c

echo "=== step 2: compile C sources ==="
$PYACES_CC -DNDEBUG -fPIC -O2 -I. -I$NUMPY_INC -I$PY_INC -I$F2PY_SRC \
    -c gen/aces2pymodule.c -o build/aces2pymodule.o
$PYACES_CC -DNDEBUG -fPIC -O2 -I. -I$NUMPY_INC -I$PY_INC -I$F2PY_SRC \
    -c $F2PY_SRC/fortranobject.c -o build/fortranobject.o

echo "=== step 3: compile Fortran sources ==="
FFLAGS="-diag-disable=10448 -i8 -fPIC -fp-model strict -O1 -assume minus0 -D__fortran -D__fortran77 -D_INTEL -D_RECL_IS_WORDS_ -I$ACES/include -I$ACES/acescore/include -I. -FI"

for f in $SOURCES; do
    base=$(basename "$f" .F)
    echo "compiling $base"
    $PYACES_FC $FFLAGS -c "$f" -o build/${base}.o
done
$PYACES_FC $FFLAGS -c gen/aces2py-f2pywrappers.f -o build/aces2py-f2pywrappers.o

# runee.o's own new CALL PROPS(...) (see vprops:props_ in tools/localize.txt)
# must be redefined to match vprops's renamed copy the SAME way, since it's
# a pyaces-side .o, not one of aces2's own archives patched in the loop
# below.
objcopy --redefine-sym props_=pyaces_vprops_props_ build/runee.o

echo "=== step 4: link ==="
# Generated fresh every build (not cached) so this can never point at a
# stale/dead location -- ACES is always this same build's own aces2 tree.
LIBDIRS=$(find $ACES -maxdepth 2 -iname "*.a" | xargs -n1 dirname | sort -u | sed 's/^/-L/' | tr '\n' ' ')
LIBNAMES=$(find $ACES -maxdepth 2 -iname "*.a" | xargs -n1 basename | sed -E 's/^lib(.*)\.a$/\1/' | sort -u | sed 's/^/-l/' | tr '\n' ' ')
# 2026-10-21: standalone aces2 programs pyaces never calls, left out of the
# link. They were re-enabled in aces2 after 08/2026 (723989ab) and carry
# their own copies of vcc routines; sorted ahead of -lvcc they won the link
# (rpa: coreng, e4seng, dmptgss, linear, ...), and UCC(4) case 025 came out
# 2.5e-5 off while standalone xvcc stayed right.
for m in a2mix aux rpa scrnc; do
  LIBNAMES=$(echo " $LIBNAMES " | sed "s/ -l$m / /")
done

# 12 aces2 modules (intprc, vmol, joda, vcc, vscf, vtran, hbar, lambda, vee,
# pccd, props, vmol2ja) each have their own Fortran PROGRAM entry (MAIN__) --
# needed for their standalone xJODA/xVSCF/etc executables (which the classic
# ACES II driver/test suite depends on), but fatal ("multiple definition of
# MAIN__") once more than one is pulled into this single shared library.
# Rather than patching aces2's own sources (that tree must stay untouched --
# it's also the classic-driver suite's shared source of truth), localize the
# MAIN__ symbol on disposable copies of just these 12 archives so it no
# longer collides at link time. aces2's own installed .a files are never
# modified.
MAIN_COLLISION_MODULES="hbar joda lambda pccd props vcc vee vmol vmol2ja vscf vtran intprc"

# Every per-module symbol rename (and the reason for it) lives in
# tools/localize.txt, which tools/collision_audit.py reads as well.
LOCALIZE=$WORK/tools/localize.txt
EXTRA_LOCALIZE_SYMBOLS=$(sed 's/#.*//' $LOCALIZE | awk '$1=="rename"{print $2":"$3}' | tr '\n' ' ')
EXTRA_LOCALIZE_MODULES=$(sed 's/#.*//' $LOCALIZE | awk '$1=="extract"{print $2}' | tr '\n' ' ')

# a2mix, aux, rpa and scrnc were .skip module directories until 723989ab /
# aff1dc3c re-enabled them. Their archives only began competing in THIS link
# afterwards -- the classic per-executable builds never link them next to the
# modules below -- and each arrived carrying same-named copies of routines
# other modules own: a2mix alone shadows setptr, zero, seeklb, makmap, header,
# factor, piksr2, prntlo, setlmn, loadint, aainer, civpt and fmc; aux shadows
# runit; rpa shadows cmpeng, t2t1aa1, t2t1ab1, sinden, srteval and the
# dmptgss pair. Two of those were already caught the expensive way this
# session: rpa:cmpeng_ (case 025 lost the singles term and came out 2.49e-05
# wrong, silently) and a2mix:setptr_ -- a2mix's takes (last,max,type,num) while
# the 658 call sites in intgrt/vksdint/sbtools/intpack/vscf/oeptools pass
# sbtools' (stack,type,num), so the KS cases died with "setptr received
# invalid data type" and a2mix does not even call its own copy.
#
# Rather than keep discovering these one failing case at a time, rename every
# symbol these four share with any other module. --redefine-sym rewrites the
# definition AND that module's own callers together, so each keeps using its
# own copy while everyone else goes back to the module they resolved to before
# these four joined the link. MAIN__ is left to the mechanism below.
NEWLY_LINKED_MODULES="a2mix aux rpa scrnc"

ALLSYMS=$(mktemp)
for a in $ACES/lib/lib*.a; do
    b=$(basename "$a" .a)
    nm --defined-only "$a" 2>/dev/null | awk -v m="${b#lib}" '$2=="T"{print $3, m}'
done | sort -u > $ALLSYMS
for mod in $NEWLY_LINKED_MODULES; do
    shared=$(awk -v m="$mod" '$2==m{mine[$1]=1} $2!=m{other[$1]=1}
                              END{for (s in mine) if (s in other) print s}' $ALLSYMS |
             grep -v '^MAIN__$' | sort)
    [ -z "$shared" ] && continue
    echo "shadowed by $mod, renaming: $(echo $shared | tr '\n' ' ')"
    for sym in $shared; do
        EXTRA_LOCALIZE_SYMBOLS="$EXTRA_LOCALIZE_SYMBOLS $mod:$sym"
    done
    case " $EXTRA_LOCALIZE_MODULES " in
        *" $mod "*) ;;
        *) EXTRA_LOCALIZE_MODULES="$EXTRA_LOCALIZE_MODULES $mod" ;;
    esac
done
rm -f $ALLSYMS

# A pair the loop below cannot apply used to vanish without a trace: a module
# outside both lists is never extracted, and a misspelt symbol is simply not
# there for objcopy to rename. Either way the collision it was meant to fix
# stays in the .so. Refuse to link instead.
bad_pairs=
for pair in $EXTRA_LOCALIZE_SYMBOLS; do
    pmod=${pair%%:*}
    psym=${pair#*:}
    case " $MAIN_COLLISION_MODULES $EXTRA_LOCALIZE_MODULES " in
        *" $pmod "*) ;;
        *) bad_pairs="$bad_pairs
    $pair: $pmod is neither in MAIN_COLLISION_MODULES nor an extract line of tools/localize.txt"
           continue ;;
    esac
    nm --defined-only $ACES/lib/lib${pmod}.a 2>/dev/null | grep -q " T ${psym}\$" ||
        bad_pairs="$bad_pairs
    $pair: lib${pmod}.a defines no T ${psym}"
done
if [ -n "$bad_pairs" ]; then
    echo "tools/localize.txt has renames that cannot be applied:$bad_pairs" >&2
    exit 1
fi

PATCHED_LIBDIR=$WORK/build/patched_libs
# Wipe it first. This directory is linked with -L ahead of aces2's own lib/,
# so anything left here from a previous build shadows the real archive -- and
# it is only ever repopulated for the modules in the lists above, so dropping a
# module from them leaves its stale patched copy in the link line forever.
# That is not hypothetical: an experiment that briefly added intgrt to the list
# left a libintgrt.a here from before intgrt_run.F existed, and the next build
# linked that instead, failing at import with "undefined symbol: intgrt_run_"
# while nm showed the symbol sitting in aces2/lib/libintgrt.a all along.
rm -rf $PATCHED_LIBDIR
mkdir -p $PATCHED_LIBDIR
for mod in $MAIN_COLLISION_MODULES $EXTRA_LOCALIZE_MODULES; do
    src_a=$ACES/lib/lib${mod}.a
    [ -f "$src_a" ] || { echo "expected $src_a not found" >&2; exit 1; }
    workdir=$(mktemp -d)
    (cd $workdir && ar x $src_a)
    for o in $workdir/*.o; do
        nm "$o" 2>/dev/null | grep -q " T MAIN__$" && objcopy --localize-symbol=MAIN__ "$o"
        for pair in $EXTRA_LOCALIZE_SYMBOLS; do
            pmod=${pair%%:*}
            psym=${pair#*:}
            [ "$pmod" = "$mod" ] || continue
            newsym="pyaces_${pmod}_${psym}"
            objcopy --redefine-sym ${psym}=${newsym} "$o"
        done
    done
    (cd $workdir && ar rcs $PATCHED_LIBDIR/lib${mod}.a *.o)
    rm -rf $workdir
done

$PYACES_FC -shared -nofor-main -o build/aces2py${EXT_SUFFIX} \
    build/*.o \
    -L$PATCHED_LIBDIR \
    $LIBDIRS \
    -Wl,--start-group $LIBNAMES -Wl,--end-group \
    -L$MKL -Wl,--start-group -lmkl_intel_ilp64 -lmkl_lapack95_ilp64 -lmkl_sequential -lmkl_core -Wl,--end-group \
    -lpthread

echo "=== step 5: test import ==="
cd build
python3 -c "import aces2py; print('IMPORT OK')"

# A call that binds to a COMMON block instead of a subroutine faults the moment
# it runs, and nothing else here notices it (see calls_into_data in the audit).
# Unlike the audit below this one has no judgement call in it, so it fails the
# build.
python3 $WORK/tools/collision_audit.py --calls

# Which copy of a multiply-defined routine won this link is not something the
# linker will ever complain about -- it just picks one and everybody lives with
# it. Every collision fixed above was found by a test failing, several of them
# by an energy quietly coming out wrong. Report the ones that are still there,
# and shout if a new one appeared (a re-enabled .skip module is the usual
# cause). Advisory only: it must not fail a build that otherwise linked.
echo "=== step 6: link-resolution audit ==="
python3 $WORK/tools/collision_audit.py || true
