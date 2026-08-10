#!/bin/bash
set -e

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

SHIM=$SCRIPT_DIR/tools/compat_bin
export PATH=$SHIM:$PATH
module purge
module load intel/2025.1.0 openmpi/5.0.7 python/3.12

WORK=$SCRIPT_DIR
MAIN=$WORK/Main
PYF=$WORK/a2_pyf
ACES=$SCRIPT_DIR/../aces2
NUMPY_INC=$(python3 -c "import numpy; print(numpy.get_include())")
PY_INC=$(python3 -c "import sysconfig; print(sysconfig.get_path('include'))")
F2PY_SRC=/apps/python/3.12/lib/python3.12/site-packages/numpy/f2py/src
MKL=$MKLROOT/lib/intel64

cd $MAIN
rm -rf gen build
mkdir -p build

# f2py needs the ACTUAL defining source alongside the .pyf to generate real
# wrapper bodies (the .pyf alone produces empty stub functions). But passing
# EVERY source file also makes f2py auto-wrap every subroutine it finds in
# them as extra top-level Python functions we never declared -- and several
# of those (all-caps legacy F77 names) hit a case-mangling mismatch against
# our lowercase-exporting ifort build. So: only pass the files that actually
# define one of the .pyf's own declared interface routines.
PYF_SOURCES="$MAIN/init.F $PYF/dgetrecpy.F $PYF/igetrecpy.F $PYF/dputrecpy.F \
    $PYF/iputrecpy.F $PYF/cgetrecpy.F $PYF/cputrecpy.F $PYF/getlstpy.F \
    $PYF/putlstpy.F $PYF/getallpy.F $PYF/putallpy.F $PYF/getreclenpy.F \
    $PYF/get1ehpy.F $PYF/get2ehpy.F $PYF/buildhbar.F $PYF/runints.F \
    $PYF/runscf.F $PYF/runcc.F $PYF/runhbar.F $PYF/runprops.F $PYF/runpccd.F \
    $PYF/runee.F $PYF/acesjafin.F $PYF/acesjainit.F"

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
python3 -m numpy.f2py aces2.pyf $PYF_SOURCES -m aces2py --build-dir gen > $WORK/f2py_gen.log 2>&1

# f2py's F_FUNC(f,F) macro uses "f" (whatever case the Fortran source used
# for the routine name) as its lowercase form, assuming the compiler
# preserves source case. ifort/gfortran always fold external symbols to
# lowercase regardless of source case, so any routine not written in
# all-lowercase source (e.g. "Subroutine Init()") ends up with a mismatched
# undefined symbol at import time (e.g. "Init_" instead of "init_"). Force
# the first F_FUNC/F_FUNC_US argument to lowercase to match reality.
sed -i -E 's/(F_FUNC(_US)?\()([A-Za-z0-9_]+)(,)/\1\L\3\E\4/g' gen/aces2pymodule.c

echo "=== step 2: compile C sources ==="
mpicc -DNDEBUG -fPIC -O2 -I. -I$NUMPY_INC -I$PY_INC -I$F2PY_SRC \
    -c gen/aces2pymodule.c -o build/aces2pymodule.o
mpicc -DNDEBUG -fPIC -O2 -I. -I$NUMPY_INC -I$PY_INC -I$F2PY_SRC \
    -c $F2PY_SRC/fortranobject.c -o build/fortranobject.o

echo "=== step 3: compile Fortran sources ==="
FFLAGS="-i8 -fPIC -fp-model strict -O1 -assume minus0 -qopenmp -D__fortran -D__fortran77 -D_INTEL -D_RECL_IS_WORDS_ -I$ACES/include -I$ACES/acescore/include -I. -FI"

for f in $SOURCES; do
    base=$(basename "$f" .F)
    echo "compiling $base"
    ifort $FFLAGS -c "$f" -o build/${base}.o
done
ifort $FFLAGS -c gen/aces2py-f2pywrappers.f -o build/aces2py-f2pywrappers.o

echo "=== step 4: link ==="
# Generated fresh every build (not cached) so this can never point at a
# stale/dead location -- ACES is always this same build's own aces2 tree.
LIBDIRS=$(find $ACES -maxdepth 2 -iname "*.a" | xargs -n1 dirname | sort -u | sed 's/^/-L/' | tr '\n' ' ')
LIBNAMES=$(find $ACES -maxdepth 2 -iname "*.a" | xargs -n1 basename | sed -E 's/^lib(.*)\.a$/\1/' | sort -u | sed 's/^/-l/' | tr '\n' ' ')

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

# aces2's ~150 modules are historically full of same-named "generic" utility
# routines (independently written per-module, never meant to share a link
# unit -- the classic driver never links more than one module together, so
# these never collided there). aces2py.so links dozens of them into ONE
# process, so any such collision becomes real. Found empirically, one at a
# time, same as every other instance of this bug class in this project:
# lambda/finish.f's own SUBROUTINE FINISH(ICYCLE,pCCD) (2 args, unrelated to
# CC convergence) silently wins over vcc/finish.f's SUBROUTINE
# FINISH(ICYCLE,pCCD,pCCDS,pCCDTS[,PCCDTSD]) (5 args at vcc's call site --
# arity mismatch is itself benign/pre-existing, PCCDTSD is unused) when both
# libvcc.a and liblambda.a are linked together -- silently drops vcc's own
# "write the converged TOTENERG to JOBARC" step, leaving Runcc's reported
# energy stuck at the pre-convergence value mbptout.f wrote early on.
# "module:symbol" pairs to localize (module's OWN copy is neutralized so the
# OTHER module's same-named symbol wins the link) -- extend this list as new
# collisions are found instead of touching aces2's own sources. molcas also
# defines its own unrelated finish_ (no MAIN__ of its own, so it's not in
# MAIN_COLLISION_MODULES otherwise) -- same collision class, 3rd copy found.
EXTRA_LOCALIZE_SYMBOLS="lambda:finish_ molcas:finish_"
EXTRA_LOCALIZE_MODULES="molcas"

PATCHED_LIBDIR=$WORK/build/patched_libs
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
            nm "$o" 2>/dev/null | grep -q " T ${psym}\$" && objcopy --localize-symbol=${psym} "$o"
        done
    done
    (cd $workdir && ar rcs $PATCHED_LIBDIR/lib${mod}.a *.o)
    rm -rf $workdir
done

ifx -shared -nofor-main -o build/aces2py.cpython-312-x86_64-linux-gnu.so \
    build/*.o \
    -L$PATCHED_LIBDIR \
    $LIBDIRS \
    -Wl,--start-group $LIBNAMES -Wl,--end-group \
    -L$MKL -Wl,--start-group -lmkl_intel_ilp64 -lmkl_lapack95_ilp64 -lmkl_sequential -lmkl_core -Wl,--end-group \
    -lpthread

echo "=== step 5: test import ==="
cd build
python3 -c "import aces2py; print('IMPORT OK')"
