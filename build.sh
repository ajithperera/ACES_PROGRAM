#!/bin/bash
set -e

SHIM=/tmp/claude-11138/-home-perera/e5b737e2-0fde-4c73-986c-761c6577ed12/scratchpad/compat_bin
export PATH=$SHIM:$PATH
module purge
module load intel/2025.1.0 openmpi/5.0.7 python/3.12

WORK=/blue/ufhpc/perera/acespy_build
MAIN=$WORK/Main
PYF=$WORK/a2_pyf
ACES=/blue/ufhpc/perera/acesii_build_pic
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
LIBDIRS=$(cat /tmp/claude-11138/-home-perera/e5b737e2-0fde-4c73-986c-761c6577ed12/scratchpad/aces2_libdirs.txt | sed 's/^/-L/' | tr '\n' ' ')
LIBNAMES=$(cat /tmp/claude-11138/-home-perera/e5b737e2-0fde-4c73-986c-761c6577ed12/scratchpad/aces2_libnames.txt | sed 's/^/-l/' | tr '\n' ' ')

ifx -shared -nofor-main -o build/aces2py.cpython-312-x86_64-linux-gnu.so \
    build/*.o \
    $LIBDIRS \
    -Wl,--start-group $LIBNAMES -Wl,--end-group \
    -L$MKL -Wl,--start-group -lmkl_intel_ilp64 -lmkl_lapack95_ilp64 -lmkl_sequential -lmkl_core -Wl,--end-group \
    -lpthread

echo "=== step 5: test import ==="
cd build
python3 -c "import aces2py; print('IMPORT OK')"
