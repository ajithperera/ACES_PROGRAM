#!/bin/bash
set -e

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

SHIM=$SCRIPT_DIR/tools/compat_bin
export PATH=$SHIM:$PATH
# Driven from ../build.sh this must use the module set ./configure recorded
# (this script purges the environment build.sh already set up, so it has to
# be told); run standalone it falls back to the HiPerGator set.
module purge
if ! module load ${ACES_PROGRAM_MODULES:-intel/2025.1.0 openmpi/5.0.7 python/3.12}; then
    echo "module load failed: ${ACES_PROGRAM_MODULES:-intel/2025.1.0 openmpi/5.0.7 python/3.12}" >&2
    echo "Set ACES_PROGRAM_MODULES, or run ../configure --modules=\"...\" and ../build.sh pyaces" >&2
    exit 1
fi
: "${MKLROOT:?is unset -- the loaded modules provide no Intel MKL}"

WORK=$SCRIPT_DIR
MAIN=$WORK/Main
PYF=$WORK/a2_pyf
ACES=$SCRIPT_DIR/../aces2
NUMPY_INC=$(python3 -c "import numpy; print(numpy.get_include())")
PY_INC=$(python3 -c "import sysconfig; print(sysconfig.get_path('include'))")
# fortranobject.c/.h ship inside numpy itself, so ask the interpreter that is
# actually being built against rather than hardcoding one site's install path.
F2PY_SRC=$(python3 -c "import numpy, os; print(os.path.join(os.path.dirname(numpy.__file__), 'f2py', 'src'))")
MKL=$MKLROOT/lib/intel64
# oneAPI's own layout has the libraries directly under lib/ on some installs.
[ -d "$MKL" ] || MKL=$MKLROOT/lib

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
    $PYF/runee.F $PYF/rundirmp2.F $PYF/acesjafin.F $PYF/acesjainit.F"

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
FFLAGS="-diag-disable=10448 -i8 -fPIC -fp-model strict -O1 -assume minus0 -qopenmp -D__fortran -D__fortran77 -D_INTEL -D_RECL_IS_WORDS_ -I$ACES/include -I$ACES/acescore/include -I. -FI"

for f in $SOURCES; do
    base=$(basename "$f" .F)
    echo "compiling $base"
    ifort $FFLAGS -c "$f" -o build/${base}.o
done
ifort $FFLAGS -c gen/aces2py-f2pywrappers.f -o build/aces2py-f2pywrappers.o

# runee.o's own new CALL PROPS(...) (see the vprops:props_ comment above)
# must be redefined to match vprops's renamed copy the SAME way, since it's
# a pyaces-side .o, not one of aces2's own archives patched in the loop
# below.
objcopy --redefine-sym props_=pyaces_vprops_props_ build/runee.o

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
# joda/gschmidt.f is a 4th instance: 6 modules (libr, vcceh, joda, vea, vee,
# a vee/librv.skip backup) each define their own GSCHMIDT with DIFFERENT
# argument counts (vee's own davidtda.f callers pass 6 args matching vee's/
# vcceh's/vea's own 6-arg signature; joda's is a 5-arg version missing the
# leading VEC arg entirely -- a straight argument shift, exactly the shape
# that corrupts an address into an integer slot and segfaults, matching the
# observed vee_->drvtda_->davidtda_->gschmidt_ SIGSEGV backtrace). Confirmed
# via `grep -rl "call gschmidt"` that joda's OWN code never calls its own
# gschmidt (dead/unused there) -- zero risk to localize. libr/vcceh/vea DO
# have real internal callers of their own gschmidt (libr/david1.f,
# vcceh/{blcklineq,truncate}.f, vea/truncate.f) -- do NOT blindly localize
# those without confirming they aren't also needed by some other pyaces
# entry point first.
# libr's own copy (7-arg, extra TOL) turned out to be the actual winner
# once joda's was suppressed (confirmed via nm -S size-matching against the
# linked .so: 0x159 bytes matches libr/gschmidt.o exactly, not vee's
# 0x17b). libr/david1.f calls its own gschmidt internally and IS present in
# the linked .so (nm shows david1_ defined). vcceh/{blcklineq,truncate}.f
# and vea/truncate.f likewise each call their own gschmidt internally.
# 5th instance, same collision class again: MODF, called next inside vee_
# itself (vee.f:406/410/481/485, CALL MODF(ICORE,MAXCOR,IUHF,ITYPE), 4 args)
# after gschmidt_ was fixed and the crash moved past davidtda_. 8 modules
# define their own MODF -- lambda's has a genuinely different (longer)
# signature (extra IRREPT,IRREPF args, an arity mismatch matching the
# FINISH bug's shape); fsip/vcceh/hcmult/vea/lcct all share vee's own
# identical 4-arg shape but are still each a DIFFERENT module's own
# implementation, not guaranteed semantically identical just because the
# call shape matches. fsip/fs01.f, vcceh/vcceh.f, vea/{exitea,initea}.f all
# call their own MODF internally too.
# 6th instance: EXPDEN, next in the chain after MODF (vee_->doeomee_david_->
# nextdav_->tdens_->expden_, crashing inside a memset -- an out-of-bounds
# write, matching a wrong-copy/wrong-array-size symptom again). props's own
# copy has a COMPLETELY different signature/purpose (DTRUNC,DFULL,NORB,...
# vs vee's/vcceh's identical DOO,DVV,DVO,DOV,F,NMO,IRREPX,ISPIN,REFCONT).
# props/props.f and vcceh/cmpdso.f both call their own EXPDEN internally.
#
# 2026-08-11: switched this whole mechanism from --localize-symbol (hide
# the losing module's copy, forcing EVERYTHING -- including that module's
# OWN internal self-calls -- onto the winning module's copy) to
# --redefine-sym (rename the losing module's copy to a unique name,
# applied uniformly to every .o in its archive, definer AND any internal
# caller alike). Discovered why this matters the hard way: lambda/vlamcc.f
# calls its own 2-arg FINISH and its own 13-arg MODF internally
# (lambda:finish_/lambda:modf_ above) -- localizing those left vlamcc_'s
# own calls resolving to vcc's 4-arg FINISH / vee's 4-arg MODF instead of
# lambda's own, corrupting the args and SIGSEGVing inside Runee (found via
# a fresh gdb bt on case 054b: runee_->vlambda_->finish_). Every prior
# "safe to localize" judgment call above (libr/vcceh/vea's own gschmidt
# self-calls, fsip/vcceh/vea's own modf self-calls, props/vcceh's own
# expden self-calls) was exactly this same latent risk, just not yet hit by
# a live call path -- redefine-sym makes the whole reachability analysis
# moot: each losing module keeps calling ITS OWN renamed copy correctly,
# only the winning module's plain symbol name remains globally visible.
#
# 7th instance, but the FIRST triggered by ADDING a call rather than fixing
# a crash on an existing one: vprops/props.F's own SUBROUTINE PROPS(S,
# MAXCOR,IS,IMAXCOR) collides with props/props.F's SUBROUTINE PROPS(ICORE,
# MAXCOR,IUHF) (needed by Runprops) -- both modules just happen to be named
# "PROPS". Runee.F (a2_pyf, pyaces's OWN driver source, not aces2's) now
# calls vprops's PROPS directly (to populate DIPOLE_X/etc before Vee's own
# EOM transition-property code needs it, matching classic energy.f's own
# "call runit('xvprops')" placement before xvscf whenever EXCITE= is set).
# Since this is a NEW pyaces-side call site, not an existing aces2-internal
# one, redefine-sym must ALSO be applied to pyaces's own compiled runee.o
# (not just vprops's archive) -- see the loop right after step 3 below.
# Linking vprops in for the first time ALSO surfaced 5 more of its own
# helper routines colliding with liboeptools/nucatt.f's bundled same-named
# routines (setrhf_/rhftce_/fmc_/aainer_/civpt_) -- oeptools's copies are
# left as the global winners, vprops's own renamed so its internal call
# chain (props.f calling rhftce.f/fmc.f/aainer.f/civpt.f) stays self-bound.
#
# 8th instance: HBARXC, called right after NEWGES in vee/doeomee_david.F's
# own EOM Davidson iteration loop (the sigma-vector build step) -- crashed
# with "@GETLST: List (1,3) does not exist" for cases 054a/054b. vcceh AND
# lcct BOTH also define their own HBARXC (confirmed meaningfully different
# implementations via diff/size, not just a naming coincidence -- vee's own
# copy has LANCZOS/DAVID/BLOCK_DAVID-gated LISTT1/LISTT2/LISTT2IN branching
# that vcceh's/lcct's copies lack). All three modules ALSO have their own
# internal callers (vcceh: blcklineq.f/lineqy*.f/phbarxc2.f/secprp.f; lcct:
# clczeta2.f/calczeta.f/lineqy.f; vee: calczeta.f/clczeta2.f/diagall.f/
# doeomee_david.f) -- redefine, not localize, same reasoning as FINISH.
# Fixing HBARXC alone didn't fix 054a/054b -- a fresh gdb bt (via a
# temporary out-of-bounds write injected at the GETLST assertion site
# itself, since this is a clean aces_exit not a real crash -- ifx's own
# ABORT() intrinsic doesn't raise a real signal/produce a core, a genuine
# invalid memory write does) showed the REAL call chain: vee_->
# doeomee_david_->hbarxc_->dt2int2_->draolad_->t2toao_->getlst_. Every one
# of dt2int2/draolad/t2toao is ALSO a multi-module collision (9th, 10th,
# 11th instances): DT2INT2 in vcceh/hcmult/vee/lcct (4 genuinely different
# implementations, confirmed via `size`), DRAOLAD in vcceh/lambda/libr2/vee
# (vcceh's and vee's happen to be byte-identical; lambda's/libr2's differ),
# T2TOAO in vcceh/libr2/vee (vcceh's and vee's again byte-identical;
# libr2's differs), GETAOINF (12th instance) in vcceh/dens/vee/libr2
# (libr2's was missed on the first pass: it declares "subroutine getaoinf"
# in lowercase, a case-sensitive grep for "SUBROUTINE GETAOINF" silently
# skipped it -- confirmed the actual winner via
# `ld -Wl,-trace-symbol=<sym>`, which prints the exact archive member
# providing/referencing a symbol, more reliable than guessing from nm -S
# size comparisons alone; use this first for any FUTURE "which copy
# actually won" question).
#
# UNLIKE every other instance so far, this one is NOT a single global
# winner -- it's a genuine two-context ambiguity. draolad_/t2toao_/
# getaoinf_ are called from BOTH vee's own EOMEE Davidson path (via
# hbarxc_->dt2int2_, needs vee's own copy) AND vcc's own ABCDTYPE=AOBASIS
# path (via dre3en_->draolad_ directly, needs libr2's copy) -- confirmed
# by checking which library the CLASSIC per-executable build actually
# links for each (vcc/GNUmakefile's ACES_LINK_LISTS): xvcc's own link
# list is "vcc:trp:librt3:libr2:crust:pccd:...", no vee at all, so libr2
# must be what it uses; xvee's is "vee:trp:librt3:libr2:crust:hbar:vcc:
# ...", vee's own archive scanned FIRST so vee's own copy wins there.
# First attempt (renaming vcceh/lambda/libr2, letting vee win globally)
# fixed the vee-originated crash but broke the vcc-originated one instead
# (same @GETLST-class bug, now via getaoinf_'s own array-index arithmetic
# reading an uninitialized COMMON block vee's own Davidson setup would
# have populated but vcc's plain CC path never does). Fix: DON'T rename
# libr2's copies (let libr2 be the natural global winner, matching what
# vcc/most other modules actually need) -- instead rename VEE's OWN
# copies (vee:draolad_/t2toao_/getaoinf_), so vee's internal chain
# (hbarxc.f->dt2int2.f->draolad.f->t2toao.f->getaoinf.f, all within vee's
# own archive) stays correctly self-bound to ITS OWN implementations,
# while everyone else (vcc included) transparently gets libr2's. vcceh's
# copies of draolad_/t2toao_ (byte-identical to vee's OLD copies, now
# irrelevant either way) and dens's getaoinf_ (genuinely different, still
# wrong for anyone) stay renamed away as before.
#
# 13th instance, and the actual explanation for 054a/054b's own remaining
# small (~1e-5) MISMATCH after all of the above: libr2's own draolad.f
# calls Z2TOMO internally (the routine that prints "W(abef) AB
# contribution", directly comparable against classic's own log) -- but
# Z2TOMO ITSELF is a 4-way collision (vcceh/lambda/libr2/vee; tdee's and
# pccd's own "z2tomo"-named files are already safely prefixed at the
# source level as TDEE_Z2TOMO/PCCD_Z2TOMO, not real collisions) and
# lambda's copy was winning globally -- so libr2's OWN internal call fell
# through to lambda's DIFFERENT implementation, exactly the same
# self-caller trap as every other instance today. Same two-context split
# as draolad/t2toao/getaoinf (vee's own iabcont.f also calls its own
# z2tomo internally) -- rename vcceh/lambda/vee's copies, let libr2 win.
# 14th instance, found by gdb on suite case 021 (CALC=CCSDT) rather than by
# reading: runcc_ -> vcc_ -> trps_ -> w45t3drv_ -> SIGSEGV. W45T3DRV exists in
# trp (CORE,MAXCOR,IUHF,EXTRNL_CCSD -- 4 args, what trp/trps.f:599 and
# ccsdtqf/trps_pccd.f pass) and in lambda (CORE,MAXCOR,IUHF,LZOOFF,LZVOFF,
# LWOFF,LT3OFF,INIT -- 8 args, what lambda/vlamcc.f passes). lambda's copy was
# winning the link, confirmed by size: nm -S gives 0x2a74 = 10868 bytes in the
# linked .so, matching liblambda.a exactly, against libtrp.a's 0x22ee = 8942.
# So trp's 4-argument call landed in an 8-argument routine and the four
# missing argument addresses came from whatever was on the stack. Rename
# lambda's so trp's -- the one two different modules call -- wins globally;
# lambda/vlamcc.f keeps reaching its own through the rename, as everywhere
# else here. (trp/grad_trps.f calls it with only THREE arguments, which is a
# separate pre-existing arity bug in a path this suite does not reach; left
# alone deliberately rather than papered over from here.)
# 15th instance, and the first one this tree created rather than inherited:
# rpa was a .skip module until 723989ab re-enabled it, so its archive only
# started competing in this link after that commit -- which is why the
# zero-point baseline never saw it. CMPENG ("compute energy") exists in libr
# (0xb37 = 2871 bytes) and rpa (0xc79 = 3193), same 9-argument signature, so
# nothing crashes: rpa's simply won and vcc/eucc.f quietly got RPA's variant.
# They are not interchangeable -- libr's has the singles block
#   IF((METHOD.GT.9.AND.SING1.AND.METHOD.NE.21.AND.METHOD.NE.23)
#  &    .OR.NONHF.AND.SING1)THEN
# that accumulates ESING from the T1 amplitudes, while rpa's drops it and
# zeroes ESING_KEEP/EDBLS_KEEP instead. Caught on suite case 025 (CALC=UCC(4)):
# the correlation energy came out -0.227165476 against the serial run's
# -0.227215217, 4.97e-05 adrift before the UCC(4) correction even started, and
# the final energy missed the archive reference by 2.49e-05. vcc, lambda, pccd,
# ccsdtqf and ccsdtqf_mm all call CMPENG expecting libr's, so libr wins
# globally and rpa's own callers (mkt1at2.f, t2int1.f) follow the rename.
# 16th instance, and the first where the classic build itself resolves the
# same call two different ways depending on which executable you are in.
# XGEMM (the BLAS wrapper) lives in both crust and tools with identical
# argument lists but different degenerate-case handling: crust's returns
# without calling DGEMM at all when lda*ldb is zero (a zero-length operand
# contributes nothing) and sanitises with max(1,ld*) in the k=0 case, while
# tools' hands whatever it was given straight to DGEMM. ACES_LINK_LISTS puts
# crust ahead of LIST_DEF_LIBS (which is where tools lives) for vcc, vee and
# friends, so xvee gets crust's; a module without crust in its list gets
# tools'. Everything lands in one .so here, tools won, and suite case 057
# (EXCITE=EOMEE with ESTATE_PROP=EXPECTATION) died inside
# vee -> doeomee_david -> hbarxc -> dt1int2b -> xgemm with
# "** On entry to DGEMM parameter number 10 had an illegal value" -- parameter
# 10 is LDB, passed as zero. Let crust's win: it is the strictly more
# defensive of the two and identical on every non-degenerate call, so the
# modules that used to get tools' see no change.
# 17th instance, and the first one the audit script called before the test did:
# CHECK_T2 is trp's (33 lines, checks the T2 amplitudes against a stored copy)
# and ringcc's (115 lines, its own ring-CC variant). vcc's ACES_LINK_LISTS entry
# is vcc:trp:librt3:..., with no ringcc anywhere in it, so xvcc has always meant
# trp's; in the .so ringcc's won and vcc -> check_t2 -> getall asked for list 34,
# which does not exist under SINGLE_STORE=ON (clmoio.F builds that layout with
# INIT_NOREDNT_LIST instead of INIPCK). Cases 116a/b/c died on
# "@GETLST: List (1,34) does not exist". ccsdtqf, ccsdtqf_mm and pccd resolve to
# trp's as well, so trp wins globally and ringcc keeps its own by the rename.
# 18th instance, and the one that showed the audit script was only looking
# half way: a module that DEFINES a routine and calls it can still lose to
# somebody else's copy here, which the classic build makes impossible (the
# linker satisfies the call from the archive member it is already pulling in).
# vcc owns INITIN (252 bytes), INITSN, DIISLST and DODIIS0; lambda has its own
# larger versions (INITIN is 688 there) and won all four, so vcc ran lambda's
# initialisation and then asked for list 1, which the FNO run has not built:
# cases 124c/d/e/g.fno (FNO with CCSD(T)) died on "@GETLST: List (1,1) does not
# exist" while the five FNO cases that stop at CCSD were fine. ccsdtqf,
# ccsdtqf_mm, copt, pccd, psi4dbg, ringcc and ringcl all resolve these to vcc's
# in the classic link too, so vcc wins and lambda follows the rename.
EXTRA_LOCALIZE_SYMBOLS="lambda:finish_ molcas:finish_ joda:gschmidt_ libr:gschmidt_ vcceh:gschmidt_ vea:gschmidt_ lambda:modf_ fsip:modf_ vcceh:modf_ hcmult:modf_ vea:modf_ lcct:modf_ props:expden_ vcceh:expden_ vprops:props_ vprops:setrhf_ vprops:rhftce_ vprops:fmc_ vprops:aainer_ vprops:civpt_ vcceh:hbarxc_ lcct:hbarxc_ vcceh:dt2int2_ hcmult:dt2int2_ lcct:dt2int2_ vcceh:draolad_ lambda:draolad_ vee:draolad_ vcceh:t2toao_ vee:t2toao_ dens:getaoinf_ vee:getaoinf_ vcceh:z2tomo_ lambda:z2tomo_ vee:z2tomo_ lambda:w45t3drv_ tools:xgemm_ ringcc:check_t2_ lambda:initin_ lambda:initsn_ lambda:diislst_ lambda:dodiis0_"
EXTRA_LOCALIZE_MODULES="molcas libr vcceh vea fsip hcmult lcct vprops libr2 dens tools ringcc"

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
            newsym="pyaces_${pmod}_${psym}"
            objcopy --redefine-sym ${psym}=${newsym} "$o" 2>/dev/null || true
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

# Which copy of a multiply-defined routine won this link is not something the
# linker will ever complain about -- it just picks one and everybody lives with
# it. Every collision fixed above was found by a test failing, several of them
# by an energy quietly coming out wrong. Report the ones that are still there,
# and shout if a new one appeared (a re-enabled .skip module is the usual
# cause). Advisory only: it must not fail a build that otherwise linked.
echo "=== step 6: link-resolution audit ==="
python3 $WORK/tools/collision_audit.py || true
