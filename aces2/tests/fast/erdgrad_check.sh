# Sourced by the erd_grad_* SCRIPTs, with a ZMAT in the working directory:
# a one-cycle optimization (OPT_MAXCYC=1, "*ACES2" alone on its line, no
# INTEGRALS keyword). Runs it twice:
#   vmol/  the default route, xvmol + xvdint
#   .      INTEGRALS=ERD, xerdints + xerdgrad
# and fails unless both wrote a GRADIENT record and the records agree to
# 1e-8 hartree/bohr (and the TOTENERG records to 1e-9 hartree). The ERD JOBARC stays in the working directory, so the
# runner's TOTENERG check is on the ERD energy. OPT_MAXCYC=1 makes xaces2
# exit nonzero after the gradient, so the exit status is not checked.
test -n "$GENBAS" || exit 1
mkdir vmol || exit 1
cp ZMAT vmol || exit 1
test -f ECPDATA && { cp ECPDATA vmol || exit 1; }
(cd vmol && ln -s $GENBAS GENBAS && xaces2 >out 2>&1)
awk '{print} /^\*ACES2 *$/{print "INTEGRALS=ERD"}' vmol/ZMAT >ZMAT || exit 1
grep -q '^INTEGRALS=ERD' ZMAT || { echo "no *ACES2 line in ZMAT"; exit 1; }
ln -s $GENBAS GENBAS || exit 1
xaces2 >out 2>&1
grep -q 'Executing "xerdgrad"' out || { echo "xerdgrad did not run"; tail -15 out; exit 1; }
# jareq can print a warning first (the one-cycle job leaves JOBARC
# unfinalized); keep only the lines that are a single number.
num='^ *[-+]?[0-9.]+([eEdD][-+]?[0-9]+)? *$'
n=$(xa2proc jareq i NATOMS 1 | grep -E "$num" | tail -1)
n=$((3*n))
(cd vmol && xa2proc jareq d GRADIENT $n) | grep -E "$num" >g.vmol
xa2proc jareq d GRADIENT $n | grep -E "$num" >g.erd
paste g.vmol g.erd | awk -v n=$n '
  NF==2 {k++; d=$1-$2; if (d<0) d=-d; if (d>m) m=d}
  END {printf "GRADIENT: %d of %d components, max |ERD-VMOL| = %.2e\n", k, n, m
       exit !(k==n && m<1e-8)}' || exit 1
e1=$( (cd vmol && xa2proc jareq d TOTENERG 1) | grep -E "$num")
e2=$(xa2proc jareq d TOTENERG 1 | grep -E "$num")
echo "$e1 $e2" | awk '{d=$1-$2; if (d<0) d=-d
  printf "TOTENERG: VMOL %s  ERD %s  diff %.1e\n", $1, $2, d; exit !(NF==2 && d<1e-9)}' || exit 1
