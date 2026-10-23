# Files and data: JOBARC and the list files

Member programs share data through two kinds of files in the working
directory. Small named records go into `JOBARC`; large arrays of integrals,
amplitudes and densities go into five *list files*.

## JOBARC records

`JOBARC` is a direct-access file of named records; `JAINDX` holds the
index (names, positions and lengths). A record name is 1 to 8 characters.
Records are written with `PUTREC` and read with `GETREC` (both in
`acescore/`):

```fortran
      call putrec(iflag, 'JOBARC', 'RECNAME ', length, source)
      call getrec(iflag, 'JOBARC', 'RECNAME ', length, dest)
```

The `length` is always in **integer words**, whatever the data type: a
double-precision array of *n* elements has length `n*IINTFP`, where
`IINTFP` is the number of integers per double (from the `MACHSP` common
block). The current build compiles Fortran with `-i8`, so `IINTFP` is 1;
with 4-byte integers it would be 2. Always write `n*IINTFP`, never a
literal factor.

`PUTREC` with `iflag = 0` writes zeros instead of `source`. For `GETREC`:

| `iflag` | Behavior |
|---|---|
| > 0 | Read the record; stop the program if it does not exist |
| < 0 | Read the record; if it does not exist, fill `dest(1:length)` with zeros |
| 0 | Do not read; return the record's length in `length`, or -1 if it does not exist |

The `a2getrec`/`a2putrec` wrappers are used in the driver. By convention
records stay small, at most about the square of the number of basis
functions; anything larger belongs in a list file.

{doc}`jobarc` lists every record the code uses with a literal name, the
directories that write and read it, and, where known, what it holds.

## The list files

Large arrays are stored as *lists*. A list is a two-dimensional array of
doubles, addressed by a pair of indices `(ILEFT, IRIGHT)`: `IRIGHT` is the
list number and `ILEFT` the sublist, usually the irrep of the columns. Lists
live in five files of 100 lists each:

| Lists | File |
|---|---|
| 1–100 | `MOINTS` |
| 101–200 | `GAMLAM` |
| 201–300 | `MOABCD` |
| 301–400 | `DERINT` |
| 401–500 | `DERGAM` |

(`acescore/include/lists.h`: `_MAX_IO_FAMS` = 500 lists, `_FAMS_PER_FILE` =
100, up to `_MAX_IO_GRPS` = 10 sublists; the files are Fortran units 50–54.)

A list is created with `aces_list_touch` (older code calls `updmoi`, which
is kept for compatibility), and read or written a range of columns at a
time:

```fortran
c     columns iFirst .. iFirst+nCols-1 of list (iLeft, iRight)
      call getlst(dest, iFirst, nCols, 1, iLeft, iRight)
      call putlst(src,  iFirst, nCols, 1, iLeft, iRight)
```

The fourth argument is obsolete and ignored. The I/O system caches
physical records in memory (`CACHE_RECS`, `FILE_RECSIZ`); `INCORE`
keeps whole lists in memory instead.

### Symmetry distribution types

Rows and columns of a list are *distributions*: pairs of orbital indices
with a given total symmetry. Each list records the type of its row and
column distributions (set up with `INIPCK`, stored in `JOBARC` record
`ISYMTYP`), using these codes. Capital letters are alpha spin and lower case
beta; *A, B* are virtual and *I, J* occupied orbitals.

| Type | Pair | Type | Pair | Type | Pair |
|---|---|---|---|---|---|
| 1 | A < B | 9 | A I | 17 | i a |
| 2 | a < b | 10 | a i | 18 | I a |
| 3 | I < J | 11 | A i | 19 | A B |
| 4 | i < j | 12 | a I | 20 | a b |
| 5 | A ≤ B | 13 | A b | 21 | I J |
| 6 | a ≤ b | 14 | I j | 22 | i j |
| 7 | I ≤ J | 15 | A b (same as 13) | | |
| 8 | i ≤ j | 16 | I A | | |

For example, list 16, ⟨Ab|Ij⟩ stored as (A,b ; I,j), is created with
`INIPCK(1,13,14,16,...)`: rows of type 13, columns of type 14.

### List numbers

The table below is the catalog of list numbers from the 1990s *ACES II
Programmer's Manual*. Most entries still hold, but it predates many
programs (pCCD, ring CC, EOM extensions, the higher-order CC codes), and the
file layout has changed since (it put a file `SECDER` at lists 301–400).
Check the code before relying on an entry. *Sublists* gives the range of
`ILEFT`; the last column marks lists that RHF calculations do not use.
Storage is given as (row indices ; column indices).

| Sublists | List | Contents (rows ; columns) | UHF/ROHF only |
|---|---|---|---|
| 1-NIRREP | 5 | &lt;AB\|IJ> A,I ; B,J (NOT ANTISYMMETRIZED) |  |
| 1-NIRREP | 6 | &lt;ab\|ij> a,i ; b,j (NOT ANTISYMMETRIZED) |  |
| 1-NIRREP | 11 | &lt;IJ\|\|KL> I&lt;J ; K&lt;L | ✓ |
| 1-NIRREP | 12 | &lt;ij\|\|kl> i&lt;j ; k&lt;l | ✓ |
| 1-NIRREP | 13 | &lt;Ij\|\|Kl> I,j ; K,l |  |
| 1-NIRREP | 7 | &lt;IJ\|\|KA> I&lt;J ; K,A |  |
| 1-NIRREP | 8 | &lt;ij\|\|ka> i&lt;j ; k,a | ✓ |
| 1-NIRREP | 9 | &lt;Ij\|Ak> I,j ; A,k | ✓ |
| 1-NIRREP | 10 | &lt;Ij\|Ka> I,j ; K,a |  |
| 1-NIRREP | 14 | &lt;AB\|\|IJ> A&lt;B ; I&lt;J |  |
| 1-NIRREP | 15 | &lt;ab\|\|ij> a&lt;b ; i&lt;j | ✓ |
| 1-NIRREP | 16 | &lt;Ab\|Ij> A,b ; I,j |  |
| 1-NIRREP | 17 | &lt;Ab\|Ij> b,j ; A,I | ✓ |
| 1-NIRREP | 18 | &lt;Ab\|Ij> A,I ; b,j |  |
| 1-NIRREP | 19 | &lt;AB\|\|IJ> A,I ; B,J |  |
| 1-NIRREP | 20 | &lt;ab\|\|ij> a,i ; b,j | ✓ |
| 1-NIRREP | 21 | &lt;Ab\|Ij> A,j ; b,I |  |
| 1-NIRREP | 22 | &lt;Ab\|Ij> b,I ; A,j | ✓ |
| 1-NIRREP | 23 | &lt;IA\|\|JB> B,I ; A,J |  |
| 1-NIRREP | 24 | &lt;ia\|\|jb> b,i ; a,j | ✓ |
| 1-NIRREP | 25 | &lt;iA\|jB> B,i ; A,j |  |
| 1-NIRREP | 26 | &lt;aI\|bJ> b,I ; a,J | ✓ |
| 1-NIRREP | 27 | &lt;AB\|\|CI> A&lt;B ; C,I |  |
| 1-NIRREP | 28 | &lt;ab\|\|ci> a&lt;b ; c,i | ✓ |
| 1-NIRREP | 29 | &lt;Ab\|Ic> A,b ; I,c | ✓ |
| 1-NIRREP | 30 | &lt;Ab\|Ci> A,b ; C,i |  |
| 1-NIRREP | 34 | T2(IJ,AB) A,J ; B,I |  |
| 1-NIRREP | 35 | T2(ij,ab) a,j ; b,i | ✓ |
| 1-NIRREP | 36 | T2(Ij,Ab) b,j ; A,I | ✓ |
| 1-NIRREP | 37 | T2(Ij,Ab) A,I ; b,j |  |
| 1-NIRREP | 38 | T2(Ij,Ab) b,I ; A,j | ✓ |
| 1-NIRREP | 39 | T2(Ij,Ab) A,j ; b,I |  |
| 1-NIRREP | 40 | T2(IJ,AB) (increment)A,J ; B,I (rings) | ✓ |
| 1-NIRREP | 41 | T2(ij,ab) (increment)a,j ; b,i (rings) | ✓ |
| 1-NIRREP | 42 | T2(Ij,Ab) (increment)A,I ; b,j (rings) |  |
| 1-NIRREP | 43 | T2(Ij,Ab) (increment)A,j ; b,I (rings) |  |
| 1-NIRREP | 44 | T2(IJ,AB) A&lt;B ; I&lt;J |  |
| 1-NIRREP | 45 | T2(ij,ab) a&lt;b ; i&lt;j | ✓ |
| 1-NIRREP | 46 | T2(Ij,Ab) A,b ; I,j |  |
| 1-NIRREP | 48 | D(IJ,AB) A&lt;B ; I&lt;J |  |
| 1-NIRREP | 49 | D(ij,ab) a&lt;b ; i&lt;j | ✓ |
| 1-NIRREP | 50 | D(Ij,Ab) A,b ; I,j |  |
| 1-NIRREP | 51 | W(MN,IJ) intermed. M&lt;N ; I&lt;J | ✓ |
| 1-NIRREP | 52 | W(mn,ij) intermed. m&lt;n ; i&lt;j | ✓ |
| 1-NIRREP | 53 | W(Mn,Ij) intermed. M,n ; I,j |  |
| 1-NIRREP | 54 | W(MB,EJ) intermed. E,M ; B,J |  |
| 1-NIRREP | 55 | W(mb,ej) intermed. e,m ; b,j | ✓ |
| 1-NIRREP | 56 | W(Mb,Ej) intermed. E,M ; b,j |  |
| 1-NIRREP | 57 | W(mB,eJ) intermed. e,m ; B,J | ✓ |
| 1-NIRREP | 58 | W(mB,Ej) intermed. E,m ; B,j |  |
| 1-NIRREP | 59 | W(Mb,eJ) intermed. e,M ; b,J | ✓ |
| 1-NIRREP | 61 | T2(IJ,AB) increment A&lt;B ; I&lt;J; also L2(IJ,AB) increment A&lt;B ; I&lt;J in Lambda iterations; also T2(IJ,AB)[3] A&lt;B ; I&lt;J iterative ROHF | ✓ |
| 1-NIRREP | 62 | T2(ij,ab) increment a&lt;b ; i&lt;j; also L2(ij,ab) increment a&lt;b ; i&lt;j in Lambda iterations; also T2(ij,ab)[3] a&lt;b ; i&lt;j iterative ROHF | ✓ |
| 1-NIRREP | 63 | T2(Ij,Ab) increment A,b ; I,j; also L2(Ij,Ab) increment A,b ; I,j in Lambda iterations; also T2(Ij,Ab)[3] A,b ; I,j iterative ROHF |  |
| 1-NIRREP | 64 | Reciprocal D(IJ,AB) A&lt;B ; I&lt;J |  |
| 9 | 64 | Reciprocal D(I,A) A,I |  |
| 1-NIRREP | 65 | Reciprocal D(ij,ab) a&lt;b ; i&lt;j | ✓ |
| 9 | 65 | Reciprocal D(i,a) a,i | ✓ |
| 1-NIRREP | 66 | Reciprocal D(Ij,Ab) A,b ; I,j |  |
| 1-3 | 70-89 | Reserved for RLE Jacobi iterates |  |
| 1 | 90 | T1AA A,I |  |
| 2 | 90 | T1BB a,i | ✓ |
| 3 | 90 | T1AA AND L1AA INCREM A,I |  |
| 4 | 90 | T1BB AND L1BB INCREM a,i | ✓ |
| 5 | 90 | AB MO OVERLAP MATRIX I,j (o-o) | ✓ |
| 6 | 90 | AB MO OVERLAP MATRIX A,j (v-o) | ✓ |
| 7 | 90 | AB MO OVERLAP MATRIX I,b (o-v) | ✓ |
| 8 | 90 | AB MO OVERLAP MATRIX A,b (v-v) | ✓ |
| 9 | 90 | T1AA[2]; iterative A,I ROHF |  |
| 10 | 90 | T1BB[2]; iterative a,i ROHF |  |
| 1 | 91 | F(MI) intermed. M,I |  |
| 2 | 91 | F(mi) intermed. m,i | ✓ |
| 3 | 91 | f(M,I) (fock matrix) M,I (only non-HF methods) |  |
| 4 | 91 | f(m,i) (fock matrix) m,i (only non-HF methods) | ✓ |
| 1 | 92 | F(EA) intermed. E,A |  |
| 2 | 92 | F(ea) intermed. e,a | ✓ |
| 3 | 92 | f(E,A) (fock matrix) E,A (only non-HF methods) |  |
| 4 | 92 | f(e,a) (fock matrix) e,a (only non-HF methods) | ✓ |
| 1 | 93 | F(AI) intermed. A,I |  |
| 2 | 93 | F(ai) intermed. a,i | ✓ |
| 3 | 93 | f(A,I) (fock matrix) A,I (only non-HF methods) |  |
| 4 | 93 | f(a,i) (fock matrix) a,i (only non-HF methods) | ✓ |
| 1-NIRREP | 94 | T2(IJ,AB)[2] iterative A&lt;B ; I&lt;J ROHF |  |
| 1-NIRREP | 95 | T2(ij,ab)[2] iterative a&lt;b ; i&lt;j ROHF |  |
| 1-NIRREP | 96 | T2(Ij,Ab)[2] iterative A,b ; I,j ROHF |  |
| 1-NIRREP | 98 | T2(Ij,Ab) cumulative A,b ; I,j (MBPT S^2) ROHF |  |
| 1-NIRREP | 107 | G(IJ,KA) I&lt;J ; K,A; also W(IJ,KA) I&lt;J ; K,A in VLAMCC | ✓ |
| 1-NIRREP | 108 | G(ij,ka) i&lt;j ; k,a; also W(ij,ka) i&lt;j ; k,a in VLAMCC | ✓ |
| 1-NIRREP | 109 | G(Ij,Ak) I,j ; A,k; also W(Ij,Ak) I,j ; A,k in VLAMCC | ✓ |
| 1-NIRREP | 110 | G(Ij,Ka) I,j ; K,a; also W(Ij,Ka) I,j ; K,a in VLAMCC |  |
| 1-NIRREP | 111 | G(IJ,KL) I&lt;J ; K&lt;L | ✓ |
| 1-NIRREP | 112 | G(ij,kl) i&lt;j ; k&lt;l | ✓ |
| 1-NIRREP | 113 | G(Ij,Kl) I,j ; K,l |  |
| 1-NIRREP | 114 | G(IJ,AB) I&lt;J ; A&lt;B |  |
| 1-NIRREP | 115 | G(ij,ab) i&lt;j ; a&lt;b | ✓ |
| 1-NIRREP | 116 | G(Ij,Ab) I,j ; A,b |  |
| 1-NIRREP | 117 | G(iA,bJ) b,j ; A,I; also resorted in VDENS to b,I ; A,j | ✓ |
| 1-NIRREP | 118 | G(Ia,Bj) B,J ; a,i; also resorted in VDENS to B,i ; a,J |  |
| 1-NIRREP | 123 | G(IA,JB) B,I ; A,J; also resorted in VDENS to B,J ; A,I |  |
| 1-NIRREP | 124 | G(ia,jb) b,i ; a,j; also resorted in VDENS to b,j ; a,i | ✓ |
| 1-NIRREP | 125 | G(iA,jB) B,i ; A,j; also resorted in VDENS to B,j ; A,i |  |
| 1-NIRREP | 126 | G(Ia,Jb) b,I ; a,J; also resorted in VDENS to b,J ; a,I | ✓ |
| 1-NIRREP | 127 | G(AB,CI) A&lt;B ; C,I; also W(AB,CI) A&lt;B ; C,I in VLAMCC | ✓ |
| 1-NIRREP | 128 | G(ab,ci) a&lt;b ; c,i; also W(ab,ci) a&lt;b ; c,i in VLAMCC | ✓ |
| 1-NIRREP | 129 | G(Ab,Ic) A,b ; I,c; also W(Ab,Ic) A,b ; I,c in VLAMCC | ✓ |
| 1-NIRREP | 130 | G(Ab,Ci) A,b ; C,i; also W(Ab,Ci) A,b ; C,i in VLAMCC |  |
| 1-NIRREP | 131 | G(AB,CD) A&lt;B ; C&lt;D | ✓ |
| 1-NIRREP | 132 | G(ab,cd) a&lt;b ; c&lt;d | ✓ |
| 1-NIRREP | 133 | G(Ab,Cd) A,b ; C,d |  |
| 1-NIRREP | 134 | L2(IJ,AB) A,J ; B,I |  |
| 1-NIRREP | 135 | L2(ij,ab) a,j ; b,i | ✓ |
| 1-NIRREP | 136 | L2(Ij,Ab) b,j ; A,I | ✓ |
| 1-NIRREP | 137 | L2(Ij,Ab) A,I ; b,j |  |
| 1-NIRREP | 138 | L2(Ij,Ab) b,I ; A,j | ✓ |
| 1-NIRREP | 139 | L2(Ij,Ab) A,j ; b,I |  |
| 1-NIRREP | 144 | L2(IJ,AB) A&lt;B ; I&lt;J |  |
| 1-NIRREP | 145 | L2(ij,ab) a&lt;b ; i&lt;j | ✓ |
| 1-NIRREP | 146 | L2(Ij,Ab) A,b ; I,j |  |
| | | **Derivatives of one-electron integrals in AO basis** | |
| 1-NIRREP | 150 | S(mu,nu)^chi ipert ; mu >= nu |  |
| 1-NIRREP | 151 | f(mu,nu^(chi) (alpha)ipert ; mu >= nu |  |
| 1-NIRREP | 152 | f(mu,nu^(chi) (beta) ipert ; mu >= nu | ✓ |
| 1-3 | 153 | dipole(mu,nu) ixyz ; mu >= nu |  |
| 1-NIRREP | 154 | A^(chi)*DREL (alpha)ipert ; mu >= nu |  |
| 1-NIRREP | 155 | A^(chi)*DREL (beta) ipert ; mu >= nu | ✓ |
| | | **Relaxed density matrix (correlation correction only)** | |
| 1 | 160 | D(I,J) occ.-occ. block of relaxed density matrix (alpha) |  |
| 2 | 160 | D(i,j) occ.-occ. block of relaxed density matrix (beta) | ✓ |
| 3 | 160 | D(A,B) virt.-virt. block of relaxed density matrix (alpha) |  |
| 4 | 160 | D(a,b) virt.-virt. block of relaxed density matrix (beta) | ✓ |
| 5 | 160 | D(A,I) occ.-virt. block of relaxed density matrix (alpha) |  |
| 6 | 160 | D(a,i) occ.-virt. block of relaxed density matrix (beta) | ✓ |
| | | **Intermediate matrix I(p,q) (correlation correction only)** | |
| 1 | 161 | I(I,J) occ.-occ. block of intermediate matrix (alpha) |  |
| 2 | 161 | I(i,j) occ.-occ. block of intermediate matrix (beta) | ✓ |
| 3 | 161 | I(A,B) virt.-virt. block of intermediate matrix (alpha) |  |
| 4 | 161 | I(a,b) virt.-virt. block of intermediate matrix (beta) | ✓ |
| 5 | 161 | I(A,I) virt.-occ. block of intermediate matrix (alpha) |  |
| 6 | 161 | I(a,i) virt.-occ. block of intermediate matrix (beta); also Gamme-intermediates created in the triple code for CCSD(T), ... | ✓ |
| 1-NIRREP | 162 | G^trip(IJ,KA) I&lt;J ; K,A | ✓ |
| 1-NIRREP | 163 | G^trip(ij,ka) i&lt;j ; k,a | ✓ |
| 1-NIRREP | 164 | G^trip(Ij,Ak) I,j ; A,k | ✓ |
| 1-NIRREP | 165 | G^trip(Ij,Ka) I,j ; K,a |  |
| 1-NIRREP | 166 | G^trip(AB,CI) A&lt;B ; C,I | ✓ |
| 1-NIRREP | 167 | G^trip(ab,ci) a&lt;b ; c,i | ✓ |
| 1-NIRREP | 168 | G^trip(Ab,Ic) A,b ; I,c | ✓ |
| 1-NIRREP | 169 | G^trip(Ab,Ci) A,b ; C,i |  |
| | | **Derivatives of Overlap and Fock matrices and CPHF coefficients in MO basis** | |
| 1-NIRREP | 170 | S(I,J)^chi |  |
| 1-NIRREP | 171 | S(i,j)^chi | ✓ |
| 1-NIRREP | 172 | S(A,B)^chi |  |
| 1-NIRREP | 173 | S(a,b)^chi | ✓ |
| 1-NIRREP | 174 | S(A,I)^chi |  |
| 1-NIRREP | 175 | S(a,i)^chi | ✓ |
| 1-NIRREP | 176 | F(I,J)^chi |  |
| 1-NIRREP | 177 | F(i,j)^chi | ✓ |
| 1-NIRREP | 178 | F(A,B)^chi |  |
| 1-NIRREP | 179 | F(a,b)^chi | ✓ |
| 1-NIRREP | 180 | F(A,I)^chi |  |
| 1-NIRREP | 181 | F(a,i)^chi | ✓ |
| 1-NIRREP | 182 | U(A,I)^chi |  |
| 1-NIRREP | 183 | U(a,i)^chi | ✓ |
| 1-NIRREP | 184 | additional lists used in CPHF (RHF and UHF) |  |
| 1-NIRREP | 185 | additional lists used in CPHF (UHF only) | ✓ |
| 1-NIRREP | 186 | A^(chi)*D (I,J) |  |
| 1-NIRREP | 187 | A^(chi)*D (i,j) | ✓ |
| 1-NIRREP | 188 | A^(chi)*D (A,I) |  |
| 1-NIRREP | 189 | A^(chi)*D (a,i) | ✓ |
| 1 | 190 | L1AA A,I |  |
| 2 | 190 | L1BB a,i | ✓ |
| 1 | 191 | G(MI) intermed. M,I |  |
| 2 | 191 | G(mi) intermed. m,i | ✓ |
| 1 | 192 | G(A,E) intermed. E,A |  |
| 2 | 192 | G(a,e) intermed. e,a | ✓ |
| 1-NIRREP | 231 | &lt;AB\|\|CD> A&lt;B ; C&lt;D; also W(AB,CD) A&lt;B ; C&lt;D DURING VLAMCC | ✓ |
| 1-NIRREP | 232 | &lt;ab\|\|cd> a&lt;b ; c&lt;d; also W(ab,cd) a&lt;b ; c&lt;d DURING VLAMCC | ✓ |
| 1-NIRREP | 233 | &lt;Ab\|Cd> A,b ; C,d; also &lt;Ab\|Cd> A&lt;=b; C,d single point RHF; also W(Ab,Cd) A,b ; C,d DURING VLAMCC |  |
| | | **Additional lists for MBPT(2) second derivatives** | |
| | | **These are** | |
| 1-NIRREP | 314 | d&lt;IJ\|\|AB>/dx A&lt;B ; I&lt;J,IPERT | ✓ |
| 1-NIRREP | 315 | d&lt;ij\|\|ab>/dx a&lt;b ; i&lt;j,IPERT | ✓ |
| 1-NIRREP | 316 | d&lt;Ij\|\|Ab>/dx A,b ; I,j,IPERT |  |
| | | **Part II: extra lists for excitation energies** | |
| 1-NIRREP | 434 | C(IJ,AB) A,J ; B,I [INCREMENT] | ✓ |
| 1-NIRREP | 435 | C(ij,ab) a,j ; b,i [INCREMENT] | ✓ |
| 1-NIRREP | 436 | C(Ij,Ab) b,j ; A,I [INCREMENT] | ✓ |
| 1-NIRREP | 437 | C(Ij,Ab) A,I ; b,j [INCREMENT] |  |
| 1-NIRREP | 438 | C(Ij,Ab) b,I ; A,j [INCREMENT] | ✓ |
| 1-NIRREP | 439 | C(Ij,Ab) A,j ; b,I [INCREMENT] |  |
| 1-NIRREP | 444 | C(AB,IJ) A&lt;B ; I&lt;J | ✓ |
| 1-NIRREP | 445 | C(ab,ij) a&lt;b ; i&lt;j | ✓ |
| 1-NIRREP | 446 | C(Ab,Ij) A,b ; I,j |  |
| 1-NIRREP | 448 | D(IJ,AB) A&lt;B ; I&lt;J |  |
| 1-NIRREP | 449 | D(ij,ab) a&lt;b ; i&lt;j | ✓ |
| 1-NIRREP | 450 | D(Ij,Ab) A,b ; I,j |  |
| 1-NIRREP | 461 | C(AB,IJ) A&lt;B ; I&lt;J [INCREMENT] | ✓ |
| 1-NIRREP | 462 | C(ab,ij) a&lt;b ; i&lt;j [INCREMENT] | ✓ |
| 1-NIRREP | 463 | C(Ab,Ij) A,b ; I,j [INCREMENT]; also 1-NRECS 470 C vectors (one/record); also 1-NRECS 471 HC vectors (one/record) |  |
| 1 | 472 | Converged right C vector (on one record) |  |
| 2 | 472 | Converged left C vector (on one record) |  |
| 1 | 490 | C1(AI) AI ; |  |
| 2 | 490 | C1(ai) ai ; |  |
| 3 | 490 | C1(AI) AI ; [INCREMENT] |  |
| 4 | 490 | C1(ai) ai ; [INCREMENT] | ✓ |
| 1 | 491 | X(IJ) IJ ; |  |
| 2 | 491 | X(ij) ij ; | ✓ |
| 1 | 492 | X(AB) AB ; |  |
| 2 | 492 | X(ab) ab ; | ✓ |
| 1 | 493 | X(AI) AB ; |  |
| 2 | 493 | X(ai) ab ; | ✓ |
| 1-NIRREP | 94 | TDA eigenvectors (one/record) |  |
| 1-NIRREP | 95 | TDA eigenvalues (one/record) |  |

In analytic gradients, `xanti` and `xbcktrn` reorganize `MOINTS` and
`GAMLAM` to hold the two-particle density, first in the MO basis (lists 1–4
and 6–9 for the spin cases AAAA, AABB, ABAB, ABCD, sublists 1 to the number
of symmetry types, with the beta-spin lists at sublist 51 and up) and then
backtransformed to the AO basis.

## Symmetry records

`xjoda` finds two point groups: the *full* point group of the molecule and
the Abelian *computational* point group the calculation runs in. Their
data are stored in parallel `JOBARC` records with the prefixes `FULL` and
`COMP`:

| Records | Contents |
|---|---|
| `FULLPTGP`, `COMPPTGP` | Point-group symbols |
| `FULLORDR`, `COMPORDR` | Group orders (read these first) |
| `FULLSYOP`, `COMPSYOP` | 3 × 3 Cartesian matrices of all operations, ordered by class |
| `FULLPERM`, `COMPPERM` | Permutation vectors: element *k* of the vector for operation *i* is the atom that atom *k* maps to (NATOMS × order) |
| `FULLNORB`, `COMPNORB` | Number of orbits (sets of symmetry-equivalent atoms) |
| `FULLSTGP`, `COMPSTGP` | Site group of each orbit |
| `FULLPOPV`, `COMPPOPV` | Number of atoms in each orbit |
| `FULLMEMB`, `COMPMEMB` | Atoms ordered by orbit (use with `*POPV`) |
