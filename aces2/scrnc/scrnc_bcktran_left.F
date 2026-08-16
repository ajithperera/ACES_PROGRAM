
      Subroutine Scrnc_bcktran_left(JJ_respn_oo,JJ_respn_xx,
     +                              Eveca,Scr,
     +                              Nij_rows,Nxj_rows,Nxx_rows,
     +                              Nao,Nmo,Irrep_l,Ispin,Iuhf)

C Transforms a full occ-occ MO pair matrix (both members of the pair,
C spin Ispin, overall pair symmetry Irrep_l) to the AO-pair basis, per
C irrep block:
C   JJ_respn_xx(X,Y) = Sum_ij C(X,i) JJ_respn_oo(i,j) C(Y,j)
C i.e. Eveca*Z*Eveca^T -- the same convention as the general ACES II
C bcktrn module's SYTRAA (aces2/bcktrn/sytraa.f: "Q(XX) = C(X,A) Z(A,I)
C C^T(X,I)"), confirmed 2026-08-16 by reading sytraa.f/bxaaaa2.f: no AO
C overlap-matrix (S) correction anywhere in that module. This is NOT
C the algebraic inverse of the AO2MO integral transform used elsewhere
C in scrnc (scrnc_i1tran.F's ICASE=3, "M_MO = C^T M_AO C"), which would
C need an S-metric sandwich to invert exactly in a non-orthogonal AO
C basis -- it is the separate, standard convention for re-expressing a
C density/response-like tensor (not a raw integral) in a different
C one-particle basis.
C
C 2026-08-16, rewritten to fix two real indexing bugs found while
C re-deriving this against scrnc_i1tran.F's ICASE=3 block (the trusted,
C already-working AO2MO reference this is the adjoint of -- same
C AOSYM/SYM/SYMPOP offset conventions):
C 1. The per-irrep offset into the input MO pair array used to be
C    "Ioff_ji = Ioff_ji + Nsum" (Nsum = Pop(Irrep_j,Ispin) alone) --
C    wrong, since each Irrep_i block of the input actually holds
C    Pop(Irrep_j,Ispin)*Pop(Irrep_i,Ispin) elements (a full rectangular
C    sub-block, not just one row's worth). Fixed to increment by the
C    full block size.
C 2. The column offsets used to slice Eveca ("Ioff_j"/"Ioff_i") were
C    set to raw MO population COUNTS (Pop(Irrep,Ispin)) instead of
C    actual MO column OFFSETS (Ioffo(Irrep,Ispin), the same array
C    scrnc_i1tran.F uses directly as a Fortran column index into
C    EVEC1(Nao,Nmo)).
C Both bugs are silently invisible on a single-irrep (sym=off) system,
C since Nirrep=1 makes every offset trivially 0/1 regardless -- they
C only bite once real point-group symmetry is turned on. Also added
C explicit Ispin/Iuhf arguments -- this used to be hardcoded to alpha
C (SCFEVECA, Pop(*,1)) unconditionally, wrong for the UHF-only systems
C scrnc requires (Scrnc_set_vars.F rejects RHF outright). No cross-spin
C (AABB/ABAB) blocks are built here since the only currently-planned
C consumer (a Green's-function-type quantity) is spin-diagonal; a
C genuine two-particle-density consumer would need those added the way
C aces2/bcktrn/bxaaaa2.f's ITYPE loop does.

      Implicit Double Precision (A-H,O-Z)

      Double Precision JJ_respn_oo(Nij_rows)
      Double Precision JJ_respn_xx(Nxx_rows)
      Double Precision Scr(Nxj_rows)
      Dimension Eveca(Nao,Nmo)
      Character*8 Labelx(2)

      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),ioffv(8,2),ioffo(8,2),
     &             IRPDPDAO(8),IRPDPDAOMO_OCCBK(8,2),
     &             IRPDPDAOMO_VRTBK(8,2),IRPDPDAOMO_OCCKB(8,2),
     &             IRPDPDAOMO_VRTKB(8,2),
     &             IRPDPDAOS(8),
     &             ISTART(8,8),ISTARTMO(8,3)
#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "sym.com"
#include "flags2.com"
#include "flags.com"
#include "parcc.com"
#include "info.com"

      Data Zilch,One /0.0D0,1.0D0/
      Data Labelx /'SCFEVECA','SCFEVECB'/

C Same RHF/UHF eigenvector-record convention as scrnc_i1tran.F: for
C RHF (Iuhf=0) both spins share SCFEVECA; for UHF, spin 2 (beta) reads
C SCFEVECB.

      Ilabel = Ispin
      If (Iuhf .EQ. 0) Ilabel = 1
      Call Getrec(20,"JOBARC",Labelx(Ilabel),Nao*Nmo*Iintfp,Eveca)

      Ioffoo = 1
      Ioffxx = 1

      Do Irrep_j = 1, Nirrep
         Irrep_i = Dirprd(Irrep_j,Irrep_l)

         Numi = Pop(Irrep_i,Ispin)
         Numj = Pop(Irrep_j,Ispin)
         Nao1 = Iaopop(Irrep_i)
         Nao2 = Iaopop(Irrep_j)

         Ioffr1 = Ioffo(Irrep_i,Ispin)
         Ioffl1 = Ioffao(Irrep_i)
         Ioffr2 = Ioffo(Irrep_j,Ispin)
         Ioffl2 = Ioffao(Irrep_j)

C Step 1: Scr(X1,J) = Sum_I Eveca(X1,I) * Z(I,J)   [X1 in irrep_i AO,
C I in irrep_i MO, J in irrep_j MO -- transforms the FIRST pair member]

         If (Min(Nao1,Numj,Numi) .NE. 0) Then
            Call Xgemm("N","N",Nao1,Numj,Numi,One,
     +                  Eveca(Ioffl1,Ioffr1),Nao,JJ_respn_oo(Ioffoo),
     +                  Numi,Zilch,Scr,Nao1)
         Else
            Call Dzero(Scr,Nao1*Numj)
         Endif

C Step 2: JJ_respn_xx(X1,X2) = Sum_J Scr(X1,J) * Eveca(X2,J)^T
C [X2 in irrep_j AO -- transforms the SECOND pair member]

         If (Min(Nao1,Nao2,Numj) .NE. 0) Then
            Call Xgemm("N","T",Nao1,Nao2,Numj,One,Scr,Nao1,
     +                  Eveca(Ioffl2,Ioffr2),Nao,Zilch,
     +                  JJ_respn_xx(Ioffxx),Nao1)
         Else
            Call Dzero(JJ_respn_xx(Ioffxx),Nao1*Nao2)
         Endif

C Self-contained correctness check (2026-08-16, kept as a permanent
C diagnostic, not exploratory cruft): independently recompute this same
C per-block Xgemm result via brute-force nested loops -- direct
C evaluation of Result(X1,X2) = Sum_I Sum_J Eveca(X1,I)*Z(I,J)*Eveca(X2,J)
C with no Xgemm, no leading-dimension arithmetic, nothing shared with
C the code above except the raw offsets -- and report the max
C difference. Catches argument-order/leading-dimension bugs in the
C Xgemm calls above, which is exactly the bug class this routine had
C twice already (see the file header). Nirrep=1 (sym=off) test systems
C can't exercise the multi-irrep offset logic (Ioffo/Ioffao collapse to
C trivial values), so this check matters most once symmetry is ever
C turned on for a scrnc run -- until then it's a cheap no-cost sanity
C gate on every run.
#ifdef _DEBUG_LVL0
         Xbmax = 0.0D0
         Do Ix2 = 1, Nao2
            Do Ix1 = 1, Nao1
               Brute = 0.0D0
               Do Ij = 1, Numj
                  Do Ii = 1, Numi
                     Brute = Brute +
     +                  Eveca(Ioffl1-1+Ix1,Ioffr1-1+Ii) *
     +                  JJ_respn_oo(Ioffoo-1+(Ij-1)*Numi+Ii) *
     +                  Eveca(Ioffl2-1+Ix2,Ioffr2-1+Ij)
                  Enddo
               Enddo
               Xgot = JJ_respn_xx(Ioffxx-1+(Ix2-1)*Nao1+Ix1)
               Xbmax = Max(Xbmax,Abs(Brute-Xgot))
            Enddo
         Enddo
         Write(6,"(a,i2,a,i2,a,D15.6)")
     +      " Scrnc_bcktran_left brute-force check, Irrep_i=",Irrep_i,
     +      " Irrep_j=",Irrep_j," max|brute-xgemm| = ",Xbmax
#endif

         Ioffoo = Ioffoo + Numi*Numj
         Ioffxx = Ioffxx + Nao1*Nao2

      Enddo

      Return
      End
