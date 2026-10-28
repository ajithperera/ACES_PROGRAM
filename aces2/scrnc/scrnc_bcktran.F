      Subroutine Scrnc_bcktran(Work,Maxmem,Iuhf,JJ_respn_ijkl,Npairs,
     +                         Irrepx,Ispin)

C Transforms a full occ-occ MO pair matrix JJ_respn_ijkl(Npairs,Npairs)
C (spin Ispin, GLOBAL occ-pair index running over ALL irreps -- see
C ordering convention below) to the AO-pair basis (XX,XX), via
C   (IJ,KL) -> (XX,KL) -> (KL,XX) -> (KL,ZW) -> (ZW,KL)transposed=(XX,XX)
C using Scrnc_bcktran_left (the fixed, verified transform kernel, see
C its own header) to transform each pair-member in turn, matching the
C same "transform left pair, transpose, transform right pair" pattern
C aces2/bcktrn/bxaaaa2.f uses for the analogous MO-density backtransform.
C
C GLOBAL occ-pair index ordering convention (2026-08-16): Npairs is
C assumed laid out as consecutive per-irrep blocks in irrep order
C 1..Nirrep, each block of size Irpdpd(Irrep,20+Ispin) (the same
C IRPDPD column convention scrnc_xxxx_pqxx.F already uses for IFLAG=3,
C "NUMIJ=IRPDPD(IRREPIJ,20+ISPIN)"). This is scrnc's own convention
C (no other scrnc routine stores a dense all-irrep occ-pair matrix to
C confirm against), so treat it as a documented assumption to verify
C once a real caller exists, not an established fact.
C
C 2026-08-16, fixed two real bugs found while rewriting this against
C the corrected Scrnc_bcktran_left.F:
C 1. "I040 = Iend" (old line 53) used Iend before it was ever assigned
C    -- undefined-variable read, likely garbage/zero depending on
C    compiler. That whole I040 buffer was never actually referenced
C    again in the routine, so it's removed rather than fixed forward.
C 2. JJ_respn_ijkl(1,Icol) always started reading from GLOBAL row 1 and
C    LOCAL column Icol=1..Nij_cols, for every Irrep_r -- correct only
C    for the first irrep block. Any irrep past the first silently
C    re-read the first irrep's rows/columns instead of its own. Fixed
C    with explicit cumulative per-irrep offsets (Nij_cumoff), computed
C    once before the main loop, indexed by irrep number directly (not
C    tied to loop order) so it stays correct even if Irrepx is ever
C    something other than the totally-symmetric irrep 1 that scrnc
C    always uses today.
C Also added the Ispin argument (threaded into Scrnc_bcktran_left and
C into the IRPDPD/IRPDPDAOMO_OCCBK column selection) -- this used to be
C implicitly alpha-only.

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Maxmem)
      Double Precision JJ_respn_ijkl(Npairs,Npairs)
      Dimension Nij_cumoff(9)

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

#ifdef _DEBUG_LVL0
      Write(6,"(1x,a)") "-----Entering scrnc_bcktran------"
      Write(6,"(a)") " The input JJ response matrx"
      call output(JJ_respn_ijkl,1,Npairs,1,Npairs,Npairs,Npairs,1)
      Write(6,*)
      Write(6,"(a)") " The screened <xx|xx> integrals"
#endif

      Nmo  = Nocco(1)+Nvrto(1)
      Ione = 1
      Call Getrec(20,"JOBARC","NBASTOT",Ione,Nao)

      Nij_cumoff(1) = 0
      Do Irrep = 1, Nirrep
         Nij_cumoff(Irrep+1) = Nij_cumoff(Irrep) +
     +                          Irpdpd(Irrep,20+Ispin)
      Enddo

      Do irrep_r = 1, Nirrep
         Irrep_l = Dirprd(Irrep_r,Irrepx)

         Nij_cols  = Irpdpd(Irrep_r,20+Ispin)
         Nij_rows  = Irpdpd(Irrep_l,20+Ispin)
         Icol_off  = Nij_cumoff(Irrep_r)
         Irow_off  = Nij_cumoff(Irrep_l)
         Nxj_rows  = Irpdpdaomo_occbk(Irrep_l,Ispin)
         Nxj_cols  = Irpdpdaomo_occbk(Irrep_r,Ispin)
         Nxx_rows  = Irpdpdao(Irrep_l)

         I000 = 1
         I010 = I000 + Nxx_rows * Nij_cols
         I020 = I010 + Nao * Nmo
         I030 = I020 + Nxj_rows
         Iend = I030
         if (Iend .GE. Maxmem) Call Insmem("Scrnc_backtrn",Iend,
     +                                      Maxmem)
         Ioff_xx_l = I000
         Do Icol = 1, Nij_cols

            Call Scrnc_bcktran_left(
     +           JJ_respn_ijkl(Irow_off+1,Icol_off+Icol),
     +           Work(Ioff_xx_l),Work(I010),
     +           Work(I020),Nij_rows,Nxj_rows,
     +           Nxx_rows,Nao,Nmo,Irrep_l,Ispin,Iuhf)
            Ioff_xx_l = Ioff_xx_l + Nxx_rows
         Enddo

#ifdef _DEBUG_LVL0
      call checksum("(xx|IJ): ",Work(I000),Nxx_rows*Nij_cols)
#endif

C After completion of the above loop, we should have (xx|pq) batches of
C integrals. Transpose it form (pq|xx) follwed by the loop to form (xx|xx).

         Ncols = Nij_cols
         Nrows = Nxx_rows
         Call Transp(Work(I000),Work(I030),Ncols,Nrows)
         Call Dcopy(Ncols*Nrows,Work(I030),1,Work(I000),1)

         Nxx_cols = Nxx_rows
         Nxx_rows = Irpdpdao(Irrep_r)
         Nxj_rows = Nxj_cols

         I020 = I010 + Nao * Nmo
         I030 = I020 + Nxx_rows * Nxx_cols
         I040 = I030 + Nxj_rows
         Iend = I040
         if (Iend .GE. Maxmem) Call Insmem("Scrnc_backtrn",Iend,
     +                                      Maxmem)
         Ioff_xx_r = I000
         Ioff_xx_l = I020

         Do Icol = 1, Nxx_cols

            Call Scrnc_bcktran_left(Work(Ioff_xx_r),Work(Ioff_xx_l),
     +                              Work(I010),Work(I030),Nij_cols,
     +                              Nxj_rows,Nxx_rows,Nao,Nmo,Irrep_r,
     +                              Ispin,Iuhf)
            Ioff_xx_r = Ioff_xx_r + Nij_cols
            Ioff_xx_l = Ioff_xx_l + Nxx_rows

         Enddo
#ifdef _DEBUG_LVL0
      call checksum("(xx|xx): ",Work(I020),Nxx_rows*Nxx_cols)
#endif

      Enddo

      Return
      End
