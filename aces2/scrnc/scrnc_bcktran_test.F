#include "flags.h"
      Subroutine Scrnc_bcktran_test(Work,Maxmem,Iuhf)

C Temporary, self-contained validation harness for Scrnc_bcktran_left
C (2026-08-16) -- NOT part of the physics pipeline, no caller feeds it
C real data yet (scrnc_bcktran/scrnc_bcktran_left have no wired-in
C consumer -- see SCRNC_DEBUG_LOG.md). Builds a synthetic, deterministic
C occ-occ MO pair "response" matrix, calls Scrnc_bcktran_left on it, and
C relies on that routine's own internal _DEBUG_LVL0 brute-force
C cross-check (Xgemm result vs direct nested-loop evaluation) to confirm
C the fixed transform kernel is self-consistent. Called from
C Scrnc_driver.F, guarded by _DEBUG_LVL0, right after Scrnc_getaoinf so
C the AOSYM/SYM/SYMPOP common blocks are already populated. Safe to
C delete once a real caller exists and exercises the same code path
C with genuine data.

      Implicit Double Precision (A-H,O-Z)
      Dimension Work(Maxmem)

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
#include "info.com"

      Nmo = Nocco(1)+Nvrto(1)
      Ione = 1
      Call Getrec(20,"JOBARC","NBASTOT",Ione,Nao)

      Write(6,*)
      Write(6,"(a)") " ---- Scrnc_bcktran_left validation harness ----"

      Do Ispin = 1, Iuhf+1
         Nij_rows = Irpdpd(1,20+Ispin)*Irpdpd(1,20+Ispin)
         Nxx_rows = Irpdpdao(1)
         Nxj_rows = Irpdpdaomo_occbk(1,Ispin)

         I000 = 1
         I010 = I000 + Nij_rows
         I020 = I010 + Nxx_rows
         I030 = I020 + Nao*Nmo
         I040 = I030 + Max(Nxj_rows,Nao*Nmo)
         If (I040 .GT. Maxmem) Call Insmem("Scrnc_bcktran_test",
     +                                      I040,Maxmem)

         Do I = 1, Nij_rows
            Work(I000-1+I) = 0.01D0*Dble(I) - 0.5D0*Dble(Ispin)
         Enddo

         Write(6,"(a,i2,a,i4,a,i4)") " Ispin=",Ispin," Nij_rows=",
     +                                Nij_rows," Nxx_rows=",Nxx_rows

         Call Scrnc_bcktran_left(Work(I000),Work(I010),Work(I020),
     +                           Work(I030),Nij_rows,Nxj_rows,
     +                           Nxx_rows,Nao,Nmo,1,Ispin,Iuhf)
      Enddo

      Write(6,"(a)") " ---- end Scrnc_bcktran_left validation ----"
      Write(6,*)

      Return
      End
