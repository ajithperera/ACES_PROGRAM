      Subroutine Runpccd(Iuhf)

      Implicit None
      Integer*8 Iuhf

#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "jobarc.com"

C Init() has already run Joda_ and aces_init, as for Runcc. Single-point
C pCCD as xrunpccd's HF/KS branch runs it: xvmol, xvmol2ja, xvscf, xvtran,
C xintprc, xvcc. Neither xrunpccd nor xaces2 runs xlambda for it (Vlambda
C here stops on a missing list), and the orbital-optimized variants
C (PROGRAM=PCCD) are xrunpccd's own loop, not reachable from here.

      Call Ints_(Icore(i0),Icrsiz,Iuhf)
      Call V2ja(icore(i0),Icrsiz,Iuhf)
      Call Vscf(Icore(i0),Icrsiz,Iuhf)
      Call aces_init_chemsys()
      Call Vtran(icore(i0),Icrsiz,Iuhf)
      Call aces_init_chemsys()
      Call Init_legacies()
      Call Intproc(icore(i0),Icrsiz,Iuhf)
      Call Vcc(icore(i0),Icrsiz,Iuhf)

C Do not call aces_fin (use might be intending to run something
C else after external manupilations).

CSSS      Call aces_fin
CSSS      Call C_free(Icore)
 
      Return
      End

