#include "flags.h"
      Program Oorot
C
C Pilot utility for the OO-pCCD exact-Newton (finite-difference Hessian)
C test. The mode is read from the first line of oorot.in:
C   SAVE  : SCFEVCA0 (AO x MO)            -> c0.dat
C   APPLY : C = C0*U, C0 from c0.dat,
C           U (MO x MO) from u.dat         -> SCFEVCA0, SCFEVCB0
C   DUMP  : OBRTGRDA (MO x MO), TOTENERG   -> grad.dat
C All text files are column-major with a "nrow ncol" header line.
C
      Implicit Double Precision (A-H,O-Z)
      Character*8 Mode
#include "icore.com"
#include "istart.com"
#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "info.com"
      Call Aces_init(icore, i0, icrsiz, iuhf, .true.)
      Open(Unit=31,File="oorot.in",Status="Old")
      Read(31,"(a)") Mode
      Close(31)
      Call Getrec(20,"JOBARC","NBASTOT ",1,Nao)
      Nmo  = Nocco(1) + Nvrto(1)
      I000 = i0
      I010 = I000 + Nao*Nmo*Iintfp
      I020 = I010 + Nmo*Nmo*Iintfp
      Iend = I020 + Nao*Nmo*Iintfp
      If (Iend-i0 .Gt. icrsiz) Call Insmem("oorot",Iend-i0,icrsiz)
      Call Oorot_do(icore(I000),icore(I010),icore(I020),Nao,Nmo,Mode)
      Call Aces_fin
      Stop
      End

