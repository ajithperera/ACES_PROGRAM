      Subroutine Bwstart(Iuhf)
C
C MR-BW-CC (a bwdata file) runs entirely in Bwcc, which sets iref
C before its own Crapsi, so that MOINTS is the per-reference file
C MOINTS.01, and stops. It must therefore come before xvcc's Crapsi.
C
      Implicit Integer (A-Z)
#include "bwcc.com"

      Call Bwread(.True.,Iuhf)
      If (Isbwcc) Call Bwcc(.False.)

      Return
      End
