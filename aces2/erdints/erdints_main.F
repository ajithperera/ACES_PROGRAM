      Program Erdints_main
C
C xerdints: ACES II integrals from Flocke's OED/ERD instead of VMOL.
C It replaces xvmol + xvmol2ja: it builds the SO basis itself, writes
C the one- and two-electron SO integrals in the VMOL file format
C (IIII, IIJJ, IJIJ, IJKL) and the JOBARC records vmol2ja writes.
C Ajith Perera, 10/2026.
C
      Implicit None
      Integer Iuhf

#include "icore.com"
#include "istart.com"
#include "machsp.com"

      Call Aces_init(Icore,I0,Icrsiz,Iuhf,.True.)
      Call Erd_driver(Icore(I0),Icrsiz)
      Call Aces_fin

      Stop
      End
