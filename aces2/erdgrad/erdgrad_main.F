      Program Erdgrad_main
C
C xerdgrad: ACES II analytic gradients from Flocke's OED/ERD
C derivative integrals instead of VMOL (xvdint). It uses the basis,
C symmetry and SO basis of xerdints (liberdints), turns the SO
C densities into densities over ERD functions and contracts them with
C the derivative integrals of the symmetry unique shell pairs and
C quartets. Only gradients; second derivatives and properties stay
C with xvdint.
C Ajith Perera, 10/2026.
C
      Implicit None
      Integer Iuhf

#include "icore.com"
#include "istart.com"
#include "machsp.com"

      Call Aces_init(Icore,I0,Icrsiz,Iuhf,.True.)
      Call Erg_driver(Icore(I0),Icrsiz,Iuhf)
      Call Aces_fin

      Stop
      End
