      Program Aux

      Implicit None
      Integer Ndim,Iuhf

      Parameter(Ndim=14)
      Double Precision G2(Ndim,Ndim,Ndim,Ndim)

#include "machsp.com"
#include "istart.com"
#include "icore.com"

      call aces_init(Icore,I0,Icrsiz,Iuhf,.false.)
      Call Get2eints("TWOELSUP",G2,Ndim)

      Stop
      End
