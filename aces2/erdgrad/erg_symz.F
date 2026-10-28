      Subroutine Erg_symz(G,Gs)
C
C Symmetrize a gradient made from the symmetry unique shell pairs and
C quartets (weighted with their orbit sizes): for an operation R with
C coordinate signs s(k) and atom map p, the image of a contribution
C on atom a is s(k) times it on atom p(a). Averaging over the group
C gives the full gradient (the operations are their own inverses).
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Dimension G(3,Natm), Gs(3,Natm)

      Do Ia = 1, Natm
         Do K = 1, 3
            Sum = 0.0D0
            Do Iop = 0, Nops-1
               Sum = Sum + Dble(Iopsgn(K,Iop))*G(K,Iatmap(Ia,Iop))
            Enddo
            Gs(K,Ia) = Sum/Dble(Nops)
         Enddo
      Enddo

      Return
      End
