      Subroutine Pccd_rsp_gett(T2,Nrow,T)
C T(i,a) = T2(aa,ii) from the irrep-1 block of list 46 (rows A-b by
C irrep of b, A fastest; columns I-j likewise), MOs in Irp order.
      Implicit Double Precision(A-H,O-Z)
#include "syminf.com"
#include "sym.com"
      Dimension T2(Nrow,*)
      Common /Pccdrsp/ Nb,No,Nv,Nov,Npr,Lh,Ljx,Lkx,Le,Lw,Lk,Lde,Ldw,
     +                 Ldk,Lg1,Lg2,Lg3,Lz,Lbm,Lt,Llm,Lmk,Lt1,Ll1,Lvw,
     +                 Lr1,Lr2,La,Lpv,Lsc,Ldv,Ldvl,Ixs,Ievs
      Dimension T(No,Nv)
      Ic0 = 0
      I0  = 0
      Do Iri = 1, Nirrep
         Npi = Pop(Iri,1)
         Ir0 = 0
         Ia0 = 0
         Do Ira = 1, Nirrep
            Nva = Vrt(Ira,1)
            Do Ia = 1, Nva
               Irow = Ir0 + (Ia-1)*Nva + Ia
               Do I = 1, Npi
                  Icol = Ic0 + (I-1)*Npi + I
                  T(I0+I,Ia0+Ia) = T2(Irow,Icol)
               Enddo
            Enddo
            Ir0 = Ir0 + Nva*Nva
            Ia0 = Ia0 + Nva
         Enddo
         Ic0 = Ic0 + Npi*Npi
         I0  = I0 + Npi
      Enddo
      Return
      End
