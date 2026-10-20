      Subroutine Pccd_pair_putt(T,T2,Nrow)
C T2(aa,ii) = T(i,a) in the irrep-1 block of list 46 (layout of
C Pccd_rsp_gett), every other element zero.
      Implicit Double Precision(A-H,O-Z)
#include "syminf.com"
#include "sym.com"
#include "sympop.com"
      Dimension T2(Nrow,*)
      Common /Pccdrsp/ Nb,No,Nv,Nov,Npr,Lh,Ljx,Lkx,Le,Lw,Lk,Lde,Ldw,
     +                 Ldk,Lg1,Lg2,Lg3,Lz,Lbm,Lt,Llm,Lmk,Lt1,Ll1,Lvw,
     +                 Lr1,Lr2,La,Lpv,Lsc,Ldv,Ldvl,Ixs,Ievs
      Dimension T(No,Nv)
      Ncol = Irpdpd(1,14)
      Call Dzero(T2,Nrow*Ncol)
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
                  T2(Irow,Icol) = T(I0+I,Ia0+Ia)
               Enddo
            Enddo
            Ir0 = Ir0 + Nva*Nva
            Ia0 = Ia0 + Nva
         Enddo
         Ic0 = Ic0 + Npi*Npi
         I0  = I0 + Npi
      Enddo
      Call Putlst(T2,1,Ncol,1,1,46)
      Return
      End
