      Subroutine Pccd_rsp_save(G,Hh,Wy,Nbas)
C After Pccd_ohess_asm: h and the two integral classes of the response.
      Implicit Double Precision(A-H,O-Z)
      Dimension G(Nbas,Nbas,Nbas,Nbas),Hh(Nbas,Nbas),
     +          Wy(Nbas,Nbas,Nbas)
#include "machsp.com"
      N3 = Nbas*Nbas*Nbas
      Do Iy = 1, Nbas
         Do Iq = 1, Nbas
            Do Ip = 1, Nbas
               Wy(Ip,Iq,Iy) = G(Ip,Iq,Iy,Iy)
            Enddo
         Enddo
      Enddo
      Call Putrec(20,"JOBARC","RSP_JX  ",N3*Iintfp,Wy)
      Do Iy = 1, Nbas
         Do Iq = 1, Nbas
            Do Ip = 1, Nbas
               Wy(Ip,Iq,Iy) = G(Ip,Iy,Iq,Iy)
            Enddo
         Enddo
      Enddo
      Call Putrec(20,"JOBARC","RSP_KX  ",N3*Iintfp,Wy)
      Call Putrec(20,"JOBARC","RSP_HMO ",Nbas*Nbas*Iintfp,Hh)
      Return
      End
