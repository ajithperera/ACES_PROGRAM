      Subroutine Erd_so1e(A,N,B,Asym)
C
C AO -> SO transformation of a symmetric one-electron matrix: A is
C over the basis functions (N = Nfnerd, Erd_ao1e), B(Nso,Nso) over
C the SOs,
C
C    B(p,q) = sum_m,n s(m,p) s(n,q) A(mu(m,p),mu(n,q))
C
C with mu(m,p) the function of SO p on member m of its orbit. The full
C B is returned so that the caller can check that it is block
C diagonal; Asym is the largest element coupling different irreps (it
C must vanish).
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Dimension A(N,N), B(Nso,Nso)

      Asym = 0.0D0
      Do Ip = 1, Nso
         Iop = Isoorb(Ip)
         Do Iq = 1, Ip
            Ioq = Isoorb(Iq)
            Sum = 0.0D0
            Do M = 1, Norbm(Iop)
               Mu = Ioffn(Iorbsh(M,Iop)) + Isofn(Ip)
               Do K = 1, Norbm(Ioq)
                  Nu = Ioffn(Iorbsh(K,Ioq)) + Isofn(Iq)
                  Sum = Sum + Dble(Isosgn(M,Ip)*Isosgn(K,Iq))*A(Mu,Nu)
               Enddo
            Enddo
            B(Ip,Iq) = Sum
            B(Iq,Ip) = Sum
            If (Isoirr(Ip) .Ne. Isoirr(Iq))
     &         Asym = Max(Asym,Dabs(Sum))
         Enddo
      Enddo

      Return
      End
