      Subroutine Erg_bso(B,N)
C
C The SOs as combinations of ERD functions: SO p = sum_i B(i,p) f_i.
C An SO is a +-1 sum of one VMOL AO on each member shell of its orbit
C (Erd_sobld, Erd_soord). A VMOL AO is FNFAC times the ERD function
C for Cartesian functions (and L <= 1), and sum_m WVSPH(m,j,L) times
C ERD spherical function m for VMOL harmonic j (Erd_vsph).
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Dimension B(N,N), Tmp(13)

      Call Zero(B,N*N)
      Do Ip = 1, Nso
         Io = Isoorb(Ip)
         Do M = 1, Norbm(Io)
            Is = Iorbsh(M,Io)
            B(Ioffn(Is)+Isofn(Ip),Ip) = Dble(Isosgn(M,Ip))
         Enddo
      Enddo

      Do Is = 1, Nshl
         L  = Lshl(Is)
         I0 = Ioffn(Is)
         If (.Not. Spher .Or. L .Le. 1) Then
            Do If = 1, Nfnsh(Is)
               Call Dscal(N,Fnfac(I0+If),B(I0+If,1),N)
            Enddo
         Else
            Ns = 2*L + 1
            Do Ic = 1, Nctsh(Is)
               Ic0 = I0 + (Ic-1)*Ns
               Do Ip = 1, N
                  Do M = 1, Ns
                     Sum = 0.0D0
                     Do J = 1, Ns
                        Sum = Sum + Wvsph(M,J,L)*B(Ic0+J,Ip)
                     Enddo
                     Tmp(M) = Sum
                  Enddo
                  Do M = 1, Ns
                     B(Ic0+M,Ip) = Tmp(M)
                  Enddo
               Enddo
            Enddo
         Endif
      Enddo

      Return
      End
