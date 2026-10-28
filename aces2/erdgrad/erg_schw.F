      Subroutine Erg_schw(Q,Fd,Zscr,Lzscr,Iscr,Liscr)
C
C Bounds for screening the derivative ERIs (Erg_two). For every shell
C pair (s1 >= s2, Q(s1*(s1-1)/2+s2)) the Schwarz factor
C
C    Q = sqrt( max |(ij|i'j')| ),  i, i' in s1 and j, j' in s2,
C
C the largest element of the (s1 s2|s1 s2) batch (Erd_eri), which is
C not smaller than the largest diagonal one; and for every shell Fd =
C 2 sqrt(largest exponent), the size of a derivative of its functions
C relative to the functions themselves.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Parameter (Mxpk = 5000)
      Dimension Q(*), Fd(Nshl), Zscr(Lzscr), Alp(Mxpk), Cc(Mxpk)
      Integer Iscr(Liscr), List(4), Icb(Mxpk), Ice(Mxpk)
      Logical Screen
      Data Screen /.False./

      Do Is = 1, Nshl
         Amax = 0.0D0
         Do Ip = 1, Nprsh(Is)
            Amax = Max(Amax,Alpha(Iofalp(Is)+Ip))
         Enddo
         Fd(Is) = 2.0D0*Dsqrt(Amax)
      Enddo

      Do Is1 = 1, Nshl
      Do Is2 = 1, Is1
         List(1) = Is1
         List(2) = Is2
         List(3) = Is1
         List(4) = Is2
         Call Erd_eri(List,Zscr,Lzscr,Iscr,Liscr,Alp,Cc,Icb,Ice,Mxpk,
     &                Screen,Nbatch,Nfirst)
         X = 0.0D0
         If (Nbatch .Gt. 0) Then
            Do I = 0, Nfnsh(Is1)**2*Nfnsh(Is2)**2-1
               X = Max(X,Dabs(Zscr(Nfirst+I)))
            Enddo
         Endif
         Q(Is1*(Is1-1)/2+Is2) = Dsqrt(X)
      Enddo
      Enddo

      Return
      End
