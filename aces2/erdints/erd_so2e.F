      Subroutine Erd_so2e(Lu,Tla,Lfull,W1,W2,Lw,Zscr,Lzscr,Iscr,Liscr,
     &                    Asym,Nwr)
C
C Two-electron SO integrals from ERD, written in VMOL's format: unit
C Lu(1) gets the (ii|ii) class (IIII), Lu(2) (ii|jj) (IIJJ), Lu(3)
C (ij|ij) (IJIJ) and Lu(4) (ij|kl) (IJKL), all four irreps different.
C Each unique integral is written once, with p <= q, r <= s in its
C label (p,q|r,s); integrals smaller than Tla are dropped.
C
C The loop is over quartets of shell orbits (K1 K2|K3 K4), K1 >= K2,
C K3 >= K4, (K1 K2) >= (K3 K4). For each, ERD gives the AO batches of
C the member shells, and four quarter transformations (Erd_qtr) turn
C every member index into an SO index. Only totally symmetric SO
C quartets are kept.
C
C Lfull = .True. sums over all members of the first orbit, and
C returns in Asym the largest SO integral that symmetry forbids (it
C must vanish: a check of the SO basis). Otherwise only the
C representative of the first orbit is used and the sum is scaled by
C the orbit size, which is exact for the totally symmetric integrals.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"
#include "machsp.com"

      Parameter (Lnbuf = 600, Mxpk = 5000)
      Integer Lu(4), Nwr(4), Iscr(Liscr), Nut(4), Ibuf(Lnbuf,4)
      Integer List(4), Icb(Mxpk), Ice(Mxpk), Kq(4), Ns(4), Nf(4),
     &        Nm(4), Isr(4)
      Dimension W1(Lw), W2(Lw), Zscr(Lzscr), Buf(Lnbuf,4),
     &          Alp(Mxpk), Cc(Mxpk)
      Logical Lfull, Screen
      Data Screen /.True./

      Ipack(I,J,K,L) = Ior(Ior(Ior(I,Ishft(J,Ibitwd)),
     &                 Ishft(K,2*Ibitwd)),Ishft(L,3*Ibitwd))

      Asym = 0.0D0
      Do I = 1, 4
         Nut(I) = 0
         Nwr(I) = 0
      Enddo

      Do K1 = 1, Norb
      Do K2 = 1, K1
      Do K3 = 1, K1
         K4top = K3
         If (K3 .Eq. K1) K4top = K2
      Do K4 = 1, K4top
         Kq(1) = K1
         Kq(2) = K2
         Kq(3) = K3
         Kq(4) = K4
         Do I = 1, 4
            Isr(I) = Iorbsh(1,Kq(I))
            Nf(I)  = Nfnsh(Isr(I))
            Nm(I)  = Norbm(Kq(I))
         Enddo
         If (Lfull) Then
            M1hi = Nm(1)
         Else
            M1hi = 1
         Endif
         Ntot = Nf(1)*M1hi*Nf(2)*Nm(2)*Nf(3)*Nm(3)*Nf(4)*Nm(4)
         If (Ntot .Gt. Lw) Then
            Write(6,"(a,2i12)") " @Erd_so2e: work space too small ",
     &                          Ntot, Lw
            Call Errex
         Endif
         Nd1 = Nf(1)*M1hi
         Nd2 = Nf(2)*Nm(2)
         Nd3 = Nf(3)*Nm(3)
         Nd4 = Nf(4)*Nm(4)

C AO integrals over the member shells, stored as
C W1(f1,m1,f2,m2,f3,m3,f4,m4).

         Do M1 = 1, M1hi
         Do M2 = 1, Nm(2)
         Do M3 = 1, Nm(3)
         Do M4 = 1, Nm(4)
            List(1) = Iorbsh(M1,K1)
            List(2) = Iorbsh(M2,K2)
            List(3) = Iorbsh(M3,K3)
            List(4) = Iorbsh(M4,K4)
            Call Erd_eri(List,Zscr,Lzscr,Iscr,Liscr,Alp,Cc,Icb,Ice,
     &                   Mxpk,Screen,Nbatch,Nfirst)
            Do I4 = 1, Nf(4)
            Do I3 = 1, Nf(3)
            Do I2 = 1, Nf(2)
            Do I1 = 1, Nf(1)
               Iw = I1 + Nf(1)*(M1-1)
     &            + Nd1*(I2-1 + Nf(2)*(M2-1)
     &            + Nd2*(I3-1 + Nf(3)*(M3-1)
     &            + Nd3*(I4-1 + Nf(4)*(M4-1))))
               If (Nbatch .Gt. 0) Then
                  W1(Iw) = Zscr(Nfirst + I1-1 + Nf(1)*(I2-1 +
     &                     Nf(2)*(I3-1 + Nf(3)*(I4-1))))
               Else
                  W1(Iw) = 0.0D0
               Endif
            Enddo
            Enddo
            Enddo
            Enddo
         Enddo
         Enddo
         Enddo
         Enddo

C Members -> SOs for indices 2, 3 and 4. The first index is done
C while the integrals are picked out below: it only needs the signs
C of the members of orbit K1 that were computed (all of them, or the
C representative scaled by the orbit size).

         Call Erd_qtr(W1,W2,Nd1,Nf(2),Nm(2),Nd3*Nd4,K2,1,1.0D0)
         Call Erd_qtr(W2,W1,Nd1*Nd2,Nf(3),Nm(3),Nd4,K3,1,1.0D0)
         Call Erd_qtr(W1,W2,Nd1*Nd2*Nd3,Nf(4),Nm(4),1,K4,1,1.0D0)
         Wt = Dble(Nm(1))/Dble(M1hi)

C Pick out the unique, totally symmetric SO integrals. W2 is now
C (f1,m1,f2,j2,f3,j3,f4,j4), f running fastest.

         Do J4 = 1, Nm(4)
         Do I4 = 1, Nf(4)
            Is = Lorbso(Iorbso(K4) + (I4-1)*Nm(4) + J4)
         Do J3 = 1, Nm(3)
         Do I3 = 1, Nf(3)
            Ir = Lorbso(Iorbso(K3) + (I3-1)*Nm(3) + J3)
         Do J2 = 1, Nm(2)
         Do I2 = 1, Nf(2)
            Iq = Lorbso(Iorbso(K2) + (I2-1)*Nm(2) + J2)
            Ib = Nd1*(I2-1 + Nf(2)*(J2-1) + Nd2*(I3-1 + Nf(3)*(J3-1)
     &         + Nd3*(I4-1 + Nf(4)*(J4-1))))
         Do J1 = 1, Nm(1)
         Do I1 = 1, Nf(1)
            Ip = Lorbso(Iorbso(K1) + (I1-1)*Nm(1) + J1)
            X = 0.0D0
            Do M1 = 1, M1hi
               X = X + Dble(Isosgn(M1,Ip))*W2(Ib + I1 + Nf(1)*(M1-1))
            Enddo
            X = X*Wt
            Ixor = Ieor(Ieor(Isoirr(Ip)-1,Isoirr(Iq)-1),
     &                  Ieor(Isoirr(Ir)-1,Isoirr(Is)-1))
            If (Ixor .Ne. 0) Then
               Asym = Max(Asym,Dabs(X))
               Goto 10
            Endif
            If (K1 .Eq. K2 .And. Ip .Lt. Iq) Goto 10
            If (K3 .Eq. K4 .And. Ir .Lt. Is) Goto 10
            Ipq = Max(Ip,Iq)*(Max(Ip,Iq)-1)/2 + Min(Ip,Iq)
            Irs = Max(Ir,Is)*(Max(Ir,Is)-1)/2 + Min(Ir,Is)
            If (K1 .Eq. K3 .And. K2 .Eq. K4 .And. Ipq .Lt. Irs) Goto 10
            If (Dabs(X) .Lt. Tla) Goto 10

            I = Min(Ip,Iq)
            J = Max(Ip,Iq)
            K = Min(Ir,Is)
            L = Max(Ir,Is)
            Ii = Isoirr(I)
            Ij = Isoirr(J)
            Ik = Isoirr(K)
            Il = Isoirr(L)
            If (Ii .Eq. Ij .And. Ik .Eq. Il) Then
               If (Ii .Eq. Ik) Then
                  Ifl = 1
               Else
                  Ifl = 2
               Endif
            Else If (Ii .Eq. Ik .And. Ij .Eq. Il) Then
               Ifl = 3
            Else
               Ifl = 4
            Endif
            Nut(Ifl) = Nut(Ifl) + 1
            Buf(Nut(Ifl),Ifl)  = X
            Ibuf(Nut(Ifl),Ifl) = Ipack(I,J,K,L)
            Nwr(Ifl) = Nwr(Ifl) + 1
            If (Nut(Ifl) .Eq. Lnbuf) Then
               Write(Lu(Ifl)) (Buf(Ii,Ifl),Ii=1,Lnbuf),
     &                        (Ibuf(Ii,Ifl),Ii=1,Lnbuf), Nut(Ifl)
               Nut(Ifl) = 0
            Endif
 10         Continue
         Enddo
         Enddo
         Enddo
         Enddo
         Enddo
         Enddo
         Enddo
         Enddo
      Enddo
      Enddo
      Enddo
      Enddo

      Do Ifl = 1, 4
         If (Lu(Ifl) .Gt. 0) Then
            If (Nut(Ifl) .Ne. 0) Write(Lu(Ifl))
     &         (Buf(Ii,Ifl),Ii=1,Lnbuf), (Ibuf(Ii,Ifl),Ii=1,Lnbuf),
     &         Nut(Ifl)
            Nut(Ifl) = -1
            Write(Lu(Ifl)) (Buf(Ii,Ifl),Ii=1,Lnbuf),
     &                     (Ibuf(Ii,Ifl),Ii=1,Lnbuf), Nut(Ifl)
         Else If (Nwr(Ifl) .Ne. 0) Then
            Write(6,"(a,i2)") " @Erd_so2e: integrals for closed file",
     &                        Ifl
            Call Errex
         Endif
      Enddo

      Return
      End
