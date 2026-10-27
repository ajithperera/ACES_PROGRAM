      Subroutine Erg_gnsp(List,Fac,Gam,Bx,N,Ng,Irg,Ncol,Icol0,Lcol,
     &                    Ipx,Npr,Ioffg,Gso,W1,W2,Lw)
C
C Add Fac times the nonseparable two-particle density of one shell
C quartet, over ERD functions, to Gam (N1 x N2 x N3 x N4, first index
C fastest):
C
C    Gam(i,j,k,l) += Fac sum Bx(i,p) Bx(j,q) Bx(k,r) Bx(l,s) Gam(pq,rs)
C
C Gam(pq,rs) from Erg_gam over the Ng SOs of the lists, Bx their
C expansion in ERD functions. The functions of shell s only enter the
C Ncol(s) SOs Lcol(Icol0(s)+1 ...) (Erg_bcol). The SO block is
C gathered into W1 and transformed one index at a time.
C
      Implicit Double Precision (A-H,O-Z)
#include "mxatms.par"
#include "erdints.com"

      Integer List(4), Irg(Ng), Ncol(*), Icol0(*), Lcol(*),
     &        Ipx(Ng,Ng), Npr(8), Ioffg(8), Nf(4), Nb(4), Io(4), I0(4)
      Dimension Gam(*), Bx(N,Ng), Gso(*), W1(Lw), W2(Lw)

      Do I = 1, 4
         Is    = List(I)
         Nf(I) = Nfnsh(Is)
         Nb(I) = Ncol(Is)
         Io(I) = Icol0(Is)
         I0(I) = Ioffn(Is)
      Enddo
      Nso4 = Nb(1)*Nb(2)*Nb(3)*Nb(4)
      If (2*Nso4 .Gt. Lw) Then
         Write(6,"(a)") " @Erg_gnsp: SO quartet buffer too small"
         Call Errex
      Endif

C Gather the SO block, W1(p,q,r,s).

      Ind = 0
      Do L = 1, Nb(4)
         Is4 = Lcol(Io(4)+L)
         Do K = 1, Nb(3)
            Is3 = Lcol(Io(3)+K)
            Igrs = Ieor(Irg(Is3)-1,Irg(Is4)-1) + 1
            Irs  = Ipx(Is3,Is4)
            Do J = 1, Nb(2)
               Is2 = Lcol(Io(2)+J)
               Do I = 1, Nb(1)
                  Is1 = Lcol(Io(1)+I)
                  Ind = Ind + 1
                  Igpq = Ieor(Irg(Is1)-1,Irg(Is2)-1) + 1
                  If (Igpq .Eq. Igrs) Then
                     W1(Ind) = Gso(Ioffg(Igpq) +
     &                         (Ipx(Is1,Is2)-1)*Npr(Igpq) + Irs)
                  Else
                     W1(Ind) = 0.0D0
                  Endif
               Enddo
            Enddo
         Enddo
      Enddo

C B blocks of the four shells go to the end of W2 (rows I0+1..,
C columns the SOs they enter); the transformed blocks alternate between
C W1 and W2.
C  1: W2(i,qrs)   = B1(i,p) W1(p,qrs)

      Lb  = Lw
      Ib1 = Lb - Nf(1)*Nb(1) + 1
      Call Erg_bblk(Bx,N,Ng,Lcol,I0(1),Nf(1),Io(1),Nb(1),W2(Ib1))
      M = Nb(2)*Nb(3)*Nb(4)
      Call Dgemm('N','N',Nf(1),M,Nb(1),1.0D0,W2(Ib1),Nf(1),W1,Nb(1),
     &           0.0D0,W2,Nf(1))

C  2: W1(i,j,rs)  = W2(i,q,rs) B2(j,q)

      Ib2 = Lw - Nf(2)*Nb(2) + 1
      Call Erg_bblk(Bx,N,Ng,Lcol,I0(2),Nf(2),Io(2),Nb(2),W1(Ib2))
      M  = Nb(3)*Nb(4)
      Do Irs = 1, M
         Call Dgemm('N','T',Nf(1),Nf(2),Nb(2),1.0D0,
     &              W2(1+(Irs-1)*Nf(1)*Nb(2)),Nf(1),W1(Ib2),Nf(2),
     &              0.0D0,W2(1+Lw/2+(Irs-1)*Nf(1)*Nf(2)),Nf(1))
      Enddo
      Nij = Nf(1)*Nf(2)
      Call Dcopy(Nij*M,W2(1+Lw/2),1,W1,1)

C  3: W2(ij,k,s)  = W1(ij,r,s) B3(k,r)

      Ib3 = Lw - Nf(3)*Nb(3) + 1
      Call Erg_bblk(Bx,N,Ng,Lcol,I0(3),Nf(3),Io(3),Nb(3),W2(Ib3))
      Do Is = 1, Nb(4)
         Call Dgemm('N','T',Nij,Nf(3),Nb(3),1.0D0,
     &              W1(1+(Is-1)*Nij*Nb(3)),Nij,W2(Ib3),Nf(3),
     &              0.0D0,W2(1+(Is-1)*Nij*Nf(3)),Nij)
      Enddo

C  4: Gam(ijk,l) += Fac W2(ijk,s) B4(l,s)

      Ib4 = Lw - Nf(4)*Nb(4) + 1
      Call Erg_bblk(Bx,N,Ng,Lcol,I0(4),Nf(4),Io(4),Nb(4),W1(Ib4))
      Nijk = Nij*Nf(3)
      Call Dgemm('N','T',Nijk,Nf(4),Nb(4),Fac,W2,Nijk,W1(Ib4),Nf(4),
     &           1.0D0,Gam,Nijk)

      Return
      End
