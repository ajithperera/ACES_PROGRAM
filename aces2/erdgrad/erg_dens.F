      Subroutine Erg_dens(Iuhf,N,B,Da,Db,W,C,Occ,Eps,Scr)
C
C SCF densities over ERD functions. For each spin the SO density and
C energy weighted density are
C
C    D(SO) = C n C(T),    W(SO) = C n e C(T)
C
C (n the occupation numbers ORB_OCC*, e the eigenvalues SCFEVL*), and
C over ERD functions D = B D(SO) B(T) (Erg_bso). Da and Db are the
C alpha and beta densities (Db = Da for RHF); W is the total energy
C weighted density.
C
      Implicit Double Precision (A-H,O-Z)

#include "machsp.com"

      Dimension B(N,N), Da(N,N), Db(N,N), W(N,N), C(N,N), Occ(N),
     &          Eps(N), Scr(N,N,2)
      Character*8 Labc(2), Labe(2), Labo(2)
      Data Labc /'SCFEVCA0','SCFEVCB0'/
      Data Labe /'SCFEVLA0','SCFEVLB0'/
      Data Labo /'ORB_OCCA','ORB_OCCB'/

      Call Zero(W,N*N)
      Do Ispin = 1, Iuhf+1
         Call Getrec(20,'JOBARC',Labc(Ispin),N*N*Iintfp,C)
         Call Getrec(20,'JOBARC',Labe(Ispin),N*Iintfp,Eps)
         Call Getrec(20,'JOBARC',Labo(Ispin),N*Iintfp,Occ)
         Fac = 1.0D0
         If (Iuhf .Eq. 0) Fac = 2.0D0

C D(SO) and W(SO) into Scr(,,1), then over ERD functions.

         Do Ipass = 1, 2
            Do J = 1, N
               If (Ipass .Eq. 1) Then
                  X = Occ(J)
               Else
                  X = Fac*Occ(J)*Eps(J)
               Endif
               Do I = 1, N
                  Scr(I,J,2) = X*C(I,J)
               Enddo
            Enddo
            Call Dgemm('N','T',N,N,N,1.0D0,Scr(1,1,2),N,C,N,0.0D0,
     &                 Scr(1,1,1),N)
            Call Dgemm('N','N',N,N,N,1.0D0,B,N,Scr(1,1,1),N,0.0D0,
     &                 Scr(1,1,2),N)
            If (Ipass .Eq. 1) Then
               If (Ispin .Eq. 1) Then
                  Call Dgemm('N','T',N,N,N,1.0D0,Scr(1,1,2),N,B,N,
     &                       0.0D0,Da,N)
               Else
                  Call Dgemm('N','T',N,N,N,1.0D0,Scr(1,1,2),N,B,N,
     &                       0.0D0,Db,N)
               Endif
            Else
               Call Dgemm('N','T',N,N,N,1.0D0,Scr(1,1,2),N,B,N,
     &                    1.0D0,W,N)
            Endif
         Enddo
      Enddo
      If (Iuhf .Eq. 0) Call Dcopy(N*N,Da,1,Db,1)

      Return
      End
